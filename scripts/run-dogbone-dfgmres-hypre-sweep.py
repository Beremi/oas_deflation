#!/usr/bin/env python3
"""Run Dogbone DFGMRES with hypre as a right preconditioner over basis sizes."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "data/cases/Dogbone"
MASTER = CASE_DIR / "master.inp"
RUNNER = ROOT / "scripts/run-oas-profile.sh"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str | None, default: float = math.nan) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def as_int(value: str | None, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except ValueError:
        return default


def md_table(headers: list[str], rows: list[list[object]]) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(format_cell(item) for item in row) + " |")
    return "\n".join(lines)


def format_cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if abs(value) >= 1000 or (0 < abs(value) < 1e-3):
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)


def read_solver_out(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip().lstrip("#").split("\t")
        rows = []
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, values)))
    return rows


def run_variant(args: argparse.Namespace, out_root: Path, name: str, nvec: int, collect: bool, prefix: bool, eps: str, diagnostic: bool = False) -> dict[str, object]:
    out_dir = out_root / name
    title = f"Dogbone DFGMRES hypre {'diagnostic ' if diagnostic else ''}{'prefix ' if prefix else ''}N={nvec} eps={eps} Linear-Solve Profile"
    overrides = [
        "linear_solver_profile_matrix_delta=1",
        "runtime_phase_profile=1",
        "runtime_phase_profile_file=runtime_profile",
        "dfgmres_tolerance=1e-6",
        "dfgmres_true_tolerance=1e-6",
        "dfgmres_max_iterations=500",
        "dfgmres_restart=80",
        f"dfgmres_deflation_vectors={nvec}",
        f"dfgmres_deflation_eps={eps}",
        f"dfgmres_collect_newton_steps={1 if collect else 0}",
        "dfgmres_preconditioner=hypre",
        "dfgmres_reorthogonalize_on_matrix_change=1",
        "amgcl_near_nullspace=1",
        "amgcl_elastic_full_lift=1",
        "amgcl_use_block_backend=1",
        "hypre_coarsen_type=8",
        "hypre_interp_type=6",
        "hypre_strong_threshold=0.25",
        "hypre_nodal=4",
        "hypre_relax_type=6",
        "hypre_p_max=4",
        "hypre_boomer_max_iterations=1",
        "hypre_use_dof_functions=0",
        "hypre_use_interp_vectors=0",
    ]
    if prefix:
        overrides.extend([
            "total_time=0.01",
            "max_iterations=5",
            "min_iterations=5",
            "limit_tolerance=10",
            "linear_solver_replay_dump=0",
        ])

    command = [
        str(RUNNER),
        str(CASE_DIR),
        str(MASTER),
        str(args.oas_bin),
        str(args.threads),
        "DeflatedFGMRES",
        title,
        "Dogbone",
        str(out_dir),
        *overrides,
    ]
    env = os.environ.copy()
    env.setdefault("MKLROOT", args.mklroot)
    env["OAS_TIMEOUT"] = args.prefix_timeout if prefix else args.full_timeout
    result = subprocess.run(command, cwd=ROOT, env=env, text=True)
    return summarize_run(name, nvec, collect, prefix, eps, diagnostic, result.returncode, out_dir)


def summarize_run(name: str, nvec: int, collect: bool, prefix: bool, eps: str, diagnostic: bool, returncode: int, out_dir: Path) -> dict[str, object]:
    events = read_tsv(out_dir / "linear_profile_events.tsv")
    solver_rows = read_solver_out(out_dir / "solver.out")
    solve_events = [row for row in events if row.get("phase") == "solve"]
    factorize_events = [row for row in events if row.get("phase") == "factorize"]
    iterations = [as_int(row.get("solver_iterations"), -1) for row in solve_events if as_int(row.get("solver_iterations"), -1) >= 0]
    residuals = [as_float(row.get("solver_error")) for row in solve_events if as_float(row.get("solver_error")) >= 0]
    successes = [as_int(row.get("success"), 0) for row in solve_events]
    max_step = max([as_int(row.get("step"), 0) for row in solver_rows], default=0)
    nonlinear_iterations = sum(as_int(row.get("iterations"), 0) for row in solver_rows)
    total_wall = max([as_float(row.get("elapsed_time"), 0.0) for row in solver_rows], default=sum(as_float(row.get("duration_seconds"), 0.0) for row in events))
    solve_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in solve_events)
    setup_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in factorize_events)
    analyze_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in events if row.get("phase") == "analyze")
    precond_apply = sum(as_float(row.get("preconditioner_apply_seconds"), 0.0) for row in solve_events)
    orth = sum(as_float(row.get("orthogonalization_seconds"), 0.0) for row in solve_events)
    least_squares = sum(as_float(row.get("least_squares_seconds"), 0.0) for row in solve_events)
    matvec = sum(as_float(row.get("matvec_seconds"), 0.0) for row in solve_events)
    deflation = sum(as_float(row.get("deflation_seconds"), 0.0) for row in solve_events)
    basis_values = [as_int(row.get("deflation_basis_size"), 0) for row in solve_events]
    log_basis_values: list[int] = []
    log_path = out_dir / "oas.log"
    if log_path.exists():
        for match in re.finditer(r"active_basis_size=(\d+)", log_path.read_text(encoding="utf-8", errors="replace")):
            log_basis_values.append(int(match.group(1)))
    max_basis = max(basis_values + log_basis_values, default=0)
    final_basis = log_basis_values[-1] if log_basis_values else (basis_values[-1] if basis_values else 0)
    discarded = max([as_int(row.get("deflation_discarded_count"), 0) for row in solve_events], default=0)
    raw_candidates = max([as_int(row.get("deflation_raw_candidate_count"), 0) for row in solve_events], default=0)
    capacity_evictions = max([as_int(row.get("deflation_capacity_eviction_count"), 0) for row in solve_events], default=0)
    low_anorm_discards = max([as_int(row.get("deflation_low_anorm_discard_count"), 0) for row in solve_events], default=0)
    max_offdiag = max([as_float(row.get("deflation_orthogonality_max_offdiag"), 0.0) for row in solve_events], default=0.0)
    max_diag_error = max([as_float(row.get("deflation_orthogonality_max_diag_error"), 0.0) for row in solve_events], default=0.0)
    last_final_anorm = next((as_float(row.get("deflation_last_final_anorm"), math.nan) for row in reversed(solve_events) if row.get("deflation_last_final_anorm")), math.nan)
    value_hashes = {row.get("value_hash", "") for row in factorize_events if row.get("value_hash")}
    pattern_hashes = {row.get("pattern_hash", "") for row in factorize_events if row.get("pattern_hash")}

    required_solves = 5 if prefix else 1
    completed = returncode == 0 and (prefix or max_step >= 100)
    linear_ok = len(solve_events) >= required_solves and all(successes) and (max(residuals, default=math.inf) <= 1e-6)
    basis_ok = (not collect) or nvec == 0 or max_basis > 1
    orthogonality_ok = max_basis <= 1 or (max_offdiag <= 1e-6 and max_diag_error <= 1e-6)
    if completed and linear_ok and basis_ok and orthogonality_ok:
        status = "passed"
    elif linear_ok and not basis_ok:
        status = "basis_failed"
    elif linear_ok and not orthogonality_ok:
        status = "orthogonality_failed"
    elif linear_ok:
        status = "linear_prefix_passed" if prefix else "linear_ok_partial"
    elif solve_events:
        status = "linear_failed"
    else:
        status = "no_profile"

    return {
        "name": name,
        "N": nvec,
        "collect": collect,
        "prefix": prefix,
        "eps": eps,
        "diagnostic": diagnostic,
        "returncode": returncode,
        "status": status,
        "out_dir": out_dir,
        "steps": max_step,
        "nonlinear_iterations": nonlinear_iterations,
        "linear_solves": len(solve_events),
        "unique_patterns": len(pattern_hashes),
        "unique_matrices": len(value_hashes),
        "total_outer_iterations": sum(iterations),
        "median_iterations": statistics.median(iterations) if iterations else math.nan,
        "max_iterations": max(iterations, default=math.nan),
        "max_true_residual": max(residuals, default=math.nan),
        "max_basis": max_basis,
        "final_basis": final_basis,
        "discarded": discarded,
        "raw_candidates": raw_candidates,
        "capacity_evictions": capacity_evictions,
        "low_anorm_discards": low_anorm_discards,
        "max_offdiag": max_offdiag,
        "max_diag_error": max_diag_error,
        "last_final_anorm": last_final_anorm,
        "wall_seconds": total_wall,
        "analyze_seconds": analyze_seconds,
        "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds,
        "preconditioner_apply_seconds": precond_apply,
        "orthogonalization_seconds": orth,
        "least_squares_seconds": least_squares,
        "matvec_seconds": matvec,
        "deflation_seconds": deflation,
        "nonlinear_overhead_seconds": max(total_wall - analyze_seconds - setup_seconds - solve_seconds, 0.0),
        "runtime_hotspots": top_runtime_hotspots(out_dir),
    }


def top_runtime_hotspots(out_dir: Path, limit: int = 8) -> list[dict[str, object]]:
    rows = read_tsv(out_dir / "runtime_profile_summary.tsv")
    parsed: list[dict[str, object]] = []
    for row in rows:
        parsed.append({
            "phase": row.get("phase", ""),
            "detail": row.get("detail", ""),
            "count": as_int(row.get("count"), 0),
            "total_seconds": as_float(row.get("total_seconds"), 0.0),
            "mean_seconds": as_float(row.get("mean_seconds"), 0.0),
            "max_seconds": as_float(row.get("max_seconds"), 0.0),
        })
    parsed.sort(key=lambda item: float(item["total_seconds"]), reverse=True)
    return parsed[:limit]


def make_plots(out_root: Path, run_summaries: list[dict[str, object]]) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []

    plots: list[str] = []
    completed = [row for row in run_summaries if not row["prefix"] and not row["diagnostic"] and Path(row["out_dir"], "linear_profile_events.tsv").exists()]
    if completed:
        fig, ax = plt.subplots(figsize=(10, 4.5))
        for row in completed:
            events = read_tsv(Path(row["out_dir"]) / "linear_profile_events.tsv")
            solve = [event for event in events if event.get("phase") == "solve"]
            ax.plot(range(1, len(solve) + 1), [as_int(event.get("solver_iterations"), 0) for event in solve], label=f"N={row['N']}")
        ax.set_xlabel("Linear solve index")
        ax.set_ylabel("Outer iterations")
        ax.set_title("DFGMRES iterations per solve")
        ax.legend()
        fig.tight_layout()
        name = "iterations-by-solve.png"
        fig.savefig(out_root / name, dpi=150)
        plt.close(fig)
        plots.append(name)

        fig, ax = plt.subplots(figsize=(10, 4.5))
        for row in completed:
            events = read_tsv(Path(row["out_dir"]) / "linear_profile_events.tsv")
            solve = [event for event in events if event.get("phase") == "solve"]
            ax.plot(range(1, len(solve) + 1), [as_int(event.get("deflation_basis_size"), 0) for event in solve], label=f"N={row['N']}")
        ax.set_xlabel("Linear solve index")
        ax.set_ylabel("Basis size")
        ax.set_title("Deflation basis size")
        ax.legend()
        fig.tight_layout()
        name = "basis-size-by-solve.png"
        fig.savefig(out_root / name, dpi=150)
        plt.close(fig)
        plots.append(name)

        labels = [f"N={row['N']}" for row in completed]
        buckets = [
            ("setup_seconds", "hypre setup/factorize"),
            ("preconditioner_apply_seconds", "preconditioner apply"),
            ("orthogonalization_seconds", "orthogonalization"),
            ("matvec_seconds", "matvec"),
            ("least_squares_seconds", "least squares"),
            ("deflation_seconds", "deflation"),
            ("nonlinear_overhead_seconds", "nonlinear/other"),
        ]
        bottoms = [0.0] * len(completed)
        fig, ax = plt.subplots(figsize=(10, 4.8))
        for key, label in buckets:
            values = [float(row[key]) for row in completed]
            ax.bar(labels, values, bottom=bottoms, label=label)
            bottoms = [a + b for a, b in zip(bottoms, values)]
        ax.set_ylabel("Seconds")
        ax.set_title("Timing breakdown by deflation size")
        ax.legend(fontsize="small")
        fig.tight_layout()
        name = "timing-breakdown-by-n.png"
        fig.savefig(out_root / name, dpi=150)
        plt.close(fig)
        plots.append(name)

        fig, ax = plt.subplots(figsize=(10, 4.5))
        for row in completed:
            events = read_tsv(Path(row["out_dir"]) / "linear_profile_events.tsv")
            factorize = [event for event in events if event.get("phase") == "factorize"]
            seen: set[str] = set()
            timeline = []
            for event in factorize:
                seen.add(event.get("value_hash", ""))
                timeline.append(len(seen))
            ax.step(range(1, len(timeline) + 1), timeline, where="post", label=f"N={row['N']}")
        ax.set_xlabel("Factorization/update index")
        ax.set_ylabel("Unique matrix value hashes")
        ax.set_title("Unique matrix timeline")
        ax.legend()
        fig.tight_layout()
        name = "unique-matrix-timeline.png"
        fig.savefig(out_root / name, dpi=150)
        plt.close(fig)
        plots.append(name)

    return plots


def write_summary(out_root: Path, summaries: list[dict[str, object]], plots: list[str], prefix_passed: bool) -> None:
    summary = out_root / "summary.md"
    run_rows = []
    timing_rows = []
    runtime_rows = []
    for row in summaries:
        run_rows.append([
            row["name"],
            "diagnostic" if row["diagnostic"] else ( "prefix" if row["prefix"] else "full" ),
            row["N"],
            row["eps"],
            1 if row["collect"] else 0,
            row["status"],
            row["steps"],
            row["unique_matrices"],
            row["nonlinear_iterations"],
            row["linear_solves"],
            row["total_outer_iterations"],
            row["median_iterations"],
            row["max_iterations"],
            row["max_true_residual"],
            row["raw_candidates"],
            row["max_basis"],
            row["final_basis"],
            row["discarded"],
            row["low_anorm_discards"],
            row["capacity_evictions"],
            row["max_offdiag"],
            row["max_diag_error"],
        ])
        timing_rows.append([
            row["name"],
            row["wall_seconds"],
            row["analyze_seconds"],
            row["setup_seconds"],
            row["solve_seconds"],
            row["preconditioner_apply_seconds"],
            row["orthogonalization_seconds"],
            row["matvec_seconds"],
            row["least_squares_seconds"],
            row["deflation_seconds"],
            row["nonlinear_overhead_seconds"],
        ])
        for hotspot in row.get("runtime_hotspots", []):
            runtime_rows.append([
                row["name"],
                hotspot["phase"],
                hotspot["detail"],
                hotspot["count"],
                hotspot["total_seconds"],
                hotspot["mean_seconds"],
                hotspot["max_seconds"],
            ])

    full_rows = [row for row in summaries if not row["prefix"] and not row["diagnostic"]]
    deflated_rows = [row for row in full_rows if int(row["N"]) > 0]
    fastest_wall = min(full_rows, key=lambda row: float(row["wall_seconds"])) if full_rows else None
    fastest_solve = min(full_rows, key=lambda row: float(row["solve_seconds"])) if full_rows else None
    basis_limited = bool(deflated_rows) and max(int(row["max_basis"]) for row in deflated_rows) <= 1

    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# Dogbone DFGMRES Hypre Deflation Sweep\n\n")
        handle.write("Outer solver is native `DeflatedFGMRES`; hypre BoomerAMG is used only as a right preconditioner apply. ")
        handle.write("The non-deflated baseline is `N=0` with Newton-step collection disabled. Corrected deflated runs use newest-first A-orthogonal basis rebuilds.\n\n")
        if not prefix_passed:
            handle.write("Full Dogbone runs were skipped because one or more prefix gates failed.\n\n")
        handle.write("## Findings\n\n")
        finding_rows = []
        if fastest_wall:
            finding_rows.append(["fastest wall time", f"{fastest_wall['name']} ({float(fastest_wall['wall_seconds']):.3f} s)"])
        if fastest_solve:
            finding_rows.append(["fastest linear solve time", f"{fastest_solve['name']} ({float(fastest_solve['solve_seconds']):.3f} s)"])
        if basis_limited:
            finding_rows.append([
                "basis growth",
                "The selected epsilon still admitted only one vector; inspect diagnostic epsilon rows before trusting N-sweep timing.",
            ])
            finding_rows.append([
                "N sweep interpretation",
                "`N=5,10,20,30` are effectively the same one-vector deflation run for this tolerance.",
            ])
        elif deflated_rows:
            finding_rows.append([
                "basis growth",
                f"Corrected deflated runs reached max basis size {max(int(row['max_basis']) for row in deflated_rows)}.",
            ])
            finding_rows.append([
                "A-orthogonality",
                f"Max offdiag {max(float(row['max_offdiag']) for row in deflated_rows):.3e}, max diagonal error {max(float(row['max_diag_error']) for row in deflated_rows):.3e}.",
            ])
        if full_rows:
            finding_rows.append([
                "hypre setup reuse",
                f"{full_rows[0]['unique_matrices']} unique matrix values over {full_rows[0]['linear_solves']} linear solves; setup/factorize is reused between OAS matrix updates.",
            ])
        handle.write(md_table(["item", "finding"], finding_rows) + "\n\n")
        handle.write("## Run Summary\n\n")
        handle.write(md_table([
            "run",
            "scope",
            "N",
            "eps",
            "collect",
            "status",
            "steps",
            "unique matrices",
            "nonlinear iters",
            "linear solves",
            "outer iters",
            "median iter",
            "max iter",
            "max true relres",
            "raw candidates",
            "max basis",
            "final basis",
            "discarded",
            "low A-norm discards",
            "capacity evictions",
            "max offdiag",
            "max diag error",
        ], run_rows) + "\n\n")
        handle.write("## Timing\n\n")
        handle.write(md_table([
            "run",
            "wall s",
            "analyze s",
            "hypre setup/factorize s",
            "linear solve s",
            "precond apply s",
            "orthogonalization s",
            "matvec s",
            "least squares s",
            "deflation s",
            "nonlinear/other s",
        ], timing_rows) + "\n\n")
        if runtime_rows:
            handle.write("## Runtime Phase Hotspots\n\n")
            handle.write(md_table([
                "run",
                "phase",
                "detail",
                "count",
                "total s",
                "mean s",
                "max s",
            ], runtime_rows) + "\n\n")
        handle.write("## Artifacts\n\n")
        artifact_rows = []
        for row in summaries:
            path = Path(row["out_dir"]).resolve()
            try:
                display_path = path.relative_to(ROOT.resolve())
            except ValueError:
                display_path = path
            artifact_rows.append([row["name"], display_path])
        handle.write(md_table(["run", "folder"], artifact_rows) + "\n\n")
        if plots:
            handle.write("## Plots\n\n")
            for plot in plots:
                handle.write(f"![{plot}]({plot})\n\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", default="16")
    parser.add_argument("--oas-bin", default=str(ROOT / "../oas_deflation-build/release/bin/OAS"))
    parser.add_argument("--mklroot", default="/opt/intel/oneapi/mkl/latest")
    parser.add_argument("--eps", default="1e-15")
    parser.add_argument("--diagnostic-eps", default="1e-3,1e-6,1e-9,1e-12,1e-15")
    parser.add_argument("--prefix-timeout", default=os.environ.get("OAS_PREFIX_TIMEOUT", "30m"))
    parser.add_argument("--full-timeout", default=os.environ.get("OAS_FULL_TIMEOUT", "12h"))
    parser.add_argument("--out-root", default="")
    args = parser.parse_args()

    args.oas_bin = Path(args.oas_bin)
    if not args.oas_bin.exists():
        raise SystemExit(f"Missing OAS binary: {args.oas_bin}")
    if not MASTER.exists():
        raise SystemExit(f"Missing Dogbone master input: {MASTER}")
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_root = Path(args.out_root) if args.out_root else ROOT / f"results/dogbone-dfgmres-hypre-deflation-fixed-{timestamp}"
    out_root.mkdir(parents=True, exist_ok=True)

    summaries: list[dict[str, object]] = []
    for eps in [item.strip() for item in args.diagnostic_eps.split(",") if item.strip()]:
        safe_eps = eps.replace("-", "m").replace("+", "p")
        summaries.append(run_variant(args, out_root, f"diag-eps-{safe_eps}", 30, True, True, eps, diagnostic=True))

    prefix_variants = [
        ("prefix-N0", 0, False),
        ("prefix-N30", 30, True),
    ]
    for name, nvec, collect in prefix_variants:
        summaries.append(run_variant(args, out_root, name, nvec, collect, prefix=True, eps=args.eps))

    prefix_passed = all(row["status"] == "passed" for row in summaries if row["prefix"] and not row["diagnostic"])
    if prefix_passed:
        for nvec in [0, 5, 10, 20, 30]:
            collect = nvec > 0
            name = f"full-N{nvec}"
            summaries.append(run_variant(args, out_root, name, nvec, collect, prefix=False, eps=args.eps))

    plots = make_plots(out_root, summaries)
    write_summary(out_root, summaries, plots, prefix_passed)
    print(out_root / "summary.md")
    return 0 if prefix_passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
