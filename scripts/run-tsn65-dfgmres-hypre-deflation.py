#!/usr/bin/env python3
"""Run TS-N_65 two-step DFGMRES+hypre deflation experiments."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "data/cases/TS-N_65"
MASTER_FILE = CASE_DIR / "master.inp"
PROFILE_DIR = CASE_DIR / "results"
DEFAULT_OAS_BIN = ROOT.parent / "oas_deflation-build/release/bin/OAS"


@dataclass
class RunSpec:
    name: str
    nvec: int


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


def format_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        if abs(value) >= 1000 or (0 < abs(value) < 1e-3):
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)


def md_table(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(format_cell(item) for item in row) + " |")
    return "\n".join(lines)


def read_solver_out(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip().lstrip("#").split("\t")
        rows = []
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def patch_key_value_file(path: Path, updates: dict[str, str], marker: str) -> str:
    original = path.read_text(encoding="utf-8", errors="replace")
    lines = original.splitlines()
    out: list[str] = []
    seen: set[str] = set()
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            out.append(line)
            continue
        key = stripped.split()[0]
        if key in updates:
            out.append(f"{key}\t{updates[key]}")
            seen.add(key)
        else:
            out.append(line)
    missing = [key for key in updates if key not in seen]
    if missing:
        if out and out[-1].strip():
            out.append("")
        out.append(f"# {marker}")
        for key in missing:
            out.append(f"{key}\t{updates[key]}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return original


def clean_case_results() -> None:
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    for pattern in [
        "linear_profile_events.tsv",
        "linear_profile_iterations.tsv",
        "runtime_profile_events.tsv",
        "runtime_profile_summary.tsv",
        "solver.out",
        "LD.out",
    ]:
        for path in PROFILE_DIR.glob(pattern):
            if path.is_file():
                path.unlink()


def solver_updates(nvec: int, eps: str, total_time: str) -> dict[str, str]:
    return {
        "total_time": total_time,
        "solver_type": "DeflatedFGMRES",
        "linear_solver_profile": "1",
        "linear_solver_profile_matrix_delta": "0",
        "linear_solver_profile_file": "linear_profile",
        "runtime_phase_profile": "1",
        "runtime_phase_profile_file": "runtime_profile",
        "stiffness_matrix_iter_update": "-5",
        "stiffness_matrix_step_update": "1",
        "max_iterations": "1000",
        "dfgmres_tolerance": "1e-6",
        "dfgmres_true_tolerance": "1e-6",
        "dfgmres_max_iterations": "500",
        "dfgmres_restart": "80",
        "dfgmres_deflation_vectors": str(nvec),
        "dfgmres_deflation_eps": eps,
        "dfgmres_collect_newton_steps": "1",
        "dfgmres_preconditioner": "hypre",
        "dfgmres_reorthogonalize_on_matrix_change": "1",
        "amgcl_near_nullspace": "1",
        "amgcl_elastic_full_lift": "1",
        "amgcl_use_block_backend": "1",
        "hypre_coarsen_type": "8",
        "hypre_interp_type": "6",
        "hypre_strong_threshold": "0.25",
        "hypre_nodal": "4",
        "hypre_relax_type": "6",
        "hypre_p_max": "4",
        "hypre_boomer_max_iterations": "1",
        "hypre_use_dof_functions": "0",
        "hypre_use_interp_vectors": "0",
    }


def copy_artifacts(run_dir: Path) -> None:
    for name in [
        "linear_profile_events.tsv",
        "linear_profile_iterations.tsv",
        "runtime_profile_events.tsv",
        "runtime_profile_summary.tsv",
        "solver.out",
        "LD.out",
        "version.txt",
    ]:
        src = PROFILE_DIR / name
        if src.exists():
            shutil.copy2(src, run_dir / name)


def analyze_linear_profile(run_dir: Path, title: str, threads: int) -> None:
    if not (run_dir / "linear_profile_events.tsv").exists():
        return
    cmd = [
        sys.executable,
        str(ROOT / "scripts/analyze-linear-profile.py"),
        "--profile-dir",
        str(run_dir),
        "--out-dir",
        str(run_dir),
        "--title",
        title,
        "--case",
        "TS-N_65",
        "--threads",
        str(threads),
    ]
    subprocess.run(cmd, cwd=str(ROOT), check=False)


def summarize_run(run_dir: Path, spec: RunSpec, returncode: int, wall_seconds: float, eps: str) -> dict[str, Any]:
    events = read_tsv(run_dir / "linear_profile_events.tsv")
    solver_rows = read_solver_out(run_dir / "solver.out")
    solve = [row for row in events if row.get("phase") == "solve"]
    factorize = [row for row in events if row.get("phase") == "factorize"]
    analyze = [row for row in events if row.get("phase") == "analyze"]
    iterations = [as_int(row.get("solver_iterations"), -1) for row in solve if as_int(row.get("solver_iterations"), -1) >= 0]
    residuals = [as_float(row.get("solver_error")) for row in solve if as_float(row.get("solver_error")) >= 0]
    successes = [as_int(row.get("success"), 0) for row in solve]
    max_step = max([as_int(row.get("step"), 0) for row in solver_rows], default=0)
    nonlinear_iterations = sum(as_int(row.get("iterations"), 0) for row in solver_rows)
    elapsed_solver = max([as_float(row.get("elapsed_time"), 0.0) for row in solver_rows], default=0.0)
    elapsed_wall = max(elapsed_solver, wall_seconds)
    solve_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in solve)
    setup_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in factorize)
    analyze_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in analyze)
    precond_apply = sum(as_float(row.get("preconditioner_apply_seconds"), 0.0) for row in solve)
    orth = sum(as_float(row.get("orthogonalization_seconds"), 0.0) for row in solve)
    matvec = sum(as_float(row.get("matvec_seconds"), 0.0) for row in solve)
    least_squares = sum(as_float(row.get("least_squares_seconds"), 0.0) for row in solve)
    deflation = sum(as_float(row.get("deflation_seconds"), 0.0) for row in solve)
    linear_seconds = analyze_seconds + setup_seconds + solve_seconds
    value_hashes = {row.get("value_hash", "") for row in factorize if row.get("value_hash")}
    basis_values = [as_int(row.get("deflation_basis_size"), 0) for row in solve]
    log_text = (run_dir / "oas.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "oas.log").exists() else ""
    log_basis_values = [int(match.group(1)) for match in re.finditer(r"active_basis_size=(\d+)", log_text)]
    max_basis = max(basis_values + log_basis_values, default=0)
    final_basis = log_basis_values[-1] if log_basis_values else (basis_values[-1] if basis_values else 0)
    max_offdiag = max([as_float(row.get("deflation_orthogonality_max_offdiag"), 0.0) for row in solve], default=0.0)
    max_diag_error = max([as_float(row.get("deflation_orthogonality_max_diag_error"), 0.0) for row in solve], default=0.0)
    discarded = max([as_int(row.get("deflation_discarded_count"), 0) for row in solve], default=0)
    capacity_evictions = max([as_int(row.get("deflation_capacity_eviction_count"), 0) for row in solve], default=0)
    status = "completed_2_steps" if returncode == 0 and max_step >= 2 else "partial"
    if solve and (not all(successes) or max(residuals, default=math.inf) > 1e-6):
        status = "linear_failed"
    if returncode != 0 and max_step == 0 and not solve:
        status = "failed"

    return {
        "name": spec.name,
        "N": spec.nvec,
        "eps": eps,
        "status": status,
        "returncode": returncode,
        "steps": max_step,
        "nonlinear_iterations": nonlinear_iterations,
        "linear_solves": len(solve),
        "factorizations": len(factorize),
        "unique_matrices": len(value_hashes),
        "outer_iterations": sum(iterations),
        "median_iterations": median(iterations) if iterations else None,
        "max_iterations": max(iterations, default=None),
        "max_true_residual": max(residuals, default=None),
        "max_basis": max_basis,
        "final_basis": final_basis,
        "discarded": discarded,
        "capacity_evictions": capacity_evictions,
        "max_offdiag": max_offdiag,
        "max_diag_error": max_diag_error,
        "wall_seconds": elapsed_wall,
        "analyze_seconds": analyze_seconds,
        "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds,
        "preconditioner_apply_seconds": precond_apply,
        "orthogonalization_seconds": orth,
        "matvec_seconds": matvec,
        "least_squares_seconds": least_squares,
        "deflation_seconds": deflation,
        "linear_share": linear_seconds / elapsed_wall if elapsed_wall > 0 else None,
        "other_seconds": max(elapsed_wall - linear_seconds, 0.0),
        "runtime_hotspots": top_runtime_hotspots(run_dir),
    }


def top_runtime_hotspots(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    rows = read_tsv(run_dir / "runtime_profile_summary.tsv")
    parsed: list[dict[str, Any]] = []
    for row in rows:
        parsed.append(
            {
                "phase": row.get("phase", ""),
                "detail": row.get("detail", ""),
                "count": as_int(row.get("count"), 0),
                "total_seconds": as_float(row.get("total_seconds"), 0.0),
                "mean_seconds": as_float(row.get("mean_seconds"), 0.0),
                "max_seconds": as_float(row.get("max_seconds"), 0.0),
            }
        )
    parsed.sort(key=lambda row: row["total_seconds"], reverse=True)
    return parsed[:limit]


def run_spec(spec: RunSpec, out_dir: Path, oas_bin: Path, threads: int, eps: str, total_time: str) -> dict[str, Any]:
    run_dir = out_dir / spec.name
    run_dir.mkdir(parents=True, exist_ok=True)
    solver_inp = CASE_DIR / "solver.inp"
    original = solver_inp.read_text(encoding="utf-8", errors="replace")
    (run_dir / "solver.inp.original").write_text(original, encoding="utf-8")
    started = datetime.now().isoformat(timespec="seconds")
    start_time = time.monotonic()
    returncode = 999
    try:
        patch_key_value_file(solver_inp, solver_updates(spec.nvec, eps, total_time), "Local TS-N_65 DFGMRES-hypre deflation settings.")
        (run_dir / "solver.inp.effective").write_text(solver_inp.read_text(encoding="utf-8"), encoding="utf-8")
        clean_case_results()
        log_path = run_dir / "oas.log"
        env = os.environ.copy()
        env.setdefault("MKLROOT", "/opt/intel/oneapi/mkl/latest")
        env["MKL_NUM_THREADS"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
        cmd = [str(oas_bin), "-j", str(threads), str(MASTER_FILE)]
        print(f"[{started}] starting {spec.name}: {' '.join(cmd)}", flush=True)
        with log_path.open("w", encoding="utf-8", errors="replace") as log:
            log.write(f"Starting {spec.name} at {started}\n")
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write(f"Threads: {threads}\n")
            log.write(f"Deflation vectors: {spec.nvec}\n")
            log.write(f"Deflation eps: {eps}\n")
            log.flush()
            process = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT, text=True)
            while True:
                try:
                    returncode = process.wait(timeout=60)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - start_time
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {spec.name} still running ({elapsed:.0f}s)", flush=True)
                    log.write(f"\n[runner] still running at {elapsed:.0f}s\n")
                    log.flush()
            log.write(f"\nOAS exit status: {returncode}\nFinished at {datetime.now().isoformat(timespec='seconds')}\n")
        print(f"[{datetime.now().isoformat(timespec='seconds')}] finished {spec.name} with status {returncode}", flush=True)
    finally:
        solver_inp.write_text(original, encoding="utf-8")

    elapsed = time.monotonic() - start_time
    copy_artifacts(run_dir)
    analyze_linear_profile(run_dir, f"TS-N_65 DeflatedFGMRES hypre N={spec.nvec} Two-Step Profile", threads)
    manifest = summarize_run(run_dir, spec, returncode, elapsed, eps)
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def make_plots(out_dir: Path, manifests: list[dict[str, Any]]) -> list[str]:
    plots: list[str] = []
    completed = [row for row in manifests if (out_dir / row["name"] / "linear_profile_events.tsv").exists()]
    if not completed:
        return plots
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for row in completed:
        events = read_tsv(out_dir / row["name"] / "linear_profile_events.tsv")
        solve = [event for event in events if event.get("phase") == "solve"]
        ax.plot(range(1, len(solve) + 1), [as_int(event.get("solver_iterations"), 0) for event in solve], label=f"N={row['N']}")
    ax.set_xlabel("Linear solve index")
    ax.set_ylabel("DFGMRES iterations")
    ax.legend()
    fig.tight_layout()
    name = "iterations-by-solve.png"
    fig.savefig(out_dir / name, dpi=150)
    plt.close(fig)
    plots.append(name)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for row in completed:
        events = read_tsv(out_dir / row["name"] / "linear_profile_events.tsv")
        solve = [event for event in events if event.get("phase") == "solve"]
        ax.plot(range(1, len(solve) + 1), [as_int(event.get("deflation_basis_size"), 0) for event in solve], label=f"N={row['N']}")
    ax.set_xlabel("Linear solve index")
    ax.set_ylabel("Basis size")
    ax.legend()
    fig.tight_layout()
    name = "basis-size-by-solve.png"
    fig.savefig(out_dir / name, dpi=150)
    plt.close(fig)
    plots.append(name)
    return plots


def write_summary(out_dir: Path, manifests: list[dict[str, Any]]) -> None:
    plots = make_plots(out_dir, manifests)
    run_rows = []
    timing_rows = []
    hotspot_rows = []
    for row in manifests:
        run_rows.append([
            row["name"],
            row["N"],
            row["eps"],
            row["status"],
            row["steps"],
            row["nonlinear_iterations"],
            row["linear_solves"],
            row["factorizations"],
            row["unique_matrices"],
            row["outer_iterations"],
            row["median_iterations"],
            row["max_iterations"],
            row["max_true_residual"],
            row["max_basis"],
            row["final_basis"],
            row["discarded"],
            row["capacity_evictions"],
            row["max_offdiag"],
        ])
        timing_rows.append([
            row["name"],
            row["wall_seconds"],
            row["setup_seconds"],
            row["solve_seconds"],
            row["preconditioner_apply_seconds"],
            row["orthogonalization_seconds"],
            row["matvec_seconds"],
            row["least_squares_seconds"],
            row["deflation_seconds"],
            row["other_seconds"],
            row["linear_share"],
        ])
        for hotspot in row.get("runtime_hotspots", [])[:5]:
            hotspot_rows.append([
                row["name"],
                hotspot["phase"],
                hotspot["detail"],
                hotspot["count"],
                hotspot["total_seconds"],
                hotspot["mean_seconds"],
                hotspot["max_seconds"],
            ])
    with (out_dir / "summary.md").open("w", encoding="utf-8") as handle:
        handle.write("# TS-N_65 DFGMRES Hypre Deflation Two-Step Runs\n\n")
        handle.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n\n")
        handle.write("Outer solver is native `DeflatedFGMRES`; hypre BoomerAMG is used only as a right preconditioner. Runs target the first two TS-N_65 load steps with `total_time=0.01`, 16 shared-memory threads, and matrix-delta profiling disabled.\n\n")
        handle.write("## Run Summary\n\n")
        handle.write(md_table([
            "run", "N", "eps", "status", "steps", "nonlinear iters", "linear solves", "factorizations", "unique matrices", "outer iters", "median iter", "max iter", "max true relres", "max basis", "final basis", "discarded", "capacity evictions", "max offdiag",
        ], run_rows) + "\n\n")
        handle.write("## Timing\n\n")
        handle.write(md_table([
            "run", "wall s", "setup s", "solve s", "precond apply s", "orth s", "matvec s", "least squares s", "deflation s", "other s", "linear share",
        ], timing_rows) + "\n\n")
        if hotspot_rows:
            handle.write("## Runtime Phase Hotspots\n\n")
            handle.write(md_table([
                "run", "phase", "detail", "count", "total s", "mean s", "max s",
            ], hotspot_rows) + "\n\n")
        handle.write("## Artifacts\n\n")
        artifact_rows = []
        for row in manifests:
            path = (out_dir / row["name"]).resolve()
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
    print(out_dir / "summary.md")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oas-bin", type=Path, default=DEFAULT_OAS_BIN)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--eps", default="1e-15")
    parser.add_argument("--total-time", default="1.000000e-02")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--report-only", type=Path, default=None)
    args = parser.parse_args()

    if args.report_only:
        manifests = []
        for run_json in sorted(args.report_only.glob("*/run.json")):
            manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
        write_summary(args.report_only, manifests)
        return 0

    if not args.oas_bin.exists() or not os.access(args.oas_bin, os.X_OK):
        raise SystemExit(f"OAS binary is missing or not executable: {args.oas_bin}")
    if not MASTER_FILE.exists():
        raise SystemExit(f"TS-N_65 master input is missing: {MASTER_FILE}")

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = args.out_dir or ROOT / f"results/tsn65-dfgmres-hypre-deflation-{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    specs = [
        RunSpec("dfgmres-hypre-N0", 0),
        RunSpec("dfgmres-hypre-N5", 5),
        RunSpec("dfgmres-hypre-N10", 10),
    ]
    if args.only:
        specs = [spec for spec in specs if any(token in spec.name for token in args.only)]
    manifests = []
    for run_json in sorted(out_dir.glob("*/run.json")):
        manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
    seen = {item["name"] for item in manifests}
    for spec in specs:
        if spec.name in seen:
            continue
        manifest = run_spec(spec, out_dir, args.oas_bin, args.threads, args.eps, args.total_time)
        manifests.append(manifest)
        write_summary(out_dir, manifests)
    write_summary(out_dir, manifests)
    (out_dir / "campaign.json").write_text(json.dumps(manifests, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
