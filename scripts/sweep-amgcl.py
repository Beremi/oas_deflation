#!/usr/bin/env python3
"""Run small local AMGCL/OAS profiling sweeps and summarize the results."""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OAS_BIN = ROOT.parent / "oas_deflation-build" / "release" / "bin" / "OAS"


def split_values(raw: str, cast=str):
    return [cast(item) for item in raw.split(",") if item != ""]


def format_value(value: object) -> str:
    text = str(value)
    return text.replace("+", "").replace(".", "p").replace("-", "m")


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def read_solver_out(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    with path.open("r", encoding="utf-8") as handle:
        columns = handle.readline().strip().lstrip("#").split("\t")
    return pd.read_csv(path, sep="\t", comment="#", names=columns)


def summarize_run(out_dir: Path, run_id: str, params: dict[str, object], status: int) -> dict[str, object]:
    out_dir = out_dir.resolve()
    events_path = out_dir / "linear_profile_events.tsv"
    solver_out_path = out_dir / "solver.out"
    report_path = out_dir / "linear-profile.md"
    try:
        report = str(report_path.relative_to(ROOT))
    except ValueError:
        report = str(report_path)
    row: dict[str, object] = {
        "run_id": run_id,
        "status": status,
        "report": report if report_path.exists() else "",
    }
    row.update(params)

    if events_path.exists():
        events = pd.read_csv(events_path, sep="\t")
        solve = events[events["phase"] == "solve"].copy()
        factorize = events[events["phase"] == "factorize"].copy()
        row["solves"] = len(solve)
        row["solve_success"] = int(solve["success"].sum()) if "success" in solve else ""
        row["solve_seconds"] = f"{solve['duration_seconds'].sum():.6g}" if not solve.empty else "0"
        row["factorize_seconds"] = f"{factorize['duration_seconds'].sum():.6g}" if not factorize.empty else "0"
        if not solve.empty and "solver_iterations" in solve:
            iterative = solve[solve["solver_iterations"] >= 0]
            if not iterative.empty:
                row["iter_median"] = f"{iterative['solver_iterations'].median():.6g}"
                row["iter_max"] = f"{iterative['solver_iterations'].max():.6g}"
                row["relres_median"] = f"{iterative['solver_error'].median():.6g}"
                row["relres_max"] = f"{iterative['solver_error'].max():.6g}"
    else:
        row["solves"] = 0
        row["solve_success"] = 0

    solver_out = read_solver_out(solver_out_path)
    if not solver_out.empty:
        row["steps"] = int(solver_out["step"].max())
        row["nonlinear_iters"] = int(solver_out["iterations"].sum())
        row["wall_seconds"] = f"{float(solver_out['elapsed_time'].iloc[-1]):.6g}"
    else:
        row["steps"] = ""
        row["nonlinear_iters"] = ""
        row["wall_seconds"] = ""

    return row


def build_combinations(args: argparse.Namespace) -> list[dict[str, object]]:
    tolerances = split_values(args.tolerances)
    max_iterations = split_values(args.max_iterations, int)
    backend = split_values(args.backend)
    use_block_backend = split_values(args.use_block_backend, int)
    elastic_full_lift = split_values(args.elastic_full_lift, int)
    block_relaxation = split_values(args.block_relaxation)
    block_size = split_values(args.block_size, int)
    eps_strong = split_values(args.eps_strong)
    relax = split_values(args.relax)
    npre = split_values(args.npre, int)
    npost = split_values(args.npost, int)
    ncycle = split_values(args.ncycle, int)
    near_nullspace = split_values(args.near_nullspace, int)
    reuse_initial_guess = split_values(args.reuse_initial_guess, int)
    coarse_enough = split_values(args.coarse_enough, int)
    estimate_spectral_radius = split_values(args.estimate_spectral_radius, int)
    check_matrix = split_values(args.check_matrix, int)

    combos = []
    for values in itertools.product(
        tolerances,
        max_iterations,
        backend,
        use_block_backend,
        elastic_full_lift,
        block_relaxation,
        block_size,
        eps_strong,
        relax,
        npre,
        npost,
        ncycle,
        near_nullspace,
        reuse_initial_guess,
        coarse_enough,
        estimate_spectral_radius,
        check_matrix,
    ):
        (
            tol,
            maxit,
            backend_value,
            use_block,
            full_lift,
            block_relax,
            block_size_value,
            eps,
            relax_value,
            npre_value,
            npost_value,
            ncycle_value,
            nullspace,
            reuse_guess,
            coarse,
            spectral_radius,
            matrix_check,
        ) = values
        combos.append(
            {
                "amgcl_tolerance": tol,
                "amgcl_max_iterations": maxit,
                "amgcl_backend": backend_value,
                "amgcl_use_block_backend": use_block,
                "amgcl_elastic_full_lift": full_lift,
                "amgcl_block_relaxation": block_relax,
                "amgcl_block_size": block_size_value,
                "amgcl_eps_strong": eps,
                "amgcl_relax": relax_value,
                "amgcl_npre": npre_value,
                "amgcl_npost": npost_value,
                "amgcl_ncycle": ncycle_value,
                "amgcl_near_nullspace": nullspace,
                "amgcl_reuse_initial_guess": reuse_guess,
                "amgcl_coarse_enough": coarse,
                "amgcl_estimate_spectral_radius": spectral_radius,
                "amgcl_check_matrix": matrix_check,
            }
        )

    if args.limit:
        combos = combos[: args.limit]
    return combos


def mode_overrides(mode: str) -> tuple[str, str, dict[str, object]]:
    if mode == "dogbone-prefix":
        return "Dogbone", "Dogbone AMGCL-CG Prefix Linear-Solve Profile", {
            "total_time": "0.01",
            "max_iterations": "5",
            "min_iterations": "5",
            "limit_tolerance": "10",
        }
    if mode == "dogbone-full":
        return "Dogbone", "Dogbone AMGCL-CG Full Linear-Solve Profile", {}
    if mode == "tsn65-prefix":
        return "TS-N_65", "TS-N_65 AMGCL-CG 8-thread First-Step Prefix Linear-Solve Profile", {
            "total_time": "5.000000e-03",
            "max_iterations": "5",
            "min_iterations": "5",
            "limit_tolerance": "10",
        }
    if mode == "tsn65-firstsolve":
        return "TS-N_65", "TS-N_65 AMGCL-CG First Linear-Solve Profile", {
            "total_time": "5.000000e-03",
            "max_iterations": "1",
            "min_iterations": "1",
            "limit_tolerance": "10",
            "linear_solver_profile_matrix_delta": "0",
            "linear_solver_replay_dump": "1",
            "linear_solver_replay_limit": "1",
            "linear_solver_replay_dir": "linear_replay",
        }
    if mode == "tsn65-10steps":
        return "TS-N_65", "TS-N_65 AMGCL-CG 8-thread Ten-Step Linear-Solve Profile", {
            "total_time": "5.000000e-02",
        }
    raise ValueError(f"Unsupported mode: {mode}")


def default_paths(mode: str, args: argparse.Namespace) -> tuple[Path, Path]:
    if args.case_dir and args.master_file:
        return Path(args.case_dir), Path(args.master_file)
    if mode.startswith("dogbone"):
        case_dir = ROOT / "data" / "cases" / "Dogbone"
        return case_dir, case_dir / "master.inp"
    case_dir = ROOT / "data" / "cases" / "TS-N_65"
    return case_dir, case_dir / "master.inp"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["dogbone-prefix", "dogbone-full", "tsn65-prefix", "tsn65-firstsolve", "tsn65-10steps"], required=True)
    parser.add_argument("--threads", default="4")
    parser.add_argument("--oas-bin", default=str(DEFAULT_OAS_BIN))
    parser.add_argument("--case-dir")
    parser.add_argument("--master-file")
    parser.add_argument("--out-root")
    parser.add_argument("--timeout", default="12h")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--tolerances", default="1e-6")
    parser.add_argument("--max-iterations", default="500")
    parser.add_argument("--backend", default="auto")
    parser.add_argument("--use-block-backend", default="1")
    parser.add_argument("--elastic-full-lift", default="1")
    parser.add_argument("--block-relaxation", default="ilu0")
    parser.add_argument("--block-size", default="6")
    parser.add_argument("--eps-strong", default="0")
    parser.add_argument("--relax", default="1")
    parser.add_argument("--npre", default="1")
    parser.add_argument("--npost", default="1")
    parser.add_argument("--ncycle", default="1")
    parser.add_argument("--near-nullspace", default="1")
    parser.add_argument("--reuse-initial-guess", default="0")
    parser.add_argument("--coarse-enough", default="0")
    parser.add_argument("--estimate-spectral-radius", default="0")
    parser.add_argument("--check-matrix", default="0")
    args = parser.parse_args()

    case_name, title, base_overrides = mode_overrides(args.mode)
    case_dir, master_file = default_paths(args.mode, args)
    oas_bin = Path(args.oas_bin)
    if not oas_bin.exists():
        raise SystemExit(f"Missing OAS executable: {oas_bin}")
    if not master_file.exists():
        raise SystemExit(f"Missing OAS master file: {master_file}")

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_root = Path(args.out_root) if args.out_root else ROOT / "results" / f"amgcl-{args.mode}-{stamp}"
    out_root.mkdir(parents=True, exist_ok=True)

    combinations = build_combinations(args)
    rows: list[dict[str, object]] = []
    env = os.environ.copy()
    env["OAS_TIMEOUT"] = args.timeout

    for index, params in enumerate(combinations, start=1):
        short = [
            f"tol{format_value(params['amgcl_tolerance'])}",
            f"it{params['amgcl_max_iterations']}",
            f"be{format_value(params['amgcl_backend'])}",
            f"blk{params['amgcl_use_block_backend']}",
            f"lift{params['amgcl_elastic_full_lift']}",
            f"br{format_value(params['amgcl_block_relaxation'])}",
            f"bs{params['amgcl_block_size']}",
            f"eps{format_value(params['amgcl_eps_strong'])}",
            f"r{format_value(params['amgcl_relax'])}",
            f"pre{params['amgcl_npre']}",
            f"post{params['amgcl_npost']}",
            f"ns{params['amgcl_near_nullspace']}",
            f"reuse{params['amgcl_reuse_initial_guess']}",
        ]
        run_id = f"{index:03d}-" + "-".join(short)
        out_dir = out_root / run_id
        overrides = {**base_overrides, **params}
        cmd = [
            "scripts/run-oas-profile.sh",
            str(case_dir),
            str(master_file),
            str(oas_bin),
            str(args.threads),
            "AmgclCG",
            title,
            case_name,
            str(out_dir),
        ] + [f"{key}={value}" for key, value in overrides.items()]

        print(f"[{index}/{len(combinations)}] {run_id}", flush=True)
        completed = subprocess.run(cmd, cwd=ROOT, env=env)
        rows.append(summarize_run(out_dir, run_id, params, completed.returncode))

    columns = [
        "run_id",
        "status",
        "amgcl_tolerance",
        "amgcl_max_iterations",
        "amgcl_backend",
        "amgcl_use_block_backend",
        "amgcl_elastic_full_lift",
        "amgcl_block_relaxation",
        "amgcl_block_size",
        "amgcl_eps_strong",
        "amgcl_relax",
        "amgcl_npre",
        "amgcl_npost",
        "amgcl_ncycle",
        "amgcl_near_nullspace",
        "amgcl_reuse_initial_guess",
        "amgcl_coarse_enough",
        "amgcl_estimate_spectral_radius",
        "amgcl_check_matrix",
        "solves",
        "solve_success",
        "iter_median",
        "iter_max",
        "relres_median",
        "relres_max",
        "solve_seconds",
        "factorize_seconds",
        "steps",
        "nonlinear_iters",
        "wall_seconds",
        "report",
    ]

    summary_tsv = out_root / "summary.tsv"
    pd.DataFrame(rows).to_csv(summary_tsv, sep="\t", index=False)
    summary_md = out_root / "summary.md"
    summary_md.write_text(
        f"# AMGCL Sweep: {args.mode}\n\n"
        f"Case: `{case_name}`  \n"
        f"Threads: `{args.threads}`  \n"
        f"Runs: `{len(rows)}`\n\n"
        + markdown_table(rows, columns)
        + "\n",
        encoding="utf-8",
    )
    print(summary_md)
    return 0 if all(row["status"] == 0 for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
