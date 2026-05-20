#!/usr/bin/env python3
"""Run the CP4 arc-length prototype checkpoint on a small benchmark."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from summarize_oas_run import format_tsv, parse_log  # noqa: E402


ARC_RE = re.compile(
    r"^ARC_LENGTH\s+step\s+(?P<step>\d+)\s+it\s+(?P<it>\d+)\s+"
    r"lambda\s+(?P<lambda>\S+)\s+dlambda\s+(?P<dlambda>\S+)\s+"
    r"radius\s+(?P<radius>\S+)"
)

PROPORTIONAL_FUNCTIONS = """PWLFunction 2 0 1 0 0
PWLFunction 2 0 1 0 17
PWLFunction 2 0 1 0 13
PWLFunction 2 0 1 0 11
PWLFunction 2 0 1 0 7
PWLFunction 2 0 1 0 5
PWLFunction 2 0 1 0 3
PWLFunction 2 0 1 0 1
"""

LOAD_SOLVER = """SteadyStateNonLinearSolver
time_step 0.1
total_time 0.3
solver_type EigenSparseLU
nonlinear_control load
tolerance 1e-10
max_iterations 20
limit_tolerance 0
stiffness_matrix_iter_update 1
"""
LEGACY_LOAD_SOLVER = LOAD_SOLVER.replace("nonlinear_control load\n", "")

ARC_SOLVER = """SteadyStateNonLinearSolver
time_step 0.1
total_time 0.3
solver_type EigenSparseLU
nonlinear_control arc_length
arc_length_radius_initial 0.1
arc_length_radius_min 1e-5
arc_length_radius_max 0.1
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 3
arc_length_max_iterations 20
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
tolerance 1e-10
max_iterations 20
limit_tolerance 0
stiffness_matrix_iter_update 1
"""

ARC_SOLVER_WITH_LINE_SEARCH = ARC_SOLVER + """nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_min_alpha 0.5
nonlinear_line_search_max_trials 2
nonlinear_line_search_cutback_on_fail 0
nonlinear_line_search_report_trials 1
"""


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def run_git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def prepare_case(source: Path, run_dir: Path, solver: str) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(source, run_dir)
    (run_dir / "solver.inp").write_text(solver, encoding="utf-8")
    (run_dir / "functions.inp").write_text(PROPORTIONAL_FUNCTIONS, encoding="utf-8")


def run_command(command: list[str], cwd: Path, env: dict[str, str], timeout_seconds: int) -> tuple[int, float, bool]:
    start = time.monotonic()
    timed_out = False
    with (cwd / "solver.out").open("w", encoding="utf-8", errors="replace") as log:
        try:
            completed = subprocess.run(
                command,
                cwd=cwd,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=timeout_seconds,
                check=False,
            )
            exit_status = completed.returncode
        except subprocess.TimeoutExpired:
            exit_status = 124
            timed_out = True
            log.write(f"\nRUNNER_TIMEOUT after {timeout_seconds} seconds\n")
    return exit_status, time.monotonic() - start, timed_out


def parse_arc_rows(log: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for raw in log.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ARC_RE.match(raw.strip())
        if not match:
            continue
        rows.append(
            {
                "step": int(match.group("step")),
                "iteration": int(match.group("it")),
                "lambda": float(match.group("lambda")),
                "dlambda": float(match.group("dlambda")),
                "radius": float(match.group("radius")),
            }
        )
    return rows


def write_arc_tsv(path: Path, rows: list[dict[str, float]]) -> None:
    lines = ["step\titeration\tlambda\tdlambda\tradius"]
    for row in rows:
        lines.append(
            f"{int(row['step'])}\t{int(row['iteration'])}\t{row['lambda']:.12e}\t"
            f"{row['dlambda']:.12e}\t{row['radius']:.12e}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_env(jobs: int) -> dict[str, str]:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": str(jobs),
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "OMP_PLACES": "cores",
        }
    )
    return env


def run_variant(root: Path, source: Path, exe: Path, jobs: int, timeout: int, name: str, solver: str) -> dict[str, Any]:
    run_dir = root / "runs" / name
    prepare_case(source, run_dir, solver)
    env = make_env(jobs)
    command = [str(exe), "-j", str(jobs), "master.inp"]
    started = now_iso()
    exit_status, elapsed, timed_out = run_command(command, run_dir, env, timeout)
    summary = parse_log(run_dir / "solver.out")
    arc_rows = parse_arc_rows(run_dir / "solver.out")
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (run_dir / "step_summary.tsv").write_text(format_tsv(summary), encoding="utf-8")
    write_arc_tsv(run_dir / "arc_length.tsv", arc_rows)
    (run_dir / "run_metadata.txt").write_text(
        f"started={started}\n"
        f"ended={now_iso()}\n"
        f"elapsed_seconds={elapsed:.3f}\n"
        f"exit_status={exit_status}\n"
        f"timed_out={int(timed_out)}\n"
        f"command={' '.join(command)}\n",
        encoding="utf-8",
    )
    return {
        "name": name,
        "run_dir": str(run_dir),
        "exit_status": exit_status,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "summary": summary,
        "arc_rows": arc_rows,
    }


def verdict(result: dict[str, Any]) -> str:
    summary = result["summary"]
    if result["exit_status"] != 0 or result["timed_out"]:
        return "failed"
    if not summary.get("end_of_calculation"):
        return "failed"
    if summary.get("fallback_acceptance_count", 0) != 0:
        return "failed"
    if summary.get("nan_count", 0) != 0:
        return "failed"
    return "passed"


def write_report(root: Path, repo: Path, exe: Path, results: list[dict[str, Any]]) -> None:
    commit = run_git(repo, "rev-parse", "HEAD")
    status = run_git(repo, "status", "--short")
    (root / "metadata.txt").write_text(
        f"created_at={now_iso()}\n"
        f"repo={repo}\n"
        f"branch={run_git(repo, 'branch', '--show-current')}\n"
        f"commit={commit}\n"
        f"executable={exe}\n"
        "source_benchmark=src/benchmark/Timoshenko_beam3D\n",
        encoding="utf-8",
    )
    with (root / "cp4_arc_length.tsv").open("w", encoding="utf-8") as handle:
        handle.write("variant\tverdict\tsteps\trows\tfinal_lambda\tarc_rows\texit_status\telapsed_seconds\tmin_alpha\tmax_ls_trials\n")
        for result in results:
            summary = result["summary"]
            steps = summary.get("steps") or []
            min_alpha_values = [step["min_alpha"] for step in steps if step.get("min_alpha") is not None]
            max_trials = max([step.get("max_ls_trials", 0) for step in steps] or [0])
            final_lambda = ""
            if result["arc_rows"]:
                final_lambda = f"{result['arc_rows'][-1]['lambda']:.12g}"
            elif steps:
                final_lambda = f"{steps[-1].get('time'):.12g}"
            handle.write(
                "\t".join(
                    [
                        result["name"],
                        verdict(result),
                        str(summary.get("step_count", 0)),
                        str(summary.get("total_iters", 0)),
                        final_lambda,
                        str(len(result["arc_rows"])),
                        str(result["exit_status"]),
                        f"{result['elapsed_seconds']:.3f}",
                        "" if not min_alpha_values else f"{min(min_alpha_values):.12g}",
                        str(max_trials),
                    ]
                )
                + "\n"
            )
    lines = [
        "# CP4 Arc-Length Prototype Report",
        "",
        f"Generated: {now_iso()}",
        "",
        "## Code State",
        "",
        f"- Repository: `{repo}`",
        f"- Commit: `{commit}`",
        f"- Git status at report time: `{status if status else 'clean'}`",
        f"- Executable: `{exe}`",
        "",
        "## Benchmark",
        "",
        "Small benchmark source: `src/benchmark/Timoshenko_beam3D`, copied into each run directory.",
        "The runner rewrites the load functions as proportional PWL loads so the prototype reference load is `f_ext(1)-f_ext(0)`.",
        "",
        "## Results",
        "",
        "| variant | verdict | steps | rows | final time/lambda | arc rows | exit | elapsed s |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        summary = result["summary"]
        steps = summary.get("step_count", 0)
        rows = summary.get("total_iters", 0)
        final_time = None
        if summary.get("steps"):
            final_time = summary["steps"][-1].get("time")
        arc_rows = result["arc_rows"]
        if arc_rows:
            final_time = arc_rows[-1]["lambda"]
        lines.append(
            f"| `{result['name']}` | {verdict(result)} | {steps} | {rows} | "
            f"{final_time if final_time is not None else 'n/a'} | {len(arc_rows)} | "
            f"{result['exit_status']} | {result['elapsed_seconds']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Checkpoint Notes",
            "",
            "- `CP4-legacy-load-default` is the no-keyword legacy load-control run.",
            "- `CP4-load-control-regression` exercises the explicit `nonlinear_control load` path.",
            "- `CP4-arc-length-prototype` exercises the spherical Schur-complement predictor/corrector path.",
            "- `CP4-arc-length-line-search` runs the same prototype with snapshot rollback and arc-length trial line-search logging enabled.",
            "- Large run folders remain under `results/*/runs/` and are not intended for commit.",
        ]
    )
    (root / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--exe", type=Path, default=Path("/tmp/oas_tsn65_full_baseline_build/bin/OAS"))
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--root", type=Path)
    args = parser.parse_args()

    repo = args.repo.resolve()
    root = args.root or repo / "results" / f"arc-length-cp4-{now_stamp()}"
    root.mkdir(parents=True, exist_ok=True)
    source = repo / "src/benchmark/Timoshenko_beam3D"

    variants = [
        ("CP4-legacy-load-default", LEGACY_LOAD_SOLVER),
        ("CP4-load-control-regression", LOAD_SOLVER),
        ("CP4-arc-length-prototype", ARC_SOLVER),
        ("CP4-arc-length-line-search", ARC_SOLVER_WITH_LINE_SEARCH),
    ]
    results = [run_variant(root, source, args.exe, args.jobs, args.timeout, name, solver) for name, solver in variants]
    write_report(root, repo, args.exe, results)
    print(root)


if __name__ == "__main__":
    main()
