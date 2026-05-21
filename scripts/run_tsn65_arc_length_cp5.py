#!/usr/bin/env python3
"""Run CP5-style TS-N65 arc-length experiments.

The script keeps the large run folders under results/*/runs and writes compact
TSV/Markdown summaries for comparing against the strict 8-quarter baseline.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_tsn65_globalization_sweep import BASELINE_SOLVER, prepare_run_dir  # noqa: E402
from summarize_oas_run import TSN65_TARGET_ITERS, TSN65_TARGET_TOTAL, format_tsv, parse_log  # noqa: E402


ARC_RE = re.compile(
    r"^ARC_LENGTH\s+step\s+(?P<step>\d+)\s+it\s+(?P<it>\d+)\s+"
    r"lambda\s+(?P<lambda>\S+)\s+dlambda\s+(?P<dlambda>\S+)\s+"
    r"radius\s+(?P<radius>\S+)"
)


VARIANTS: dict[str, str] = {
    "CP5-arc-proportional-smoke": """
nonlinear_control arc_length
arc_length_radius_initial 1.25e-3
arc_length_radius_min 1.25e-5
arc_length_radius_max 1.25e-3
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 40
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
""",
    "CP5-arc-fd-r1p25e-3": """
nonlinear_control arc_length
arc_length_reference finite_difference
arc_length_reference_delta 1.25e-3
arc_length_radius_initial 1.25e-3
arc_length_radius_min 1.25e-5
arc_length_radius_max 1.25e-3
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 300
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
nonlinear_material_snapshot_rollback 1
""",
    "CP5-arc-fd-r2p5e-3": """
nonlinear_control arc_length
arc_length_reference finite_difference
arc_length_reference_delta 1.25e-3
arc_length_radius_initial 2.5e-3
arc_length_radius_min 1.25e-5
arc_length_radius_max 2.5e-3
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 300
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
nonlinear_material_snapshot_rollback 1
""",
    "CP5-arc-fd-r1p35e-2": """
nonlinear_control arc_length
arc_length_reference finite_difference
arc_length_reference_delta 1.25e-3
arc_length_radius_initial 1.35e-2
arc_length_radius_min 1.25e-5
arc_length_radius_max 1.35e-2
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 300
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
nonlinear_material_snapshot_rollback 1
""",
    "CP5-arc-fd-r1p25e-3-line-search": """
nonlinear_control arc_length
arc_length_reference finite_difference
arc_length_reference_delta 1.25e-3
arc_length_radius_initial 1.25e-3
arc_length_radius_min 1.25e-5
arc_length_radius_max 1.25e-3
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 300
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 1
""",
    "CP5-arc-fd-r1p35e-2-line-search": """
nonlinear_control arc_length
arc_length_reference finite_difference
arc_length_reference_delta 1.25e-3
arc_length_radius_initial 1.35e-2
arc_length_radius_min 1.25e-5
arc_length_radius_max 1.35e-2
arc_length_psi 1
arc_length_shrink 0.5
arc_length_expand 1
arc_length_target_iterations 8
arc_length_max_iterations 300
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 1
""",
}


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def run_git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def selected_variants(selection: str | None) -> list[str]:
    if not selection:
        return list(VARIANTS)
    names = [name.strip() for name in selection.split(",") if name.strip()]
    unknown = [name for name in names if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")
    return names


def solver_for(variant: str, total_time: float | None) -> str:
    solver = BASELINE_SOLVER
    if total_time is not None:
        solver += f"\ntotal_time {total_time:.12g}\n"
    return solver + "\n" + VARIANTS[variant].strip() + "\n"


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
        if match:
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


def run_variant(root: Path, deck: Path, exe: Path, jobs: int, timeout: int, variant: str, total_time: float | None) -> dict[str, Any]:
    run_dir = root / "runs" / variant
    prepare_run_dir(deck, run_dir, solver_for(variant, total_time), copy_deck=False)
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
        "variant": variant,
        "run_dir": str(run_dir),
        "exit_status": exit_status,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "summary": summary,
        "arc_rows": arc_rows,
    }


def verdict(result: dict[str, Any], target_time: float) -> str:
    summary = result["summary"]
    arc_rows = result["arc_rows"]
    warnings = summary.get("warnings", [])
    if any("nonzero proportional reference load" in warning for warning in warnings):
        return "unsupported_zero_reference"
    if any("Error:" in warning for warning in warnings):
        return "failed"
    if result["exit_status"] != 0 or result["timed_out"]:
        return "failed"
    if not summary.get("end_of_calculation"):
        return "failed"
    if summary.get("fallback_acceptance_count", 0) != 0 or summary.get("nan_count", 0) != 0:
        return "failed"
    final_lambda = arc_rows[-1]["lambda"] if arc_rows else (summary.get("steps") or [{}])[-1].get("time", 0.0)
    if final_lambda is None or final_lambda + 1e-10 < target_time:
        return "partial"
    rows = summary.get("total_iters", 0)
    if rows < TSN65_TARGET_TOTAL:
        return "improves"
    if rows == TSN65_TARGET_TOTAL:
        return "neutral"
    return "worsens"


def aggregate_rows(results: list[dict[str, Any]], target_time: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        summary = result["summary"]
        arc_rows = result["arc_rows"]
        steps = summary.get("steps") or []
        final_lambda = arc_rows[-1]["lambda"] if arc_rows else ""
        final_radius = arc_rows[-1]["radius"] if arc_rows else ""
        min_radius = min((row["radius"] for row in arc_rows), default="")
        min_dlambda = min((abs(row["dlambda"]) for row in arc_rows if row["dlambda"] != 0), default="")
        rows.append(
            {
                "variant": result["variant"],
                "exit_status": result["exit_status"],
                "verdict": verdict(result, target_time),
                "steps": summary.get("step_count", 0),
                "rows": summary.get("total_iters", 0),
                "duration": summary.get("total_duration") or "",
                "elapsed_seconds": f"{result['elapsed_seconds']:.3f}",
                "warnings": summary.get("warning_count", 0),
                "nans": summary.get("nan_count", 0),
                "fallback": summary.get("fallback_acceptance_count", 0),
                "cutbacks": summary.get("cutback_count", 0),
                "arc_rows": len(arc_rows),
                "final_lambda": final_lambda,
                "final_radius": final_radius,
                "min_radius": min_radius,
                "min_dlambda": min_dlambda,
                "last_step": steps[-1]["step"] if steps else "",
                "last_iter": steps[-1]["last_iter"] if steps else "",
                "last_residual": steps[-1]["residual"] if steps else "",
                "last_displacement": steps[-1]["displacement"] if steps else "",
                "last_energy": steps[-1]["energy"] if steps else "",
            }
        )
    return rows


def write_report(root: Path, repo: Path, deck: Path, exe: Path, rows: list[dict[str, Any]], variants: list[str], target_time: float) -> None:
    (root / "metadata.txt").write_text(
        f"created_at={now_iso()}\n"
        f"repo={repo}\n"
        f"branch={run_git(repo, 'branch', '--show-current')}\n"
        f"commit={run_git(repo, 'rev-parse', 'HEAD')}\n"
        f"git_status_short<<EOF\n{run_git(repo, 'status', '--short')}\nEOF\n"
        f"executable={exe}\n"
        f"deck={deck}\n"
        f"target_time={target_time:.12g}\n"
        f"variants={','.join(variants)}\n",
        encoding="utf-8",
    )

    headers = [
        "variant",
        "exit_status",
        "verdict",
        "steps",
        "rows",
        "duration",
        "elapsed_seconds",
        "warnings",
        "nans",
        "fallback",
        "cutbacks",
        "arc_rows",
        "final_lambda",
        "final_radius",
        "min_radius",
        "min_dlambda",
        "last_step",
        "last_iter",
        "last_residual",
        "last_displacement",
        "last_energy",
    ]
    with (root / "cp5_arc_length.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[h]) for h in headers) + "\n")

    lines = [
        "# CP5 TS-N65 Arc-Length Experiments",
        "",
        f"Generated: {now_iso()}",
        "",
        "## Baseline Target",
        "",
        f"- Final pseudo-time / load parameter: `{target_time:.12g}`",
        f"- Strict baseline iteration sequence: `{','.join(str(i) for i in TSN65_TARGET_ITERS)}`",
        f"- Strict baseline total rows: `{TSN65_TARGET_TOTAL}`",
        "",
        "## Result Table",
        "",
        "| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | final lambda | min radius | min dlambda |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {variant} | {exit_status} | {verdict} | {steps} | {rows} | {duration} | {warnings} | {nans} | {fallback} | {cutbacks} | {final_lambda} | {min_radius} | {min_dlambda} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `CP5-arc-proportional-smoke` checks the CP4 proportional-load arc-length path directly on TS-N65.",
            "- `CP5-arc-fd-*` uses `arc_length_reference finite_difference`, which derives the reference direction from residual changes caused by the prescribed displacement/load parameter.",
            "- A variant is only promotable if it reaches the same final parameter as the strict `0.01` baseline with no fallback acceptance or NaNs and fewer rows or less wall time.",
            "",
            "## Files",
            "",
            "- `cp5_arc_length.tsv`: compact aggregate table.",
            "- `runs/<variant>/solver.out`: ignored full solver log.",
            "- `runs/<variant>/arc_length.tsv`: ignored arc-length iteration trace.",
        ]
    )
    (root / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--deck", type=Path, default=Path("results/tsn65-zip-runs-20260519-190000/input/TS-N_65"))
    parser.add_argument("--exe", type=Path, default=Path("/tmp/oas_tsn65_full_baseline_build/bin/OAS"))
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--variants", default=None, help="Comma-separated variant names")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--total-time", type=float, default=1e-2)
    args = parser.parse_args()

    repo = args.repo.resolve()
    deck = (repo / args.deck).resolve() if not args.deck.is_absolute() else args.deck.resolve()
    exe = args.exe.resolve()
    variants = selected_variants(args.variants)
    root = args.root or repo / "results" / f"tsn65-arc-length-cp5-{now_stamp()}"
    root.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, Any]] = []
    for variant in variants:
        print(f"[cp5] running {variant}", flush=True)
        results.append(run_variant(root, deck, exe, args.jobs, args.timeout, variant, args.total_time))
        rows = aggregate_rows(results, args.total_time)
        write_report(root, repo, deck, exe, rows, variants, args.total_time)

    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
