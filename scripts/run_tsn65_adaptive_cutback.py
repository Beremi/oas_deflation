#!/usr/bin/env python3
"""Run the CP2 TS-N65 adaptive cutback/stagnation experiments."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_tsn65_globalization_sweep import BASELINE_SOLVER, prepare_run_dir, run_git  # noqa: E402
from summarize_oas_run import TSN65_TARGET_TOTAL, format_tsv, parse_log  # noqa: E402


ADAPTIVE_CONTROLS = """max_iterations 40
min_time_step 1.25e-5
max_time_step 1.25e-3
critical_step_decrease 0.5
step_decrease 0.5
step_increase 1.2
enlargeIt 6
shortenIt 15
limit_tolerance 0
nonlinear_line_search_cutback_on_fail 1
nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25
"""


VARIANTS: dict[str, str] = {
    "G2-stagnation-only": "",
    "G2-backtracking-actual": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
""",
    "G2-backtracking-frozen-actual": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
""",
    "G2-backtracking-actual-adaptive-K": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_adaptive_matrix_update 1
nonlinear_rebuild_on_small_alpha 0.5
nonlinear_rebuild_on_merit_growth 1
""",
    "G2-fixed-damping-05": """
nonlinear_damping_type fixed
nonlinear_damping_factor 0.5
""",
}


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def selected_variants(selection: str | None) -> list[str]:
    if not selection:
        return list(VARIANTS)
    names = [name.strip() for name in selection.split(",") if name.strip()]
    unknown = [name for name in names if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")
    return names


def adaptive_solver_text(extra: str) -> str:
    lines: list[str] = []
    for raw in BASELINE_SOLVER.splitlines():
        key = raw.split(None, 1)[0] if raw.strip() else ""
        if key in {
            "max_iterations",
            "min_time_step",
            "max_time_step",
            "limit_tolerance",
        }:
            continue
        lines.append(raw)
    solver = "\n".join(lines).rstrip() + "\n\n" + ADAPTIVE_CONTROLS.rstrip() + "\n"
    extra = extra.strip()
    if extra:
        solver += "\n" + extra + "\n"
    return solver


def tail_rows(summary: dict[str, Any]) -> int:
    return sum(step["iters"] for step in summary["steps"][5:8])


def accepted_steps(summary: dict[str, Any]) -> int:
    return sum(1 for step in summary["steps"] if step.get("converged_by_tolerance"))


def min_dt(summary: dict[str, Any]) -> float | None:
    reported = summary.get("min_reported_dt")
    if reported is not None:
        return float(reported)
    values = [step.get("dt") for step in summary["steps"] if step.get("dt") is not None]
    return min(values) if values else None


def verdict(summary: dict[str, Any], exit_status: int | None) -> str:
    if exit_status is None:
        return "not_run"
    if exit_status != 0:
        return "failed"
    if summary["nan_count"] or summary["fallback_acceptance_count"]:
        return "failed"
    if not summary["end_of_calculation"]:
        return "failed"
    final_time = summary["steps"][-1]["time"] if summary["steps"] else 0.0
    if final_time < 9.999e-3:
        return "failed"
    if summary["total_iters"] < TSN65_TARGET_TOTAL:
        return "improves"
    if summary["total_iters"] == TSN65_TARGET_TOTAL:
        return "neutral"
    return "worsens"


def read_exit_status(run_dir: Path) -> int | None:
    timing = run_dir / "wrapper_timing.txt"
    if not timing.exists():
        return None
    for line in timing.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("exit_status="):
            value = line.split("=", 1)[1].strip()
            try:
                return int(value)
            except ValueError:
                return None
    return None


def ordered_rows(rows_by_variant: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [rows_by_variant[name] for name in VARIANTS if name in rows_by_variant]


def collect_existing_rows(run_root: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for variant in VARIANTS:
        run_dir = run_root / variant
        log = run_dir / "solver.out"
        if not log.exists():
            continue
        summary = parse_log(log, "tsn65")
        exit_status = read_exit_status(run_dir)
        rows[variant] = {
            "variant": variant,
            "run_dir": str(run_dir),
            "exit_status": exit_status,
            "summary": summary,
            "verdict": verdict(summary, exit_status),
        }
    return rows


def write_report(root: Path, rows: list[dict[str, Any]], metadata: dict[str, str]) -> None:
    with (root / "cp2_adaptive.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "variant\texit_status\tverdict\tsteps\taccepted_steps\ttotal_iters\ttail_6_8_iters\tmin_dt\tduration\twarnings\tnans\tcutbacks\trestarts\tshortens\tfallback_accepts\tmin_alpha\tmax_ls_trials\n"
        )
        for row in rows:
            summary = row.get("summary") or {}
            steps = summary.get("steps") or []
            min_alpha_values = [step["min_alpha"] for step in steps if step.get("min_alpha") is not None]
            max_trials = max([step.get("max_ls_trials", 0) for step in steps] or [0])
            value_min_dt = min_dt(summary) if steps else None
            handle.write(
                "\t".join(
                    [
                        row["variant"],
                        "" if row.get("exit_status") is None else str(row["exit_status"]),
                        row.get("verdict", "not_run"),
                        str(summary.get("step_count", 0)),
                        str(accepted_steps(summary)) if steps else "0",
                        str(summary.get("total_iters", 0)),
                        str(tail_rows(summary)) if steps else "",
                        "" if value_min_dt is None else f"{value_min_dt:.12g}",
                        str(summary.get("total_duration") or ""),
                        str(summary.get("warning_count", 0)),
                        str(summary.get("nan_count", 0)),
                        str(summary.get("cutback_count", 0)),
                        str(summary.get("restart_count", 0)),
                        str(summary.get("shortening_count", 0)),
                        str(summary.get("fallback_acceptance_count", 0)),
                        "" if not min_alpha_values else f"{min(min_alpha_values):.12g}",
                        str(max_trials),
                    ]
                )
                + "\n"
            )

    lines = [
        "# CP2 TS-N65 Adaptive Cutback/Stagnation",
        "",
        f"Generated: {now_iso()}",
        "",
        "## Metadata",
        "",
        "```text",
    ]
    lines.extend(f"{key}={value}" for key, value in metadata.items())
    lines.extend(
        [
            "```",
            "",
            "## Adaptive Deck Override",
            "",
            "```text",
            ADAPTIVE_CONTROLS.rstrip(),
            "```",
            "",
            "## Result Table",
            "",
            "| variant | exit | verdict | steps | accepted | rows | tail 6-8 | min dt | duration | warnings | NaNs | cutbacks | restarts | shortens | fallback | min alpha | max ls trials |",
            "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        summary = row.get("summary") or {}
        steps = summary.get("steps") or []
        min_alpha_values = [step["min_alpha"] for step in steps if step.get("min_alpha") is not None]
        max_trials = max([step.get("max_ls_trials", 0) for step in steps] or [0])
        value_min_dt = min_dt(summary) if steps else None
        lines.append(
            "| {variant} | {exit_status} | {verdict} | {steps} | {accepted} | {rows} | {tail} | {min_dt} | {duration} | {warnings} | {nans} | {cutbacks} | {restarts} | {shortens} | {fallback} | {alpha} | {trials} |".format(
                variant=row["variant"],
                exit_status="" if row.get("exit_status") is None else row["exit_status"],
                verdict=row.get("verdict", "not_run"),
                steps=summary.get("step_count", 0),
                accepted=accepted_steps(summary) if steps else 0,
                rows=summary.get("total_iters", 0),
                tail=tail_rows(summary) if steps else "",
                min_dt="" if value_min_dt is None else f"{value_min_dt:.6g}",
                duration=summary.get("total_duration") or "",
                warnings=summary.get("warning_count", 0),
                nans=summary.get("nan_count", 0),
                cutbacks=summary.get("cutback_count", 0),
                restarts=summary.get("restart_count", 0),
                shortens=summary.get("shortening_count", 0),
                fallback=summary.get("fallback_acceptance_count", 0),
                alpha="" if not min_alpha_values else f"{min(min_alpha_values):.6g}",
                trials=max_trials,
            )
        )

    executed = [row for row in rows if row.get("summary")]
    passed = [row for row in executed if row.get("verdict") in {"improves", "neutral", "worsens"}]
    failed = [row for row in executed if row.get("verdict") == "failed"]
    furthest_time = 0.0
    for row in executed:
        for step in row["summary"].get("steps", []):
            if step.get("converged_by_tolerance") and step.get("time") is not None:
                furthest_time = max(furthest_time, float(step["time"]))

    lines.extend(["", "## Checkpoint Verdict", ""])
    if passed:
        lines.append(
            "At least one executed CP2 variant reached the target time without NaNs or fallback acceptance. "
            "Promote only rows marked `improves` after reviewing the load/displacement curve."
        )
    elif failed:
        lines.append(
            "CP2 did not pass for the executed variants. The runs avoided NaNs and fallback acceptance, "
            "but reached `min_time_step` and failed before `total_time=1e-2`."
        )
        lines.append("")
        lines.append(
            f"Furthest accepted physical time across executed variants: `{furthest_time:.8g}`. "
            "This is a useful negative checkpoint and points the next experiment toward CP3 "
            "indirect displacement control before implementing arc-length."
        )
    else:
        lines.append("No CP2 variants have been executed yet.")

    if executed:
        lines.extend(["", "## Per-Variant Step Tables", ""])
    for row in executed:
        summary = row["summary"]
        lines.extend(
            [
                f"### {row['variant']}",
                "",
                "| step | time | dt | rows | converged | residual | displacement | energy | min alpha | ls trials | cutback reasons |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
            ]
        )
        for step in summary.get("steps", []):
            lines.append(
                "| {step} | {time} | {dt} | {iters} | {conv} | {res} | {disp} | {ene} | {alpha} | {trials} | {reasons} |".format(
                    step=step.get("step", ""),
                    time="" if step.get("time") is None else f"{step['time']:.8g}",
                    dt="" if step.get("dt") is None else f"{step['dt']:.8g}",
                    iters=step.get("iters", ""),
                    conv=int(bool(step.get("converged_by_tolerance"))),
                    res="" if step.get("residual") is None else f"{step['residual']:.6g}",
                    disp="" if step.get("displacement") is None else f"{step['displacement']:.6g}",
                    ene="" if step.get("energy") is None else f"{step['energy']:.6g}",
                    alpha="" if step.get("min_alpha") is None else f"{step['min_alpha']:.6g}",
                    trials=step.get("max_ls_trials", ""),
                    reasons=", ".join(step.get("cutback_reasons") or []),
                )
            )
        lines.append("")

    lines.extend(
        [
            "",
            "## Baseline Comparison",
            "",
            "- Strict target sequence: `6,6,10,13,17,183,187,163`.",
            "- Strict target rows: `585`.",
            "- CP2 pass requires reaching `total_time=1e-2` without fallback acceptance or NaNs.",
            "",
            "## Files",
            "",
            "- `cp2_adaptive.tsv`: machine-readable aggregate table.",
            "- `runs/<variant>/solver.out`: ignored full logs.",
            "- `runs/<variant>/summary.json`: ignored parsed full summary.",
            "- `runs/<variant>/step_summary.tsv`: ignored parsed step table.",
        ]
    )
    (root / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--deck", type=Path, default=Path("results/tsn65-zip-runs-20260519-190000/input/TS-N_65"))
    parser.add_argument("--exe", type=Path, default=Path("/tmp/oas_tsn65_full_baseline_build/bin/OAS"))
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--variants", help="Comma-separated variant names. Default: all CP2 variants.")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    parser.add_argument("--copy-deck", action="store_true", help="Copy deck files instead of symlinking them.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--summarize-existing", action="store_true", help="Rewrite aggregate reports from existing run logs only.")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    deck = (repo / args.deck).resolve() if not args.deck.is_absolute() else args.deck
    exe = args.exe.resolve()
    root = args.root or repo / f"results/tsn65-adaptive-cutback-{dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
    root = root.resolve()
    run_root = root / "runs"
    run_root.mkdir(parents=True, exist_ok=True)

    variants = selected_variants(args.variants)
    metadata = {
        "created_at": now_iso(),
        "repo": str(repo),
        "branch": run_git(["branch", "--show-current"], repo),
        "commit": run_git(["rev-parse", "HEAD"], repo),
        "executable": str(exe),
        "deck": str(deck),
        "variants": ",".join(variants),
    }
    (root / "metadata.txt").write_text("\n".join(f"{k}={v}" for k, v in metadata.items()) + "\n", encoding="utf-8")

    rows_by_variant = collect_existing_rows(run_root)
    if args.summarize_existing:
        write_report(root, ordered_rows(rows_by_variant), metadata)
        print(root)
        return 0

    for variant in variants:
        run_dir = run_root / variant
        prepare_run_dir(deck, run_dir, adaptive_solver_text(VARIANTS[variant]), args.copy_deck)

        command = [str(exe), "-j", str(args.jobs), "master.inp"]
        env = os.environ.copy()
        env.update(
            {
                "OMP_NUM_THREADS": str(args.jobs),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            }
        )
        (run_dir / "run_command.txt").write_text(
            " ".join(f"{key}={env[key]}" for key in ["OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES"])
            + " "
            + " ".join(command)
            + "\n",
            encoding="utf-8",
        )

        row: dict[str, Any] = {"variant": variant, "run_dir": str(run_dir)}
        if args.dry_run:
            row["exit_status"] = None
            row["verdict"] = "not_run"
            rows_by_variant[variant] = row
            write_report(root, ordered_rows(rows_by_variant), metadata)
            continue

        start = time.monotonic()
        start_iso = now_iso()
        timed_out = False
        with (run_dir / "solver.out").open("w", encoding="utf-8", errors="replace") as log:
            try:
                completed = subprocess.run(
                    command,
                    cwd=run_dir,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=args.timeout_seconds,
                    check=False,
                )
                exit_status = completed.returncode
            except subprocess.TimeoutExpired:
                exit_status = 124
                timed_out = True
                log.write(f"\nRUNNER_TIMEOUT after {args.timeout_seconds} seconds\n")
        elapsed = time.monotonic() - start
        (run_dir / "wrapper_timing.txt").write_text(
            f"start={start_iso}\nend={now_iso()}\nelapsed_seconds={elapsed:.3f}\nexit_status={exit_status}\ntimed_out={int(timed_out)}\n",
            encoding="utf-8",
        )

        summary = parse_log(run_dir / "solver.out", "tsn65")
        (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (run_dir / "step_summary.tsv").write_text(format_tsv(summary), encoding="utf-8")
        row["exit_status"] = exit_status
        row["summary"] = summary
        row["verdict"] = verdict(summary, exit_status)
        rows_by_variant[variant] = row
        write_report(root, ordered_rows(rows_by_variant), metadata)

    write_report(root, ordered_rows(rows_by_variant), metadata)
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
