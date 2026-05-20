#!/usr/bin/env python3
"""Run the CP3 TS-N65 indirect-displacement-control checkpoint."""

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
from summarize_oas_run import format_tsv, parse_log  # noqa: E402


TARGET_TOTAL_TIME = 1.0e-2
TARGET_FIRST_TIME = 1.25e-3
TARGET_V01_FINAL = 1.519697e-5
TARGET_SCALE_FROM_UNSCALED_STRICT = 1.2321920851549628
BASELINE_FIRST = {
    "time": 1.25e-3,
    "displ[m]": -6.0e-5,
    "load[N]": -3.740127e4,
    "v01": 2.613423e-6,
}
BASELINE_FINAL = {
    "time": 1.0e-2,
    "displ[m]": -4.8e-4,
    "load[N]": -2.169625e5,
    "v01": TARGET_V01_FINAL,
}

IDC_BLOCK_TEMPLATE = """indirect_control 2
ic_xcoords 2.1 2.1
ic_ycoords 0 0
ic_zcoords 0 -0.4
ic_directions 2 2
ic_displ_weights {weights}
ic_force_weights 0 0
ic_function 4
"""

VARIANTS: dict[str, dict[str, Any]] = {
    "CP3-idc-minus-calibration": {
        "weights": "-1 1",
        "total_time": TARGET_FIRST_TIME,
        "purpose": "sign calibration, roadmap weight order",
    },
    "CP3-idc-plus-calibration": {
        "weights": "1 -1",
        "total_time": TARGET_FIRST_TIME,
        "purpose": "sign calibration, flipped weight order",
    },
    "CP3-idc-minus-strict": {
        "weights": "-1 1",
        "total_time": TARGET_TOTAL_TIME,
        "target_scale": 1.0,
        "purpose": "full strict IDC run, roadmap weight order",
    },
    "CP3-idc-plus-strict": {
        "weights": "1 -1",
        "total_time": TARGET_TOTAL_TIME,
        "target_scale": 1.0,
        "purpose": "full strict IDC run, flipped weight order",
    },
    "CP3-idc-minus-final-v01-calibrated-strict": {
        "weights": "-1 1",
        "total_time": TARGET_TOTAL_TIME,
        "target_scale": TARGET_SCALE_FROM_UNSCALED_STRICT,
        "purpose": "full strict IDC run with target function scaled from unscaled strict final-v01 miss",
    },
}


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def selected_variants(selection: str | None) -> list[str]:
    if not selection:
        return ["benchmark", "calibration", "selected-strict"]
    names = [name.strip() for name in selection.split(",") if name.strip()]
    allowed = set(VARIANTS) | {"benchmark", "calibration", "selected-strict"}
    unknown = [name for name in names if name not in allowed]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")
    return names


def solver_with_total_time(total_time: float, idc_block: str) -> str:
    lines: list[str] = []
    for raw in BASELINE_SOLVER.splitlines():
        key = raw.split(None, 1)[0] if raw.strip() else ""
        if key == "total_time":
            lines.append(f"total_time {total_time:.12g}")
        else:
            lines.append(raw)
    return "\n".join(lines).rstrip() + "\n\n" + idc_block.rstrip() + "\n"


def append_control_function(run_dir: Path, target_scale: float = 1.0) -> int:
    original = run_dir / "functions.inp"
    if original.is_symlink():
        source = original.resolve()
        original.unlink()
        text = source.read_text(encoding="utf-8", errors="replace")
    else:
        text = original.read_text(encoding="utf-8", errors="replace")
    if not text.endswith("\n"):
        text += "\n"
    target = TARGET_V01_FINAL * target_scale
    text += f"PWLFunction 2 0 {TARGET_TOTAL_TIME:.12g} 0 {target:.12g} # CP3 IDC v01 target, scale={target_scale:.12g}\n"
    original.write_text(text, encoding="utf-8")
    return len([line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]) - 1


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


def write_run_metadata(run_dir: Path, command: list[str], env: dict[str, str], start_iso: str, elapsed: float, exit_status: int, timed_out: bool) -> None:
    (run_dir / "run_command.txt").write_text(
        " ".join(f"{key}={env[key]}" for key in ["OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES"])
        + " "
        + " ".join(command)
        + "\n",
        encoding="utf-8",
    )
    (run_dir / "wrapper_timing.txt").write_text(
        f"start={start_iso}\nend={now_iso()}\nelapsed_seconds={elapsed:.3f}\nexit_status={exit_status}\ntimed_out={int(timed_out)}\n",
        encoding="utf-8",
    )


def read_ld(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"header": [], "rows": [], "first": {}, "last": {}}
    header: list[str] = []
    rows: list[dict[str, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.strip():
            continue
        if raw.startswith("#"):
            header = raw.lstrip("#").split()
            continue
        parts = raw.split()
        if not header or len(parts) < len(header):
            continue
        row: dict[str, float] = {}
        for key, value in zip(header, parts):
            try:
                row[key] = float(value)
            except ValueError:
                pass
        if row:
            rows.append(row)
    return {"header": header, "rows": rows, "first": rows[0] if rows else {}, "last": rows[-1] if rows else {}}


def sign(value: float | None) -> int:
    if value is None:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def sign_score(first: dict[str, float]) -> int:
    score = 0
    for key in ["displ[m]", "load[N]", "v01"]:
        if sign(first.get(key)) == sign(BASELINE_FIRST[key]):
            score += 1
    return score


def v01_relative_error(last: dict[str, float]) -> float | None:
    if "v01" not in last:
        return None
    return abs(last["v01"] - TARGET_V01_FINAL) / max(abs(TARGET_V01_FINAL), 1e-30)


def run_benchmark(root: Path, repo: Path, exe: Path, jobs: int, timeout_seconds: int) -> dict[str, Any]:
    run_dir = root / "runs" / "CP3-indirect-control-benchmark"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(repo / "src/benchmark/Indirect_Control", run_dir)
    env = make_env(jobs)
    command = [str(exe), "-j", str(jobs), "master.inp"]
    start_iso = now_iso()
    exit_status, elapsed, timed_out = run_command(command, run_dir, env, timeout_seconds)
    write_run_metadata(run_dir, command, env, start_iso, elapsed, exit_status, timed_out)
    summary = parse_log(run_dir / "solver.out")
    ld = read_ld(run_dir / "results/LD.out")
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (run_dir / "step_summary.tsv").write_text(format_tsv(summary), encoding="utf-8")
    return {
        "variant": "CP3-indirect-control-benchmark",
        "kind": "benchmark",
        "run_dir": str(run_dir),
        "exit_status": exit_status,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "summary": summary,
        "ld": ld,
        "verdict": "passed" if exit_status == 0 and summary.get("end_of_calculation") else "failed",
    }


def run_tsn65_variant(root: Path, deck: Path, exe: Path, jobs: int, timeout_seconds: int, variant: str) -> dict[str, Any]:
    spec = VARIANTS[variant]
    run_dir = root / "runs" / variant
    if run_dir.exists():
        shutil.rmtree(run_dir)
    idc_block = IDC_BLOCK_TEMPLATE.format(weights=spec["weights"])
    prepare_run_dir(deck, run_dir, solver_with_total_time(spec["total_time"], idc_block), copy_deck=False)
    target_scale = float(spec.get("target_scale", 1.0))
    function_index = append_control_function(run_dir, target_scale)
    env = make_env(jobs)
    command = [str(exe), "-j", str(jobs), "master.inp"]
    start_iso = now_iso()
    exit_status, elapsed, timed_out = run_command(command, run_dir, env, timeout_seconds)
    write_run_metadata(run_dir, command, env, start_iso, elapsed, exit_status, timed_out)
    summary = parse_log(run_dir / "solver.out", "tsn65")
    ld = read_ld(run_dir / "results/LD.out")
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (run_dir / "step_summary.tsv").write_text(format_tsv(summary), encoding="utf-8")
    first = ld["first"]
    last = ld["last"]
    final_time = last.get("time")
    reached_target_time = bool(final_time is not None and final_time >= spec["total_time"] * 0.999)
    rel_err = v01_relative_error(last)
    return {
        "variant": variant,
        "kind": "tsn65",
        "purpose": spec["purpose"],
        "weights": spec["weights"],
        "total_time_target": spec["total_time"],
        "function_index": function_index,
        "target_scale": target_scale,
        "run_dir": str(run_dir),
        "exit_status": exit_status,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "summary": summary,
        "ld": ld,
        "sign_score": sign_score(first),
        "v01_relative_error": rel_err,
    }
    if variant.endswith("-calibration"):
        result["verdict"] = (
            "passed"
            if exit_status == 0
            and summary.get("end_of_calculation")
            and result["sign_score"] > 0
            and summary.get("nan_count", 0) == 0
            and summary.get("fallback_acceptance_count", 0) == 0
            else "failed"
        )
    else:
        result["verdict"] = (
            "passed"
            if exit_status == 0
            and summary.get("end_of_calculation")
            and reached_target_time
            and summary.get("nan_count", 0) == 0
            and summary.get("fallback_acceptance_count", 0) == 0
            and (rel_err is None or rel_err < 5e-2)
            else "failed"
        )
    return result


def collect_existing_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for run_dir in sorted((root / "runs").glob("CP3-*")):
        log = run_dir / "solver.out"
        if not log.exists():
            continue
        timing = run_dir / "wrapper_timing.txt"
        exit_status = None
        elapsed = 0.0
        timed_out = False
        if timing.exists():
            for line in timing.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("exit_status="):
                    try:
                        exit_status = int(line.split("=", 1)[1])
                    except ValueError:
                        pass
                elif line.startswith("elapsed_seconds="):
                    try:
                        elapsed = float(line.split("=", 1)[1])
                    except ValueError:
                        pass
                elif line.startswith("timed_out="):
                    timed_out = line.split("=", 1)[1].strip() == "1"
        is_tsn65 = run_dir.name in VARIANTS
        summary = parse_log(log, "tsn65" if is_tsn65 else None)
        ld = read_ld(run_dir / "results/LD.out")
        row: dict[str, Any] = {
            "variant": run_dir.name,
            "kind": "tsn65" if is_tsn65 else "benchmark",
            "run_dir": str(run_dir),
            "exit_status": exit_status,
            "timed_out": timed_out,
            "elapsed_seconds": elapsed,
            "summary": summary,
            "ld": ld,
            "verdict": "passed" if exit_status == 0 and summary.get("end_of_calculation") else "failed",
        }
        if is_tsn65:
            target_scale = float(VARIANTS[run_dir.name].get("target_scale", 1.0))
            final_time = ld["last"].get("time")
            reached_target_time = bool(final_time is not None and final_time >= VARIANTS[run_dir.name]["total_time"] * 0.999)
            rel_err = v01_relative_error(ld["last"])
            row.update(
                {
                    "purpose": VARIANTS[run_dir.name]["purpose"],
                    "weights": VARIANTS[run_dir.name]["weights"],
                    "total_time_target": VARIANTS[run_dir.name]["total_time"],
                    "target_scale": target_scale,
                    "sign_score": sign_score(ld["first"]),
                    "v01_relative_error": rel_err,
                }
            )
            if run_dir.name.endswith("-calibration"):
                row["verdict"] = (
                    "passed"
                    if exit_status == 0
                    and summary.get("end_of_calculation")
                    and row["sign_score"] > 0
                    and summary.get("nan_count", 0) == 0
                    and summary.get("fallback_acceptance_count", 0) == 0
                    else "failed"
                )
            else:
                row["verdict"] = (
                    "passed"
                    if exit_status == 0
                    and summary.get("end_of_calculation")
                    and reached_target_time
                    and summary.get("nan_count", 0) == 0
                    and summary.get("fallback_acceptance_count", 0) == 0
                    and (rel_err is None or rel_err < 5e-2)
                    else "failed"
                )
        rows.append(row)
    return rows


def choose_sign(rows: list[dict[str, Any]]) -> str | None:
    calibration = [row for row in rows if row["variant"].endswith("-calibration") and row.get("exit_status") == 0]
    if not calibration:
        return None
    calibration.sort(key=lambda row: (row.get("sign_score", 0), -abs(row["ld"]["first"].get("v01", 0.0))), reverse=True)
    return "minus" if "minus" in calibration[0]["variant"] else "plus"


def write_report(root: Path, rows: list[dict[str, Any]], metadata: dict[str, str]) -> None:
    rows_by_name = {row["variant"]: row for row in rows}
    ordered_names = [
        "CP3-indirect-control-benchmark",
        "CP3-idc-minus-calibration",
        "CP3-idc-plus-calibration",
        "CP3-idc-minus-strict",
        "CP3-idc-plus-strict",
        "CP3-idc-minus-final-v01-calibrated-strict",
    ]
    rows = [rows_by_name[name] for name in ordered_names if name in rows_by_name]
    selected = choose_sign(rows)

    with (root / "cp3_indirect_control.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "variant\texit_status\tverdict\tsteps\trows\tduration\twarnings\tnans\tfallback\tcutbacks\tmin_dt\tidc_weights\tsign_score\tfinal_time\tfinal_idc_time\tfinal_displ\tfinal_load\tfinal_v01\tv01_rel_err\tlast_step\tlast_iter\tlast_residual\tlast_displacement\tlast_energy\n"
        )
        for row in rows:
            summary = row["summary"]
            last = row["ld"]["last"]
            last_step = summary.get("steps", [])[-1] if summary.get("steps") else {}
            handle.write(
                "\t".join(
                    [
                        row["variant"],
                        "" if row.get("exit_status") is None else str(row["exit_status"]),
                        row.get("verdict", ""),
                        str(summary.get("step_count", 0)),
                        str(summary.get("total_iters", 0)),
                        str(summary.get("total_duration") or ""),
                        str(summary.get("warning_count", 0)),
                        str(summary.get("nan_count", 0)),
                        str(summary.get("fallback_acceptance_count", 0)),
                        str(summary.get("cutback_count", 0)),
                        "" if summary.get("min_reported_dt") is None else f"{summary['min_reported_dt']:.12g}",
                        row.get("weights", ""),
                        "" if row.get("sign_score") is None else str(row.get("sign_score")),
                        "" if "time" not in last else f"{last['time']:.12g}",
                        "" if "idc_time" not in last else f"{last['idc_time']:.12g}",
                        "" if "displ[m]" not in last else f"{last['displ[m]']:.12g}",
                        "" if "load[N]" not in last else f"{last['load[N]']:.12g}",
                        "" if "v01" not in last else f"{last['v01']:.12g}",
                        "" if row.get("v01_relative_error") is None else f"{row['v01_relative_error']:.12g}",
                        "" if "step" not in last_step else str(last_step["step"]),
                        "" if "last_iter" not in last_step else str(last_step["last_iter"]),
                        "" if last_step.get("residual") is None else f"{last_step['residual']:.12g}",
                        "" if last_step.get("displacement") is None else f"{last_step['displacement']:.12g}",
                        "" if last_step.get("energy") is None else f"{last_step['energy']:.12g}",
                    ]
                )
                + "\n"
            )

    lines = [
        "# CP3 TS-N65 Indirect Displacement Control",
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
            "## Function Index Note",
            "",
            "OAS function indices are zero-based. The handout block listed `ic_function 2`, "
            "but in the TS-N65 deck function index `2` is the steel material envelope, not the "
            "displacement ramp. CP3 appends a local IDC target function and uses `ic_function 4` "
            "inside each ignored run directory.",
            "",
            "The appended target is:",
            "",
            "```text",
            f"PWLFunction 2 0 {TARGET_TOTAL_TIME:.12g} 0 {TARGET_V01_FINAL:.12g} # CP3 IDC v01 target",
            "```",
            "",
            "The final-v01 calibrated strict run scales that target by "
            f"`{TARGET_SCALE_FROM_UNSCALED_STRICT:.12g}`, derived from "
            "`target_v01 / CP3-idc-minus-strict_final_v01`.",
            "",
            "## Baseline Gauges",
            "",
            "| point | time | displ[m] | load[N] | v01 |",
            "| --- | ---: | ---: | ---: | ---: |",
            f"| first strict step | {BASELINE_FIRST['time']:.8g} | {BASELINE_FIRST['displ[m]']:.8g} | {BASELINE_FIRST['load[N]']:.8g} | {BASELINE_FIRST['v01']:.8g} |",
            f"| strict target | {BASELINE_FINAL['time']:.8g} | {BASELINE_FINAL['displ[m]']:.8g} | {BASELINE_FINAL['load[N]']:.8g} | {BASELINE_FINAL['v01']:.8g} |",
            "",
            "## Result Table",
            "",
            "| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | min dt | weights | sign score | final time | idc time | displ[m] | load[N] | v01 | v01 rel err |",
            "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        summary = row["summary"]
        last = row["ld"]["last"]
        lines.append(
            "| {variant} | {exit_status} | {verdict} | {steps} | {rows} | {duration} | {warnings} | {nans} | {fallback} | {cutbacks} | {min_dt} | {weights} | {score} | {time} | {idc_time} | {displ} | {load} | {v01} | {err} |".format(
                variant=row["variant"],
                exit_status="" if row.get("exit_status") is None else row["exit_status"],
                verdict=row.get("verdict", ""),
                steps=summary.get("step_count", 0),
                rows=summary.get("total_iters", 0),
                duration=summary.get("total_duration") or "",
                warnings=summary.get("warning_count", 0),
                nans=summary.get("nan_count", 0),
                fallback=summary.get("fallback_acceptance_count", 0),
                cutbacks=summary.get("cutback_count", 0),
                min_dt="" if summary.get("min_reported_dt") is None else f"{summary['min_reported_dt']:.6g}",
                weights=row.get("weights", ""),
                score="" if row.get("sign_score") is None else row.get("sign_score"),
                time="" if "time" not in last else f"{last['time']:.8g}",
                idc_time="" if "idc_time" not in last else f"{last['idc_time']:.8g}",
                displ="" if "displ[m]" not in last else f"{last['displ[m]']:.8g}",
                load="" if "load[N]" not in last else f"{last['load[N]']:.8g}",
                v01="" if "v01" not in last else f"{last['v01']:.8g}",
                err="" if row.get("v01_relative_error") is None else f"{row['v01_relative_error']:.6g}",
            )
        )

    failed_details: list[str] = []
    for row in rows:
        summary = row["summary"]
        last_step = summary.get("steps", [])[-1] if summary.get("steps") else {}
        if row.get("exit_status") == 0 and summary.get("end_of_calculation"):
            continue
        if not last_step:
            continue
        failed_details.append(
            "| {variant} | {exit_status} | {step} | {it} | {res:.6e} | {disp:.6e} | {energy:.6e} |".format(
                variant=row["variant"],
                exit_status="" if row.get("exit_status") is None else row["exit_status"],
                step=last_step.get("step", ""),
                it=last_step.get("last_iter", ""),
                res=last_step.get("residual") or 0.0,
                disp=last_step.get("displacement") or 0.0,
                energy=last_step.get("energy") or 0.0,
            )
        )

    if failed_details:
        lines.extend(
            [
                "",
                "## Failed Or Incomplete Last Rows",
                "",
                "| variant | exit | last step | last iter | residual | displacement | energy |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        lines.extend(failed_details)

    strict_rows = [row for row in rows if row["variant"].endswith("-strict")]
    passed_strict = [row for row in strict_rows if row.get("verdict") == "passed"]
    lines.extend(["", "## Checkpoint Verdict", ""])
    if selected:
        lines.append(f"Sign calibration selected `{selected}` weights based on the first-step direction score.")
    if passed_strict:
        best = min(passed_strict, key=lambda row: row.get("v01_relative_error") or 0.0)
        lines.append(
            f"CP3 passed for `{best['variant']}`: it reached the strict target time with no fallback acceptance or NaNs."
        )
    elif strict_rows:
        lines.append(
            "CP3 did not pass for the strict TS-N65 IDC run(s). Preserve the run directories and use the compact table to decide whether to try adaptive IDC or CP4 arc-length."
        )
    else:
        lines.append("Strict TS-N65 IDC has not been run yet.")

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `cp3_indirect_control.tsv`: machine-readable aggregate table.",
            "- `runs/<variant>/solver.out`: ignored full logs.",
            "- `runs/<variant>/results/LD.out`: ignored gauge output.",
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
    parser.add_argument("--variants", help="Comma-separated variant names, or benchmark/calibration/selected-strict.")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    deck = (repo / args.deck).resolve() if not args.deck.is_absolute() else args.deck
    exe = args.exe.resolve()
    root = args.root or repo / f"results/tsn65-indirect-control-{dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
    root = root.resolve()
    (root / "runs").mkdir(parents=True, exist_ok=True)

    metadata = {
        "created_at": now_iso(),
        "repo": str(repo),
        "branch": run_git(["branch", "--show-current"], repo),
        "commit": run_git(["rev-parse", "HEAD"], repo),
        "executable": str(exe),
        "deck": str(deck),
        "target_v01_final": f"{TARGET_V01_FINAL:.12g}",
    }
    (root / "metadata.txt").write_text("\n".join(f"{k}={v}" for k, v in metadata.items()) + "\n", encoding="utf-8")

    if args.summarize_existing:
        rows = collect_existing_rows(root)
        write_report(root, rows, metadata)
        print(root)
        return 0

    rows = collect_existing_rows(root)
    if args.dry_run:
        write_report(root, rows, metadata)
        print(root)
        return 0

    selection = selected_variants(args.variants)

    if "benchmark" in selection:
        rows.append(run_benchmark(root, repo, exe, args.jobs, min(args.timeout_seconds, 1200)))
        write_report(root, rows, metadata)

    if "calibration" in selection:
        for variant in ["CP3-idc-minus-calibration", "CP3-idc-plus-calibration"]:
            rows.append(run_tsn65_variant(root, deck, exe, args.jobs, min(args.timeout_seconds, 1200), variant))
            write_report(root, rows, metadata)

    explicit_variants = [name for name in selection if name in VARIANTS]
    for variant in explicit_variants:
        rows.append(run_tsn65_variant(root, deck, exe, args.jobs, args.timeout_seconds, variant))
        write_report(root, rows, metadata)

    if "selected-strict" in selection:
        selected = choose_sign(rows)
        if selected is None:
            raise SystemExit("cannot run selected-strict without a successful calibration run")
        variant = f"CP3-idc-{selected}-strict"
        rows.append(run_tsn65_variant(root, deck, exe, args.jobs, args.timeout_seconds, variant))
        write_report(root, rows, metadata)

    write_report(root, rows, metadata)
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
