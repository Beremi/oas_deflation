#!/usr/bin/env python3
"""Run TS-N65 gauge-constrained arc-length experiments.

CP5 showed that a spherical load-factor arc-length continuation is not directly
usable for TS-N65 because the benchmark is driven by prescribed displacement
history.  This runner tests the next usable variant: keep the arc-length load
factor as the equilibrium load parameter, but constrain the corrector with the
same displacement gauge used by the CP3 indirect-control experiment.
"""

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

from run_tsn65_globalization_sweep import BASELINE_SOLVER, prepare_run_dir  # noqa: E402
from summarize_oas_run import TSN65_TARGET_ITERS, TSN65_TARGET_TOTAL, format_tsv, parse_log  # noqa: E402


TARGET_TOTAL_TIME = 1.0e-2
TARGET_FIRST_TIME = 1.25e-3
TARGET_V01_FIRST = 2.613423e-6
TARGET_V01_FINAL = 1.519697e-5
TARGET_CONTROL_SCALE_FIRST = 1.205765978081034
TARGET_CONTROL_SCALE_FINAL = 1.2321920851549628

IDC_TEMPLATE = """indirect_control 2
ic_xcoords 2.1 2.1
ic_ycoords 0 0
ic_zcoords 0 -0.4
ic_directions 2 2
ic_displ_weights {weights}
ic_force_weights 0 0
ic_function {function_index}
"""

ARC_RE = re.compile(
    r"^ARC_LENGTH\s+step\s+(?P<step>\d+)\s+it\s+(?P<it>\d+)\s+"
    r"lambda\s+(?P<lambda>\S+)\s+dlambda\s+(?P<dlambda>\S+)\s+"
    r"radius\s+(?P<radius>\S+)"
    r"(?:.*?\bgauge\s+(?P<gauge>\S+)\s+gauge_target\s+(?P<target>\S+)\s+"
    r"gauge_error\s+(?P<gauge_error>\S+)\s+control_time\s+(?P<control_time>\S+))?"
)

VARIANTS: dict[str, dict[str, Any]] = {
    "CP5b-gauge-first-quarter": {
        "total_time": TARGET_FIRST_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FIRST,
        "weights": "-1 1",
        "line_search": False,
        "purpose": "first-quarter smoke test against the strict baseline first-step gauge",
    },
    "CP5b-gauge-first-quarter-line-search": {
        "total_time": TARGET_FIRST_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FIRST,
        "weights": "-1 1",
        "line_search": True,
        "purpose": "first-quarter smoke test with existing actual-residual backtracking enabled",
    },
    "CP5b-gauge-first-quarter-calibrated": {
        "total_time": TARGET_FIRST_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FIRST,
        "control_scale": TARGET_CONTROL_SCALE_FIRST,
        "weights": "-1 1",
        "line_search": False,
        "purpose": "first-quarter target with the internal IDC gauge scaled to match exported LD v01",
    },
    "CP5b-gauge-first-quarter-calibrated-line-search": {
        "total_time": TARGET_FIRST_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FIRST,
        "control_scale": TARGET_CONTROL_SCALE_FIRST,
        "weights": "-1 1",
        "line_search": True,
        "purpose": "calibrated first-quarter target with existing actual-residual backtracking enabled",
    },
    "CP5b-gauge-full": {
        "total_time": TARGET_TOTAL_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FINAL,
        "weights": "-1 1",
        "line_search": False,
        "purpose": "full eight-quarter target with gauge-constrained arc-length",
    },
    "CP5b-gauge-full-line-search": {
        "total_time": TARGET_TOTAL_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FINAL,
        "weights": "-1 1",
        "line_search": True,
        "purpose": "full eight-quarter target with gauge arc-length plus existing backtracking",
    },
    "CP5b-gauge-full-calibrated": {
        "total_time": TARGET_TOTAL_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FINAL,
        "control_scale": TARGET_CONTROL_SCALE_FINAL,
        "weights": "-1 1",
        "line_search": False,
        "purpose": "full eight-quarter target with the CP3-calibrated internal IDC gauge scale",
    },
    "CP5b-gauge-full-calibrated-line-search": {
        "total_time": TARGET_TOTAL_TIME,
        "time_step": TARGET_FIRST_TIME,
        "target_v01": TARGET_V01_FINAL,
        "control_scale": TARGET_CONTROL_SCALE_FINAL,
        "weights": "-1 1",
        "line_search": True,
        "purpose": "full calibrated gauge arc-length plus existing actual-residual backtracking",
    },
}


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def run_git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def selected_variants(selection: str | None) -> list[str]:
    if not selection:
        return ["CP5b-gauge-first-quarter", "CP5b-gauge-first-quarter-line-search"]
    names = [name.strip() for name in selection.split(",") if name.strip()]
    unknown = [name for name in names if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")
    return names


def replace_solver_value(solver: str, key: str, value: str) -> str:
    lines: list[str] = []
    replaced = False
    for raw in solver.splitlines():
        raw_key = raw.split(None, 1)[0] if raw.strip() else ""
        if raw_key == key:
            lines.append(f"{key} {value}")
            replaced = True
        else:
            lines.append(raw)
    if not replaced:
        lines.append(f"{key} {value}")
    return "\n".join(lines).rstrip() + "\n"


def append_control_function(run_dir: Path, total_time: float, target_control: float) -> int:
    path = run_dir / "functions.inp"
    if path.is_symlink():
        source = path.resolve()
        path.unlink()
        text = source.read_text(encoding="utf-8", errors="replace")
    else:
        text = path.read_text(encoding="utf-8", errors="replace")
    if not text.endswith("\n"):
        text += "\n"
    non_comment_count = len([line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")])
    text += (
        f"PWLFunction 2 0 {total_time:.12g} 0 {target_control:.12g} "
        f"# CP5b gauge arc-length internal control target\n"
    )
    path.write_text(text, encoding="utf-8")
    return non_comment_count


def solver_for(spec: dict[str, Any], function_index: int) -> str:
    solver = BASELINE_SOLVER
    solver = replace_solver_value(solver, "time_step", f"{spec['time_step']:.12g}")
    solver = replace_solver_value(solver, "total_time", f"{spec['total_time']:.12g}")
    solver = replace_solver_value(solver, "min_time_step", f"{spec['time_step']:.12g}")
    solver = replace_solver_value(solver, "max_time_step", f"{spec['time_step']:.12g}")
    controls = [
        "",
        "nonlinear_control arc_length",
        "arc_length_constraint gauge",
        "arc_length_reference finite_difference",
        "arc_length_reference_delta 1.25e-3",
        "arc_length_sign_strategy monotone_lambda",
        "arc_length_auto_radius 1",
        "arc_length_gauge_tolerance 1e-3",
        "arc_length_radius_initial 1e-3",
        "arc_length_radius_min 1e-8",
        "arc_length_radius_max 1",
        "arc_length_psi 1",
        "arc_length_shrink 0.5",
        "arc_length_expand 1",
        "arc_length_target_iterations 8",
        "arc_length_max_iterations 300",
        "nonlinear_material_snapshot_rollback 1",
    ]
    if spec.get("line_search"):
        controls.extend(
            [
                "nonlinear_line_search backtracking",
                "nonlinear_line_search_evaluation actual",
                "nonlinear_line_search_merit mixed",
                "nonlinear_line_search_min_alpha 0.03125",
                "nonlinear_line_search_max_trials 6",
                "nonlinear_line_search_cutback_on_fail 1",
            ]
        )
    controls.append(IDC_TEMPLATE.format(weights=spec["weights"], function_index=function_index).rstrip())
    return solver.rstrip() + "\n" + "\n".join(controls).rstrip() + "\n"


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
        if not match:
            continue
        row: dict[str, float] = {
            "step": float(match.group("step")),
            "iteration": float(match.group("it")),
            "lambda": float(match.group("lambda")),
            "dlambda": float(match.group("dlambda")),
            "radius": float(match.group("radius")),
        }
        for key in ["gauge", "target", "gauge_error", "control_time"]:
            value = match.group(key)
            row[key] = float(value) if value is not None else float("nan")
        rows.append(row)
    return rows


def write_arc_tsv(path: Path, rows: list[dict[str, float]]) -> None:
    headers = ["step", "iteration", "lambda", "dlambda", "radius", "gauge", "target", "gauge_error", "control_time"]
    lines = ["\t".join(headers)]
    for row in rows:
        lines.append("\t".join(f"{row[h]:.12e}" for h in headers))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_ld(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"header": [], "rows": [], "first": {}, "last": {}}
    header: list[str] = []
    rows: list[dict[str, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("#"):
            header = raw.lstrip("#").split()
            continue
        if not raw.strip() or not header:
            continue
        parts = raw.split()
        if len(parts) < len(header):
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


def run_variant(root: Path, deck: Path, exe: Path, jobs: int, timeout: int, variant: str) -> dict[str, Any]:
    spec = VARIANTS[variant]
    run_dir = root / "runs" / variant
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    prepare_run_dir(deck, run_dir, "# temporary solver, rewritten after CP5b function append\n", copy_deck=False)
    target_control = float(spec["target_v01"]) * float(spec.get("control_scale", 1.0))
    function_index = append_control_function(run_dir, float(spec["total_time"]), target_control)
    (run_dir / "solver.inp").write_text(solver_for(spec, function_index), encoding="utf-8")
    env = make_env(jobs)
    command = [str(exe), "-j", str(jobs), "master.inp"]
    started = now_iso()
    exit_status, elapsed, timed_out = run_command(command, run_dir, env, timeout)
    summary = parse_log(run_dir / "solver.out", "tsn65")
    arc_rows = parse_arc_rows(run_dir / "solver.out")
    ld = read_ld(run_dir / "results" / "LD.out")
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (run_dir / "step_summary.tsv").write_text(format_tsv(summary), encoding="utf-8")
    write_arc_tsv(run_dir / "arc_length.tsv", arc_rows)
    (run_dir / "run_metadata.txt").write_text(
        f"started={started}\n"
        f"ended={now_iso()}\n"
        f"elapsed_seconds={elapsed:.3f}\n"
        f"exit_status={exit_status}\n"
        f"timed_out={int(timed_out)}\n"
        f"command={' '.join(command)}\n"
        f"function_index={function_index}\n"
        f"target_control={target_control:.12g}\n",
        encoding="utf-8",
    )
    return {
        "variant": variant,
        "spec": spec,
        "run_dir": str(run_dir),
        "exit_status": exit_status,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "summary": summary,
        "arc_rows": arc_rows,
        "ld": ld,
        "target_control": target_control,
    }


def verdict(result: dict[str, Any]) -> str:
    summary = result["summary"]
    spec = result["spec"]
    if result["exit_status"] != 0 or result["timed_out"]:
        return "failed"
    if not summary.get("end_of_calculation"):
        return "failed"
    if summary.get("nan_count", 0) != 0 or summary.get("fallback_acceptance_count", 0) != 0:
        return "failed"
    last = result["ld"]["last"]
    final_time = last.get("time")
    if final_time is None or final_time + 1e-12 < float(spec["total_time"]):
        return "partial"
    if "v01" in last:
        rel_err = abs(last["v01"] - float(spec["target_v01"])) / max(abs(float(spec["target_v01"])), 1e-30)
        if rel_err > 5e-2:
            return "wrong_gauge"
    rows = summary.get("total_iters", 0)
    target_rows = TSN65_TARGET_ITERS[0] if float(spec["total_time"]) <= TARGET_FIRST_TIME + 1e-15 else TSN65_TARGET_TOTAL
    if rows < target_rows:
        return "improves"
    if rows == target_rows:
        return "neutral"
    return "worsens"


def aggregate_rows(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        summary = result["summary"]
        arc_rows = result["arc_rows"]
        last_arc = arc_rows[-1] if arc_rows else {}
        last_ld = result["ld"]["last"]
        final_v01 = last_ld.get("v01")
        target_v01 = float(result["spec"]["target_v01"])
        target_control = float(result["target_control"])
        rel_err = "" if final_v01 is None else abs(final_v01 - target_v01) / max(abs(target_v01), 1e-30)
        rows.append(
            {
                "variant": result["variant"],
                "exit_status": result["exit_status"],
                "verdict": verdict(result),
                "steps": summary.get("step_count", 0),
                "rows": summary.get("total_iters", 0),
                "duration": summary.get("total_duration") or "",
                "elapsed_seconds": f"{result['elapsed_seconds']:.3f}",
                "warnings": summary.get("warning_count", 0),
                "nans": summary.get("nan_count", 0),
                "fallback": summary.get("fallback_acceptance_count", 0),
                "cutbacks": summary.get("cutback_count", 0),
                "arc_rows": len(arc_rows),
                "final_time": last_ld.get("time", ""),
                "final_lambda": last_arc.get("lambda", ""),
                "final_dlambda": last_arc.get("dlambda", ""),
                "final_radius": last_arc.get("radius", ""),
                "final_gauge": last_arc.get("gauge", ""),
                "target_gauge": target_control,
                "gauge_rel_err": rel_err,
                "final_v01": final_v01 if final_v01 is not None else "",
                "target_v01": target_v01,
                "control_scale": result["spec"].get("control_scale", 1.0),
                "purpose": result["spec"]["purpose"],
            }
        )
    return rows


def collect_existing_results(root: Path) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for run_dir in sorted((root / "runs").glob("CP5b-*")):
        if run_dir.name not in VARIANTS or not (run_dir / "solver.out").exists():
            continue
        exit_status = 0
        timed_out = False
        elapsed = 0.0
        target_control = float(VARIANTS[run_dir.name]["target_v01"]) * float(VARIANTS[run_dir.name].get("control_scale", 1.0))
        metadata = run_dir / "run_metadata.txt"
        if metadata.exists():
            for raw in metadata.read_text(encoding="utf-8", errors="replace").splitlines():
                key, _, value = raw.partition("=")
                if key == "exit_status":
                    try:
                        exit_status = int(value)
                    except ValueError:
                        pass
                elif key == "timed_out":
                    timed_out = value.strip() == "1"
                elif key == "elapsed_seconds":
                    try:
                        elapsed = float(value)
                    except ValueError:
                        pass
                elif key == "target_control":
                    try:
                        target_control = float(value)
                    except ValueError:
                        pass
        summary = parse_log(run_dir / "solver.out", "tsn65")
        arc_rows = parse_arc_rows(run_dir / "solver.out")
        ld = read_ld(run_dir / "results" / "LD.out")
        results.append(
            {
                "variant": run_dir.name,
                "spec": VARIANTS[run_dir.name],
                "run_dir": str(run_dir),
                "exit_status": exit_status,
                "timed_out": timed_out,
                "elapsed_seconds": elapsed,
                "summary": summary,
                "arc_rows": arc_rows,
                "ld": ld,
                "target_control": target_control,
            }
        )
    return results


def write_report(root: Path, repo: Path, deck: Path, exe: Path, variants: list[str], results: list[dict[str, Any]]) -> None:
    rows = aggregate_rows(results)
    (root / "metadata.txt").write_text(
        f"created_at={now_iso()}\n"
        f"repo={repo}\n"
        f"branch={run_git(repo, 'branch', '--show-current')}\n"
        f"commit={run_git(repo, 'rev-parse', 'HEAD')}\n"
        f"git_status_short<<EOF\n{run_git(repo, 'status', '--short')}\nEOF\n"
        f"executable={exe}\n"
        f"deck={deck}\n"
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
        "final_time",
        "final_lambda",
        "final_dlambda",
        "final_radius",
        "final_gauge",
        "target_gauge",
        "gauge_rel_err",
        "final_v01",
        "target_v01",
    ]
    with (root / "cp5b_gauge_arc_length.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[h]) for h in headers) + "\n")

    lines = [
        "# CP5b TS-N65 Gauge Arc-Length Experiments",
        "",
        f"Generated: {now_iso()}",
        "",
        "## Baseline Target",
        "",
        f"- Strict eight-quarter iteration sequence: `{','.join(str(i) for i in TSN65_TARGET_ITERS)}`",
        f"- Strict total nonlinear rows: `{TSN65_TARGET_TOTAL}`",
        f"- First-quarter gauge target: `{TARGET_V01_FIRST:.12g}`",
        f"- Full gauge target: `{TARGET_V01_FINAL:.12g}`",
        f"- First-quarter internal-control calibration scale: `{TARGET_CONTROL_SCALE_FIRST:.12g}`",
        f"- Full internal-control calibration scale from CP3: `{TARGET_CONTROL_SCALE_FINAL:.12g}`",
        "",
        "## Result Table",
        "",
        "| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | final time | lambda | radius | gauge | target | gauge rel err |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {variant} | {exit_status} | {verdict} | {steps} | {rows} | {duration} | {warnings} | {nans} | {fallback} | {cutbacks} | {final_time} | {final_lambda} | {final_radius} | {final_gauge} | {target_gauge} | {gauge_rel_err} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `arc_length_constraint gauge` uses the CP3 `v01` displacement gauge as the continuation constraint.",
            "- The runner appends a local target function to `functions.inp` and writes the actual zero-based `ic_function` index into each run's `solver.inp`.",
            "- Uncalibrated runs target the internal IDC gauge directly. Calibrated runs scale that internal target because the exported `LD.out` interpolation gauge differs from the internal IDC coordinate selection.",
            "- A result is promotable only if it reaches the target gauge with no fallback acceptance or NaNs and fewer rows or less wall time than the matching strict baseline target.",
            "",
            "## Files",
            "",
            "- `cp5b_gauge_arc_length.tsv`: compact aggregate table.",
            "- `runs/<variant>/solver.out`: ignored full solver log.",
            "- `runs/<variant>/arc_length.tsv`: ignored arc-length trace with gauge columns.",
            "- `runs/<variant>/results/LD.out`: ignored load/displacement gauge output.",
        ]
    )
    (root / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--deck", type=Path, default=Path("results/tsn65-zip-runs-20260519-190000/input/TS-N_65"))
    parser.add_argument("--exe", type=Path, default=Path("/tmp/oas_tsn65_full_baseline_build/bin/OAS"))
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--variants", default=None, help="Comma-separated variant names")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    deck = (repo / args.deck).resolve() if not args.deck.is_absolute() else args.deck.resolve()
    exe = args.exe.resolve()
    root = args.root or repo / "results" / f"tsn65-gauge-arc-cp5b-{now_stamp()}"
    root.mkdir(parents=True, exist_ok=True)
    variants = selected_variants(args.variants)
    results = collect_existing_results(root)
    if args.summarize_existing:
        write_report(root, repo, deck, exe, variants, results)
        print(root)
        return 0
    for variant in variants:
        print(f"[cp5b] running {variant}", flush=True)
        results = [result for result in results if result["variant"] != variant]
        results.append(run_variant(root, deck, exe, args.jobs, args.timeout, variant))
        write_report(root, repo, deck, exe, variants, results)
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
