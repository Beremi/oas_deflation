#!/usr/bin/env python3
"""Summarize OAS nonlinear solver logs.

The parser is intentionally log-format oriented and dependency free so it can be
used from benchmark run folders without importing OAS Python modules.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


TSN65_TARGET_ITERS = [6, 6, 10, 13, 17, 183, 187, 163]
TSN65_TARGET_TOTAL = sum(TSN65_TARGET_ITERS)


STEP_RE = re.compile(
    r"^#+ Solving step\s+(?P<step>\d+)\s+at time\s+(?P<time>[^;]+);\s+time step\s+(?P<dt>\S+)"
)
ITER_RE = re.compile(
    r"^\s*(?P<iter>\d+)\s+"
    r"(?P<res>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"(?P<disp>---|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"(?P<energy>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?P<extra>.*)$"
)
DURATION_RE = re.compile(r"^step duration:\s*(?P<duration>\S+)")
TOTAL_DURATION_RE = re.compile(r"^#+ total duration:\s*(?P<duration>\S+)")
SETUP_RE = re.compile(r"DeflatedFGMRES setup complete: (?P<body>.*)$")


def as_float(value: str | None) -> float | None:
    if value is None or value == "---":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_key_values(line: str) -> dict[str, str]:
    tokens = line.strip().split()
    out: dict[str, str] = {}
    for i in range(0, len(tokens) - 1, 2):
        out[tokens[i]] = tokens[i + 1]
    return out


def parse_setup(line: str) -> dict[str, str]:
    match = SETUP_RE.search(line)
    if not match:
        return {}
    out: dict[str, str] = {}
    for part in match.group("body").split(","):
        if "=" in part:
            key, value = part.strip().split("=", 1)
            out[key.strip()] = value.strip()
    return out


def new_step(step: int, time: float | None, dt: float | None) -> dict[str, Any]:
    return {
        "step": step,
        "time": time,
        "dt": dt,
        "iters": 0,
        "last_iter": None,
        "residual": None,
        "displacement": None,
        "energy": None,
        "duration": None,
        "min_alpha": None,
        "max_ls_trials": 0,
        "matrix_rebuild_rows": 0,
        "cutback_reasons": [],
        "converged_by_tolerance": False,
    }


def update_min(current: float | None, candidate: float | None) -> float | None:
    if candidate is None or not math.isfinite(candidate):
        return current
    if current is None:
        return candidate
    return min(current, candidate)


def parse_log(path: Path, baseline: str | None = None) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    setup_lines: list[dict[str, str]] = []
    warnings: list[str] = []
    nan_lines: list[str] = []
    line_search = {"trial_lines": 0, "actual_lines": 0, "fail_lines": 0, "min_alpha": None}
    cutbacks: list[str] = []
    matrix_rebuild_reason_count: dict[str, int] = {}
    total_duration: str | None = None
    end_of_calculation = False
    fallback_acceptance_count = 0
    factorization_count = 0
    restart_count = 0
    shortening_count = 0
    enlargement_count = 0

    with path.open("r", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            lower = stripped.lower()

            step_match = STEP_RE.match(stripped)
            if step_match:
                current = new_step(
                    int(step_match.group("step")),
                    as_float(step_match.group("time")),
                    as_float(step_match.group("dt")),
                )
                steps.append(current)
                continue

            if "end of calculation" in lower:
                end_of_calculation = True

            total_match = TOTAL_DURATION_RE.match(stripped)
            if total_match:
                total_duration = total_match.group("duration")
                continue

            duration_match = DURATION_RE.match(stripped)
            if duration_match and current is not None:
                current["duration"] = duration_match.group("duration")
                continue

            if "factorizing system matrix" in lower:
                factorization_count += 1

            if "enlarging step" in lower:
                enlargement_count += 1
            if "shortening step" in lower:
                shortening_count += 1
            if "restarting step" in lower:
                restart_count += 1

            setup = parse_setup(stripped)
            if setup:
                setup_lines.append(setup)

            if lower.startswith("ls_trial"):
                line_search["trial_lines"] += 1
                fields = parse_key_values(stripped)
                line_search["min_alpha"] = update_min(line_search["min_alpha"], as_float(fields.get("alpha")))
            elif lower.startswith("ls_actual"):
                line_search["actual_lines"] += 1
                fields = parse_key_values(stripped)
                line_search["min_alpha"] = update_min(line_search["min_alpha"], as_float(fields.get("alpha")))
            elif lower.startswith("ls_fail"):
                line_search["fail_lines"] += 1

            if lower.startswith("nonlinear cutback reason:"):
                reason = stripped.split(":", 1)[1].strip()
                cutbacks.append(reason)
                if current is not None:
                    current["cutback_reasons"].append(reason)
                matrix_rebuild_reason_count[reason] = matrix_rebuild_reason_count.get(reason, 0) + 1

            if "fallback acceptance" in lower or "tolerance increased in this step" in lower:
                fallback_acceptance_count += 1

            if (
                "did not converge" in lower
                or "warning: performed" in lower
                or "error:" in lower
                or "failed" in lower
            ):
                warnings.append(f"{line_number}: {stripped}")

            if re.search(r"\bnan\b", lower):
                nan_lines.append(f"{line_number}: {stripped}")

            iter_match = ITER_RE.match(line)
            if iter_match and current is not None:
                current["iters"] += 1
                current["last_iter"] = int(iter_match.group("iter"))
                current["residual"] = as_float(iter_match.group("res"))
                current["displacement"] = as_float(iter_match.group("disp"))
                current["energy"] = as_float(iter_match.group("energy"))

                extra = iter_match.group("extra").split()
                if len(extra) >= 5:
                    alpha = as_float(extra[0])
                    ls_trials = as_float(extra[1])
                    cutback_reason = extra[4]
                    current["min_alpha"] = update_min(current["min_alpha"], alpha)
                    if ls_trials is not None:
                        current["max_ls_trials"] = max(current["max_ls_trials"], int(ls_trials))
                    if cutback_reason != "-":
                        current["cutback_reasons"].append(cutback_reason)
                        matrix_rebuild_reason_count[cutback_reason] = matrix_rebuild_reason_count.get(cutback_reason, 0) + 1
                if len(extra) >= 6 and extra[5] == "1":
                    current["matrix_rebuild_rows"] += 1

    for step in steps:
        res_ok = step["residual"] is not None and step["residual"] <= 1e-3
        disp_ok = step["displacement"] is None or step["displacement"] <= 1e-3
        ene_ok = step["energy"] is not None and step["energy"] <= 1e-3
        step["converged_by_tolerance"] = bool(res_ok and disp_ok and ene_ok)

    target_iters: list[int] = []
    if baseline == "tsn65":
        target_iters = TSN65_TARGET_ITERS
    for idx, step in enumerate(steps):
        step["target_iters"] = target_iters[idx] if idx < len(target_iters) else None
        step["matches_target"] = step["target_iters"] == step["iters"] if step["target_iters"] is not None else None

    total_iters = sum(step["iters"] for step in steps)
    result = {
        "log": str(path),
        "baseline": baseline,
        "steps": steps,
        "step_count": len(steps),
        "total_iters": total_iters,
        "target_total_iters": TSN65_TARGET_TOTAL if baseline == "tsn65" else None,
        "matches_target_total": total_iters == TSN65_TARGET_TOTAL if baseline == "tsn65" else None,
        "matches_target_sequence": [step["iters"] for step in steps] == target_iters if baseline == "tsn65" else None,
        "total_duration": total_duration,
        "end_of_calculation": end_of_calculation,
        "fallback_acceptance_count": fallback_acceptance_count,
        "warning_count": len(warnings),
        "warnings": warnings,
        "nan_count": len(nan_lines),
        "nan_lines": nan_lines,
        "line_search": line_search,
        "cutback_count": len(cutbacks),
        "cutbacks": cutbacks,
        "matrix_rebuild_reason_count": matrix_rebuild_reason_count,
        "factorization_count": factorization_count,
        "restart_count": restart_count,
        "shortening_count": shortening_count,
        "enlargement_count": enlargement_count,
        "first_setup": setup_lines[0] if setup_lines else {},
        "last_setup": setup_lines[-1] if setup_lines else {},
    }
    return result


def format_tsv(summary: dict[str, Any]) -> str:
    lines = [
        "\t".join(
            [
                "step",
                "time",
                "dt",
                "iters",
                "target_iters",
                "matches_target",
                "duration",
                "converged_by_tolerance",
                "residual",
                "displacement",
                "energy",
                "min_alpha",
                "max_ls_trials",
                "matrix_rebuild_rows",
                "cutback_reasons",
            ]
        )
    ]
    for step in summary["steps"]:
        lines.append(
            "\t".join(
                [
                    str(step["step"]),
                    "" if step["time"] is None else f"{step['time']:.12g}",
                    "" if step["dt"] is None else f"{step['dt']:.12g}",
                    str(step["iters"]),
                    "" if step["target_iters"] is None else str(step["target_iters"]),
                    "" if step["matches_target"] is None else str(int(step["matches_target"])),
                    step["duration"] or "",
                    str(int(step["converged_by_tolerance"])),
                    "" if step["residual"] is None else f"{step['residual']:.12g}",
                    "" if step["displacement"] is None else f"{step['displacement']:.12g}",
                    "" if step["energy"] is None else f"{step['energy']:.12g}",
                    "" if step["min_alpha"] is None else f"{step['min_alpha']:.12g}",
                    str(step["max_ls_trials"]),
                    str(step["matrix_rebuild_rows"]),
                    "|".join(step["cutback_reasons"]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def format_markdown(summary: dict[str, Any]) -> str:
    status = "pass" if summary.get("matches_target_sequence") else "check"
    lines = [
        f"# OAS Run Summary: {Path(summary['log']).name}",
        "",
        f"- Status: `{status}`",
        f"- Steps: `{summary['step_count']}`",
        f"- Total nonlinear rows: `{summary['total_iters']}`",
        f"- Target rows: `{summary['target_total_iters']}`",
        f"- Matches target sequence: `{summary['matches_target_sequence']}`",
        f"- OAS total duration: `{summary['total_duration']}`",
        f"- End of calculation: `{summary['end_of_calculation']}`",
        f"- Warnings: `{summary['warning_count']}`",
        f"- NaN lines: `{summary['nan_count']}`",
        f"- Cutbacks: `{summary['cutback_count']}`",
        f"- Fallback accepts: `{summary['fallback_acceptance_count']}`",
        "",
        "| step | time | dt | rows | target | match | duration | conv | residual | displacement | energy | alpha | ls trials | K rows | cutbacks |",
        "| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for step in summary["steps"]:
        lines.append(
            "| {step} | {time} | {dt} | {iters} | {target} | {match} | {duration} | {conv} | {residual} | {disp} | {energy} | {alpha} | {trials} | {krows} | {cutbacks} |".format(
                step=step["step"],
                time="" if step["time"] is None else f"{step['time']:.6g}",
                dt="" if step["dt"] is None else f"{step['dt']:.6g}",
                iters=step["iters"],
                target="" if step["target_iters"] is None else step["target_iters"],
                match="" if step["matches_target"] is None else ("yes" if step["matches_target"] else "no"),
                duration=step["duration"] or "",
                conv=int(step["converged_by_tolerance"]),
                residual="" if step["residual"] is None else f"{step['residual']:.6g}",
                disp="" if step["displacement"] is None else f"{step['displacement']:.6g}",
                energy="" if step["energy"] is None else f"{step['energy']:.6g}",
                alpha="" if step["min_alpha"] is None else f"{step['min_alpha']:.6g}",
                trials=step["max_ls_trials"],
                krows=step["matrix_rebuild_rows"],
                cutbacks="<br>".join(step["cutback_reasons"]),
            )
        )
    if summary["warnings"]:
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- `{warning}`" for warning in summary["warnings"][:20])
    if summary["nan_lines"]:
        lines.extend(["", "## NaNs", ""])
        lines.extend(f"- `{line}`" for line in summary["nan_lines"][:20])
    return "\n".join(lines) + "\n"


def format_text(summary: dict[str, Any]) -> str:
    sequence = ",".join(str(step["iters"]) for step in summary["steps"])
    lines = [
        f"log: {summary['log']}",
        f"steps: {summary['step_count']}",
        f"iteration_sequence: {sequence}",
        f"total_iters: {summary['total_iters']}",
        f"target_total_iters: {summary['target_total_iters']}",
        f"matches_target_sequence: {summary['matches_target_sequence']}",
        f"total_duration: {summary['total_duration']}",
        f"end_of_calculation: {summary['end_of_calculation']}",
        f"warnings: {summary['warning_count']}",
        f"nan_lines: {summary['nan_count']}",
        f"cutbacks: {summary['cutback_count']}",
        f"fallback_acceptance_count: {summary['fallback_acceptance_count']}",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="Path to solver.out or another OAS log")
    parser.add_argument("--baseline", choices=["tsn65"], help="Compare against a known baseline")
    parser.add_argument("--format", choices=["text", "json", "tsv", "markdown"], default="text")
    args = parser.parse_args(argv)

    summary = parse_log(args.log, args.baseline)
    if args.format == "json":
        print(json.dumps(summary, indent=2, sort_keys=True))
    elif args.format == "tsv":
        print(format_tsv(summary), end="")
    elif args.format == "markdown":
        print(format_markdown(summary), end="")
    else:
        print(format_text(summary), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
