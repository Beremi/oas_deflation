#!/usr/bin/env python3
"""Compare OAS runtime profile artifacts from two saved runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


DEFAULT_PHASES = [
    "forces.total:active",
    "forces.integrate_internal_forces:active",
    "forces.simplex_volumetric_strains:active",
    "forces.external_reactions:active",
    "matrix.stiffness_update:secant",
    "matrix.factorize_linear_system_total:-",
    "matrix.constraint_transform:-",
    "step_finalize.material_status_update:-",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: Any, default: float = math.nan) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except (TypeError, ValueError):
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
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(format_cell(item) for item in row) + " |")
    return "\n".join(lines)


def phase_key(spec: str) -> tuple[str, str]:
    if ":" not in spec:
        return spec, ""
    phase, detail = spec.split(":", 1)
    return phase, detail


def phase_map(run_dir: Path) -> dict[tuple[str, str], dict[str, Any]]:
    out: dict[tuple[str, str], dict[str, Any]] = {}
    for row in read_tsv(run_dir / "runtime_profile_summary.tsv"):
        key = (row.get("phase", ""), row.get("detail", ""))
        out[key] = {
            "count": as_int(row.get("count")),
            "total_seconds": as_float(row.get("total_seconds")),
            "mean_seconds": as_float(row.get("mean_seconds")),
            "max_seconds": as_float(row.get("max_seconds")),
        }
    return out


def step_elapsed(run_dir: Path) -> dict[int, float]:
    elapsed: dict[int, float] = {}
    for row in read_tsv(run_dir / "solver.out"):
        step = as_int(row.get("#step") or row.get("step"))
        value = as_float(row.get("elapsed_time"))
        if step > 0 and math.isfinite(value):
            elapsed[step] = value
    return elapsed


def linear_step_map(run_dir: Path) -> dict[int, dict[str, float]]:
    out: dict[int, dict[str, float]] = {}
    for row in read_tsv(run_dir / "linear_profile_events.tsv"):
        if row.get("phase") != "solve":
            continue
        step = as_int(row.get("step"))
        slot = out.setdefault(step, {"solves": 0, "outer_iterations": 0, "solve_seconds": 0.0})
        slot["solves"] += 1
        slot["outer_iterations"] += as_int(row.get("solver_iterations"))
        slot["solve_seconds"] += as_float(row.get("duration_seconds"), 0.0)
    return out


def ratio(old: float, new: float) -> float:
    if not math.isfinite(old) or not math.isfinite(new) or new == 0:
        return math.nan
    return old / new


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--baseline-label", default="baseline")
    parser.add_argument("--candidate-label", default="candidate")
    parser.add_argument("--phase", action="append", dest="phases", default=[])
    args = parser.parse_args()

    phases = args.phases or DEFAULT_PHASES
    baseline_json = read_json(args.baseline / "run.json")
    candidate_json = read_json(args.candidate / "run.json")

    metric_rows = []
    for key, label in [
        ("status", "status"),
        ("steps", "steps"),
        ("wall_seconds", "wall s"),
        ("nonlinear_iterations", "nonlinear iters"),
        ("linear_solves", "linear solves"),
        ("outer_iterations", "DFGMRES outer iters"),
        ("max_iterations", "max linear iters"),
        ("final_true_residual", "final true relres"),
        ("final_nonlinear_residual", "final nonlinear residual"),
    ]:
        old = baseline_json.get(key)
        new = candidate_json.get(key)
        speed = ratio(as_float(old), as_float(new)) if key not in {"status", "steps"} else None
        metric_rows.append([label, old, new, speed])

    print(md_table(["metric", args.baseline_label, args.candidate_label, "baseline/candidate"], metric_rows))

    baseline_phases = phase_map(args.baseline)
    candidate_phases = phase_map(args.candidate)
    phase_rows = []
    for spec in phases:
        key = phase_key(spec)
        old = baseline_phases.get(key, {})
        new = candidate_phases.get(key, {})
        old_total = as_float(old.get("total_seconds"))
        new_total = as_float(new.get("total_seconds"))
        phase_rows.append([
            spec,
            old.get("count", ""),
            old_total,
            old.get("mean_seconds", ""),
            new.get("count", ""),
            new_total,
            new.get("mean_seconds", ""),
            ratio(old_total, new_total),
        ])

    print()
    print(md_table([
        "phase",
        f"{args.baseline_label} count",
        f"{args.baseline_label} total s",
        f"{args.baseline_label} mean s",
        f"{args.candidate_label} count",
        f"{args.candidate_label} total s",
        f"{args.candidate_label} mean s",
        "speedup",
    ], phase_rows))

    base_steps = step_elapsed(args.baseline)
    cand_steps = step_elapsed(args.candidate)
    if base_steps or cand_steps:
        rows = []
        previous_base = 0.0
        previous_cand = 0.0
        for step in sorted(set(base_steps) | set(cand_steps)):
            base_elapsed = base_steps.get(step, math.nan)
            cand_elapsed = cand_steps.get(step, math.nan)
            base_step = base_elapsed - previous_base if math.isfinite(base_elapsed) else math.nan
            cand_step = cand_elapsed - previous_cand if math.isfinite(cand_elapsed) else math.nan
            previous_base = base_elapsed if math.isfinite(base_elapsed) else previous_base
            previous_cand = cand_elapsed if math.isfinite(cand_elapsed) else previous_cand
            rows.append([step, base_step, cand_step, ratio(base_step, cand_step)])
        print()
        print(md_table(["step", f"{args.baseline_label} step s", f"{args.candidate_label} step s", "speedup"], rows))

    base_linear = linear_step_map(args.baseline)
    cand_linear = linear_step_map(args.candidate)
    if base_linear or cand_linear:
        rows = []
        for step in sorted(set(base_linear) | set(cand_linear)):
            old = base_linear.get(step, {})
            new = cand_linear.get(step, {})
            rows.append([
                step,
                old.get("solves", ""),
                old.get("outer_iterations", ""),
                old.get("solve_seconds", ""),
                new.get("solves", ""),
                new.get("outer_iterations", ""),
                new.get("solve_seconds", ""),
            ])
        print()
        print(md_table([
            "step",
            f"{args.baseline_label} solves",
            f"{args.baseline_label} outer",
            f"{args.baseline_label} solve s",
            f"{args.candidate_label} solves",
            f"{args.candidate_label} outer",
            f"{args.candidate_label} solve s",
        ], rows))


if __name__ == "__main__":
    main()
