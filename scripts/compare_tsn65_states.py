#!/usr/bin/env python3
"""Compare compact TS-N65 nonlinear material-state dump files."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


SUMMARY_KEYS = [
    "time",
    "dt",
    "nonlinear_iteration_count",
    "residual_error",
    "displacement_error",
    "energy_error",
    "global_displacement_norm",
    "step_increment_norm",
    "number_of_CSL_statuses",
    "number_of_CSL_statuses_with_damage_growth",
    "number_of_CSL_statuses_with_damage_gt_0",
    "number_of_CSL_statuses_with_temp_damage_gt_damage",
    "damage_max",
    "damage_mean",
    "damage_p90",
    "damage_p95",
    "damage_p99",
    "damage_p999",
    "damage_increment_max",
    "damage_increment_mean",
    "damage_increment_p90",
    "damage_increment_p95",
    "damage_increment_p99",
    "damage_increment_p999",
]


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def rel_diff(a: Any, b: Any) -> str:
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        return ""
    denom = max(abs(float(a)), abs(float(b)), 1e-300)
    value = abs(float(a) - float(b)) / denom
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def top_file(summary_path: Path) -> Path | None:
    candidate = Path(str(summary_path).replace("_summary.json", "_top_damage.csv"))
    return candidate if candidate.exists() else None


def load_top_ids(path: Path | None) -> set[tuple[str, str]]:
    if path is None:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {(row["element_id"], row["ip"]) for row in reader}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path, help="Baseline step_NNN_accepted_summary.json")
    parser.add_argument("candidate", type=Path, help="Candidate step_NNN_accepted_summary.json")
    parser.add_argument("--top", type=int, default=1000, help="Top-damage overlap denominator")
    args = parser.parse_args()

    baseline = load_json(args.baseline)
    candidate = load_json(args.candidate)

    print("# TS-N65 State Comparison")
    print()
    print(f"baseline: `{args.baseline}`")
    print(f"candidate: `{args.candidate}`")
    print()

    hash_equal = baseline.get("material_status_hash_global") == candidate.get("material_status_hash_global")
    print(f"global status hash equal: `{str(hash_equal).lower()}`")
    print()
    print("| metric | baseline | candidate | relative difference |")
    print("| --- | ---: | ---: | ---: |")
    for key in SUMMARY_KEYS:
        print(f"| `{key}` | {baseline.get(key, '')} | {candidate.get(key, '')} | {rel_diff(baseline.get(key), candidate.get(key))} |")

    baseline_top = load_top_ids(top_file(args.baseline))
    candidate_top = load_top_ids(top_file(args.candidate))
    if baseline_top or candidate_top:
        overlap = len(baseline_top & candidate_top)
        denom = max(min(len(baseline_top), len(candidate_top), args.top), 1)
        print()
        print(f"top-damage set overlap: `{overlap}/{denom}` (`{overlap / denom:.6f}`)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
