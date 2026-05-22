#!/usr/bin/env python3
"""Run TS-N65 nonlinear solver experiment sweeps.

The runner creates one ignored run directory per variant, writes the full
baseline solver deck plus variant controls, runs OAS, summarizes solver.out with
summarize_oas_run.py, and writes compact TSV/Markdown reports.
"""

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

from summarize_oas_run import TSN65_TARGET_TOTAL, format_tsv, parse_log  # noqa: E402


BASELINE_SOLVER = """SteadyStateNonLinearSolver
time_step 1.25e-3
total_time 1e-2

solver_type DeflatedFGMRES
max_iterations 300
stiffness_matrix_iter_update 10
stiffness_matrix_step_update 1
limit_tolerance 0
tolerance 1e-3

dfgmres_tolerance 1e-1
dfgmres_true_tolerance 1e-1
dfgmres_max_iterations 500
dfgmres_restart 80
dfgmres_deflation_vectors 20
dfgmres_deflation_eps 1e-15
dfgmres_collect_newton_steps 1
dfgmres_preconditioner hypre
dfgmres_reorthogonalize_on_matrix_change 1
dfgmres_elastic_reorder 2

hypre_coarsen_type 8
hypre_interp_type 6
hypre_strong_threshold 0.25
hypre_nodal 4
hypre_relax_type 16
hypre_num_sweeps 3
hypre_p_max 4
hypre_boomer_max_iterations 1
hypre_cheby_order 3
hypre_cheby_fraction -1
hypre_elastic_reorder 2
hypre_threads 0
hypre_use_dof_functions 0
hypre_use_interp_vectors 0

max_time_step 1.25e-3
min_time_step 1.25e-3
stiff_matrix_type tangent
"""


VARIANTS: dict[str, str] = {
    "G1-baseline": "",
    "G1-backtracking-frozen-actual": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "G1-backtracking-actual": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "G1-bisection-frozen-actual": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search bisection
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "G1-fixed-damping-05": """
nonlinear_damping_type fixed
nonlinear_damping_factor 0.5
""",
    "G1-adaptive-damping": """
nonlinear_damping_type adaptive
nonlinear_damping_factor 0.5
""",
    "G1-trust-stepnorm": """
nonlinear_trust_region step_norm
nonlinear_trust_radius_initial 1e-3
nonlinear_trust_radius_min 1e-8
nonlinear_trust_radius_max 1
nonlinear_trust_shrink 0.5
nonlinear_trust_expand 2
""",
    "G1-line-search-adaptive-K": """
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
nonlinear_adaptive_matrix_update 1
nonlinear_rebuild_on_small_alpha 0.5
nonlinear_rebuild_on_merit_growth 1
""",
    "G1-firstit-elastic": """
first_iteration_stiff_matrix_type elastic
""",
    "G1-firstit-secant": """
first_iteration_stiff_matrix_type secant
""",
    "G1-firstit-archived-csl": """
first_iteration_stiff_matrix_type archived_csl_damage_tangent
""",
    "P1-laststep-a025": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 0.25
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
""",
    "P1-laststep-a050": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 0.5
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
""",
    "P1-laststep-a075": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 0.75
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
""",
    "P1-laststep-a100": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 1.0
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
""",
    "P2-laststep-a050-plus-backtracking": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 0.5
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "P2-laststep-a050-plus-adaptive-K": """
nonlinear_initial_guess last_step
nonlinear_initial_guess_alpha 0.5
nonlinear_initial_guess_guard 1
nonlinear_initial_guess_accept_ratio 1.0
nonlinear_initial_guess_frozen_eval 1
nonlinear_initial_guess_max_norm_ratio 1.0
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation frozen_then_actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
nonlinear_adaptive_matrix_update 1
nonlinear_rebuild_on_small_alpha 0.5
nonlinear_rebuild_on_merit_growth 1
""",
    "D0-baseline-state-dump": """
nonlinear_state_dump 1
nonlinear_state_dump_steps 6 7
nonlinear_state_dump_top_damage 1000
nonlinear_state_dump_include_coordinates 1
""",
    "D1-backtracking-state-dump": """
nonlinear_state_dump 1
nonlinear_state_dump_steps 6 7
nonlinear_state_dump_top_damage 1000
nonlinear_state_dump_include_coordinates 1
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "D2-archived-csl-backtracking-state-dump": """
stiff_matrix_type archived_csl_damage_tangent
stiffness_matrix_iter_update -1
nonlinear_state_dump 1
nonlinear_state_dump_steps 6 7
nonlinear_state_dump_top_damage 1000
nonlinear_state_dump_include_coordinates 1
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual
nonlinear_line_search_merit mixed
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_cutback_on_fail 0
""",
    "S1-cslbeta005-gamma005-stepstart": """
stiff_matrix_type csl_stabilized_tangent
stiffness_matrix_iter_update -1
csl_tangent_beta 0.05
csl_tangent_softening_limit 0.05
csl_tangent_active_only 1
csl_tangent_log_stats 1
""",
    "S2-cslbeta010-gamma005-stepstart": """
stiff_matrix_type csl_stabilized_tangent
stiffness_matrix_iter_update -1
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.05
csl_tangent_active_only 1
csl_tangent_log_stats 1
""",
    "S3-cslbeta020-gamma005-stepstart": """
stiff_matrix_type csl_stabilized_tangent
stiffness_matrix_iter_update -1
csl_tangent_beta 0.20
csl_tangent_softening_limit 0.05
csl_tangent_active_only 1
csl_tangent_log_stats 1
""",
    "S4-cslbeta010-gamma010-stepstart": """
stiff_matrix_type csl_stabilized_tangent
stiffness_matrix_iter_update -1
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
""",
    "S5-cslbeta020-gamma010-stepstart": """
stiff_matrix_type csl_stabilized_tangent
stiffness_matrix_iter_update -1
csl_tangent_beta 0.20
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
""",
    "LM1-legacy-mu1e-4": """
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept merit
nonlinear_lm_max_trials 4
""",
    "LM2-legacy-mu1e-3": """
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-3
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept merit
nonlinear_lm_max_trials 4
""",
    "LM3-cslbeta005-gamma005-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.05
csl_tangent_softening_limit 0.05
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept merit
nonlinear_lm_max_trials 4
""",
    "LM4-cslbeta010-gamma005-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.05
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept merit
nonlinear_lm_max_trials 4
""",
    "LM5-cslbeta010-gamma010-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept merit
nonlinear_lm_max_trials 4
""",
    "LM6-legacy-relaxed102-mu1e-4": """
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.02
nonlinear_lm_max_trials 4
""",
    "LM7-cslbeta010-gamma010-relaxed102-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.02
nonlinear_lm_max_trials 4
""",
    "LM8-legacy-relaxed105-mu1e-4": """
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.05
nonlinear_lm_max_trials 4
""",
    "LM9-legacy-relaxed110-mu1e-4": """
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.10
nonlinear_lm_max_trials 4
""",
    "LM10-cslbeta010-gamma010-relaxed105-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.05
nonlinear_lm_max_trials 4
""",
    "LM11-cslbeta010-gamma010-relaxed110-mu1e-4": """
stiff_matrix_type csl_stabilized_tangent
csl_tangent_beta 0.10
csl_tangent_softening_limit 0.10
csl_tangent_active_only 1
csl_tangent_log_stats 1
nonlinear_lm_regularization 1
nonlinear_lm_mu_initial 1e-4
nonlinear_lm_mu_growth 10
nonlinear_lm_mu_shrink 0.25
nonlinear_lm_diag abs_diag
nonlinear_lm_accept relaxed_merit
nonlinear_lm_accept_growth 1.10
nonlinear_lm_max_trials 4
""",
}


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def run_git(args: list[str], cwd: Path) -> str:
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def selected_variants(selection: str | None) -> list[str]:
    if not selection:
        return list(VARIANTS)
    names = [name.strip() for name in selection.split(",") if name.strip()]
    unknown = [name for name in names if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")
    return names


def prepare_run_dir(deck: Path, run_dir: Path, solver_text: str, copy_deck: bool) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    for item in deck.iterdir():
        target = run_dir / item.name
        if item.name == "solver.inp":
            continue
        if target.exists() or target.is_symlink():
            continue
        if copy_deck:
            if item.is_dir():
                shutil.copytree(item, target)
            else:
                shutil.copy2(item, target)
        else:
            target.symlink_to(item.resolve())
    (run_dir / "solver.inp").write_text(solver_text, encoding="utf-8")


def verdict(summary: dict[str, Any], exit_status: int | None) -> str:
    if exit_status is None:
        return "not_run"
    if exit_status != 0:
        return "failed"
    if not summary["end_of_calculation"] or summary["nan_count"] or summary["fallback_acceptance_count"]:
        return "failed"
    if summary["step_count"] != 8:
        return "failed"
    if summary["total_iters"] < TSN65_TARGET_TOTAL:
        return "improves"
    if summary["total_iters"] == TSN65_TARGET_TOTAL:
        return "neutral"
    return "worsens"


def tail_rows(summary: dict[str, Any]) -> int:
    return sum(step["iters"] for step in summary["steps"][5:8])


def write_aggregate(root: Path, rows: list[dict[str, Any]], metadata: dict[str, str]) -> None:
    tsv = root / "cp1_sweep.tsv"
    with tsv.open("w", encoding="utf-8") as handle:
        handle.write(
            "variant\texit_status\tverdict\tsteps\ttotal_iters\ttail_6_8_iters\tduration\twarnings\tnans\tcutbacks\tfallback_accepts\tmin_alpha\tmax_ls_trials\tinitial_guess_accepted\tinitial_guess_rejected\tinitial_guess_min_ratio\n"
        )
        for row in rows:
            summary = row.get("summary") or {}
            steps = summary.get("steps") or []
            min_alpha_values = [step["min_alpha"] for step in steps if step.get("min_alpha") is not None]
            max_trials = max([step.get("max_ls_trials", 0) for step in steps] or [0])
            initial_guess = summary.get("initial_guess") or {}
            handle.write(
                "\t".join(
                    [
                        row["variant"],
                        "" if row.get("exit_status") is None else str(row["exit_status"]),
                        row.get("verdict", "not_run"),
                        str(summary.get("step_count", 0)),
                        str(summary.get("total_iters", 0)),
                        str(tail_rows(summary)) if steps else "",
                        str(summary.get("total_duration") or ""),
                        str(summary.get("warning_count", 0)),
                        str(summary.get("nan_count", 0)),
                        str(summary.get("cutback_count", 0)),
                        str(summary.get("fallback_acceptance_count", 0)),
                        "" if not min_alpha_values else f"{min(min_alpha_values):.12g}",
                        str(max_trials),
                        str(initial_guess.get("accepted", 0)),
                        str(initial_guess.get("rejected", 0)),
                        "" if initial_guess.get("min_predictor_merit_ratio") is None else f"{initial_guess['min_predictor_merit_ratio']:.12g}",
                    ]
                )
                + "\n"
            )

    report = root / "report.md"
    lines = [
        "# TS-N65 Nonlinear Solver Sweep",
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
            "## Result Table",
            "",
            "| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials | IG accept | IG reject | IG min ratio |",
            "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        summary = row.get("summary") or {}
        steps = summary.get("steps") or []
        min_alpha_values = [step["min_alpha"] for step in steps if step.get("min_alpha") is not None]
        max_trials = max([step.get("max_ls_trials", 0) for step in steps] or [0])
        initial_guess = summary.get("initial_guess") or {}
        lines.append(
            "| {variant} | {exit_status} | {verdict} | {steps} | {rows} | {tail} | {duration} | {warnings} | {nans} | {cutbacks} | {fallback} | {alpha} | {trials} | {ig_accept} | {ig_reject} | {ig_ratio} |".format(
                variant=row["variant"],
                exit_status="" if row.get("exit_status") is None else row["exit_status"],
                verdict=row.get("verdict", "not_run"),
                steps=summary.get("step_count", 0),
                rows=summary.get("total_iters", 0),
                tail=tail_rows(summary) if steps else "",
                duration=summary.get("total_duration") or "",
                warnings=summary.get("warning_count", 0),
                nans=summary.get("nan_count", 0),
                cutbacks=summary.get("cutback_count", 0),
                fallback=summary.get("fallback_acceptance_count", 0),
                alpha="" if not min_alpha_values else f"{min(min_alpha_values):.6g}",
                trials=max_trials,
                ig_accept=initial_guess.get("accepted", 0),
                ig_reject=initial_guess.get("rejected", 0),
                ig_ratio="" if initial_guess.get("min_predictor_merit_ratio") is None else f"{initial_guess['min_predictor_merit_ratio']:.6g}",
            )
        )
    lines.extend(
        [
            "",
            "## Baseline Gate",
            "",
            "- Target sequence: `6,6,10,13,17,183,187,163`.",
            "- Target total rows: `585`.",
            "- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.",
            "",
            "## Files",
            "",
            "- `cp1_sweep.tsv`: machine-readable aggregate table.",
            "- `runs/<variant>/solver.out`: ignored full logs.",
            "- `runs/<variant>/summary.json`: ignored parsed full summary.",
            "- `runs/<variant>/step_summary.tsv`: ignored parsed step table.",
        ]
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--deck", type=Path, default=Path("results/tsn65-zip-runs-20260519-190000/input/TS-N_65"))
    parser.add_argument("--exe", type=Path, default=Path("/tmp/oas_tsn65_full_baseline_build/bin/OAS"))
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--variants", help="Comma-separated variant names. Default: all CP1 variants.")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    parser.add_argument("--copy-deck", action="store_true", help="Copy deck files instead of symlinking them.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    deck = (repo / args.deck).resolve() if not args.deck.is_absolute() else args.deck
    exe = args.exe.resolve()
    root = args.root or repo / f"results/tsn65-globalization-sweep-{dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
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

    rows: list[dict[str, Any]] = []
    for variant in variants:
        run_dir = run_root / variant
        solver_text = BASELINE_SOLVER.rstrip() + "\n"
        extra = VARIANTS[variant].strip()
        if extra:
            solver_text += "\n" + extra + "\n"
        prepare_run_dir(deck, run_dir, solver_text, args.copy_deck)

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
            rows.append(row)
            write_aggregate(root, rows, metadata)
            continue

        start = time.monotonic()
        start_iso = now_iso()
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
                timed_out = False
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
        rows.append(row)
        write_aggregate(root, rows, metadata)

    write_aggregate(root, rows, metadata)
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
