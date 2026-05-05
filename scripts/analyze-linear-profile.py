#!/usr/bin/env python3
"""Generate a compact Markdown/PNG report from OAS linear solver profile TSVs."""

from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def read_solver_out(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
    columns = header.lstrip("#").split("\t")
    return pd.read_csv(path, sep="\t", comment="#", names=columns)


def read_solver_inp(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    params: dict[str, str] = {}
    solver_class = ""
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 1 and not solver_class:
            solver_class = parts[0]
            continue
        if len(parts) >= 2:
            params[parts[0]] = " ".join(parts[1:])
    if solver_class:
        params["oas_solver_class"] = solver_class
    return params


def git_revision() -> str:
    root = Path(__file__).resolve().parents[1]
    try:
        sha = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--short=12", "HEAD"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        ).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return ""


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    headers = [str(col) for col in df.columns]
    rows = [[format_cell(value) for value in row] for row in df.itertuples(index=False, name=None)]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def format_cell(value) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        if abs(value) >= 1000 or (abs(value) < 0.001 and value != 0):
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)


def save_stacked_step_time(per_step: pd.DataFrame, out_dir: Path) -> str | None:
    if per_step.empty:
        return None
    columns = [col for col in ["analyze", "factorize", "solve", "nonlinear_overhead"] if col in per_step]
    ax = per_step.set_index("step")[columns].plot(kind="bar", stacked=True, figsize=(13, 5), width=0.9)
    ax.set_xlabel("Step")
    ax.set_ylabel("Seconds")
    ax.set_title("Per-step wall time split")
    ax.legend(loc="upper right")
    every = max(1, len(per_step) // 20)
    for label in ax.get_xticklabels():
        visible = int(label.get_text()) % every == 0
        label.set_visible(visible)
    plt.tight_layout()
    name = "step-time-stacked.png"
    plt.savefig(out_dir / name, dpi=150)
    plt.close()
    return name


def save_iterations_solves(per_step: pd.DataFrame, out_dir: Path) -> str | None:
    if per_step.empty:
        return None
    fig, ax1 = plt.subplots(figsize=(12, 4.5))
    ax1.plot(per_step["step"], per_step["iterations"], color="#1f77b4", label="Nonlinear iterations")
    ax1.set_xlabel("Step")
    ax1.set_ylabel("Nonlinear iterations", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax2 = ax1.twinx()
    ax2.plot(per_step["step"], per_step["solve_count"], color="#d62728", label="Linear solves")
    ax2.set_ylabel("Linear solves", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    fig.suptitle("Nonlinear iterations and linear solves per step")
    fig.tight_layout()
    name = "iterations-solves.png"
    plt.savefig(out_dir / name, dpi=150)
    plt.close()
    return name


def save_matrix_delta(factorize: pd.DataFrame, out_dir: Path) -> str | None:
    data = factorize[factorize["matrix_relative_delta"] >= 0].copy()
    if data.empty:
        return None
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(range(1, len(data) + 1), data["matrix_relative_delta"], marker=".", linestyle="-")
    ax.set_xlabel("Factorization index with previous comparable matrix")
    ax.set_ylabel("Relative matrix value delta")
    ax.set_yscale("log")
    ax.set_title("Matrix change between factorizations")
    fig.tight_layout()
    name = "matrix-delta.png"
    plt.savefig(out_dir / name, dpi=150)
    plt.close()
    return name


def save_linear_share(per_step: pd.DataFrame, out_dir: Path) -> str | None:
    if per_step.empty:
        return None
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(per_step["step"], 100.0 * per_step["linear_share"], marker=".", linestyle="-")
    ax.set_xlabel("Step")
    ax.set_ylabel("Linear share of step wall time [%]")
    ax.set_title("Linear solver share per step")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    name = "linear-share.png"
    plt.savefig(out_dir / name, dpi=150)
    plt.close()
    return name


def save_solve_duration_distribution(solve: pd.DataFrame, out_dir: Path) -> str | None:
    if solve.empty:
        return None
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for rhs_kind, group in solve.groupby("rhs_kind"):
        ax.hist(group["duration_seconds"], bins=50, alpha=0.55, label=rhs_kind)
    ax.set_xlabel("Solve duration [s]")
    ax.set_ylabel("Count")
    ax.set_title("Linear solve duration distribution")
    ax.legend()
    fig.tight_layout()
    name = "solve-duration-distribution.png"
    plt.savefig(out_dir / name, dpi=150)
    plt.close()
    return name


def compute_per_step(events: pd.DataFrame, solver_out: pd.DataFrame) -> pd.DataFrame:
    if solver_out.empty:
        return pd.DataFrame()
    per_step = solver_out[["step", "time", "iterations", "elapsed_time"]].copy()
    per_step["step_wall_seconds"] = per_step["elapsed_time"].diff()
    per_step.loc[per_step.index[0], "step_wall_seconds"] = per_step.loc[per_step.index[0], "elapsed_time"]
    step_events = events[events["step"] > 0]
    phase_seconds = step_events.pivot_table(index="step", columns="phase", values="duration_seconds", aggfunc="sum", fill_value=0.0)
    phase_counts = step_events.pivot_table(index="step", columns="phase", values="event_index", aggfunc="count", fill_value=0)
    for phase in ["analyze", "factorize", "solve"]:
        if phase not in phase_seconds:
            phase_seconds[phase] = 0.0
        if phase not in phase_counts:
            phase_counts[phase] = 0
    per_step = per_step.merge(phase_seconds[["analyze", "factorize", "solve"]], left_on="step", right_index=True, how="left")
    per_step = per_step.merge(
        phase_counts[["analyze", "factorize", "solve"]].rename(columns={c: f"{c}_count" for c in ["analyze", "factorize", "solve"]}),
        left_on="step",
        right_index=True,
        how="left",
    )
    for col in ["analyze", "factorize", "solve", "analyze_count", "factorize_count", "solve_count"]:
        per_step[col] = per_step[col].fillna(0)
    per_step["linear_seconds"] = per_step["analyze"] + per_step["factorize"] + per_step["solve"]
    per_step["nonlinear_overhead"] = (per_step["step_wall_seconds"] - per_step["linear_seconds"]).clip(lower=0)
    per_step["linear_share"] = per_step["linear_seconds"] / per_step["step_wall_seconds"].replace(0, math.nan)
    return per_step


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-dir", default="data/cases/Dogbone/results")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--title", default="Dogbone EigenLDLT Linear-Solve Profile")
    parser.add_argument("--case", default="Dogbone")
    parser.add_argument("--threads", default="")
    args = parser.parse_args()

    profile_dir = Path(args.profile_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    events_path = profile_dir / "linear_profile_events.tsv"
    iterations_path = profile_dir / "linear_profile_iterations.tsv"
    solver_out_path = profile_dir / "solver.out"
    runtime_summary_path = profile_dir / "runtime_profile_summary.tsv"
    solver_inp = {}
    for solver_inp_path in [
        profile_dir / "solver.inp.effective",
        profile_dir / "solver.inp",
        profile_dir.parent / "solver.inp",
    ]:
        solver_inp = read_solver_inp(solver_inp_path)
        if solver_inp:
            break
    if not events_path.exists():
        raise SystemExit(f"Missing profile events file: {events_path}")

    events = pd.read_csv(events_path, sep="\t")
    iterations = pd.read_csv(iterations_path, sep="\t") if iterations_path.exists() else pd.DataFrame()
    solver_out = read_solver_out(solver_out_path)

    solve = events[events["phase"] == "solve"].copy()
    factorize = events[events["phase"] == "factorize"].copy()
    analyze = events[events["phase"] == "analyze"].copy()
    per_step = compute_per_step(events, solver_out)

    total_wall = float(solver_out["elapsed_time"].iloc[-1]) if not solver_out.empty else float(events["duration_seconds"].sum())
    analyze_seconds = float(analyze["duration_seconds"].sum())
    factorize_seconds = float(factorize["duration_seconds"].sum())
    solve_seconds = float(solve["duration_seconds"].sum())
    linear_seconds = analyze_seconds + factorize_seconds + solve_seconds
    overhead_seconds = max(total_wall - linear_seconds, 0.0)

    runtime_share = pd.DataFrame(
        [
            ["analyzePattern", len(analyze), analyze_seconds, analyze_seconds / total_wall],
            ["factorize", len(factorize), factorize_seconds, factorize_seconds / total_wall],
            ["solve", len(solve), solve_seconds, solve_seconds / total_wall],
            ["nonlinear/assembly/export/other", "", overhead_seconds, overhead_seconds / total_wall],
            ["total wall", "", total_wall, 1.0],
        ],
        columns=["bucket", "count", "seconds", "share_of_wall"],
    )
    runtime_phase_hotspots = pd.DataFrame()
    if runtime_summary_path.exists():
        runtime_phase_hotspots = pd.read_csv(runtime_summary_path, sep="\t")
        if not runtime_phase_hotspots.empty and "total_seconds" in runtime_phase_hotspots:
            runtime_phase_hotspots = runtime_phase_hotspots.sort_values("total_seconds", ascending=False).head(20)

    internal_timing = pd.DataFrame()
    timing_columns = [
        ("preconditioner_apply_seconds", "preconditioner apply"),
        ("orthogonalization_seconds", "orthogonalization"),
        ("least_squares_seconds", "least squares"),
        ("matvec_seconds", "matrix-vector"),
        ("deflation_seconds", "deflation projection"),
    ]
    available_timing = [(col, label) for col, label in timing_columns if col in solve.columns]
    if available_timing:
        rows = []
        for col, label in available_timing:
            seconds = float(solve[col].sum())
            rows.append([label, seconds, seconds / solve_seconds if solve_seconds > 0 else math.nan, seconds / total_wall if total_wall > 0 else math.nan])
        accounted = sum(row[1] for row in rows)
        rows.append(["unclassified solve time", max(solve_seconds - accounted, 0.0), max(solve_seconds - accounted, 0.0) / solve_seconds if solve_seconds > 0 else math.nan, max(solve_seconds - accounted, 0.0) / total_wall if total_wall > 0 else math.nan])
        internal_timing = pd.DataFrame(rows, columns=["bucket", "seconds", "share_of_solve_time", "share_of_wall"])

    system_rows = []
    for (system_kind, rhs_kind, solver_type, solver_name), group in solve.groupby(["system_kind", "rhs_kind", "solver_type", "solver_name"], dropna=False):
        matching_factorize = factorize[factorize["system_kind"] == system_kind]
        solver_iterations = group[group["solver_iterations"] >= 0]["solver_iterations"]
        solver_error = group[group["solver_error"] >= 0]["solver_error"]
        system_rows.append(
            {
                "system_kind": system_kind,
                "rhs_kind": rhs_kind,
                "solver_type": solver_type,
                "solver_name": solver_name,
                "rows": int(group["rows"].max()),
                "nnz": int(group["nnz"].max()),
                "solve_count": len(group),
                "factorize_count": len(matching_factorize),
                "success_count": int(group["success"].sum()) if "success" in group else "",
                "iter_median": solver_iterations.median() if not solver_iterations.empty else math.nan,
                "iter_max": solver_iterations.max() if not solver_iterations.empty else math.nan,
                "relres_median": solver_error.median() if not solver_error.empty else math.nan,
                "relres_max": solver_error.max() if not solver_error.empty else math.nan,
                "solve_seconds": group["duration_seconds"].sum(),
                "factorize_seconds": matching_factorize["duration_seconds"].sum(),
            }
        )
    systems = pd.DataFrame(system_rows)

    solver_stats = pd.DataFrame()
    if not solve.empty and "solver_iterations" in solve:
        iterative = solve[solve["solver_iterations"] >= 0].copy()
        if not iterative.empty:
            solver_stats = iterative.groupby(["system_kind", "rhs_kind", "solver_type"], dropna=False).agg(
                solves=("event_index", "count"),
                success=("success", "sum"),
                iter_min=("solver_iterations", "min"),
                iter_median=("solver_iterations", "median"),
                iter_max=("solver_iterations", "max"),
                relres_min=("solver_error", "min"),
                relres_median=("solver_error", "median"),
                relres_max=("solver_error", "max"),
            ).reset_index()

    delta = factorize[factorize["matrix_relative_delta"] >= 0]["matrix_relative_delta"]
    matrix_evolution = pd.DataFrame(
        [
            ["unique_patterns", int(factorize["pattern_hash"].nunique()) if not factorize.empty else 0],
            ["pattern_changes", int((factorize["pattern_hash"] != factorize["pattern_hash"].shift()).sum() - (0 if factorize.empty else 1))],
            ["factorizations", len(factorize)],
            ["solves_per_factorization", len(solve) / len(factorize) if len(factorize) else math.nan],
            ["relative_delta_min", delta.min() if not delta.empty else math.nan],
            ["relative_delta_median", delta.median() if not delta.empty else math.nan],
            ["relative_delta_max", delta.max() if not delta.empty else math.nan],
        ],
        columns=["metric", "value"],
    )

    metadata = pd.DataFrame(
        [
            ["case", args.case],
            ["commit", git_revision()],
            ["threads", args.threads],
            ["profile_dir", str(profile_dir)],
            ["total_wall_seconds", total_wall],
            ["current_solver_inp_type", solver_inp.get("solver_type", "")],
            ["profiled_solver_type", ", ".join(sorted(str(x) for x in solve["solver_type"].dropna().unique()))],
            ["linear_solver_class", ", ".join(sorted(str(x) for x in solve["solver_name"].dropna().unique()))],
            ["oas_solver_class", solver_inp.get("oas_solver_class", "")],
            ["stiffness_matrix_iter_update", solver_inp.get("stiffness_matrix_iter_update", "")],
            ["stiffness_matrix_step_update", solver_inp.get("stiffness_matrix_step_update", "")],
            ["events", len(events)],
            ["steps", int(solver_out["step"].max()) if not solver_out.empty else ""],
            ["nonlinear_iterations", int(solver_out["iterations"].sum()) if not solver_out.empty else ""],
            ["linear_solves", len(solve)],
            ["factorizations", len(factorize)],
        ],
        columns=["field", "value"],
    )

    plot_names = [
        save_stacked_step_time(per_step, out_dir),
        save_iterations_solves(per_step, out_dir),
        save_matrix_delta(factorize, out_dir),
        save_linear_share(per_step, out_dir),
        save_solve_duration_distribution(solve, out_dir),
    ]
    plot_names = [name for name in plot_names if name]

    top_steps = pd.DataFrame()
    if not per_step.empty:
        top_steps = per_step.sort_values("linear_seconds", ascending=False).head(12)[
            ["step", "iterations", "solve_count", "factorize_count", "step_wall_seconds", "linear_seconds", "linear_share"]
        ]

    report = out_dir / "linear-profile.md"
    with report.open("w", encoding="utf-8") as handle:
        handle.write(f"# {args.title}\n\n")
        handle.write("## Metadata\n\n")
        handle.write(md_table(metadata) + "\n\n")
        handle.write("## Runtime Share\n\n")
        handle.write(md_table(runtime_share) + "\n\n")
        if not runtime_phase_hotspots.empty:
            handle.write("## Runtime Phase Hotspots\n\n")
            handle.write(md_table(runtime_phase_hotspots) + "\n\n")
        if not internal_timing.empty and internal_timing["seconds"].sum() > 0:
            handle.write("## DFGMRES Internal Timings\n\n")
            handle.write(md_table(internal_timing) + "\n\n")
        handle.write("## Linear System Types\n\n")
        handle.write(md_table(systems) + "\n\n")
        if not solver_stats.empty:
            handle.write("## Iterative Solver Stats\n\n")
            handle.write(md_table(solver_stats) + "\n\n")
        handle.write("## Matrix Evolution\n\n")
        handle.write(md_table(matrix_evolution) + "\n\n")
        handle.write("## Highest Linear-Time Steps\n\n")
        handle.write(md_table(top_steps) + "\n\n")
        if not iterations.empty:
            handle.write("## Iteration Summary Sample\n\n")
            sample = iterations.sort_values(["step", "iteration"]).head(12)
            handle.write(md_table(sample) + "\n\n")
        handle.write("## Graphs\n\n")
        for name in plot_names:
            handle.write(f"![{name}]({name})\n\n")

    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
