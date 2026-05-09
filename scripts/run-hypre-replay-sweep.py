#!/usr/bin/env python3
"""Run standalone hypre replay sweeps on a dumped OAS TS-N65 linear system."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import itertools
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIR = ROOT.parent / "oas_deflation-build" / "release"
DEFAULT_OAS_BIN = DEFAULT_BUILD_DIR / "bin" / "OAS"
CASE_DIR = ROOT / "data" / "cases" / "TS-N_65"
MASTER_FILE = CASE_DIR / "master.inp"


def split_csv(text: str, cast=str):
    return [cast(item.strip()) for item in text.split(",") if item.strip()]


def reset_affinity() -> str:
    cpu_count = os.cpu_count() or 1
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(range(cpu_count)))
        return ",".join(str(cpu) for cpu in sorted(os.sched_getaffinity(0)))
    return "unavailable"


def run_command(command: list[str], *, env: dict[str, str] | None = None, cwd: Path = ROOT, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True, check=check)


def benchmark_binary(build_dir: Path) -> Path:
    candidates = [
        build_dir / "bin" / "benchmark_hypre_replay",
        build_dir / "benchmark_hypre_replay",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def build_benchmark(build_dir: Path, jobs: int) -> Path:
    command = ["cmake", "--build", str(build_dir), "--target", "benchmark_hypre_replay", f"-j{jobs}"]
    print("Building benchmark:", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)
    binary = benchmark_binary(build_dir)
    if not binary.exists():
        raise RuntimeError(f"benchmark target did not produce {binary}")
    return binary


def generate_replay(args: argparse.Namespace, out_dir: Path) -> Path:
    replay_run = out_dir / "fresh-replay-source"
    replay_run.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.replay_threads)
    env["MKL_NUM_THREADS"] = str(args.replay_threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    env.setdefault("OAS_TIMEOUT", args.oas_timeout)
    affinity = reset_affinity()

    overrides = [
        "total_time=5.000000e-03",
        "time_step=5.000000e-03",
        "stiffness_matrix_iter_update=-5",
        "stiffness_matrix_step_update=1",
        "linear_solver_profile=1",
        "linear_solver_profile_matrix_delta=0",
        "linear_solver_replay_dump=1",
        "linear_solver_replay_limit=1",
        "linear_solver_replay_dir=linear_replay",
        "runtime_phase_profile=1",
        "runtime_phase_profile_file=runtime_profile",
        "max_iterations=1000",
        "dfgmres_tolerance=1e-1",
        "dfgmres_true_tolerance=1e-1",
        "dfgmres_max_iterations=500",
        "dfgmres_restart=80",
        "dfgmres_deflation_vectors=20",
        "dfgmres_deflation_eps=1e-15",
        "dfgmres_collect_newton_steps=1",
        "dfgmres_preconditioner=hypre",
        "dfgmres_reorthogonalize_on_matrix_change=1",
        "dfgmres_elastic_reorder=2",
        "amgcl_near_nullspace=1",
        "amgcl_elastic_full_lift=1",
        "amgcl_use_block_backend=1",
        "hypre_coarsen_type=8",
        "hypre_interp_type=6",
        "hypre_strong_threshold=0.25",
        "hypre_nodal=4",
        "hypre_relax_type=16",
        "hypre_num_sweeps=3",
        "hypre_p_max=4",
        "hypre_boomer_max_iterations=1",
        "hypre_cheby_order=3",
        "hypre_cheby_fraction=-1",
        "hypre_elastic_reorder=2",
        "hypre_threads=0",
        "hypre_use_dof_functions=0",
        "hypre_use_interp_vectors=0",
    ]
    command = [
        str(ROOT / "scripts/run-oas-profile.sh"),
        str(CASE_DIR),
        str(MASTER_FILE),
        str(args.oas_bin),
        str(args.replay_threads),
        "DeflatedFGMRES",
        "TS-N65 DFGMRES hypre replay dump",
        "TS-N_65",
        str(replay_run),
        *overrides,
    ]
    (replay_run / "affinity-before-launch.txt").write_text(affinity + "\n", encoding="utf-8")
    print("Generating fresh replay dump:", replay_run, flush=True)
    completed = run_command(command, env=env, check=False)
    (replay_run / "profile-wrapper.stdout").write_text(completed.stdout, encoding="utf-8")
    (replay_run / "profile-wrapper.stderr").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"fresh replay generation failed with status {completed.returncode}; see {replay_run}")
    replay_dir = replay_run / "linear_replay" / "solve_1"
    required = ["matrix.mtx", "rhs.tsv", "solution.tsv", "metadata.tsv", "near_nullspace.tsv", "summary.tsv"]
    missing = [name for name in required if not (replay_dir / name).exists()]
    if missing:
        raise RuntimeError(f"fresh replay dump is missing {missing} in {replay_dir}")
    return replay_dir


def run_benchmark(
    binary: Path,
    replay_dir: Path,
    out_dir: Path,
    params: dict[str, Any],
    repeat_index: int = 0,
    timeout_seconds: int = 0,
) -> dict[str, str]:
    threads = int(params["threads"])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    affinity = reset_affinity()

    command = [
        str(binary),
        "--replay-dir",
        str(replay_dir),
        "--representation",
        str(params["representation"]),
        "--nullspace",
        str(params["nullspace"]),
        "--solver",
        str(params.get("solver", "fgmres")),
        "--threads",
        str(threads),
        "--nodal",
        str(params.get("nodal", 4)),
        "--interp-type",
        str(params.get("interp_type", 6)),
        "--strong-threshold",
        str(params.get("strong_threshold", 0.25)),
        "--relax-type",
        str(params.get("relax_type", 16)),
        "--num-sweeps",
        str(params.get("num_sweeps", 3)),
        "--non-galerkin-tol",
        str(params.get("non_galerkin_tol", -1.0)),
        "--interp-vec-variant",
        str(params.get("interp_vec_variant", 2)),
        "--tolerance",
        str(params.get("tolerance", 1e-1)),
        "--max-iterations",
        str(params.get("max_iterations", 500)),
        "--apply-repeat",
        str(params.get("apply_repeat", 5)),
        "--apply-warmup",
        str(params.get("apply_warmup", 1)),
        "--csv",
    ]
    label = (
        f"{params['representation']}__ns-{params['nullspace']}__solver-{params.get('solver', 'fgmres')}"
        f"__t{threads}__nodal-{params.get('nodal', 4)}__interp-{params.get('interp_type', 6)}"
        f"__strong-{params.get('strong_threshold', 0.25)}__relax-{params.get('relax_type', 16)}"
        f"__sweeps-{params.get('num_sweeps', 3)}__ng-{params.get('non_galerkin_tol', -1.0)}__r{repeat_index}"
    ).replace("/", "_")
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    print("Running", label, flush=True)
    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout_seconds if timeout_seconds > 0 else None,
        )
        stdout = completed.stdout
        stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = (exc.stderr or "") + f"\nTimed out after {timeout_seconds} seconds.\n"
        returncode = 124

    (log_dir / f"{label}.stdout").write_text(stdout, encoding="utf-8", errors="replace")
    (log_dir / f"{label}.stderr").write_text(stderr, encoding="utf-8", errors="replace")

    rows = list(csv.DictReader(line for line in stdout.splitlines() if line.strip()))
    if rows:
        row = rows[-1]
    else:
        row = {"success": "0", "message": f"no csv output; returncode={returncode}"}
        for key, value in params.items():
            csv_key = {
                "interp_type": "interp_type",
                "strong_threshold": "strong_threshold",
                "relax_type": "relax_type",
                "num_sweeps": "num_sweeps",
                "non_galerkin_tol": "non_galerkin_tol",
                "interp_vec_variant": "interp_vec_variant",
            }.get(key, key)
            row[csv_key] = str(value)
    if returncode == 124:
        row["success"] = "0"
        row["message"] = f"timeout_after_{timeout_seconds}_seconds"
    row["returncode"] = str(returncode)
    row["repeat_index"] = str(repeat_index)
    row["affinity"] = affinity
    row["stage"] = str(params.get("stage", ""))
    return row


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def as_float(row: dict[str, str], key: str, default: float = math.nan) -> float:
    try:
        return float(row.get(key, "") or default)
    except ValueError:
        return default


def as_int(row: dict[str, str], key: str, default: int = 0) -> int:
    try:
        return int(float(row.get(key, "") or default))
    except ValueError:
        return default


def total_solve_time(row: dict[str, str]) -> float:
    return as_float(row, "preconditioner_setup_seconds", 0.0) + as_float(row, "solver_setup_seconds", 0.0) + as_float(row, "solve_seconds", 0.0)


def stage1_params(threads: list[int], apply_repeat: int) -> list[dict[str, Any]]:
    representations = [
        "reduced-original",
        "full-node-major",
        "full-coordinate-node-major",
        "full-rcm-node-major",
        "full-component-major",
        "compressed-free-node-major",
    ]
    nullspaces = ["off", "rbm6"]
    params = []
    for representation, nullspace, thread in itertools.product(representations, nullspaces, threads):
        params.append(
            {
                "stage": "stage1",
                "representation": representation,
                "nullspace": nullspace,
                "solver": "fgmres",
                "threads": thread,
                "nodal": 4,
                "interp_type": 6,
                "strong_threshold": 0.25,
                "relax_type": 16,
                "num_sweeps": 3,
                "non_galerkin_tol": -1.0,
                "interp_vec_variant": 2,
                "apply_repeat": apply_repeat,
            }
        )
    return params


def choose_stage2_bases(rows: list[dict[str, str]], limit: int = 2) -> list[tuple[str, str]]:
    candidates = [
        row for row in rows
        if row.get("stage") == "stage1"
        and row.get("success") == "1"
        and row.get("solver") == "fgmres"
        and row.get("threads") == "16"
        and as_int(row, "iterations", 999999) > 0
    ]
    candidates.sort(key=total_solve_time)
    bases: list[tuple[str, str]] = []
    for row in candidates:
        key = (row.get("representation", ""), row.get("nullspace", ""))
        if key not in bases:
            bases.append(key)
        if len(bases) >= limit:
            break
    return bases


def stage2_params(bases: list[tuple[str, str]], threads: list[int], apply_repeat: int, max_rows: int, full: bool) -> list[dict[str, Any]]:
    nodal_values = [1, 4]
    interp_values = [6, 10, 11, 17]
    strong_values = [0.25, 0.5, 0.7]
    relax_values = [16, 18, 6]
    sweep_values = [1, 3]
    nongal_values = [-1.0, 0.0, 0.01, 0.05]
    all_params = []
    for (representation, nullspace), nodal, interp, strong, relax, sweeps, nongal, thread in itertools.product(
        bases, nodal_values, interp_values, strong_values, relax_values, sweep_values, nongal_values, threads
    ):
        all_params.append(
            {
                "stage": "stage2",
                "representation": representation,
                "nullspace": nullspace,
                "solver": "fgmres",
                "threads": thread,
                "nodal": nodal,
                "interp_type": interp,
                "strong_threshold": strong,
                "relax_type": relax,
                "num_sweeps": sweeps,
                "non_galerkin_tol": nongal,
                "interp_vec_variant": 2,
                "apply_repeat": apply_repeat,
            }
        )
    if full or max_rows <= 0:
        return all_params
    # Keep a representative bounded slice: current settings, high-potential relaxers,
    # and a spread across strength/interpolation options.
    preferred = []
    for params in all_params:
        score = 0
        if params["threads"] == 16:
            score -= 3
        if params["nodal"] == 4:
            score -= 2
        if params["interp_type"] in (6, 17):
            score -= 2
        if params["strong_threshold"] in (0.25, 0.5):
            score -= 1
        if params["relax_type"] in (16, 18):
            score -= 2
        if params["num_sweeps"] == 3:
            score -= 1
        if params["non_galerkin_tol"] in (-1.0, 0.01):
            score -= 1
        preferred.append((score, params))
    preferred.sort(key=lambda item: (item[0], str(item[1])))
    selected = [params for _, params in preferred[:max_rows]]
    # Include matching 1-thread rows for selected 16-thread variants.
    selected_keys = {
        tuple((k, v) for k, v in params.items() if k != "threads")
        for params in selected if params["threads"] == 16
    }
    for params in all_params:
        key = tuple((k, v) for k, v in params.items() if k != "threads")
        if params["threads"] == 1 and key in selected_keys and params not in selected:
            selected.append(params)
    return selected


def markdown_table(rows: list[list[str]], headers: list[str]) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def format_seconds(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.3g}"


def generate_report(path: Path, rows: list[dict[str, str]], replay_dir: Path) -> None:
    successful = [row for row in rows if row.get("success") == "1"]
    failed = [row for row in rows if row.get("success") != "1"]
    by_variant: dict[tuple[str, str, str, str, str, str, str, str, str, str], dict[str, dict[str, str]]] = {}
    for row in successful:
        key = (
            row.get("representation", ""),
            row.get("nullspace", ""),
            row.get("solver", ""),
            row.get("nodal", ""),
            row.get("interp_type", ""),
            row.get("strong_threshold", ""),
            row.get("relax_type", ""),
            row.get("num_sweeps", ""),
            row.get("non_galerkin_tol", ""),
            row.get("interp_vec_variant", ""),
        )
        by_variant.setdefault(key, {})[row.get("threads", "")] = row

    paired_rows = []
    for key, pair in by_variant.items():
        j16 = pair.get("16")
        if not j16:
            continue
        j1 = pair.get("1")
        speedup = ""
        if j1:
            baseline = total_solve_time(j1)
            candidate = total_solve_time(j16)
            if baseline > 0 and candidate > 0:
                speedup = f"{baseline / candidate:.2f}x"
        paired_rows.append((total_solve_time(j16), key, j1, j16, speedup))
    paired_rows.sort(key=lambda item: item[0])

    fastest = paired_rows[:15]
    all_table = []
    for _, key, j1, j16, speedup in paired_rows:
        setup1 = format_seconds(as_float(j1, "preconditioner_setup_seconds")) if j1 else ""
        setup16 = format_seconds(as_float(j16, "preconditioner_setup_seconds"))
        apply1 = format_seconds(as_float(j1, "apply_mean_seconds")) if j1 else ""
        apply16 = format_seconds(as_float(j16, "apply_mean_seconds"))
        solve1 = format_seconds(as_float(j1, "solve_seconds")) if j1 else ""
        solve16 = format_seconds(as_float(j16, "solve_seconds"))
        all_table.append(
            [
                key[0],
                key[1],
                key[2],
                key[3],
                key[4],
                key[5],
                key[6],
                key[7],
                key[8],
                setup1,
                setup16,
                apply1,
                apply16,
                solve1,
                solve16,
                str(as_int(j16, "iterations", -1)),
                f"{as_float(j16, 'true_residual'):.3e}",
                f"{as_float(j16, 'solution_error'):.3e}",
                speedup,
            ]
        )

    with path.open("w", encoding="utf-8") as file:
        file.write("# TS-N65 standalone hypre replay sweep\n\n")
        file.write(f"- Replay dump: `{replay_dir}`\n")
        file.write(f"- Rows tested: `{len(rows)}`; successful: `{len(successful)}`; failed/skipped: `{len(failed)}`\n")
        if successful:
            first = successful[0]
            file.write(f"- Matrix rows: `{first.get('rows', '')}`; nnz: `{first.get('nnz', '')}`\n")
            file.write(f"- Full elastic rows: `{first.get('full_rows', '')}`; block size: `{first.get('block_size', '')}`; near-nullspace columns: `{first.get('near_nullspace_columns', '')}`\n")
        file.write("\n## Fastest 16-thread variants\n\n")
        file.write(markdown_table(
            [
                [
                    key[0],
                    key[1],
                    key[3],
                    key[4],
                    key[5],
                    key[6],
                    key[7],
                    key[8],
                    format_seconds(total_solve_time(j16)),
                    format_seconds(as_float(j16, "preconditioner_setup_seconds")),
                    format_seconds(as_float(j16, "apply_mean_seconds")),
                    format_seconds(as_float(j16, "solve_seconds")),
                    str(as_int(j16, "iterations", -1)),
                    f"{as_float(j16, 'true_residual'):.3e}",
                    speedup,
                ]
                for _, key, _, j16, speedup in fastest
            ],
            ["representation", "nullspace", "nodal", "interp", "strong", "relax", "sweeps", "non-gal", "total16 s", "setup16 s", "apply16 s", "solve16 s", "iters16", "true residual16", "speedup"],
        ))
        if successful:
            def top_rows(metric: str, count: int = 8) -> list[list[str]]:
                candidates = [row for row in successful if row.get("threads") == "16"]
                candidates.sort(key=lambda row: as_float(row, metric, math.inf))
                return [
                    [
                        row.get("representation", ""),
                        row.get("nullspace", ""),
                        row.get("nodal", ""),
                        row.get("interp_type", ""),
                        row.get("strong_threshold", ""),
                        row.get("relax_type", ""),
                        row.get("num_sweeps", ""),
                        row.get("non_galerkin_tol", ""),
                        format_seconds(as_float(row, metric)),
                        str(as_int(row, "iterations", -1)),
                        f"{as_float(row, 'true_residual'):.3e}",
                    ]
                    for row in candidates[:count]
                ]

            file.write("\n\n## Best 16-thread setup-only variants\n\n")
            file.write(markdown_table(
                top_rows("preconditioner_setup_seconds"),
                ["representation", "nullspace", "nodal", "interp", "strong", "relax", "sweeps", "non-gal", "setup16 s", "iters16", "res16"],
            ))
            file.write("\n\n## Best 16-thread apply-only variants\n\n")
            file.write(markdown_table(
                top_rows("apply_mean_seconds"),
                ["representation", "nullspace", "nodal", "interp", "strong", "relax", "sweeps", "non-gal", "apply16 s", "iters16", "res16"],
            ))
            file.write("\n\n## Best 16-thread Krylov-solve variants\n\n")
            file.write(markdown_table(
                top_rows("solve_seconds"),
                ["representation", "nullspace", "nodal", "interp", "strong", "relax", "sweeps", "non-gal", "solve16 s", "iters16", "res16"],
            ))
        file.write("\n\n## All paired successful variants\n\n")
        file.write(markdown_table(
            all_table,
            ["representation", "nullspace", "solver", "nodal", "interp", "strong", "relax", "sweeps", "non-gal", "setup1", "setup16", "apply1", "apply16", "solve1", "solve16", "iters16", "res16", "solerr16", "speedup"],
        ))
        if failed:
            file.write("\n\n## Failed or skipped variants\n\n")
            file.write(markdown_table(
                [
                    [
                        row.get("stage", ""),
                        row.get("representation", ""),
                        row.get("nullspace", ""),
                        row.get("solver", ""),
                        row.get("threads", ""),
                        row.get("message", ""),
                    ]
                    for row in failed[:80]
                ],
                ["stage", "representation", "nullspace", "solver", "threads", "message"],
            ))
        if fastest:
            best = fastest[0]
            _, key, _, j16, speedup = best
            file.write("\n\n## Recommendation from this standalone sweep\n\n")
            file.write(
                f"The fastest successful 16-thread standalone variant is `{key[0]}` with nullspace `{key[1]}`, "
                f"`nodal={key[3]}`, `interp_type={key[4]}`, `strong_threshold={key[5]}`, "
                f"`relax_type={key[6]}`, `num_sweeps={key[7]}`, and `non_galerkin_tol={key[8]}`. "
                f"Its standalone 16-thread total setup+solve time is `{format_seconds(total_solve_time(j16))} s`"
            )
            if speedup:
                file.write(f" with `{speedup}` speedup versus the matching 1-thread replay.\n")
            else:
                file.write(".\n")
            file.write("Production DFGMRES should only adopt this after a full TS-N65 first-step run confirms unchanged outer iterations and lower wall time.\n")


def main() -> int:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument("--oas-bin", type=Path, default=DEFAULT_OAS_BIN)
    parser.add_argument("--out-dir", type=Path, default=ROOT / f"results/tsn65-hypre-replay-sweep-{stamp}")
    parser.add_argument("--replay-dir", type=Path)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-generate-replay", action="store_true")
    parser.add_argument("--threads", default="1,16")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--replay-threads", type=int, default=16)
    parser.add_argument("--oas-timeout", default="2h")
    parser.add_argument("--apply-repeat", type=int, default=5)
    parser.add_argument("--stage", choices=["smoke", "stage1", "stage2", "full"], default="stage1")
    parser.add_argument("--stage2-limit", type=int, default=48)
    parser.add_argument("--full-stage2", action="store_true")
    parser.add_argument("--row-timeout", type=int, default=300)
    parser.add_argument("--seed-results", type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    binary = benchmark_binary(args.build_dir) if args.skip_build else build_benchmark(args.build_dir, args.jobs)
    if args.replay_dir:
        replay_dir = args.replay_dir
    elif args.skip_generate_replay:
        raise SystemExit("--skip-generate-replay requires --replay-dir")
    else:
        replay_dir = generate_replay(args, args.out_dir)

    threads = split_csv(args.threads, int)
    rows: list[dict[str, str]] = []
    if args.seed_results:
        with args.seed_results.open("r", encoding="utf-8") as file:
            rows.extend(csv.DictReader(file))

    if args.stage == "smoke":
        params = [{
            "stage": "smoke",
            "representation": "full-coordinate-node-major",
            "nullspace": "off",
            "solver": "fgmres",
            "threads": threads[0],
            "apply_repeat": max(1, min(args.apply_repeat, 2)),
        }]
    elif args.seed_results and args.stage in ("stage2", "full"):
        params = []
    else:
        params = stage1_params(threads, args.apply_repeat)

    for item in params:
        rows.append(run_benchmark(binary, replay_dir, args.out_dir, item, timeout_seconds=args.row_timeout))
        write_csv(args.out_dir / "results.csv", rows)

    if args.stage in ("stage2", "full"):
        bases = choose_stage2_bases(rows)
        stage2 = stage2_params(bases, threads, args.apply_repeat, args.stage2_limit, args.full_stage2 or args.stage == "full")
        for item in stage2:
            rows.append(run_benchmark(binary, replay_dir, args.out_dir, item, timeout_seconds=args.row_timeout))
            write_csv(args.out_dir / "results.csv", rows)

    generate_report(args.out_dir / "summary.md", rows, replay_dir)
    print(f"Wrote {args.out_dir / 'summary.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
