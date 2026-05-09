#!/usr/bin/env python3
"""Run standalone AMGCL replay sweeps on a dumped OAS TS-N65 linear system."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIR = ROOT.parent / "oas_deflation-build" / "release"
DEFAULT_REPLAY_DIR = ROOT / "results/tsn65-hypre-replay-sweep-20260508-smoke/fresh-replay-source/linear_replay/solve_1"
DEFAULT_HYPRE_CSV = ROOT / "results/tsn65-hypre-replay-sweep-20260508-stage2-paired/results.csv"


def split_csv(text: str, cast=str):
    return [cast(item.strip()) for item in text.split(",") if item.strip()]


def benchmark_binary(build_dir: Path) -> Path:
    candidates = [build_dir / "bin" / "benchmark_amgcl_replay", build_dir / "benchmark_amgcl_replay"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def reset_affinity() -> str:
    cpu_count = os.cpu_count() or 1
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(range(cpu_count)))
        return ",".join(str(cpu) for cpu in sorted(os.sched_getaffinity(0)))
    return "unavailable"


def build_benchmark(build_dir: Path, jobs: int) -> Path:
    subprocess.run(["cmake", "--build", str(build_dir), "--target", "benchmark_amgcl_replay", f"-j{jobs}"], cwd=ROOT, check=True)
    binary = benchmark_binary(build_dir)
    if not binary.exists():
        raise RuntimeError(f"benchmark target did not produce {binary}")
    return binary


def base_params(threads: list[int], apply_repeat: int) -> list[dict[str, Any]]:
    params: list[dict[str, Any]] = []

    scalar_reps = ["reduced-original", "compressed-free-node-major"]
    full_reps = ["full-coordinate-node-major", "full-rcm-node-major"]
    eps_values = [0.0, 0.02, 0.08]

    for rep in scalar_reps:
        for eps in eps_values:
            for thread in threads:
                params.append({
                    "stage": "scalar-stage1",
                    "representation": rep,
                    "nullspace": "off",
                    "solver": "cg",
                    "backend": "scalar",
                    "relaxation": "spai0",
                    "eps_strong": eps,
                    "relax": 1.0,
                    "block_size": 1,
                    "npre": 1,
                    "npost": 1,
                    "ncycle": 1,
                    "threads": thread,
                    "apply_repeat": apply_repeat,
                })

    for rep in full_reps:
        for backend in ["scalar", "hybrid6"]:
            for nullspace in ["off", "rbm6"]:
                for thread in threads:
                    params.append({
                        "stage": "full-stage1",
                        "representation": rep,
                        "nullspace": nullspace,
                        "solver": "cg",
                        "backend": backend,
                        "relaxation": "spai0",
                        "eps_strong": 0.0,
                        "relax": 1.0,
                        "block_size": 6,
                        "npre": 1,
                        "npost": 1,
                        "ncycle": 1,
                        "threads": thread,
                        "apply_repeat": apply_repeat,
                    })

    # Focused parameter sweep around the scalar forms that historically survive TS-N65.
    for rep in scalar_reps:
        for relaxation in ["spai0", "damped_jacobi", "ilu0"]:
            for solver in ["cg", "fgmres"]:
                for eps in [0.0, 0.02, 0.08, 0.16]:
                    if relaxation == "ilu0" and solver == "fgmres":
                        continue
                    for npre, npost in [(1, 1), (2, 2)]:
                        for thread in threads:
                            params.append({
                                "stage": "scalar-stage2",
                                "representation": rep,
                                "nullspace": "off",
                                "solver": solver,
                                "backend": "scalar",
                                "relaxation": relaxation,
                                "eps_strong": eps,
                                "relax": 1.0,
                                "block_size": 1,
                                "npre": npre,
                                "npost": npost,
                                "ncycle": 1,
                                "threads": thread,
                                "apply_repeat": apply_repeat,
                            })
    # Remove duplicate dictionaries while preserving order.
    out: list[dict[str, Any]] = []
    seen = set()
    for item in params:
        key = tuple(sorted(item.items()))
        if key not in seen:
            seen.add(key)
            out.append(item)
    return out


def run_one(binary: Path, replay_dir: Path, out_dir: Path, params: dict[str, Any], timeout: int) -> dict[str, str]:
    threads = int(params["threads"])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    affinity = reset_affinity()

    command = [
        str(binary),
        "--replay-dir", str(replay_dir),
        "--representation", str(params["representation"]),
        "--nullspace", str(params["nullspace"]),
        "--solver", str(params["solver"]),
        "--backend", str(params["backend"]),
        "--amgcl-relaxation", str(params["relaxation"]),
        "--threads", str(threads),
        "--eps-strong", str(params["eps_strong"]),
        "--relax", str(params["relax"]),
        "--block-size", str(params["block_size"]),
        "--npre", str(params["npre"]),
        "--npost", str(params["npost"]),
        "--ncycle", str(params["ncycle"]),
        "--tolerance", "1e-1",
        "--max-iterations", "500",
        "--apply-repeat", str(params["apply_repeat"]),
        "--apply-warmup", "1",
        "--csv",
    ]
    label = (
        f"{params['representation']}__{params['nullspace']}__{params['solver']}__{params['backend']}"
        f"__{params['relaxation']}__eps-{params['eps_strong']}__pre-{params['npre']}__post-{params['npost']}__t{threads}"
    ).replace("/", "_")
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    print("Running", label, flush=True)
    try:
        completed = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True, timeout=timeout, check=False)
        stdout = completed.stdout
        stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = (exc.stderr or "") + f"\nTimed out after {timeout} seconds.\n"
        returncode = 124
    (log_dir / f"{label}.stdout").write_text(stdout, encoding="utf-8", errors="replace")
    (log_dir / f"{label}.stderr").write_text(stderr, encoding="utf-8", errors="replace")
    rows = list(csv.DictReader(line for line in stdout.splitlines() if line.strip()))
    if rows:
        row = rows[-1]
    else:
        row = {key: str(value) for key, value in params.items()}
        row["success"] = "0"
        row["message"] = f"no csv output; returncode={returncode}"
    if returncode == 124:
        row["success"] = "0"
        row["message"] = f"timeout_after_{timeout}_seconds"
    row["returncode"] = str(returncode)
    row["affinity"] = affinity
    row["stage"] = str(params["stage"])
    return row


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str] | None, key: str, default: float = math.nan) -> float:
    if row is None:
        return default
    try:
        return float(row.get(key, "") or default)
    except ValueError:
        return default


def total(row: dict[str, str] | None) -> float:
    return f(row, "preconditioner_setup_seconds", 0.0) + f(row, "solver_setup_seconds", 0.0) + f(row, "solve_seconds", 0.0)


def markdown_table(rows: list[list[str]], headers: list[str]) -> str:
    if not rows:
        return "_No rows._"
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
        *("| " + " | ".join(row) + " |" for row in rows),
    ])


def variant_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(k, "") for k in [
        "representation", "nullspace", "solver", "backend", "relaxation",
        "eps_strong", "relax", "block_size_option", "npre", "npost", "ncycle",
    ])


def load_hypre_best(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as file:
        rows = [row for row in csv.DictReader(file) if row.get("success") == "1" and row.get("threads") == "16"]
    rows.sort(key=total)
    return rows[:8]


def generate_report(path: Path, rows: list[dict[str, str]], replay_dir: Path, hypre_csv: Path) -> None:
    successful = [row for row in rows if row.get("success") == "1"]
    failed = [row for row in rows if row.get("success") != "1"]
    by_key: dict[tuple[str, ...], dict[str, dict[str, str]]] = {}
    for row in successful:
        by_key.setdefault(variant_key(row), {})[row.get("threads", "")] = row

    paired = []
    for key, pair in by_key.items():
        j16 = pair.get("16")
        if not j16:
            continue
        j1 = pair.get("1")
        speedup = total(j1) / total(j16) if j1 and total(j16) > 0 else math.nan
        paired.append((total(j16), speedup, key, j1, j16))
    paired.sort(key=lambda item: item[0])

    fastest = [
        [
            key[0], key[1], key[2], key[3], key[4], key[5], key[8], key[9],
            f"{t16:.3g}",
            f"{f(j16, 'preconditioner_setup_seconds'):.3g}",
            f"{f(j16, 'apply_mean_seconds'):.3g}",
            f"{f(j16, 'solve_seconds'):.3g}",
            str(int(f(j16, "iterations", -1))),
            f"{f(j16, 'true_residual'):.3e}",
            f"{speedup:.2f}x" if math.isfinite(speedup) else "",
        ]
        for t16, speedup, key, _, j16 in paired[:20]
    ]

    solve_sorted = sorted(paired, key=lambda item: f(item[4], "solve_seconds"))
    solve_table = [
        [
            key[0], key[1], key[2], key[3], key[4], key[5], key[8], key[9],
            f"{f(j16, 'solve_seconds'):.3g}",
            f"{f(j16, 'apply_mean_seconds'):.3g}",
            str(int(f(j16, "iterations", -1))),
            f"{f(j16, 'true_residual'):.3e}",
        ]
        for _, _, key, _, j16 in solve_sorted[:15]
    ]

    hypre_rows = load_hypre_best(hypre_csv)
    hypre_table = [
        [
            row.get("representation", ""),
            row.get("nullspace", ""),
            row.get("interp_type", ""),
            row.get("strong_threshold", ""),
            row.get("relax_type", ""),
            row.get("num_sweeps", ""),
            f"{total(row):.3g}",
            f"{f(row, 'preconditioner_setup_seconds'):.3g}",
            f"{f(row, 'apply_mean_seconds'):.3g}",
            f"{f(row, 'solve_seconds'):.3g}",
            row.get("iterations", ""),
            f"{f(row, 'true_residual'):.3e}",
        ]
        for row in hypre_rows
    ]

    with path.open("w", encoding="utf-8") as file:
        file.write("# TS-N65 standalone AMGCL replay sweep\n\n")
        file.write(f"- Replay dump: `{replay_dir}`\n")
        file.write(f"- Rows tested: `{len(rows)}`; successful: `{len(successful)}`; failed/skipped: `{len(failed)}`\n")
        if successful:
            first = successful[0]
            file.write(f"- Matrix rows: `{first.get('rows', '')}`; nnz: `{first.get('nnz', '')}`\n")
        file.write("\n## Fastest AMGCL 16-thread variants\n\n")
        file.write(markdown_table(
            fastest,
            ["representation", "nullspace", "solver", "backend", "relaxation", "eps", "npre", "npost", "total16", "setup16", "apply16", "solve16", "iters16", "res16", "speedup"],
        ))
        file.write("\n\n## Fastest AMGCL Krylov solve only\n\n")
        file.write(markdown_table(
            solve_table,
            ["representation", "nullspace", "solver", "backend", "relaxation", "eps", "npre", "npost", "solve16", "apply16", "iters16", "res16"],
        ))
        file.write("\n\n## Fastest hypre rows from previous replay sweep\n\n")
        file.write(markdown_table(
            hypre_table,
            ["representation", "nullspace", "interp", "strong", "relax", "sweeps", "total16", "setup16", "apply16", "solve16", "iters16", "res16"],
        ))
        if failed:
            file.write("\n\n## Failed or skipped AMGCL variants\n\n")
            file.write(markdown_table(
                [[row.get("stage", ""), row.get("representation", ""), row.get("nullspace", ""), row.get("solver", ""), row.get("backend", ""), row.get("relaxation", ""), row.get("threads", ""), row.get("message", "")] for row in failed[:80]],
                ["stage", "representation", "nullspace", "solver", "backend", "relaxation", "threads", "message"],
            ))


def main() -> int:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument("--replay-dir", type=Path, default=DEFAULT_REPLAY_DIR)
    parser.add_argument("--hypre-results", type=Path, default=DEFAULT_HYPRE_CSV)
    parser.add_argument("--out-dir", type=Path, default=ROOT / f"results/tsn65-amgcl-replay-sweep-{stamp}")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--threads", default="1,16")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--apply-repeat", type=int, default=2)
    parser.add_argument("--row-timeout", type=int, default=180)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    binary = benchmark_binary(args.build_dir) if args.skip_build else build_benchmark(args.build_dir, args.jobs)
    params = base_params(split_csv(args.threads, int), args.apply_repeat)
    if args.limit > 0:
        params = params[:args.limit]

    rows: list[dict[str, str]] = []
    for item in params:
        rows.append(run_one(binary, args.replay_dir, args.out_dir, item, args.row_timeout))
        write_csv(args.out_dir / "results.csv", rows)
    generate_report(args.out_dir / "summary.md", rows, args.replay_dir, args.hypre_results)
    print(f"Wrote {args.out_dir / 'summary.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
