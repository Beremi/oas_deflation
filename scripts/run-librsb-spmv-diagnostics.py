#!/usr/bin/env python3
"""Run pure SpMV scaling diagnostics on an OAS linear replay dump."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPLAY_DIR = (
    ROOT
    / "results/amgcl-tsn65-first-solve-20260503-201037"
    / "002-hybrid-aggr-blocksize1/linear_replay/solve_1"
)
DEFAULT_PKG_CONFIG_PATH = Path("/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig")


def parse_args() -> argparse.Namespace:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-dir", type=Path, default=DEFAULT_REPLAY_DIR)
    parser.add_argument("--out-dir", type=Path, default=ROOT / f"results/librsb-spmv-scaling-diagnostics-{stamp}")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--bindings", default="close,spread")
    parser.add_argument("--rsb-flags", default="noflags,default")
    parser.add_argument("--subdivisions", default="1")
    parser.add_argument("--full-subdivision-sweep", action="store_true")
    parser.add_argument("--backend", default="all", choices=["eigen", "manual", "librsb", "all"])
    parser.add_argument("--vector", default="solution", choices=["solution", "rhs", "ones"])
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--repeat", type=int, default=20)
    parser.add_argument("--validate-tol", default="1e-10")
    parser.add_argument("--tune", action="store_true")
    parser.add_argument("--tune-rounds", type=int, default=2)
    parser.add_argument("--tune-seconds", default="2.0")
    parser.add_argument("--pkg-config-path", type=Path, default=DEFAULT_PKG_CONFIG_PATH)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def split_csv(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def run_checked(command: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True, check=True)


def pkg_config(args: list[str], env: dict[str, str]) -> list[str]:
    result = run_checked(["pkg-config", *args], env=env)
    return shlex.split(result.stdout.strip())


def build_benchmark(args: argparse.Namespace, env: dict[str, str]) -> Path:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    binary = args.out_dir / "benchmark-librsb-spmv-replay"
    source = ROOT / "scripts/benchmark-librsb-spmv-replay.cpp"

    cflags = pkg_config(["--cflags", "eigen3"], env)
    rsb_flags = pkg_config(["--cflags", "--libs", "librsb"], env)
    libdir_result = run_checked(["pkg-config", "--variable=libdir", "librsb"], env=env)
    libdir = libdir_result.stdout.strip()

    command = [
        args.cxx,
        "-O3",
        "-DNDEBUG",
        "-march=native",
        "-std=c++17",
        "-fopenmp",
        str(source),
        *cflags,
        *rsb_flags,
    ]
    if libdir:
        command.append(f"-Wl,-rpath,{libdir}")
    command.extend(["-o", str(binary)])

    with (args.out_dir / "build-command.txt").open("w") as file:
        file.write(" ".join(shlex.quote(part) for part in command) + "\n")

    print("Compiling benchmark:", " ".join(shlex.quote(part) for part in command))
    if not args.dry_run:
        subprocess.run(command, cwd=ROOT, env=env, check=True)
    return binary


def benchmark_command(
    binary: Path,
    args: argparse.Namespace,
    thread_count: str,
    rsb_flags: str,
    subdivisions: str,
) -> list[str]:
    command = [
        str(binary),
        "--replay-dir",
        str(args.replay_dir),
        "--backend",
        args.backend,
        "--vector",
        args.vector,
        "--threads",
        thread_count,
        "--rsb-flags",
        rsb_flags,
        "--subdivisions",
        subdivisions,
        "--warmup",
        str(args.warmup),
        "--repeat",
        str(args.repeat),
        "--validate-tol",
        str(args.validate_tol),
        "--csv",
    ]
    if args.tune:
        command.extend(["--tune", "--tune-rounds", str(args.tune_rounds), "--tune-seconds", str(args.tune_seconds)])
    return command


def rows_from_stdout(stdout: str) -> list[dict[str, str]]:
    lines = [line for line in stdout.splitlines() if line.strip()]
    if not lines:
        return []
    reader = csv.DictReader(lines)
    return list(reader)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def as_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, "") or default)
    except ValueError:
        return default


def speedups(rows: Iterable[dict[str, str]]) -> dict[tuple[str, str, str, str, str], float]:
    baselines: dict[tuple[str, str, str, str, str], float] = {}
    for row in rows:
        if row.get("threads_requested") != "1":
            continue
        key = (
            row.get("omp_proc_bind", ""),
            row.get("backend", ""),
            row.get("rsb_flags_name", ""),
            row.get("subdivision", ""),
            row.get("tuned", ""),
        )
        baselines[key] = as_float(row, "seconds_per_spmv")
    return baselines


def format_speedup(row: dict[str, str], baselines: dict[tuple[str, str, str, str, str], float]) -> str:
    key = (
        row.get("omp_proc_bind", ""),
        row.get("backend", ""),
        row.get("rsb_flags_name", ""),
        row.get("subdivision", ""),
        row.get("tuned", ""),
    )
    baseline = baselines.get(key, 0.0)
    current = as_float(row, "seconds_per_spmv")
    if baseline <= 0.0 or current <= 0.0:
        return ""
    return f"{baseline / current:.2f}x"


def generate_report(path: Path, args: argparse.Namespace, rows: list[dict[str, str]]) -> None:
    baselines = speedups(rows)
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            row.get("omp_proc_bind", ""),
            row.get("backend", ""),
            row.get("rsb_flags_name", ""),
            as_float(row, "subdivision"),
            int(row.get("threads_requested", "0") or 0),
        ),
    )

    successful = [row for row in rows if row.get("success") == "1"]
    failed = [row for row in rows if row.get("success") != "1"]
    best_16 = sorted(
        [row for row in successful if row.get("threads_requested") == "16"],
        key=lambda row: as_float(row, "seconds_per_spmv"),
    )[:10]

    with path.open("w") as file:
        file.write("# libRSB SpMV Scaling Diagnostics\n\n")
        file.write(f"- Replay dump: `{args.replay_dir}`\n")
        file.write(f"- Vector: `{args.vector}`\n")
        file.write(f"- Warmup/repeat: `{args.warmup}/{args.repeat}`\n")
        file.write(f"- Backends: `{args.backend}`\n")
        file.write(f"- Threads: `{args.threads}`\n")
        file.write(f"- Bindings: `{args.bindings}`\n")
        file.write(f"- RSB flags: `{args.rsb_flags}`\n")
        file.write(f"- Subdivisions: `{args.subdivisions}`\n")
        file.write(f"- Tuning enabled: `{args.tune}`\n")
        file.write(f"- Rows recorded: `{len(rows)}`; failed rows: `{len(failed)}`\n\n")

        if failed:
            file.write("## Failed Rows\n\n")
            file.write("| binding | backend | threads | flags | subdivision | rel error | message |\n")
            file.write("|---|---:|---:|---:|---:|---:|---|\n")
            for row in failed:
                file.write(
                    f"| {row.get('omp_proc_bind', '')} | {row.get('backend', '')} | "
                    f"{row.get('threads_requested', '')} | {row.get('rsb_flags_name', '')} | "
                    f"{row.get('subdivision', '')} | {row.get('relative_error', '')} | "
                    f"{row.get('message', '')} |\n"
                )
            file.write("\n")

        file.write("## Fastest 16 Thread Rows\n\n")
        file.write("| binding | backend | flags | subdivision | RSB threads | leaves | s/SpMV | speedup | GF/s | GiB/s | rel err |\n")
        file.write("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in best_16:
            file.write(
                f"| {row.get('omp_proc_bind', '')} | {row.get('backend', '')} | "
                f"{row.get('rsb_flags_name', '')} | {row.get('subdivision', '')} | "
                f"{row.get('rsb_actual_threads', '') or row.get('threads_requested', '')} | "
                f"{row.get('rsb_leaves', '')} | {as_float(row, 'seconds_per_spmv'):.6g} | "
                f"{format_speedup(row, baselines)} | {as_float(row, 'gflops'):.3f} | "
                f"{as_float(row, 'gib_per_s'):.3f} | {as_float(row, 'relative_error'):.3e} |\n"
            )
        file.write("\n")

        file.write("## Full Results\n\n")
        file.write("| binding | backend | flags | subdivision | threads | actual RSB threads | leaves | flags hex | s/SpMV | speedup | GF/s | GiB/s | rel err |\n")
        file.write("|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|\n")
        for row in sorted_rows:
            file.write(
                f"| {row.get('omp_proc_bind', '')} | {row.get('backend', '')} | "
                f"{row.get('rsb_flags_name', '')} | {row.get('subdivision', '')} | "
                f"{row.get('threads_requested', '')} | {row.get('rsb_actual_threads', '')} | "
                f"{row.get('rsb_leaves', '')} | {row.get('rsb_matrix_flags_hex', '')} | "
                f"{as_float(row, 'seconds_per_spmv'):.6g} | {format_speedup(row, baselines)} | "
                f"{as_float(row, 'gflops'):.3f} | {as_float(row, 'gib_per_s'):.3f} | "
                f"{as_float(row, 'relative_error'):.3e} |\n"
            )


def main() -> int:
    args = parse_args()
    if args.full_subdivision_sweep:
        args.subdivisions = "0.5,1,2,4"

    if not (args.replay_dir / "matrix.mtx").exists():
        print(f"Replay matrix not found: {args.replay_dir / 'matrix.mtx'}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    existing_pkg_config = env.get("PKG_CONFIG_PATH", "")
    env["PKG_CONFIG_PATH"] = (
        str(args.pkg_config_path)
        if not existing_pkg_config
        else f"{args.pkg_config_path}:{existing_pkg_config}"
    )

    binary = build_benchmark(args, env)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    threads = split_csv(args.threads)
    bindings = split_csv(args.bindings)
    rsb_flags = ",".join(split_csv(args.rsb_flags))
    subdivisions = ",".join(split_csv(args.subdivisions))

    rows: list[dict[str, str]] = []
    failures = 0
    run_log = args.out_dir / "run.log"
    with run_log.open("w") as log:
        for binding in bindings:
            for thread_count in threads:
                run_env = env.copy()
                run_env["OMP_NUM_THREADS"] = thread_count
                run_env["RSB_NUM_THREADS"] = thread_count
                run_env["OMP_DYNAMIC"] = "FALSE"
                run_env["OMP_PROC_BIND"] = binding
                run_env["OMP_PLACES"] = "cores"

                command = benchmark_command(binary, args, thread_count, rsb_flags, subdivisions)
                print(f"Running bind={binding} threads={thread_count}")
                log.write("$ " + " ".join(shlex.quote(part) for part in command) + "\n")
                log.write(f"OMP_NUM_THREADS={thread_count} RSB_NUM_THREADS={thread_count} OMP_PROC_BIND={binding}\n")
                log.flush()

                if args.dry_run:
                    continue

                completed = subprocess.run(command, cwd=ROOT, env=run_env, text=True, capture_output=True)
                raw_name = f"bind-{binding}-threads-{thread_count}"
                (args.out_dir / f"{raw_name}.stdout.csv").write_text(completed.stdout)
                (args.out_dir / f"{raw_name}.stderr.log").write_text(completed.stderr)
                log.write(completed.stderr)
                log.write(completed.stdout)
                log.write(f"\nreturncode={completed.returncode}\n\n")
                log.flush()

                rows.extend(rows_from_stdout(completed.stdout))
                if completed.returncode != 0:
                    failures += 1

    if not args.dry_run:
        write_csv(args.out_dir / "results.csv", rows)
        generate_report(args.out_dir / "summary.md", args, rows)
        print(f"Wrote {args.out_dir / 'results.csv'}")
        print(f"Wrote {args.out_dir / 'summary.md'}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
