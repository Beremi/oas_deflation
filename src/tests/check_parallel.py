import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def copy_case(src, dst):
    shutil.copytree(src, dst)


def replace_or_append_solver_line(solver_file, key, value):
    lines = solver_file.read_text(encoding="utf-8").splitlines()
    replaced = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.split()[0] == key:
            lines[i] = f"{key}\t{value}"
            replaced = True
            break
    if not replaced:
        lines.append(f"{key}\t{value}")
    solver_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


def configure_dfgmres_none(case_dir):
    solver_file = case_dir / "solver.inp"
    replace_or_append_solver_line(solver_file, "solver_type", "DeflatedFGMRES")
    replace_or_append_solver_line(solver_file, "dfgmres_preconditioner", "none")
    replace_or_append_solver_line(solver_file, "dfgmres_tolerance", "1e-12")
    replace_or_append_solver_line(solver_file, "dfgmres_true_tolerance", "1e-12")
    replace_or_append_solver_line(solver_file, "dfgmres_max_iterations", "200")
    replace_or_append_solver_line(solver_file, "dfgmres_restart", "40")
    replace_or_append_solver_line(solver_file, "dfgmres_deflation_vectors", "0")
    replace_or_append_solver_line(solver_file, "amgcl_check_matrix", "1")
    replace_or_append_solver_line(solver_file, "amgcl_near_nullspace", "0")
    replace_or_append_solver_line(solver_file, "amgcl_elastic_full_lift", "0")
    replace_or_append_solver_line(solver_file, "amgcl_use_block_backend", "0")


def run_case(executable, case_dir, threads):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    proc = subprocess.run(
        [str(executable), "-j", str(threads), "master.inp"],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if proc.returncode != 0:
        sys.stdout.write(proc.stdout)
        raise RuntimeError(f"OAS failed in {case_dir} with {threads} threads")


def numeric_rows(path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.replace(",", " ").split()
        try:
            rows.append([float(token) for token in tokens])
        except ValueError:
            continue
    if not rows:
        return None
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        return None
    return np.array(rows, dtype=float)


def comparable_result_files(results_dir):
    files = {}
    for path in results_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.name == "version.txt" or path.suffix.lower() in {".vtk", ".vtu", ".pvtu"}:
            continue
        data = numeric_rows(path)
        if data is not None:
            files[path.relative_to(results_dir)] = data
    return files


def compare_results(serial_dir, parallel_dir, rtol, atol):
    serial = comparable_result_files(serial_dir / "results")
    parallel = comparable_result_files(parallel_dir / "results")
    if not serial:
        raise RuntimeError("No numeric serial result files were produced")
    if set(serial) != set(parallel):
        missing = sorted(str(path) for path in set(serial) ^ set(parallel))
        raise RuntimeError(f"Result file mismatch: {missing}")
    for relpath, serial_data in serial.items():
        parallel_data = parallel[relpath]
        if serial_data.shape != parallel_data.shape:
            raise RuntimeError(f"Shape mismatch for {relpath}: {serial_data.shape} != {parallel_data.shape}")
        if not np.allclose(serial_data, parallel_data, rtol=rtol, atol=atol):
            diff = np.max(np.abs(serial_data - parallel_data))
            raise RuntimeError(f"Numeric mismatch for {relpath}; max_abs_diff={diff}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument("case_dir", type=pathlib.Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--dfgmres-none", action="store_true")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="oas-parallel-") as tmp:
        tmpdir = pathlib.Path(tmp)
        serial_dir = tmpdir / "serial"
        parallel_dir = tmpdir / "parallel"
        copy_case(args.case_dir, serial_dir)
        copy_case(args.case_dir, parallel_dir)
        if args.dfgmres_none:
            configure_dfgmres_none(serial_dir)
            configure_dfgmres_none(parallel_dir)
        run_case(args.executable, serial_dir, 1)
        run_case(args.executable, parallel_dir, args.threads)
        compare_results(serial_dir, parallel_dir, args.rtol, args.atol)


if __name__ == "__main__":
    main()
