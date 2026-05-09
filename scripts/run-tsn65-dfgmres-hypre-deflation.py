#!/usr/bin/env python3
"""Run TS-N_65 two-step DFGMRES+hypre deflation experiments."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "data/cases/TS-N_65"
MASTER_FILE = CASE_DIR / "master.inp"
PROFILE_DIR = CASE_DIR / "results"
DEFAULT_OAS_BIN = ROOT.parent / "oas_deflation-build/release/bin/OAS"
DEFAULT_REFERENCE_DIR = ROOT / "results/tsn65-two-step-comparison-20260503-093031/pardisoldlt-reference"
VTK_EXPORT_LINE = "VTKElementExporter state saveSteps 2 1 2 ascii pointData 1 displacements"


@dataclass
class RunSpec:
    name: str
    nvec: int
    tolerance: str


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str | None, default: float = math.nan) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def as_int(value: str | None, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except ValueError:
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
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(format_cell(item) for item in row) + " |")
    return "\n".join(lines)


def read_solver_out(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    header: list[str] = []
    rows: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            maybe_header = line.lstrip("#").strip()
            if maybe_header:
                header = maybe_header.split()
            continue
        if not header:
            continue
        parts = line.split()
        rows.append({key: parts[i] if i < len(parts) else "" for i, key in enumerate(header)})
    return rows


def slug_tolerance(tolerance: str) -> str:
    return tolerance.replace("-", "m").replace("+", "p").replace(".", "p")


def patch_key_value_file(path: Path, updates: dict[str, str], marker: str) -> str:
    original = path.read_text(encoding="utf-8", errors="replace")
    lines = original.splitlines()
    out: list[str] = []
    seen: set[str] = set()
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            out.append(line)
            continue
        key = stripped.split()[0]
        if key in updates:
            out.append(f"{key}\t{updates[key]}")
            seen.add(key)
        else:
            out.append(line)
    missing = [key for key in updates if key not in seen]
    if missing:
        if out and out[-1].strip():
            out.append("")
        out.append(f"# {marker}")
        for key in missing:
            out.append(f"{key}\t{updates[key]}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return original


def patch_exporters(path: Path) -> str:
    original = path.read_text(encoding="utf-8", errors="replace")
    lines = original.splitlines()
    if not any(line.strip().startswith("VTKElementExporter state ") for line in lines):
        if lines and lines[-1].strip():
            lines.append("")
        lines.append("# Local TS-N_65 DFGMRES adaptive-tolerance VTK export.")
        lines.append(VTK_EXPORT_LINE)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return original


def clean_case_results() -> None:
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    for pattern in [
        "linear_profile_events.tsv",
        "linear_profile_iterations.tsv",
        "runtime_profile_events.tsv",
        "runtime_profile_summary.tsv",
        "solver.out",
        "LD.out",
        "state_*.vtu",
        "state_*.pvtu",
        "state_*.pvd",
        "*.visit",
    ]:
        for path in PROFILE_DIR.glob(pattern):
            if path.is_file():
                path.unlink()


def solver_updates(
    nvec: int,
    eps: str,
    time_step: str | None,
    total_time: str,
    linear_tol: str,
    true_tol: str,
    stiffness_matrix_iter_update: str,
    stiffness_matrix_step_update: str,
    hypre_relax_type: str,
    hypre_num_sweeps: str,
    hypre_boomer_max_iterations: str,
    hypre_cheby_order: str,
    hypre_cheby_fraction: str,
    hypre_elastic_reorder: str,
    hypre_threads: str,
    hypre_nodal_diag: str,
    adaptive: bool,
    loose_tol: str,
    tight_tol: str,
    trigger_ratio: str,
    require_tight_convergence: bool,
) -> dict[str, str]:
    updates = {
        "total_time": total_time,
        "solver_type": "DeflatedFGMRES",
        "linear_solver_profile": "1",
        "linear_solver_profile_matrix_delta": "0",
        "linear_solver_profile_file": "linear_profile",
        "runtime_phase_profile": "1",
        "runtime_phase_profile_file": "runtime_profile",
        "stiffness_matrix_iter_update": stiffness_matrix_iter_update,
        "stiffness_matrix_step_update": stiffness_matrix_step_update,
        "max_iterations": "1000",
        "dfgmres_tolerance": linear_tol,
        "dfgmres_true_tolerance": true_tol,
        "dfgmres_max_iterations": "500",
        "dfgmres_restart": "80",
        "dfgmres_deflation_vectors": str(nvec),
        "dfgmres_deflation_eps": eps,
        "dfgmres_collect_newton_steps": "1",
        "dfgmres_preconditioner": "hypre",
        "dfgmres_reorthogonalize_on_matrix_change": "1",
        "dfgmres_elastic_reorder": hypre_elastic_reorder,
        "amgcl_near_nullspace": "1",
        "amgcl_elastic_full_lift": "1",
        "amgcl_use_block_backend": "1",
        "hypre_coarsen_type": "8",
        "hypre_interp_type": "6",
        "hypre_strong_threshold": "0.25",
        "hypre_nodal": "4",
        "hypre_nodal_diag": hypre_nodal_diag,
        "hypre_relax_type": hypre_relax_type,
        "hypre_num_sweeps": hypre_num_sweeps,
        "hypre_p_max": "4",
        "hypre_boomer_max_iterations": hypre_boomer_max_iterations,
        "hypre_cheby_order": hypre_cheby_order,
        "hypre_cheby_fraction": hypre_cheby_fraction,
        "hypre_elastic_reorder": hypre_elastic_reorder,
        "hypre_threads": hypre_threads,
        "hypre_use_dof_functions": "0",
        "hypre_use_interp_vectors": "0",
    }
    if adaptive:
        updates.update(
            {
                "linear_solver_adaptive_tolerance": "1",
                "linear_solver_adaptive_loose_tolerance": loose_tol,
                "linear_solver_adaptive_tight_tolerance": tight_tol,
                "linear_solver_adaptive_trigger_ratio": trigger_ratio,
                "linear_solver_adaptive_require_tight_convergence": "1" if require_tight_convergence else "0",
            }
        )
    if time_step:
        updates["time_step"] = time_step
    return updates


def copy_artifacts(run_dir: Path) -> None:
    for name in [
        "linear_profile_events.tsv",
        "linear_profile_iterations.tsv",
        "runtime_profile_events.tsv",
        "runtime_profile_summary.tsv",
        "solver.out",
        "LD.out",
        "state_00001.vtu",
        "state_00002.vtu",
        "version.txt",
    ]:
        src = PROFILE_DIR / name
        if src.exists():
            shutil.copy2(src, run_dir / name)


def analyze_linear_profile(run_dir: Path, title: str, threads: int) -> None:
    if not (run_dir / "linear_profile_events.tsv").exists():
        return
    cmd = [
        sys.executable,
        str(ROOT / "scripts/analyze-linear-profile.py"),
        "--profile-dir",
        str(run_dir),
        "--out-dir",
        str(run_dir),
        "--title",
        title,
        "--case",
        "TS-N_65",
        "--threads",
        str(threads),
    ]
    subprocess.run(cmd, cwd=str(ROOT), check=False)


def summarize_run(run_dir: Path, spec: RunSpec, returncode: int, wall_seconds: float, eps: str) -> dict[str, Any]:
    events = read_tsv(run_dir / "linear_profile_events.tsv")
    solver_rows = read_solver_out(run_dir / "solver.out")
    solve = [row for row in events if row.get("phase") == "solve"]
    factorize = [row for row in events if row.get("phase") == "factorize"]
    analyze = [row for row in events if row.get("phase") == "analyze"]
    iterations = [as_int(row.get("solver_iterations"), -1) for row in solve if as_int(row.get("solver_iterations"), -1) >= 0]
    residuals = [as_float(row.get("solver_error")) for row in solve if as_float(row.get("solver_error")) >= 0]
    final_solver_row = solve[-1] if solve else {}
    final_solver_error = as_float(final_solver_row.get("solver_error"), None)
    requested_tolerances: dict[str, int] = {}
    for row in solve:
        tolerance = row.get("requested_tolerance") or ""
        if tolerance:
            requested_tolerances[tolerance] = requested_tolerances.get(tolerance, 0) + 1
    successes = [as_int(row.get("success"), 0) for row in solve]
    max_step = max([as_int(row.get("step"), 0) for row in solver_rows], default=0)
    nonlinear_iterations = sum(as_int(row.get("iterations"), 0) for row in solver_rows)
    final_nonlinear_row = solver_rows[-1] if solver_rows else {}
    final_nonlinear_residual = as_float(final_nonlinear_row.get("error_residuals"), None)
    final_displacement_error = as_float(final_nonlinear_row.get("error_displacements"), None)
    final_energy_error = as_float(final_nonlinear_row.get("error_energy"), None)
    elapsed_solver = max([as_float(row.get("elapsed_time"), 0.0) for row in solver_rows], default=0.0)
    elapsed_wall = max(elapsed_solver, wall_seconds)
    solve_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in solve)
    setup_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in factorize)
    analyze_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in analyze)
    precond_apply = sum(as_float(row.get("preconditioner_apply_seconds"), 0.0) for row in solve)
    orth = sum(as_float(row.get("orthogonalization_seconds"), 0.0) for row in solve)
    matvec = sum(as_float(row.get("matvec_seconds"), 0.0) for row in solve)
    least_squares = sum(as_float(row.get("least_squares_seconds"), 0.0) for row in solve)
    deflation = sum(as_float(row.get("deflation_seconds"), 0.0) for row in solve)
    linear_seconds = analyze_seconds + setup_seconds + solve_seconds
    value_hashes = {row.get("value_hash", "") for row in factorize if row.get("value_hash")}
    basis_values = [as_int(row.get("deflation_basis_size"), 0) for row in solve]
    log_text = (run_dir / "oas.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "oas.log").exists() else ""
    log_basis_values = [int(match.group(1)) for match in re.finditer(r"active_basis_size=(\d+)", log_text)]
    max_basis = max(basis_values + log_basis_values, default=0)
    final_basis = log_basis_values[-1] if log_basis_values else (basis_values[-1] if basis_values else 0)
    max_offdiag = max([as_float(row.get("deflation_orthogonality_max_offdiag"), 0.0) for row in solve], default=0.0)
    max_diag_error = max([as_float(row.get("deflation_orthogonality_max_diag_error"), 0.0) for row in solve], default=0.0)
    discarded = max([as_int(row.get("deflation_discarded_count"), 0) for row in solve], default=0)
    capacity_evictions = max([as_int(row.get("deflation_capacity_eviction_count"), 0) for row in solve], default=0)
    status = "completed_2_steps" if returncode == 0 and max_step >= 2 else "partial"
    residual_failed = False
    for row in solve:
        residual = as_float(row.get("solver_error"))
        requested = as_float(row.get("requested_tolerance"), 1e-6)
        if math.isfinite(residual) and math.isfinite(requested) and residual > requested * (1. + 1e-8):
            residual_failed = True
            break
    if solve and (not all(successes) or residual_failed):
        status = "linear_failed"
    if returncode != 0 and max_step == 0 and not solve:
        status = "failed"

    return {
        "name": spec.name,
        "N": spec.nvec,
        "tolerance": spec.tolerance,
        "eps": eps,
        "status": status,
        "returncode": returncode,
        "steps": max_step,
        "nonlinear_iterations": nonlinear_iterations,
        "linear_solves": len(solve),
        "factorizations": len(factorize),
        "unique_matrices": len(value_hashes),
        "outer_iterations": sum(iterations),
        "median_iterations": median(iterations) if iterations else None,
        "max_iterations": max(iterations, default=None),
        "max_true_residual": max(residuals, default=None),
        "final_true_residual": final_solver_error,
        "final_nonlinear_residual": final_nonlinear_residual,
        "final_displacement_error": final_displacement_error,
        "final_energy_error": final_energy_error,
        "requested_tolerances": requested_tolerances,
        "max_basis": max_basis,
        "final_basis": final_basis,
        "discarded": discarded,
        "capacity_evictions": capacity_evictions,
        "max_offdiag": max_offdiag,
        "max_diag_error": max_diag_error,
        "wall_seconds": elapsed_wall,
        "analyze_seconds": analyze_seconds,
        "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds,
        "preconditioner_apply_seconds": precond_apply,
        "orthogonalization_seconds": orth,
        "matvec_seconds": matvec,
        "least_squares_seconds": least_squares,
        "deflation_seconds": deflation,
        "linear_share": linear_seconds / elapsed_wall if elapsed_wall > 0 else None,
        "other_seconds": max(elapsed_wall - linear_seconds, 0.0),
        "runtime_hotspots": top_runtime_hotspots(run_dir),
    }


def read_vtu_displacements(path: Path, array_name: str = "displacements") -> np.ndarray | None:
    if not path.exists():
        return None
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        return None
    arr = grid.GetPointData().GetArray(array_name)
    if arr is None:
        return None
    data = vtk_to_numpy(arr).astype(float, copy=False)
    if data.ndim == 1:
        data = data[:, None]
    return data


def compare_vtu(candidate: Path, reference: Path) -> dict[str, Any]:
    cand = read_vtu_displacements(candidate)
    ref = read_vtu_displacements(reference)
    if cand is None or ref is None:
        return {"available": False}
    if cand.shape != ref.shape:
        return {"available": False, "reason": f"shape {cand.shape} != {ref.shape}"}
    diff = cand - ref
    node_norm = np.linalg.norm(diff, axis=1)
    denom = float(np.linalg.norm(ref))
    return {
        "available": True,
        "relative_l2": float(np.linalg.norm(diff) / denom) if denom > 0 else None,
        "rms_nodal": float(np.sqrt(np.mean(node_norm * node_norm))),
        "max_nodal": float(np.max(node_norm)),
        "max_component": float(np.max(np.abs(diff))),
    }


def add_reference_comparisons(manifest: dict[str, Any], run_dir: Path, reference_dir: Path | None) -> None:
    comparisons: dict[str, Any] = {}
    if reference_dir and reference_dir.exists():
        for step in (1, 2):
            comparisons[f"step_{step}"] = compare_vtu(run_dir / f"state_{step:05d}.vtu", reference_dir / f"state_{step:05d}.vtu")
    manifest["pardiso_vtu_comparison"] = comparisons


def top_runtime_hotspots(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    rows = read_tsv(run_dir / "runtime_profile_summary.tsv")
    parsed: list[dict[str, Any]] = []
    for row in rows:
        parsed.append(
            {
                "phase": row.get("phase", ""),
                "detail": row.get("detail", ""),
                "count": as_int(row.get("count"), 0),
                "total_seconds": as_float(row.get("total_seconds"), 0.0),
                "mean_seconds": as_float(row.get("mean_seconds"), 0.0),
                "max_seconds": as_float(row.get("max_seconds"), 0.0),
            }
        )
    parsed.sort(key=lambda row: row["total_seconds"], reverse=True)
    return parsed[:limit]


def run_spec(
    spec: RunSpec,
    out_dir: Path,
    oas_bin: Path,
    threads: int,
    eps: str,
    time_step: str | None,
    total_time: str,
    stiffness_matrix_iter_update: str,
    stiffness_matrix_step_update: str,
    hypre_relax_type: str,
    hypre_num_sweeps: str,
    hypre_boomer_max_iterations: str,
    hypre_cheby_order: str,
    hypre_cheby_fraction: str,
    hypre_elastic_reorder: str,
    hypre_threads: str,
    hypre_nodal_diag: str,
    adaptive: bool,
    loose_tol: str,
    tight_tol: str,
    trigger_ratio: str,
    require_tight_convergence: bool,
    reference_dir: Path | None,
    omp_proc_bind: str,
    omp_places: str,
) -> dict[str, Any]:
    run_dir = out_dir / spec.name
    run_dir.mkdir(parents=True, exist_ok=True)
    solver_inp = CASE_DIR / "solver.inp"
    exporters_inp = CASE_DIR / "exporters.inp"
    original = solver_inp.read_text(encoding="utf-8", errors="replace")
    exporters_original = exporters_inp.read_text(encoding="utf-8", errors="replace")
    (run_dir / "solver.inp.original").write_text(original, encoding="utf-8")
    (run_dir / "exporters.inp.original").write_text(exporters_original, encoding="utf-8")
    started = datetime.now().isoformat(timespec="seconds")
    start_time = time.monotonic()
    returncode = 999
    try:
        patch_key_value_file(
            solver_inp,
            solver_updates(
                spec.nvec,
                eps,
                time_step,
                total_time,
                spec.tolerance,
                spec.tolerance,
                stiffness_matrix_iter_update,
                stiffness_matrix_step_update,
                hypre_relax_type,
                hypre_num_sweeps,
                hypre_boomer_max_iterations,
                hypre_cheby_order,
                hypre_cheby_fraction,
                hypre_elastic_reorder,
                hypre_threads,
                hypre_nodal_diag,
                adaptive,
                loose_tol,
                tight_tol,
                trigger_ratio,
                require_tight_convergence,
            ),
            "Local TS-N_65 DFGMRES-hypre deflation settings.",
        )
        patch_exporters(exporters_inp)
        (run_dir / "solver.inp.effective").write_text(solver_inp.read_text(encoding="utf-8"), encoding="utf-8")
        (run_dir / "exporters.inp.effective").write_text(exporters_inp.read_text(encoding="utf-8"), encoding="utf-8")
        clean_case_results()
        log_path = run_dir / "oas.log"
        env = os.environ.copy()
        env.setdefault("MKLROOT", "/opt/intel/oneapi/mkl/latest")
        env["MKL_NUM_THREADS"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_DYNAMIC"] = "FALSE"
        if omp_proc_bind == "unset":
            env.pop("OMP_PROC_BIND", None)
        else:
            env["OMP_PROC_BIND"] = omp_proc_bind
        if omp_places == "unset":
            env.pop("OMP_PLACES", None)
        else:
            env["OMP_PLACES"] = omp_places
        affinity_before = ""
        affinity_after = ""
        if hasattr(os, "sched_getaffinity"):
            try:
                affinity_before = ",".join(str(cpu) for cpu in sorted(os.sched_getaffinity(0)))
            except OSError:
                affinity_before = "<unavailable>"
        if omp_proc_bind == "unset" and hasattr(os, "sched_setaffinity"):
            try:
                os.sched_setaffinity(0, set(range(os.cpu_count() or 1)))
            except OSError:
                pass
        if hasattr(os, "sched_getaffinity"):
            try:
                affinity_after = ",".join(str(cpu) for cpu in sorted(os.sched_getaffinity(0)))
            except OSError:
                affinity_after = "<unavailable>"
        cmd = [str(oas_bin), "-j", str(threads), str(MASTER_FILE)]
        print(f"[{started}] starting {spec.name}: {' '.join(cmd)}", flush=True)
        with log_path.open("w", encoding="utf-8", errors="replace") as log:
            log.write(f"Starting {spec.name} at {started}\n")
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write(f"Threads: {threads}\n")
            log.write(f"OMP_DYNAMIC: {env.get('OMP_DYNAMIC', '')}\n")
            log.write(f"OMP_PROC_BIND: {env.get('OMP_PROC_BIND', '<unset>')}\n")
            log.write(f"OMP_PLACES: {env.get('OMP_PLACES', '<unset>')}\n")
            if affinity_before or affinity_after:
                log.write(f"CPU affinity before runner reset: {affinity_before}\n")
                log.write(f"CPU affinity after runner reset: {affinity_after}\n")
            log.write(f"Deflation vectors: {spec.nvec}\n")
            log.write(f"Deflation eps: {eps}\n")
            log.write(f"Linear tolerance: {spec.tolerance}\n")
            log.write(f"True tolerance: {spec.tolerance}\n")
            if adaptive:
                log.write(f"Adaptive tolerance: loose={loose_tol}, tight={tight_tol}, trigger_ratio={trigger_ratio}\n")
            log.flush()
            process = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT, text=True)
            while True:
                try:
                    returncode = process.wait(timeout=60)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - start_time
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {spec.name} still running ({elapsed:.0f}s)", flush=True)
                    log.write(f"\n[runner] still running at {elapsed:.0f}s\n")
                    log.flush()
            log.write(f"\nOAS exit status: {returncode}\nFinished at {datetime.now().isoformat(timespec='seconds')}\n")
        print(f"[{datetime.now().isoformat(timespec='seconds')}] finished {spec.name} with status {returncode}", flush=True)
    finally:
        solver_inp.write_text(original, encoding="utf-8")
        exporters_inp.write_text(exporters_original, encoding="utf-8")

    elapsed = time.monotonic() - start_time
    copy_artifacts(run_dir)
    analyze_linear_profile(run_dir, f"TS-N_65 DeflatedFGMRES hypre N={spec.nvec} Two-Step Profile", threads)
    manifest = summarize_run(run_dir, spec, returncode, elapsed, eps)
    add_reference_comparisons(manifest, run_dir, reference_dir)
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def make_plots(out_dir: Path, manifests: list[dict[str, Any]]) -> list[str]:
    plots: list[str] = []
    completed = [row for row in manifests if (out_dir / row["name"] / "linear_profile_events.tsv").exists()]
    if not completed:
        return plots
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for row in completed:
        events = read_tsv(out_dir / row["name"] / "linear_profile_events.tsv")
        solve = [event for event in events if event.get("phase") == "solve"]
        ax.plot(range(1, len(solve) + 1), [as_int(event.get("solver_iterations"), 0) for event in solve], label=f"tol={row.get('tolerance', '')}, N={row['N']}")
    ax.set_xlabel("Linear solve index")
    ax.set_ylabel("DFGMRES iterations")
    ax.legend()
    fig.tight_layout()
    name = "iterations-by-solve.png"
    fig.savefig(out_dir / name, dpi=150)
    plt.close(fig)
    plots.append(name)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for row in completed:
        events = read_tsv(out_dir / row["name"] / "linear_profile_events.tsv")
        solve = [event for event in events if event.get("phase") == "solve"]
        ax.plot(range(1, len(solve) + 1), [as_int(event.get("deflation_basis_size"), 0) for event in solve], label=f"tol={row.get('tolerance', '')}, N={row['N']}")
    ax.set_xlabel("Linear solve index")
    ax.set_ylabel("Basis size")
    ax.legend()
    fig.tight_layout()
    name = "basis-size-by-solve.png"
    fig.savefig(out_dir / name, dpi=150)
    plt.close(fig)
    plots.append(name)
    return plots


def pardiso_reference_aggregate(reference_dir: Path | None) -> list[Any] | None:
    if reference_dir is None:
        return None
    path = reference_dir / "run.json"
    if not path.exists():
        return None
    row = json.loads(path.read_text(encoding="utf-8"))
    setup_seconds = as_float(str(row.get("analyze_seconds", 0.0)), 0.0) + as_float(str(row.get("factorize_seconds", 0.0)), 0.0)
    linear_total = as_float(str(row.get("linear_seconds", 0.0)), 0.0)
    wall = max(
        as_float(str(row.get("wall_seconds_solver", 0.0)), 0.0),
        as_float(str(row.get("wall_seconds_runner", 0.0)), 0.0),
    )
    return [
        "PardisoLDLT reference",
        "direct",
        "",
        row.get("status", ""),
        row.get("completed_steps", ""),
        row.get("nonlinear_iterations", ""),
        row.get("linear_solves", ""),
        row.get("factorizations", ""),
        row.get("krylov_iter_median", ""),
        row.get("krylov_iter_max", ""),
        row.get("true_relres_max", ""),
        row.get("true_relres_max", ""),
        "",
        "",
        "",
        "",
        "",
        setup_seconds,
        row.get("solve_seconds", ""),
        linear_total,
        row.get("other_seconds", ""),
        wall,
        row.get("linear_share", ""),
    ]


def write_summary(out_dir: Path, manifests: list[dict[str, Any]], reference_dir: Path | None = DEFAULT_REFERENCE_DIR) -> None:
    plots = make_plots(out_dir, manifests)
    manifests = sorted(manifests, key=lambda row: (as_float(row.get("tolerance")), as_int(str(row.get("N", 0)))))
    aggregate_rows = []
    pardiso_row = pardiso_reference_aggregate(reference_dir)
    if pardiso_row:
        aggregate_rows.append(pardiso_row)
    run_rows = []
    timing_rows = []
    comparison_rows = []
    hotspot_rows = []
    for row in manifests:
        aggregate_rows.append([
            row["name"],
            row.get("tolerance", ""),
            row["N"],
            row["status"],
            row["steps"],
            row["nonlinear_iterations"],
            row["linear_solves"],
            row["factorizations"],
            row["median_iterations"],
            row["max_iterations"],
            row["max_true_residual"],
            row.get("final_true_residual", ""),
            row.get("final_nonlinear_residual", ""),
            row["max_basis"],
            row["final_basis"],
            row["preconditioner_apply_seconds"],
            row["deflation_seconds"],
            row["analyze_seconds"] + row["setup_seconds"],
            row["solve_seconds"],
            row["analyze_seconds"] + row["setup_seconds"] + row["solve_seconds"],
            row["other_seconds"],
            row["wall_seconds"],
            row["linear_share"],
        ])
        run_rows.append([
            row["name"],
            row.get("tolerance", ""),
            row["N"],
            row["eps"],
            row["status"],
            row["steps"],
            row["nonlinear_iterations"],
            row["linear_solves"],
            row["factorizations"],
            row["unique_matrices"],
            row["outer_iterations"],
            row["median_iterations"],
            row["max_iterations"],
            row["max_true_residual"],
            row.get("final_true_residual", ""),
            row.get("final_nonlinear_residual", ""),
            ", ".join(f"{format_cell(float(k))}:{v}" for k, v in sorted(row.get("requested_tolerances", {}).items())),
            row["max_basis"],
            row["final_basis"],
            row["discarded"],
            row["capacity_evictions"],
            row["max_offdiag"],
        ])
        for key, cmp_row in row.get("pardiso_vtu_comparison", {}).items():
            comparison_rows.append([
                row["name"],
                key.replace("_", " "),
                cmp_row.get("available"),
                cmp_row.get("relative_l2"),
                cmp_row.get("rms_nodal"),
                cmp_row.get("max_nodal"),
                cmp_row.get("max_component"),
                cmp_row.get("reason", ""),
            ])
        timing_rows.append([
            row["name"],
            row["wall_seconds"],
            row["setup_seconds"],
            row["solve_seconds"],
            row["preconditioner_apply_seconds"],
            row["orthogonalization_seconds"],
            row["matvec_seconds"],
            row["least_squares_seconds"],
            row["deflation_seconds"],
            row["other_seconds"],
            row["linear_share"],
        ])
        for hotspot in row.get("runtime_hotspots", [])[:5]:
            hotspot_rows.append([
                row["name"],
                hotspot["phase"],
                hotspot["detail"],
                hotspot["count"],
                hotspot["total_seconds"],
                hotspot["mean_seconds"],
                hotspot["max_seconds"],
            ])
    with (out_dir / "summary.md").open("w", encoding="utf-8") as handle:
        handle.write("# TS-N_65 DFGMRES Hypre Deflation Two-Step Runs\n\n")
        handle.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n\n")
        handle.write("Outer solver is native `DeflatedFGMRES`; hypre BoomerAMG is used only as a right preconditioner. Runs target the first two TS-N_65 load steps with `total_time=0.01`, 16 shared-memory threads, and matrix-delta profiling disabled.\n\n")
        handle.write("## Aggregate Comparison\n\n")
        handle.write(md_table([
            "run", "tol", "N", "status", "steps", "nonlinear iters", "linear solves", "factorizations", "median iter", "max iter", "max true relres", "final true relres", "final nonlinear residual", "max basis", "final basis", "precond apply s", "deflation s", "setup/analyze s", "solve s", "linear total s", "other s", "wall s", "linear share",
        ], aggregate_rows) + "\n\n")
        handle.write("## Run Summary\n\n")
        handle.write(md_table([
            "run", "tol", "N", "eps", "status", "steps", "nonlinear iters", "linear solves", "factorizations", "unique matrices", "outer iters", "median iter", "max iter", "max true relres", "final true relres", "final nonlinear residual", "requested tol counts", "max basis", "final basis", "discarded", "capacity evictions", "max offdiag",
        ], run_rows) + "\n\n")
        if comparison_rows:
            handle.write("## VTU Displacement Closeness vs Pardiso\n\n")
            handle.write(md_table([
                "run", "step", "available", "relative L2", "RMS nodal", "max nodal", "max component", "reason",
            ], comparison_rows) + "\n\n")
        handle.write("## Timing\n\n")
        handle.write(md_table([
            "run", "wall s", "setup s", "solve s", "precond apply s", "orth s", "matvec s", "least squares s", "deflation s", "other s", "linear share",
        ], timing_rows) + "\n\n")
        if hotspot_rows:
            handle.write("## Runtime Phase Hotspots\n\n")
            handle.write(md_table([
                "run", "phase", "detail", "count", "total s", "mean s", "max s",
            ], hotspot_rows) + "\n\n")
        handle.write("## Artifacts\n\n")
        artifact_rows = []
        for row in manifests:
            path = (out_dir / row["name"]).resolve()
            try:
                display_path = path.relative_to(ROOT.resolve())
            except ValueError:
                display_path = path
            artifact_rows.append([row["name"], display_path])
        handle.write(md_table(["run", "folder"], artifact_rows) + "\n\n")
        if plots:
            handle.write("## Plots\n\n")
            for plot in plots:
                handle.write(f"![{plot}]({plot})\n\n")
    print(out_dir / "summary.md")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oas-bin", type=Path, default=DEFAULT_OAS_BIN)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--eps", default="1e-15")
    parser.add_argument("--time-step", default=None)
    parser.add_argument("--total-time", default="1.000000e-02")
    parser.add_argument("--stiffness-matrix-iter-update", default="-5")
    parser.add_argument("--stiffness-matrix-step-update", default="1")
    parser.add_argument("--hypre-relax-type", default="16")
    parser.add_argument("--hypre-num-sweeps", default="3")
    parser.add_argument("--hypre-boomer-max-iterations", default="1")
    parser.add_argument("--hypre-cheby-order", default="3")
    parser.add_argument("--hypre-cheby-fraction", default="-1")
    parser.add_argument("--hypre-elastic-reorder", default="2")
    parser.add_argument("--hypre-threads", default="0")
    parser.add_argument("--hypre-nodal-diag", default="0")
    parser.add_argument(
        "--omp-proc-bind",
        default="unset",
        help="OpenMP binding policy for OAS runs. Default unsets OMP_PROC_BIND because TS-N_65 force kernels are slower when pinned to cores on this host.",
    )
    parser.add_argument(
        "--omp-places",
        default="unset",
        help="OpenMP places for OAS runs. Default unsets OMP_PLACES; use values such as cores only for affinity diagnostics.",
    )
    parser.add_argument("--linear-tol", default="1e-6")
    parser.add_argument("--linear-tols", default=None)
    parser.add_argument("--true-tol", default=None)
    parser.add_argument("--adaptive-tolerance", action="store_true")
    parser.add_argument("--loose-tol", default="1e-1")
    parser.add_argument("--tight-tol", default="1e-6")
    parser.add_argument("--trigger-ratio", default="10")
    parser.add_argument("--allow-loose-step-acceptance", action="store_true")
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--nvecs", default="0,5,10")
    parser.add_argument("--report-only", type=Path, default=None)
    args = parser.parse_args()

    if args.report_only:
        manifests = []
        for run_json in sorted(args.report_only.glob("*/run.json")):
            manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
        write_summary(args.report_only, manifests, args.reference_dir)
        return 0

    if not args.oas_bin.exists() or not os.access(args.oas_bin, os.X_OK):
        raise SystemExit(f"OAS binary is missing or not executable: {args.oas_bin}")
    if not MASTER_FILE.exists():
        raise SystemExit(f"TS-N_65 master input is missing: {MASTER_FILE}")

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = args.out_dir or ROOT / f"results/tsn65-dfgmres-hypre-deflation-{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    tolerances = [item.strip() for item in (args.linear_tols or args.linear_tol).split(",") if item.strip()]
    nvecs = [int(item.strip()) for item in args.nvecs.split(",") if item.strip()]
    specs = [
        RunSpec(f"dfgmres-hypre-tol{slug_tolerance(tolerance)}-N{nvec}", nvec, tolerance)
        for tolerance in tolerances
        for nvec in nvecs
    ]
    if args.only:
        specs = [spec for spec in specs if any(token in spec.name for token in args.only)]
    manifests = []
    for run_json in sorted(out_dir.glob("*/run.json")):
        manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
    seen = {item["name"] for item in manifests}
    for spec in specs:
        if spec.name in seen:
            continue
        manifest = run_spec(
            spec,
            out_dir,
            args.oas_bin,
            args.threads,
            args.eps,
            args.time_step,
            args.total_time,
            args.stiffness_matrix_iter_update,
            args.stiffness_matrix_step_update,
            args.hypre_relax_type,
            args.hypre_num_sweeps,
            args.hypre_boomer_max_iterations,
            args.hypre_cheby_order,
            args.hypre_cheby_fraction,
            args.hypre_elastic_reorder,
            args.hypre_threads,
            args.hypre_nodal_diag,
            args.adaptive_tolerance,
            args.loose_tol,
            args.tight_tol,
            args.trigger_ratio,
            not args.allow_loose_step_acceptance,
            args.reference_dir,
            args.omp_proc_bind,
            args.omp_places,
        )
        manifests.append(manifest)
        write_summary(out_dir, manifests, args.reference_dir)
    write_summary(out_dir, manifests, args.reference_dir)
    (out_dir / "campaign.json").write_text(json.dumps(manifests, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
