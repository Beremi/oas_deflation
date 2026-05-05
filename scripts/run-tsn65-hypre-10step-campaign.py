#!/usr/bin/env python3
"""Run TS-N_65 ten-step hypre-CG tolerance validation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "data/cases/TS-N_65"
MASTER_FILE = CASE_DIR / "master.inp"
PROFILE_DIR = CASE_DIR / "results"
DEFAULT_OAS_BIN = ROOT.parent / "oas_deflation-build/vtk-release/bin/OAS"
DEFAULT_REFERENCE_DIR = ROOT / "results/tsn65-two-step-comparison-20260503-093031/pardisoldlt-reference"


@dataclass
class RunSpec:
    name: str
    tolerance: str


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_solver_out(path: Path) -> list[dict[str, str]]:
    """Read OAS solver.out, which is whitespace-separated with a # header."""
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


def patch_exporters(path: Path, save_steps: list[int]) -> str:
    original = path.read_text(encoding="utf-8", errors="replace")
    lines = [
        line for line in original.splitlines()
        if not line.strip().startswith("VTKElementExporter state ")
    ]
    if lines and lines[-1].strip():
        lines.append("")
    steps = " ".join(str(step) for step in save_steps)
    lines.append("# Local TS-N_65 ten-step hypre validation VTK export.")
    lines.append(f"VTKElementExporter state saveSteps {len(save_steps)} {steps} ascii pointData 1 displacements")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return original


def solver_updates(tolerance: str, total_time: str) -> dict[str, str]:
    return {
        "total_time": total_time,
        "solver_type": "HypreBoomerAMGCG",
        "linear_solver_profile": "1",
        "linear_solver_profile_matrix_delta": "0",
        "linear_solver_profile_file": "linear_profile",
        "runtime_phase_profile": "1",
        "runtime_phase_profile_file": "runtime_profile",
        "stiffness_matrix_iter_update": "-5",
        "stiffness_matrix_step_update": "1",
        "max_iterations": "1000",
        "hypre_tolerance": tolerance,
        "hypre_max_iterations": "500",
        "hypre_coarsen_type": "8",
        "hypre_strong_threshold": "0.25",
        "hypre_nodal": "4",
        "hypre_interp_type": "6",
        "hypre_relax_type": "6",
        "hypre_p_max": "4",
        "hypre_flex": "1",
        "hypre_skip_break": "3",
        "hypre_recompute_residual": "1",
        "hypre_use_dof_functions": "0",
        "hypre_use_interp_vectors": "0",
        "hypre_interp_vec_variant": "2",
    }


def clean_case_results() -> None:
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    patterns = [
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
    ]
    for pattern in patterns:
        for path in PROFILE_DIR.glob(pattern):
            if path.is_file():
                path.unlink()


def copy_artifacts(run_dir: Path) -> None:
    for name in [
        "linear_profile_events.tsv",
        "linear_profile_iterations.tsv",
        "runtime_profile_events.tsv",
        "runtime_profile_summary.tsv",
        "solver.out",
        "LD.out",
        "version.txt",
    ]:
        src = PROFILE_DIR / name
        if src.exists():
            shutil.copy2(src, run_dir / name)
    for src in sorted(PROFILE_DIR.glob("state_*.vtu")):
        shutil.copy2(src, run_dir / src.name)


def analyze_linear_profile(run_dir: Path, spec: RunSpec, threads: int) -> None:
    if not (run_dir / "linear_profile_events.tsv").exists():
        return
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/analyze-linear-profile.py"),
            "--profile-dir",
            str(run_dir),
            "--out-dir",
            str(run_dir),
            "--title",
            f"TS-N_65 hypre BoomerAMG-CG tol {spec.tolerance} Ten-Step Profile",
            "--case",
            "TS-N_65",
            "--threads",
            str(threads),
        ],
        cwd=str(ROOT),
        check=False,
    )


def render_step(run_dir: Path, step: int, title: str) -> None:
    vtu = run_dir / f"state_{step:05d}.vtu"
    if not vtu.exists():
        return
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/render-vtu-state.py"),
            str(vtu),
            str(run_dir / f"step{step}.png"),
            "--title",
            title,
        ],
        cwd=str(ROOT),
        check=False,
    )


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


def top_runtime_hotspots(run_dir: Path, limit: int = 6) -> list[dict[str, Any]]:
    rows = read_tsv(run_dir / "runtime_profile_summary.tsv")
    parsed = []
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


def summarize_run(run_dir: Path, spec: RunSpec, returncode: int, wall_seconds: float, target_steps: int) -> dict[str, Any]:
    solver_out = read_solver_out(run_dir / "solver.out")
    events = read_tsv(run_dir / "linear_profile_events.tsv")
    solve = [row for row in events if row.get("phase") == "solve"]
    factorize = [row for row in events if row.get("phase") == "factorize"]
    analyze = [row for row in events if row.get("phase") == "analyze"]
    steps = max([as_int(row.get("step"), 0) for row in solver_out], default=0)
    nonlinear_iterations = sum(as_int(row.get("iterations"), 0) for row in solver_out)
    iterations = [as_int(row.get("solver_iterations"), -1) for row in solve if as_int(row.get("solver_iterations"), -1) >= 0]
    residuals = [as_float(row.get("solver_error")) for row in solve if as_float(row.get("solver_error")) >= 0]
    solve_success = sum(as_int(row.get("success"), 0) for row in solve)
    solve_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in solve)
    setup_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in factorize)
    analyze_seconds = sum(as_float(row.get("duration_seconds"), 0.0) for row in analyze)
    wall_solver = max([as_float(row.get("elapsed_time"), 0.0) for row in solver_out], default=0.0)
    wall = max(wall_seconds, wall_solver)
    linear_seconds = analyze_seconds + setup_seconds + solve_seconds
    status = f"completed_{target_steps}_steps" if returncode == 0 and steps >= target_steps else ("partial" if steps > 0 else "failed")
    if solve and solve_success != len(solve):
        status = "linear_failed"
    return {
        "name": spec.name,
        "solver": "HypreBoomerAMGCG",
        "tolerance": spec.tolerance,
        "status": status,
        "returncode": returncode,
        "completed_steps": steps,
        "nonlinear_iterations": nonlinear_iterations,
        "linear_solves": len(solve),
        "factorizations": len(factorize),
        "outer_iterations": sum(iterations),
        "krylov_iter_median": median(iterations) if iterations else None,
        "krylov_iter_max": max(iterations, default=None),
        "true_relres_median": median(residuals) if residuals else None,
        "true_relres_max": max(residuals, default=None),
        "wall_seconds": wall,
        "analyze_seconds": analyze_seconds,
        "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds,
        "linear_seconds": linear_seconds,
        "other_seconds": max(wall - linear_seconds, 0.0),
        "linear_share": linear_seconds / wall if wall > 0 else None,
        "vtu_steps": [path.name for path in sorted(run_dir.glob("state_*.vtu"))],
        "runtime_hotspots": top_runtime_hotspots(run_dir),
    }


def refresh_manifest(run_dir: Path, target_steps: int) -> dict[str, Any]:
    prior_path = run_dir / "run.json"
    prior = json.loads(prior_path.read_text(encoding="utf-8")) if prior_path.exists() else {}
    spec = RunSpec(prior.get("name", run_dir.name), prior.get("tolerance", ""))
    manifest = summarize_run(
        run_dir,
        spec,
        as_int(str(prior.get("returncode", 999)), 999),
        as_float(str(prior.get("wall_seconds", 0.0)), 0.0),
        target_steps,
    )
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def run_spec(spec: RunSpec, out_dir: Path, oas_bin: Path, threads: int, target_steps: int, save_steps: list[int]) -> dict[str, Any]:
    run_dir = out_dir / spec.name
    run_dir.mkdir(parents=True, exist_ok=True)
    solver_inp = CASE_DIR / "solver.inp"
    exporters_inp = CASE_DIR / "exporters.inp"
    solver_original = solver_inp.read_text(encoding="utf-8", errors="replace")
    exporters_original = exporters_inp.read_text(encoding="utf-8", errors="replace")
    (run_dir / "solver.inp.original").write_text(solver_original, encoding="utf-8")
    (run_dir / "exporters.inp.original").write_text(exporters_original, encoding="utf-8")
    total_time = f"{target_steps * 5.0e-3:.6e}"
    start = time.monotonic()
    returncode = 999
    try:
        patch_key_value_file(solver_inp, solver_updates(spec.tolerance, total_time), "Local TS-N_65 ten-step hypre validation settings.")
        patch_exporters(exporters_inp, save_steps)
        (run_dir / "solver.inp.effective").write_text(solver_inp.read_text(encoding="utf-8"), encoding="utf-8")
        (run_dir / "exporters.inp.effective").write_text(exporters_inp.read_text(encoding="utf-8"), encoding="utf-8")
        clean_case_results()
        env = os.environ.copy()
        env.setdefault("MKLROOT", "/opt/intel/oneapi/mkl/latest")
        env["MKL_NUM_THREADS"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
        cmd = [str(oas_bin), "-j", str(threads), str(MASTER_FILE)]
        print(f"[{datetime.now().isoformat(timespec='seconds')}] starting {spec.name}: {' '.join(cmd)}", flush=True)
        with (run_dir / "oas.log").open("w", encoding="utf-8", errors="replace") as log:
            log.write(f"Starting {spec.name} at {datetime.now().isoformat(timespec='seconds')}\n")
            log.write(f"Command: {' '.join(cmd)}\nThreads: {threads}\nTarget steps: {target_steps}\nTolerance: {spec.tolerance}\n")
            log.flush()
            process = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT, text=True)
            while True:
                try:
                    returncode = process.wait(timeout=60)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - start
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {spec.name} still running ({elapsed:.0f}s)", flush=True)
                    log.write(f"\n[runner] still running at {elapsed:.0f}s\n")
                    log.flush()
            log.write(f"\nOAS exit status: {returncode}\nFinished at {datetime.now().isoformat(timespec='seconds')}\n")
        print(f"[{datetime.now().isoformat(timespec='seconds')}] finished {spec.name} with status {returncode}", flush=True)
    finally:
        solver_inp.write_text(solver_original, encoding="utf-8")
        exporters_inp.write_text(exporters_original, encoding="utf-8")

    elapsed = time.monotonic() - start
    copy_artifacts(run_dir)
    analyze_linear_profile(run_dir, spec, threads)
    render_step(run_dir, target_steps, f"{spec.name} step {target_steps}")
    manifest = summarize_run(run_dir, spec, returncode, elapsed, target_steps)
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def add_comparisons(out_dir: Path, manifests: list[dict[str, Any]], reference_dir: Path, save_steps: list[int]) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for manifest in manifests:
        run_dir = out_dir / manifest["name"]
        for step in [step for step in save_steps if step <= 2]:
            cmp = compare_vtu(run_dir / f"state_{step:05d}.vtu", reference_dir / f"state_{step:05d}.vtu")
            rows.append([manifest["name"], f"vs Pardiso", step, cmp.get("available"), cmp.get("relative_l2"), cmp.get("rms_nodal"), cmp.get("max_nodal"), cmp.get("max_component")])
    if len(manifests) >= 2:
        by_name = {m["name"]: m for m in manifests}
        loose = by_name.get("hypre-boomeramg-tol1em1")
        strict = by_name.get("hypre-boomeramg-tol1em3")
        if loose and strict:
            for step in save_steps:
                cmp = compare_vtu(out_dir / loose["name"] / f"state_{step:05d}.vtu", out_dir / strict["name"] / f"state_{step:05d}.vtu")
                rows.append([loose["name"], "vs hypre tol1e-3", step, cmp.get("available"), cmp.get("relative_l2"), cmp.get("rms_nodal"), cmp.get("max_nodal"), cmp.get("max_component")])
    return rows


def write_summary(out_dir: Path, manifests: list[dict[str, Any]], reference_dir: Path, save_steps: list[int], target_steps: int) -> None:
    runtime_rows = []
    convergence_rows = []
    hotspot_rows = []
    for m in manifests:
        runtime_rows.append([m["name"], m["tolerance"], m["status"], m["completed_steps"], m["wall_seconds"], m["setup_seconds"], m["solve_seconds"], m["other_seconds"], m["linear_share"]])
        convergence_rows.append([m["name"], m["linear_solves"], m["factorizations"], m["nonlinear_iterations"], m["outer_iterations"], m["krylov_iter_median"], m["krylov_iter_max"], m["true_relres_median"], m["true_relres_max"]])
        for h in m.get("runtime_hotspots", [])[:4]:
            hotspot_rows.append([m["name"], h["phase"], h["detail"], h["count"], h["total_seconds"], h["mean_seconds"], h["max_seconds"]])
    comparison_rows = add_comparisons(out_dir, manifests, reference_dir, save_steps) if reference_dir.exists() else []
    with (out_dir / "summary.md").open("w", encoding="utf-8") as handle:
        handle.write("# TS-N_65 Hypre-CG Ten-Step Tolerance Validation\n\n")
        handle.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n\n")
        handle.write(f"Target: first `{target_steps}` load steps, 16 shared-memory threads, VTK snapshots for steps `{', '.join(map(str, save_steps))}`.\n\n")
        handle.write("## Runtime\n\n")
        handle.write(md_table(["run", "tol", "status", "steps", "wall s", "setup s", "solve s", "other s", "linear share"], runtime_rows) + "\n\n")
        handle.write("## Convergence\n\n")
        handle.write(md_table(["run", "linear solves", "factorizations", "nonlinear iters", "Krylov total", "Krylov median", "Krylov max", "true relres median", "true relres max"], convergence_rows) + "\n\n")
        if comparison_rows:
            handle.write("## VTU Displacement Closeness\n\n")
            handle.write(md_table(["run", "reference", "step", "available", "relative L2", "RMS nodal", "max nodal", "max component"], comparison_rows) + "\n\n")
        if hotspot_rows:
            handle.write("## Runtime Hotspots\n\n")
            handle.write(md_table(["run", "phase", "detail", "count", "total s", "mean s", "max s"], hotspot_rows) + "\n\n")
        handle.write("## Artifacts\n\n")
        handle.write(md_table(["run", "folder"], [[m["name"], f"./{m['name']}"] for m in manifests]) + "\n")
    print(out_dir / "summary.md", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oas-bin", type=Path, default=DEFAULT_OAS_BIN)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--target-steps", type=int, default=10)
    parser.add_argument("--save-steps", default="1,2,3,4,5,6,7,8,9,10")
    parser.add_argument("--tolerances", default="1e-1,1e-3")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--report-only", type=Path, default=None)
    args = parser.parse_args()

    save_steps = [int(item) for item in args.save_steps.split(",") if item.strip()]
    specs = [RunSpec(f"hypre-boomeramg-tol{tol.replace('-', 'm')}", tol) for tol in args.tolerances.split(",") if tol.strip()]
    if args.only:
        specs = [spec for spec in specs if any(token in spec.name for token in args.only)]

    if args.report_only:
        manifests = [refresh_manifest(path.parent, args.target_steps) for path in sorted(args.report_only.glob("*/run.json"))]
        write_summary(args.report_only, manifests, args.reference_dir, save_steps, args.target_steps)
        return 0

    if not args.oas_bin.exists() or not os.access(args.oas_bin, os.X_OK):
        raise SystemExit(f"OAS binary is missing or not executable: {args.oas_bin}")
    if not MASTER_FILE.exists():
        raise SystemExit(f"TS-N_65 master input is missing: {MASTER_FILE}")

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = args.out_dir or ROOT / f"results/tsn65-hypre-10step-validation-{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifests: list[dict[str, Any]] = []
    for run_json in sorted(out_dir.glob("*/run.json")):
        manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
    seen = {m["name"] for m in manifests}

    for spec in specs:
        if spec.name in seen:
            continue
        manifest = run_spec(spec, out_dir, args.oas_bin, args.threads, args.target_steps, save_steps)
        manifests.append(manifest)
        write_summary(out_dir, manifests, args.reference_dir, save_steps, args.target_steps)

    write_summary(out_dir, manifests, args.reference_dir, save_steps, args.target_steps)
    (out_dir / "campaign.json").write_text(json.dumps(manifests, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if all(m["status"] == f"completed_{args.target_steps}_steps" for m in manifests) else 1


if __name__ == "__main__":
    raise SystemExit(main())
