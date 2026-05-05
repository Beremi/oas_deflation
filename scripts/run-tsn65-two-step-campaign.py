#!/usr/bin/env python3
"""Run and report the TS-N_65 two-step solver comparison campaign."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "data/cases/TS-N_65"
MASTER_FILE = CASE_DIR / "master.inp"
PROFILE_DIR = CASE_DIR / "results"
DEFAULT_OAS_BIN = ROOT.parent / "oas_deflation-build/vtk-release/bin/OAS"
VTK_EXPORT_LINE = "VTKElementExporter state saveSteps 2 1 2 ascii pointData 1 displacements"


@dataclass
class RunSpec:
    name: str
    solver: str
    title: str
    variant: str
    tolerance: str = ""
    updates: dict[str, str] = field(default_factory=dict)


def campaign_specs() -> list[RunSpec]:
    base_hypre = {
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
    amgcl_common = {
        "amgcl_max_iterations": "500",
        "amgcl_near_nullspace": "1",
        "amgcl_elastic_full_lift": "1",
        "amgcl_reuse_initial_guess": "0",
        "amgcl_eps_strong": "0",
        "amgcl_npre": "1",
        "amgcl_npost": "1",
        "amgcl_ncycle": "1",
        "amgcl_coarse_enough": "0",
    }
    amgcl_variants = [
        (
            "amgcl-dogbone-best",
            "Dogbone-best block ILU0",
            {
                **amgcl_common,
                "amgcl_use_block_backend": "1",
                "amgcl_block_size": "1",
                "amgcl_block_relaxation": "ilu0",
                "amgcl_relax": "1.2",
            },
        ),
        (
            "amgcl-scalar-lift",
            "TS-N_65 scalar-lift fallback",
            {
                **amgcl_common,
                "amgcl_use_block_backend": "0",
                "amgcl_block_size": "6",
                "amgcl_relax": "1.0",
            },
        ),
        (
            "amgcl-block-spai0",
            "TS-N_65 block SPAI0 fallback",
            {
                **amgcl_common,
                "amgcl_use_block_backend": "1",
                "amgcl_block_size": "1",
                "amgcl_block_relaxation": "spai0",
                "amgcl_relax": "1.0",
            },
        ),
    ]

    specs = [
        RunSpec(
            name="pardisoldlt-reference",
            solver="PardisoLDLT",
            title="TS-N_65 PardisoLDLT 16-thread Two-Step Reference",
            variant="reference",
        )
    ]
    for tol in ["1e-1", "1e-3", "1e-6"]:
        specs.append(
            RunSpec(
                name=f"hypre-boomeramg-tol{tol.replace('-', 'm')}",
                solver="HypreBoomerAMGCG",
                title=f"TS-N_65 hypre BoomerAMG-CG tol {tol} Two-Step Profile",
                variant="hypre tuned persistent flexible-PCG",
                tolerance=tol,
                updates={**base_hypre, "hypre_tolerance": tol},
            )
        )
    for prefix, variant, updates in amgcl_variants:
        for tol in ["1e-1", "1e-3", "1e-6"]:
            specs.append(
                RunSpec(
                    name=f"{prefix}-tol{tol.replace('-', 'm')}",
                    solver="AmgclCGElastic",
                    title=f"TS-N_65 {variant} tol {tol} Two-Step Profile",
                    variant=variant,
                    tolerance=tol,
                    updates={**updates, "amgcl_tolerance": tol},
                )
            )
    return specs


def smoke_spec() -> RunSpec:
    return RunSpec(
        name="smoke-pardisoldlt-one-iter",
        solver="PardisoLDLT",
        title="TS-N_65 PardisoLDLT VTK Smoke",
        variant="smoke",
        updates={
            "total_time": "5.000000e-03",
            "max_iterations": "1",
            "min_iterations": "1",
        },
    )


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip()
    names = header.lstrip("#").split("\t")
    return pd.read_csv(path, sep="\t", comment="#", names=names)


def md_table(rows: list[dict[str, Any]] | pd.DataFrame) -> str:
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    if df.empty:
        return "_No rows._"
    headers = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for values in df.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(format_cell(v) for v in values) + " |")
    return "\n".join(lines)


def format_cell(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    if isinstance(value, float):
        if math.isinf(value) or math.isnan(value):
            return ""
        if abs(value) >= 1000 or (0 < abs(value) < 1e-3):
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)


def fmt(value: Any) -> str:
    return format_cell(value)


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
        lines.append("# Local TS-N_65 two-step comparison VTK export.")
        lines.append(VTK_EXPORT_LINE)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return original


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


def common_solver_updates(spec: RunSpec, total_time: str = "1.000000e-02") -> dict[str, str]:
    updates = {
        "total_time": total_time,
        "solver_type": spec.solver,
        "linear_solver_profile": "1",
        "linear_solver_profile_matrix_delta": "0",
        "linear_solver_profile_file": "linear_profile",
        "runtime_phase_profile": "1",
        "runtime_phase_profile_file": "runtime_profile",
        "stiffness_matrix_iter_update": "-5",
        "stiffness_matrix_step_update": "1",
        "max_iterations": "1000",
    }
    updates.update(spec.updates)
    return updates


def run_oas(spec: RunSpec, run_dir: Path, oas_bin: Path, threads: int, total_time: str) -> dict[str, Any]:
    solver_inp = CASE_DIR / "solver.inp"
    exporters_inp = CASE_DIR / "exporters.inp"
    run_dir.mkdir(parents=True, exist_ok=True)
    solver_original = solver_inp.read_text(encoding="utf-8", errors="replace")
    exporters_original = exporters_inp.read_text(encoding="utf-8", errors="replace")
    (run_dir / "solver.inp.original").write_text(solver_original, encoding="utf-8")
    (run_dir / "exporters.inp.original").write_text(exporters_original, encoding="utf-8")

    started = datetime.now().isoformat(timespec="seconds")
    start_time = time.monotonic()
    returncode = 999
    try:
        patch_key_value_file(
            solver_inp,
            common_solver_updates(spec, total_time=total_time),
            "Local TS-N_65 two-step comparison settings.",
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
        cmd = [str(oas_bin), "-j", str(threads), str(MASTER_FILE)]
        print(f"[{datetime.now().isoformat(timespec='seconds')}] starting {spec.name}: {' '.join(cmd)}", flush=True)
        with log_path.open("w", encoding="utf-8", errors="replace") as log:
            log.write(f"Starting {spec.name} at {started}\n")
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write(f"Threads: {threads}\n")
            log.write(f"Solver: {spec.solver}\n")
            if spec.tolerance:
                log.write(f"Tolerance: {spec.tolerance}\n")
            log.flush()
            process = subprocess.Popen(
                cmd,
                cwd=str(ROOT),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
            while True:
                try:
                    returncode = process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - start_time
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {spec.name} still running ({elapsed:.0f}s)", flush=True)
                    log.write(f"\n[runner] still running at {elapsed:.0f}s\n")
                    log.flush()
            log.write(f"\nOAS exit status: {returncode}\nFinished at {datetime.now().isoformat(timespec='seconds')}\n")
        print(f"[{datetime.now().isoformat(timespec='seconds')}] finished {spec.name} with status {returncode}", flush=True)
    finally:
        solver_inp.write_text(solver_original, encoding="utf-8")
        exporters_inp.write_text(exporters_original, encoding="utf-8")

    elapsed = time.monotonic() - start_time
    copy_artifacts(run_dir)
    analyze_linear_profile(run_dir, spec, threads)
    render_step2(run_dir, spec)
    manifest = summarize_run(run_dir, spec, returncode, elapsed, threads, started)
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


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


def analyze_linear_profile(run_dir: Path, spec: RunSpec, threads: int) -> None:
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
        spec.title,
        "--case",
        "TS-N_65",
        "--threads",
        str(threads),
    ]
    subprocess.run(cmd, cwd=str(ROOT), check=False)


def render_step2(run_dir: Path, spec: RunSpec) -> None:
    vtu = run_dir / "state_00002.vtu"
    if not vtu.exists():
        return
    cmd = [
        sys.executable,
        str(ROOT / "scripts/render-vtu-state.py"),
        str(vtu),
        str(run_dir / "step2.png"),
        "--title",
        spec.name,
    ]
    subprocess.run(cmd, cwd=str(ROOT), check=False)


def top_runtime_hotspots(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    path = run_dir / "runtime_profile_summary.tsv"
    if not path.exists():
        return []
    data = pd.read_csv(path, sep="\t")
    if data.empty:
        return []
    data = data.sort_values("total_seconds", ascending=False).head(limit)
    rows: list[dict[str, Any]] = []
    for _, row in data.iterrows():
        rows.append(
            {
                "phase": row.get("phase", ""),
                "detail": row.get("detail", ""),
                "count": int(row.get("count", 0)),
                "total_seconds": float(row.get("total_seconds", 0.0)),
                "mean_seconds": float(row.get("mean_seconds", 0.0)),
                "max_seconds": float(row.get("max_seconds", 0.0)),
            }
        )
    return rows


def summarize_run(run_dir: Path, spec: RunSpec, returncode: int, wall_seconds: float, threads: int, started: str) -> dict[str, Any]:
    solver_out = read_table(run_dir / "solver.out")
    ld_out = read_table(run_dir / "LD.out")
    events = pd.read_csv(run_dir / "linear_profile_events.tsv", sep="\t") if (run_dir / "linear_profile_events.tsv").exists() else pd.DataFrame()
    solve = events[events["phase"] == "solve"].copy() if not events.empty else pd.DataFrame()
    factorize = events[events["phase"] == "factorize"].copy() if not events.empty else pd.DataFrame()
    analyze = events[events["phase"] == "analyze"].copy() if not events.empty else pd.DataFrame()

    max_step = int(solver_out["step"].max()) if not solver_out.empty else 0
    completed = max_step >= 2
    log_text = (run_dir / "oas.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "oas.log").exists() else ""
    has_nan = "nan" in log_text.lower()
    if not solver_out.empty:
        has_nan = has_nan or solver_out.apply(lambda col: pd.to_numeric(col, errors="coerce")).isna().all(axis=None)
    any_linear_fail = False
    if not solve.empty and "success" in solve:
        any_linear_fail = bool((solve["success"] == 0).any())

    if completed:
        status = "completed_2_steps"
    elif has_nan:
        status = "nan_failed"
    elif any_linear_fail:
        status = "linear_failed"
    elif max_step > 0:
        status = "partial"
    elif returncode != 0:
        status = "nonlinear_failed"
    else:
        status = "partial"

    linear_seconds = float(events["duration_seconds"].sum()) if not events.empty else 0.0
    analyze_seconds = float(analyze["duration_seconds"].sum()) if not analyze.empty else 0.0
    factorize_seconds = float(factorize["duration_seconds"].sum()) if not factorize.empty else 0.0
    solve_seconds = float(solve["duration_seconds"].sum()) if not solve.empty else 0.0
    elapsed_wall = float(solver_out["elapsed_time"].iloc[-1]) if not solver_out.empty else 0.0
    elapsed_wall = max(elapsed_wall, wall_seconds)
    other_seconds = max(elapsed_wall - linear_seconds, 0.0)
    solver_iterations = solve[solve["solver_iterations"] >= 0]["solver_iterations"] if not solve.empty and "solver_iterations" in solve else pd.Series(dtype=float)
    solver_error = solve[solve["solver_error"] >= 0]["solver_error"] if not solve.empty and "solver_error" in solve else pd.Series(dtype=float)

    return {
        "name": spec.name,
        "solver": spec.solver,
        "variant": spec.variant,
        "tolerance": spec.tolerance,
        "status": status,
        "returncode": returncode,
        "threads": threads,
        "started": started,
        "wall_seconds_runner": wall_seconds,
        "wall_seconds_solver": elapsed_wall,
        "completed_steps": max_step,
        "nonlinear_iterations": int(solver_out["iterations"].sum()) if not solver_out.empty else 0,
        "linear_events": int(len(events)),
        "linear_solves": int(len(solve)),
        "factorizations": int(len(factorize)),
        "analyze_seconds": analyze_seconds,
        "factorize_seconds": factorize_seconds,
        "solve_seconds": solve_seconds,
        "linear_seconds": linear_seconds,
        "other_seconds": other_seconds,
        "linear_share": linear_seconds / elapsed_wall if elapsed_wall > 0 else None,
        "krylov_iter_median": float(solver_iterations.median()) if not solver_iterations.empty else None,
        "krylov_iter_max": float(solver_iterations.max()) if not solver_iterations.empty else None,
        "true_relres_median": float(solver_error.median()) if not solver_error.empty else None,
        "true_relres_max": float(solver_error.max()) if not solver_error.empty else None,
        "solve_success_count": int(solve["success"].sum()) if not solve.empty and "success" in solve else None,
        "ld_rows": int(len(ld_out)),
        "has_step1_vtu": (run_dir / "state_00001.vtu").exists(),
        "has_step2_vtu": (run_dir / "state_00002.vtu").exists(),
        "has_step2_png": (run_dir / "step2.png").exists(),
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


def compare_table(candidate: Path, reference: Path, step: int) -> dict[str, Any]:
    cand = read_table(candidate)
    ref = read_table(reference)
    if cand.empty or ref.empty or "step" not in cand or "step" not in ref:
        return {"available": False}
    cand_row = cand[pd.to_numeric(cand["step"], errors="coerce") == step]
    ref_row = ref[pd.to_numeric(ref["step"], errors="coerce") == step]
    if cand_row.empty or ref_row.empty:
        return {"available": False}
    c = cand_row.iloc[-1]
    r = ref_row.iloc[-1]
    cols = [col for col in cand.columns if col in ref.columns and col not in {"step", "time"}]
    diffs: list[float] = []
    refs: list[float] = []
    max_col = ""
    max_abs = -1.0
    for col in cols:
        cv = pd.to_numeric(pd.Series([c[col]]), errors="coerce").iloc[0]
        rv = pd.to_numeric(pd.Series([r[col]]), errors="coerce").iloc[0]
        if pd.isna(cv) or pd.isna(rv):
            continue
        diff = float(cv - rv)
        diffs.append(diff)
        refs.append(float(rv))
        if abs(diff) > max_abs:
            max_abs = abs(diff)
            max_col = col
    if not diffs:
        return {"available": False}
    diff_vec = np.asarray(diffs)
    ref_vec = np.asarray(refs)
    denom = float(np.linalg.norm(ref_vec))
    return {
        "available": True,
        "relative_l2": float(np.linalg.norm(diff_vec) / denom) if denom > 0 else None,
        "rms": float(np.sqrt(np.mean(diff_vec * diff_vec))),
        "max_abs": float(max_abs),
        "max_column": max_col,
    }


def write_gallery(out_dir: Path, manifests: list[dict[str, Any]]) -> str | None:
    images: list[tuple[str, Image.Image]] = []
    for manifest in manifests:
        png = out_dir / manifest["name"] / "step2.png"
        if png.exists():
            im = Image.open(png).convert("RGB")
            im.thumbnail((420, 320))
            images.append((manifest["name"], im.copy()))
    if not images:
        return None
    width = 840
    cell_h = 370
    rows = math.ceil(len(images) / 2)
    gallery = Image.new("RGB", (width, rows * cell_h), "white")
    draw = ImageDraw.Draw(gallery)
    for idx, (name, im) in enumerate(images):
        x = (idx % 2) * 420
        y = (idx // 2) * cell_h
        draw.text((x + 10, y + 8), name, fill=(20, 20, 20))
        gallery.paste(im, (x + 10, y + 35))
    name = "gallery-step2.png"
    gallery.save(out_dir / name)
    return name


def write_summary(out_dir: Path, manifests: list[dict[str, Any]]) -> None:
    if not manifests:
        (out_dir / "summary.md").write_text("# TS-N_65 Two-Step Solver Comparison\n\nNo runs executed.\n", encoding="utf-8")
        return
    ref = next((m for m in manifests if m["name"] == "pardisoldlt-reference"), None)
    ref_dir = out_dir / "pardisoldlt-reference" if ref else None

    runtime_rows = []
    convergence_rows = []
    hotspot_rows = []
    vtu_rows = []
    gauge_rows = []
    for manifest in manifests:
        runtime_rows.append(
            {
                "run": manifest["name"],
                "solver": manifest["solver"],
                "tol": manifest["tolerance"],
                "status": manifest["status"],
                "steps": manifest["completed_steps"],
                "wall_s": manifest["wall_seconds_solver"],
                "analyze_s": manifest["analyze_seconds"],
                "setup/factor_s": manifest["factorize_seconds"],
                "solve_s": manifest["solve_seconds"],
                "other_s": manifest["other_seconds"],
                "linear_share": manifest["linear_share"],
            }
        )
        convergence_rows.append(
            {
                "run": manifest["name"],
                "solves": manifest["linear_solves"],
                "factors": manifest["factorizations"],
                "nonlinear_iters": manifest["nonlinear_iterations"],
                "krylov_med": manifest["krylov_iter_median"],
                "krylov_max": manifest["krylov_iter_max"],
                "relres_med": manifest["true_relres_median"],
                "relres_max": manifest["true_relres_max"],
                "solve_success": manifest["solve_success_count"],
            }
        )
        for hotspot in manifest.get("runtime_hotspots", [])[:5]:
            hotspot_rows.append(
                {
                    "run": manifest["name"],
                    "phase": hotspot["phase"],
                    "detail": hotspot["detail"],
                    "count": hotspot["count"],
                    "total_s": hotspot["total_seconds"],
                    "mean_s": hotspot["mean_seconds"],
                    "max_s": hotspot["max_seconds"],
                }
            )
        if ref_dir and manifest["name"] != "pardisoldlt-reference":
            run_dir = out_dir / manifest["name"]
            for step in [1, 2]:
                vtu_cmp = compare_vtu(run_dir / f"state_{step:05d}.vtu", ref_dir / f"state_{step:05d}.vtu")
                vtu_rows.append(
                    {
                        "run": manifest["name"],
                        "step": step,
                        "available": vtu_cmp.get("available"),
                        "rel_l2_u": vtu_cmp.get("relative_l2"),
                        "rms_nodal": vtu_cmp.get("rms_nodal"),
                        "max_nodal": vtu_cmp.get("max_nodal"),
                        "max_component": vtu_cmp.get("max_component"),
                    }
                )
                for fname, label in [("LD.out", "LD"), ("solver.out", "solver")]:
                    gauge_cmp = compare_table(run_dir / fname, ref_dir / fname, step)
                    gauge_rows.append(
                        {
                            "run": manifest["name"],
                            "file": label,
                            "step": step,
                            "available": gauge_cmp.get("available"),
                            "rel_l2": gauge_cmp.get("relative_l2"),
                            "rms": gauge_cmp.get("rms"),
                            "max_abs": gauge_cmp.get("max_abs"),
                            "max_column": gauge_cmp.get("max_column"),
                        }
                    )

    gallery = write_gallery(out_dir, manifests)
    run_matrix = pd.DataFrame(
        [
            {
                "run": m["name"],
                "solver": m["solver"],
                "variant": m["variant"],
                "tol": m["tolerance"],
                "status": m["status"],
                "returncode": m["returncode"],
                "step2_vtu": m["has_step2_vtu"],
                "step2_png": m["has_step2_png"],
            }
            for m in manifests
        ]
    )
    reference = next((m for m in manifests if m["name"] == "pardisoldlt-reference"), None)
    completed_candidates = [m for m in manifests if m["name"] != "pardisoldlt-reference" and m["status"] == "completed_2_steps"]
    partial_candidates = [m for m in manifests if m["status"] == "partial"]
    failed_candidates = [m for m in manifests if m["status"] not in {"completed_2_steps", "partial"} and m["name"] != "pardisoldlt-reference"]

    report = out_dir / "summary.md"
    with report.open("w", encoding="utf-8") as handle:
        handle.write("# TS-N_65 Two-Step Solver Comparison\n\n")
        handle.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n\n")
        handle.write("This local campaign compares PardisoLDLT, hypre BoomerAMG-CG, and AMGCL elastic variants on the first two TS-N_65 load steps with 16 shared-memory threads. Matrix-delta profiling is disabled for this large case.\n\n")
        handle.write("## Findings\n\n")
        if reference:
            handle.write(
                f"- PardisoLDLT is the reference and reached {reference['completed_steps']} completed steps in "
                f"{fmt(reference['wall_seconds_solver'])} s wall time, with {fmt(reference['solve_seconds'])} s in linear solves "
                f"and {fmt(reference['factorize_seconds'])} s in factorization/setup.\n"
            )
        for candidate in completed_candidates:
            ratio = None
            if reference and reference["wall_seconds_solver"]:
                ratio = candidate["wall_seconds_solver"] / reference["wall_seconds_solver"]
            ratio_text = f", {fmt(ratio)}x Pardiso wall time" if ratio else ""
            handle.write(
                f"- {candidate['name']} completed both steps{ratio_text}; linear solves took "
                f"{fmt(candidate['solve_seconds'])} s with Krylov median/max "
                f"{fmt(candidate['krylov_iter_median'])}/{fmt(candidate['krylov_iter_max'])}.\n"
            )
        if partial_candidates:
            names = ", ".join(m["name"] for m in partial_candidates)
            handle.write(
                f"- Partial runs preserved useful step-1 artifacts but do not provide step-2 state comparison: {names}.\n"
            )
        if failed_candidates:
            names = ", ".join(f"{m['name']} ({m['status']})" for m in failed_candidates)
            handle.write(
                f"- Failed AMG/iterative variants were stopped after diagnostics showed NaNs, inner-solve failure, or nonlinear stagnation: {names}.\n"
            )
        handle.write("- A return code of `-15` means the OAS child process was terminated intentionally after the run had become diagnostic rather than competitive; copied artifacts and the effective input files remain in that run directory.\n\n")
        handle.write("## Run Matrix\n\n")
        handle.write(md_table(run_matrix) + "\n\n")
        handle.write("## Runtime Breakdown\n\n")
        handle.write(md_table(pd.DataFrame(runtime_rows)) + "\n\n")
        if hotspot_rows:
            handle.write("## Runtime Phase Hotspots\n\n")
            handle.write(md_table(pd.DataFrame(hotspot_rows)) + "\n\n")
        handle.write("## Krylov and Nonlinear Convergence\n\n")
        handle.write(md_table(pd.DataFrame(convergence_rows)) + "\n\n")
        handle.write("## VTU Displacement Closeness vs Pardiso\n\n")
        handle.write(md_table(pd.DataFrame(vtu_rows)) + "\n\n")
        handle.write("## Gauge Closeness vs Pardiso\n\n")
        handle.write(md_table(pd.DataFrame(gauge_rows)) + "\n\n")
        handle.write("## Step-2 Visuals\n\n")
        if gallery:
            handle.write(f"![Step-2 gallery]({gallery})\n\n")
        for manifest in manifests:
            if manifest["has_step2_png"]:
                handle.write(f"- [{manifest['name']}](./{manifest['name']}/step2.png)\n")
        handle.write("\n")
        handle.write("## Per-Run Linear Profile Reports\n\n")
        for manifest in manifests:
            if (out_dir / manifest["name"] / "linear-profile.md").exists():
                handle.write(f"- [{manifest['name']}](./{manifest['name']}/linear-profile.md)\n")
    print(report)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oas-bin", type=Path, default=DEFAULT_OAS_BIN)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--only", action="append", default=[], help="Run only specs whose names contain this substring. Repeatable.")
    parser.add_argument("--smoke", action="store_true", help="Run only the one-iteration Pardiso VTK smoke.")
    parser.add_argument("--report-only", type=Path, default=None, help="Regenerate summary from an existing campaign directory.")
    args = parser.parse_args()

    if args.report_only:
        manifests = []
        for run_json in sorted(args.report_only.glob("*/run.json")):
            manifests.append(json.loads(run_json.read_text(encoding="utf-8")))
        write_summary(args.report_only, manifests)
        return 0

    if not args.oas_bin.exists() or not os.access(args.oas_bin, os.X_OK):
        raise SystemExit(f"OAS binary is missing or not executable: {args.oas_bin}")
    if not MASTER_FILE.exists():
        raise SystemExit(f"TS-N_65 master input is missing: {MASTER_FILE}")

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = args.out_dir or ROOT / f"results/tsn65-two-step-comparison-{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    specs = [smoke_spec()] if args.smoke else campaign_specs()
    if args.only:
        specs = [spec for spec in specs if any(token in spec.name for token in args.only)]
    if not specs:
        raise SystemExit("No runs selected.")

    order = {spec.name: idx for idx, spec in enumerate(campaign_specs() + [smoke_spec()])}
    manifest_by_name: dict[str, dict[str, Any]] = {}
    for run_json in sorted(out_dir.glob("*/run.json")):
        try:
            manifest = json.loads(run_json.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if "name" in manifest:
            manifest_by_name[manifest["name"]] = manifest

    def ordered_manifests() -> list[dict[str, Any]]:
        return sorted(manifest_by_name.values(), key=lambda m: order.get(m.get("name", ""), 10_000))

    for idx, spec in enumerate(specs, start=1):
        run_dir = out_dir / spec.name
        total_time = "5.000000e-03" if args.smoke else "1.000000e-02"
        print(f"Run {idx}/{len(specs)}: {spec.name}", flush=True)
        manifest = run_oas(spec, run_dir, args.oas_bin, args.threads, total_time=total_time)
        manifest_by_name[spec.name] = manifest
        write_summary(out_dir, ordered_manifests())
        if spec.name == "pardisoldlt-reference" and manifest["status"] != "completed_2_steps" and not args.only:
            print("Pardiso reference did not complete two steps; stopping before iterative comparisons.", flush=True)
            break

    manifests = ordered_manifests()
    write_summary(out_dir, manifests)
    (out_dir / "campaign.json").write_text(json.dumps(manifests, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
