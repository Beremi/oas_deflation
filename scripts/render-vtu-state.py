#!/usr/bin/env python3
"""Render an OAS VTU state as a compact displacement-magnitude PNG.

This intentionally uses VTK for file IO/surface extraction and Matplotlib for
the final image. That avoids depending on a working OpenGL display stack on
headless benchmark machines while still producing a stable visual comparison.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def normalize_vectors(vectors: np.ndarray) -> np.ndarray:
    if vectors.ndim == 1:
        vectors = vectors[:, None]
    if vectors.shape[1] < 3:
        padded = np.zeros((vectors.shape[0], 3), dtype=float)
        padded[:, : vectors.shape[1]] = vectors
        return padded
    if vectors.shape[1] > 3:
        return vectors[:, :3]
    return vectors


def read_surface_mesh(path: Path, vector_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[list[int]]]:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        raise SystemExit(f"No points found in {path}")

    vectors_vtk = grid.GetPointData().GetArray(vector_name)
    if vectors_vtk is None:
        available = [
            grid.GetPointData().GetArrayName(i)
            for i in range(grid.GetPointData().GetNumberOfArrays())
        ]
        raise SystemExit(f"Point-data array {vector_name!r} not found in {path}; available: {available}")

    points = vtk_to_numpy(grid.GetPoints().GetData()).astype(float, copy=False)
    vectors = normalize_vectors(vtk_to_numpy(vectors_vtk).astype(float, copy=False))

    # Extract the boundary first, then map original point ids to displacements.
    geom = vtkGeometryFilter()
    geom.SetInputData(grid)
    geom.PassThroughPointIdsOn()
    geom.Update()
    surface = geom.GetOutput()
    surface_points = vtk_to_numpy(surface.GetPoints().GetData()).astype(float, copy=False)
    ids_vtk = surface.GetPointData().GetArray("vtkOriginalPointIds")
    if ids_vtk is not None:
        ids = vtk_to_numpy(ids_vtk).astype(int, copy=False)
        surface_vectors = vectors[ids]
    else:
        # Fallback should only trigger on older VTK builds; sizes commonly match
        # for point-only exports, but use a conservative slice if needed.
        surface_vectors = vectors[: surface_points.shape[0]]

    faces: list[list[int]] = []
    for cell_index in range(surface.GetNumberOfCells()):
        cell = surface.GetCell(cell_index)
        ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        if len(ids) >= 3:
            faces.append(ids)
    return surface_points, surface_vectors, np.linalg.norm(surface_vectors, axis=1), faces


def choose_scale(points: np.ndarray, vectors: np.ndarray, user_scale: float | None) -> float:
    if user_scale is not None:
        return user_scale
    span = np.ptp(points, axis=0)
    diag = float(np.linalg.norm(span))
    max_u = float(np.linalg.norm(vectors, axis=1).max(initial=0.0))
    if not math.isfinite(diag) or not math.isfinite(max_u) or max_u <= 0:
        return 1.0
    return 0.08 * diag / max_u


def downsample_points(points: np.ndarray, values: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    if points.shape[0] <= max_points:
        return points, values
    stride = max(1, math.ceil(points.shape[0] / max_points))
    return points[::stride], values[::stride]


def downsample_faces(faces: list[list[int]], max_faces: int) -> list[list[int]]:
    if len(faces) <= max_faces:
        return faces
    stride = max(1, math.ceil(len(faces) / max_faces))
    return faces[::stride]


def set_equal_axes(ax, points: np.ndarray) -> None:
    ax.set_axis_off()
    ax.view_init(elev=18, azim=-55)
    ax.set_box_aspect(np.ptp(points, axis=0) + 1e-12)


def render_mesh(
    ax,
    points: np.ndarray,
    values: np.ndarray,
    faces: list[list[int]],
    max_faces: int,
    cmap_name: str = "viridis",
):
    plot_faces = downsample_faces(faces, max_faces)
    polygons = [points[face] for face in plot_faces]
    face_values = np.array([values[face].mean() for face in plot_faces], dtype=float)
    norm = colors.Normalize(vmin=float(values.min(initial=0.0)), vmax=float(values.max(initial=1.0)))
    cmap = plt.get_cmap(cmap_name)
    collection = Poly3DCollection(
        polygons,
        facecolors=cmap(norm(face_values)),
        edgecolors=(0.05, 0.06, 0.07, 0.22),
        linewidths=0.12,
        antialiased=True,
    )
    ax.add_collection3d(collection)
    return plt.cm.ScalarMappable(norm=norm, cmap=cmap)


def render_points(ax, points: np.ndarray, values: np.ndarray, max_points: int, cmap_name: str = "viridis"):
    plot_points, plot_values = downsample_points(points, values, max_points)
    return ax.scatter(
        plot_points[:, 0],
        plot_points[:, 1],
        plot_points[:, 2],
        c=plot_values,
        s=0.6,
        cmap=cmap_name,
        linewidths=0,
        alpha=0.9,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtu", type=Path)
    parser.add_argument("png", type=Path)
    parser.add_argument("--array", default="displacements")
    parser.add_argument("--scale", type=float, default=None)
    parser.add_argument("--max-points", type=int, default=100_000)
    parser.add_argument("--max-faces", type=int, default=150_000)
    parser.add_argument("--point-cloud", action="store_true", help="Force legacy point-cloud rendering.")
    parser.add_argument("--title", default="")
    args = parser.parse_args()

    points, vectors, magnitudes, faces = read_surface_mesh(args.vtu, args.array)
    scale = choose_scale(points, vectors, args.scale)
    deformed = points + scale * vectors

    fig = plt.figure(figsize=(8.5, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    if faces and not args.point_cloud:
        mappable = render_mesh(ax, deformed, magnitudes, faces, args.max_faces)
    else:
        mappable = render_points(ax, deformed, magnitudes, args.max_points)
    set_equal_axes(ax, deformed)
    title = args.title or args.vtu.stem
    ax.set_title(f"{title}\nwarp scale {scale:.3g}", fontsize=10)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.72, pad=0.02)
    cbar.set_label("|u|", rotation=90)
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.png, dpi=160)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
