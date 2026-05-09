#!/usr/bin/env python3
"""Plot spy views for a TS-N_65 replay matrix in DFGMRES active ordering."""

from __future__ import annotations

import argparse
import csv
import math
from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.io import mmread


def read_summary(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["key"]] = row["value"]
    return rows


def read_metadata(path: Path, reduced_rows: int, full_rows: int, block_size: int, dimension: int):
    node_count = full_rows // block_size
    reduced_to_full = np.empty(reduced_rows, dtype=np.int64)
    coordinates = np.full((node_count, dimension), np.nan, dtype=np.float64)

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reduced_row = int(row["reduced_row"])
            full_row = int(row["full_elastic_row"])
            node = full_row // block_size
            reduced_to_full[reduced_row] = full_row
            coordinates[node, 0] = float(row["x"])
            coordinates[node, 1] = float(row["y"])
            if dimension == 3:
                coordinates[node, 2] = float(row["z"])

    return reduced_to_full, coordinates


def make_node_permutation(coordinates: np.ndarray) -> tuple[np.ndarray, int]:
    node_count, dimension = coordinates.shape
    known = np.where(np.isfinite(coordinates[:, 0]))[0]
    missing_count = node_count - known.size

    if dimension == 3:
        order = np.lexsort((known, coordinates[known, 2], coordinates[known, 1], coordinates[known, 0]))
    else:
        order = np.lexsort((known, coordinates[known, 1], coordinates[known, 0]))
    node_order = known[order]

    if missing_count:
        missing = np.where(~np.isfinite(coordinates[:, 0]))[0]
        node_order = np.concatenate([node_order, missing])

    old_to_new = np.empty(node_count, dtype=np.int64)
    old_to_new[node_order] = np.arange(node_count, dtype=np.int64)
    return old_to_new, missing_count


def make_rcm_node_permutation(row_nodes: np.ndarray, col_nodes: np.ndarray, node_count: int) -> np.ndarray:
    mask = row_nodes != col_nodes
    adjacency = coo_matrix(
        (
            np.ones(int(mask.sum()), dtype=np.int16),
            (row_nodes[mask], col_nodes[mask]),
        ),
        shape=(node_count, node_count),
    ).tocsr()
    adjacency = adjacency.maximum(adjacency.T)
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()

    degrees = np.diff(adjacency.indptr)
    candidates = np.lexsort((np.arange(node_count, dtype=np.int64), degrees))
    visited = np.zeros(node_count, dtype=bool)
    cm_order: list[int] = []
    queue: deque[int] = deque()

    for raw_start in candidates:
        start = int(raw_start)
        if visited[start]:
            continue
        visited[start] = True
        queue.append(start)
        while queue:
            node = queue.popleft()
            cm_order.append(node)

            begin = adjacency.indptr[node]
            end = adjacency.indptr[node + 1]
            neighbors = [int(neighbor) for neighbor in adjacency.indices[begin:end] if not visited[int(neighbor)]]
            neighbors.sort(key=lambda neighbor: (degrees[neighbor], neighbor))
            for neighbor in neighbors:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    queue.append(neighbor)

    node_order = np.array(cm_order[::-1], dtype=np.int64)
    old_to_new = np.empty(node_count, dtype=np.int64)
    old_to_new[node_order] = np.arange(node_count, dtype=np.int64)
    return old_to_new


def save_density_image(counts: np.ndarray, title: str, path: Path, extent_size: int) -> None:
    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True)
    image = ax.imshow(
        np.log1p(counts),
        origin="upper",
        interpolation="nearest",
        cmap="magma",
        extent=[0, extent_size, extent_size, 0],
        aspect="equal",
    )
    ax.set_title(title)
    ax.set_xlabel("column index")
    ax.set_ylabel("row index")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("log(1 + nnz per bin)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_zoom(new_rows: np.ndarray, new_cols: np.ndarray, zoom: int, out_path: Path, max_scatter: int) -> tuple[int, str]:
    mask = (new_rows < zoom) & (new_cols < zoom)
    count = int(mask.sum())
    if count <= max_scatter:
        fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)
        ax.scatter(new_cols[mask], new_rows[mask], s=0.03, c="black", linewidths=0)
        ax.set_xlim(0, zoom)
        ax.set_ylim(zoom, 0)
        ax.set_aspect("equal")
        ax.set_xlabel("column index")
        ax.set_ylabel("row index")
        ax.set_title(f"DFGMRES active reordered matrix, first {zoom} rows/cols")
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        return count, "scatter"

    bins = min(4096, zoom)
    counts = np.zeros((bins, bins), dtype=np.uint32)
    rows = new_rows[mask]
    cols = new_cols[mask]
    chunk = 5_000_000
    for start in range(0, count, chunk):
        stop = min(start + chunk, count)
        br = (rows[start:stop] * bins) // zoom
        bc = (cols[start:stop] * bins) // zoom
        flat = br * bins + bc
        counts += np.bincount(flat, minlength=bins * bins).reshape(bins, bins).astype(np.uint32)
    save_density_image(
        counts,
        f"DFGMRES active reordered matrix, first {zoom} rows/cols ({bins}x{bins} bins)",
        out_path,
        zoom,
    )
    return count, f"{bins}x{bins} density"


def save_index_density(row_counts: np.ndarray, col_counts: np.ndarray, path: Path, bins: int) -> None:
    size = row_counts.size
    edges = np.linspace(0, size, bins + 1, dtype=np.int64)
    x = 0.5 * (edges[:-1] + edges[1:])
    row_binned = np.array([row_counts[edges[i] : edges[i + 1]].sum() for i in range(bins)])
    col_binned = np.array([col_counts[edges[i] : edges[i + 1]].sum() for i in range(bins)])

    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    ax.plot(x, row_binned, label="row nnz per bin", linewidth=1.0)
    ax.plot(x, col_binned, label="column nnz per bin", linewidth=1.0)
    ax.set_xlabel("active reordered index")
    ax.set_ylabel("nnz per index bin")
    ax.set_title("DFGMRES active reordered matrix nnz density by index")
    ax.legend()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--ordering", choices=["coordinate", "rcm"], default="coordinate")
    parser.add_argument("--overview-bins", type=int, default=2048)
    parser.add_argument("--index-bins", type=int, default=4096)
    parser.add_argument("--zoom", type=int, action="append", default=[3000, 12000])
    parser.add_argument("--max-scatter", type=int, default=2_000_000)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = read_summary(args.replay_dir / "summary.tsv")
    reduced_rows = int(summary["rows"])
    full_rows = int(summary["elastic_full_rows"])
    block_size = int(summary["elastic_block_size"])
    dimension = int(summary["elastic_dimension"])
    node_count = full_rows // block_size

    reduced_to_full, coordinates = read_metadata(
        args.replay_dir / "metadata.tsv",
        reduced_rows,
        full_rows,
        block_size,
        dimension,
    )
    matrix = mmread(args.replay_dir / "matrix.mtx").tocoo()
    reduced_row_nodes = reduced_to_full[matrix.row] // block_size
    reduced_col_nodes = reduced_to_full[matrix.col] // block_size

    if args.ordering == "rcm":
        node_old_to_new = make_rcm_node_permutation(reduced_row_nodes, reduced_col_nodes, node_count)
        node_entry_counts = np.bincount(np.concatenate([reduced_row_nodes, reduced_col_nodes]), minlength=node_count)
        missing_nodes = int((node_entry_counts == 0).sum())
        ordering_description = "DFGMRES active RCM node-major permutation, mode `3`."
        prefix = "tsn65_rcm_active"
        title_prefix = "TS-N_65 DFGMRES active RCM matrix"
        known_node_label = "nodes with matrix entries"
        missing_node_label = "isolated/empty nodes in node graph"
        caveat = (
            "RCM is computed on the active node graph derived from the replay matrix pattern. "
            "Each node's local DOF order is preserved after the node permutation."
        )
    else:
        node_old_to_new, missing_nodes = make_node_permutation(coordinates)
        ordering_description = "DFGMRES active coordinate-sorted node-major permutation, mode `2`."
        prefix = "tsn65_reordered_active"
        title_prefix = "TS-N_65 DFGMRES active coordinate-reordered matrix"
        known_node_label = "nodes with coordinates in replay metadata"
        missing_node_label = "empty/prescribed nodes appended after coordinate sort"
        caveat = (
            "The replay metadata does not include coordinates for fully prescribed nodes with no active matrix entries. "
            "Those nodes are appended after the coordinate-sorted nodes; since their active rows and columns are empty, this does not change plotted nonzero entries except for a sub-percent index offset relative to an exact in-solver dump."
        )

    reduced_to_new = node_old_to_new[reduced_to_full // block_size] * block_size + reduced_to_full % block_size
    new_rows = reduced_to_new[matrix.row]
    new_cols = reduced_to_new[matrix.col]
    nnz = int(matrix.nnz)

    bins = args.overview_bins
    counts = np.zeros((bins, bins), dtype=np.uint32)
    chunk = 5_000_000
    for start in range(0, nnz, chunk):
        stop = min(start + chunk, nnz)
        br = (new_rows[start:stop] * bins) // full_rows
        bc = (new_cols[start:stop] * bins) // full_rows
        flat = br * bins + bc
        counts += np.bincount(flat, minlength=bins * bins).reshape(bins, bins).astype(np.uint32)

    overview_name = f"{prefix}_spy_overview_log_density.png"
    save_density_image(
        counts,
        f"{title_prefix} ({bins}x{bins} bins)",
        args.out_dir / overview_name,
        full_rows,
    )

    row_counts = np.bincount(new_rows, minlength=full_rows)
    col_counts = np.bincount(new_cols, minlength=full_rows)
    density_name = f"{prefix}_nnz_density_by_index.png"
    save_index_density(row_counts, col_counts, args.out_dir / density_name, args.index_bins)

    zoom_rows = []
    for zoom in args.zoom:
        zoom_name = f"{prefix}_spy_zoom_first_{zoom}.png"
        zoom_count, zoom_mode = save_zoom(new_rows, new_cols, zoom, args.out_dir / zoom_name, args.max_scatter)
        zoom_rows.append((zoom, zoom_name, zoom_count, zoom_mode))

    density = nnz / float(full_rows * full_rows)
    known_nodes = node_count - missing_nodes
    with (args.out_dir / "README.md").open("w", encoding="utf-8") as handle:
        handle.write("# TS-N_65 Reordered Active Matrix Spy Plot - 2026-05-07\n\n")
        handle.write(f"Source replay: `{args.replay_dir}`\n\n")
        handle.write(f"Ordering: {ordering_description}\n\n")
        handle.write("| metric | value |\n")
        handle.write("| --- | ---: |\n")
        handle.write(f"| reduced rows | {reduced_rows:,} |\n")
        handle.write(f"| active rows | {full_rows:,} |\n")
        handle.write(f"| active block size | {block_size} |\n")
        handle.write(f"| active nodes | {node_count:,} |\n")
        handle.write(f"| {known_node_label} | {known_nodes:,} |\n")
        handle.write(f"| {missing_node_label} | {missing_nodes:,} |\n")
        handle.write(f"| nnz | {nnz:,} |\n")
        handle.write(f"| active density | {density:.6e} |\n")
        handle.write(f"| overview bins | {bins} x {bins} |\n\n")
        handle.write(caveat + "\n\n")
        handle.write("Images:\n\n")
        handle.write(f"![overview]({overview_name})\n\n")
        for zoom, zoom_name, zoom_count, zoom_mode in zoom_rows:
            handle.write(f"![zoom {zoom}]({zoom_name})\n\n")
            handle.write(f"- Zoom {zoom}: {zoom_count:,} nonzeros rendered as `{zoom_mode}`.\n\n")
        handle.write(f"![density by index]({density_name})\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
