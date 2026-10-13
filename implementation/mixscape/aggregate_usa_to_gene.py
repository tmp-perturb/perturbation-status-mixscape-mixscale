#!/usr/bin/env python3
"""Aggregate a splice-aware USA count matrix to gene-level S+A counts."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from adapt_mixscape_input import (
    as_bool,
    choose_cells,
    read_h5_array,
    read_obs_column,
    read_index,
)


def prepare_metadata(
    handle: h5py.File,
    row_indices: np.ndarray,
    obs_names: list[str],
    selection: dict[str, object],
):
    target_label = read_obs_column(handle["obs"]["target_label"])
    is_ntc = [as_bool(value) for value in read_obs_column(handle["obs"]["is_ntc"])]
    assignment_structure = read_obs_column(handle["obs"]["assignment_structure"])
    selected = choose_cells(target_label, is_ntc, assignment_structure, selection)
    if not np.array_equal(selected, row_indices):
        raise ValueError("Selected cells do not match the standard Mixscape selection")

    required_metadata = {"guide_id", "target_label", "is_ntc", "batch_id"}
    metadata_columns = [
        name
        for name in handle["obs"].keys()
        if name in required_metadata
        or not str(name).startswith("guide_")
    ]
    obs = {}
    for name in metadata_columns:
        values = read_obs_column(handle["obs"][name])
        obs[name] = [values[index] for index in row_indices]
    metadata = pd.DataFrame(obs, index=[obs_names[i] for i in row_indices])
    metadata.index.name = "cell_id"
    metadata["gene"] = metadata["target_label"]
    return metadata, metadata_columns


def aggregate_contiguous_block(
    matrix: h5py.Group,
    indptr: np.ndarray,
    source_start: int,
    source_end: int,
    selected_rows: np.ndarray,
    n_genes: int,
) -> sp.csr_matrix:
    """Read one contiguous source-cell block and return selected S+A rows."""
    data_start = int(indptr[source_start])
    data_end = int(indptr[source_end])
    raw_indptr = indptr[source_start : source_end + 1] - data_start
    raw_indices = np.asarray(matrix["indices"][data_start:data_end], dtype=np.int64)
    raw_data = np.asarray(matrix["data"][data_start:data_end])
    if not np.all(np.equal(raw_data, np.rint(raw_data))):
        raise ValueError("The counts layer contains non-integer values")
    raw = sp.csr_matrix(
        (raw_data, raw_indices, raw_indptr),
        shape=(source_end - source_start, 3 * n_genes),
    )
    selected = raw[selected_rows - source_start]
    aggregated = selected[:, :n_genes] + selected[:, 2 * n_genes :]
    aggregated = aggregated.tocsr()
    aggregated.sum_duplicates()
    aggregated.eliminate_zeros()
    return aggregated


def count_aggregated_nnz(
    source_h5ad: Path,
    row_indices: np.ndarray,
    n_genes: int,
    row_chunk_size: int,
) -> int:
    total = 0
    with h5py.File(source_h5ad, "r") as handle:
        matrix = handle["layers/counts"]
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        n_source_cells = int(indptr.shape[0] - 1)
        n_chunks = (n_source_cells + row_chunk_size - 1) // row_chunk_size
        for chunk_number, source_start in enumerate(
            range(0, n_source_cells, row_chunk_size), start=1
        ):
            source_end = min(source_start + row_chunk_size, n_source_cells)
            selected_rows = row_indices[
                (row_indices >= source_start) & (row_indices < source_end)
            ]
            if selected_rows.size == 0:
                continue
            block = aggregate_contiguous_block(
                matrix, indptr, source_start, source_end, selected_rows, n_genes
            )
            total += int(np.diff(block.indptr).sum())
            if chunk_number == 1 or chunk_number % 256 == 0 or chunk_number == n_chunks:
                print(
                    f"count pass: {chunk_number}/{n_chunks} chunks, nnz={total}",
                    flush=True,
                )
    return total


def write_aggregated_matrix(
    source_h5ad: Path,
    output_path: Path,
    row_indices: np.ndarray,
    n_genes: int,
    row_chunk_size: int,
    expected_nnz: int,
) -> int:
    written_nnz = 0
    with h5py.File(source_h5ad, "r") as handle:
        matrix = handle["layers/counts"]
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        with gzip.open(output_path, "wt", encoding="utf-8", compresslevel=1) as stream:
            stream.write("%%MatrixMarket matrix coordinate integer general\n")
            stream.write("% Aggregated from USA counts using gene-level S+A.\n")
            stream.write(f"{n_genes} {len(row_indices)} {expected_nnz}\n")
            n_source_cells = int(indptr.shape[0] - 1)
            n_chunks = (n_source_cells + row_chunk_size - 1) // row_chunk_size
            written_cells = 0
            for chunk_number, source_start in enumerate(
                range(0, n_source_cells, row_chunk_size), start=1
            ):
                source_end = min(source_start + row_chunk_size, n_source_cells)
                selected_rows = row_indices[
                    (row_indices >= source_start) & (row_indices < source_end)
                ]
                if selected_rows.size == 0:
                    continue
                block = aggregate_contiguous_block(
                    matrix, indptr, source_start, source_end, selected_rows, n_genes
                )
                block_nnz = int(block.nnz)
                if block_nnz:
                    row_indices_1based = block.indices.astype(np.int64, copy=False) + 1
                    col_indices_1based = np.repeat(
                        np.arange(written_cells, written_cells + block.shape[0], dtype=np.int64)
                        + 1,
                        np.diff(block.indptr),
                    )
                    values = block.data.astype(np.int64, copy=False)
                    entries = np.column_stack(
                        (row_indices_1based, col_indices_1based, values)
                    )
                    np.savetxt(stream, entries, fmt="%d %d %d")
                    written_nnz += block_nnz
                written_cells += block.shape[0]
                if chunk_number == 1 or chunk_number % 256 == 0 or chunk_number == n_chunks:
                    print(
                        f"write pass: {chunk_number}/{n_chunks} chunks, nnz={written_nnz}",
                        flush=True,
                    )
    if written_nnz != expected_nnz:
        raise RuntimeError(f"Aggregated nnz mismatch: {written_nnz} != {expected_nnz}")
    return written_nnz


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source_h5ad = Path(config["source_h5ad"])
    output_dir = Path(config["input_dir"])
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    row_chunk_size = int(config.get("row_chunk_size", 256))
    selection = config.get(
        "selection",
        {
            "mode": "eligible",
            "eligible_assignment_structures": [
                "single_guide",
                "concordant_construct",
            ],
        },
    )

    with h5py.File(source_h5ad, "r") as handle:
        n_cells, n_features = map(int, handle["X"].attrs["shape"])
        obs_names = read_index(handle["obs"], "_index")
        feature_names = read_index(handle["var"], "_index")
        target_label = read_obs_column(handle["obs"]["target_label"])
        is_ntc = [as_bool(value) for value in read_obs_column(handle["obs"]["is_ntc"])]
        assignment_structure = read_obs_column(handle["obs"]["assignment_structure"])
        row_indices = choose_cells(
            target_label,
            is_ntc,
            assignment_structure,
            selection,
        )
        if len(feature_names) % 3 != 0:
            raise ValueError("USA feature count is not divisible by three")
        n_genes = len(feature_names) // 3
        suffixes = [name.rsplit("_", 1)[-1] for name in feature_names]
        if suffixes != ["S"] * n_genes + ["U"] * n_genes + ["A"] * n_genes:
            raise ValueError("Feature order is not the expected S/U/A block order")
        base_names = [name.rsplit("_", 1)[0] for name in feature_names]
        if not (
            base_names[:n_genes]
            == base_names[n_genes : 2 * n_genes]
            == base_names[2 * n_genes :]
        ):
            raise ValueError("S/U/A blocks do not share the same gene order")
        if tuple(handle["layers/counts"].attrs["shape"]) != (n_cells, n_features):
            raise ValueError("Counts layer shape does not match X")
        metadata, metadata_columns = prepare_metadata(
            handle, row_indices, obs_names, selection
        )

    expected_nnz = count_aggregated_nnz(
        source_h5ad, row_indices, n_genes, row_chunk_size
    )
    counts_path = output_dir / "counts.mtx.gz"
    features_path = output_dir / "features.tsv.gz"
    metadata_path = output_dir / "metadata.tsv.gz"
    manifest_path = output_dir / "adapter_manifest.json"
    partial_counts_path = output_dir / "counts.mtx.gz.partial"
    write_aggregated_matrix(
        source_h5ad,
        partial_counts_path,
        row_indices,
        n_genes,
        row_chunk_size,
        expected_nnz,
    )
    features = pd.DataFrame(
        {
            "gene_id": base_names[:n_genes],
            "gene_name": base_names[:n_genes],
            "feature_type": "Gene Expression",
        }
    )
    partial_features_path = output_dir / "features.tsv.gz.partial"
    partial_metadata_path = output_dir / "metadata.tsv.gz.partial"
    features.to_csv(
        partial_features_path,
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    metadata.to_csv(partial_metadata_path, sep="\t", compression="gzip")
    manifest = {
        "module": "M03",
        "method": "Mixscape",
        "dataset": config["dataset"],
        "source_h5ad": str(source_h5ad.resolve()),
        "source_shape": [n_cells, n_features],
        "output_shape": [n_genes, int(len(row_indices))],
        "counts_nnz": expected_nnz,
        "counts_source": "layers/counts, S+A aggregation",
        "usa_aggregation": "S+A; U excluded",
        "feature_blocks": {"S": n_genes, "U": n_genes, "A": n_genes},
        "selection": selection,
        "selected_cells": int(len(row_indices)),
        "selected_ntc_cells": int(sum(is_ntc[index] is True for index in row_indices)),
        "metadata_columns": metadata_columns + ["gene"],
        "ordering": "sorted original H5AD row order",
    }
    partial_manifest_path = output_dir / "adapter_manifest.json.partial"
    partial_manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    partial_counts_path.replace(counts_path)
    partial_features_path.replace(features_path)
    partial_metadata_path.replace(metadata_path)
    partial_manifest_path.replace(manifest_path)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
