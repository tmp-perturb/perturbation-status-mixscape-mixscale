#!/usr/bin/env python3
"""Write gene-level S+A features, metadata, and manifest after count conversion."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import h5py
import pandas as pd
import yaml

from adapt_mixscape_input import as_bool, choose_cells, read_obs_column, read_index
from aggregate_usa_to_gene import prepare_metadata


def read_matrix_header(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        if stream.readline().strip() != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Unexpected MatrixMarket header")
        line = stream.readline()
        while line.startswith("%"):
            line = stream.readline()
        return tuple(map(int, line.split()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source_h5ad = Path(config["source_h5ad"])
    output_dir = Path(config["input_dir"])
    counts_path = output_dir / "counts.mtx.gz"
    if not counts_path.is_file():
        raise FileNotFoundError(counts_path)
    rows, columns, nnz = read_matrix_header(counts_path)
    selection = config["selection"]

    with h5py.File(source_h5ad, "r") as handle:
        n_cells, n_features = map(int, handle["X"].attrs["shape"])
        obs_names = read_index(handle["obs"], "_index")
        feature_names = read_index(handle["var"], "_index")
        target_label = read_obs_column(handle["obs"]["target_label"])
        is_ntc = [as_bool(value) for value in read_obs_column(handle["obs"]["is_ntc"])]
        assignment_structure = read_obs_column(handle["obs"]["assignment_structure"])
        row_indices = choose_cells(target_label, is_ntc, assignment_structure, selection)
        metadata, metadata_columns = prepare_metadata(
            handle, row_indices, obs_names, selection
        )

    if len(feature_names) % 3 != 0:
        raise ValueError("USA feature count is not divisible by three")
    n_genes = len(feature_names) // 3
    base_names = [name.rsplit("_", 1)[0] for name in feature_names]
    if rows != n_genes or columns != len(row_indices):
        raise ValueError("Counts dimensions do not match source metadata")

    output_dir.mkdir(parents=True, exist_ok=True)
    features = pd.DataFrame(
        {
            "gene_id": base_names[:n_genes],
            "gene_name": base_names[:n_genes],
            "feature_type": "Gene Expression",
        }
    )
    features.to_csv(
        output_dir / "features.tsv.gz.partial",
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    metadata.to_csv(output_dir / "metadata.tsv.gz.partial", sep="\t", compression="gzip")
    manifest = {
        "module": "M03",
        "method": "Mixscape",
        "dataset": config["dataset"],
        "source_h5ad": str(source_h5ad.resolve()),
        "source_shape": [n_cells, n_features],
        "output_shape": [rows, columns],
        "counts_nnz": nnz,
        "counts_source": "M03 eligible MatrixMarket counts, S+A aggregation",
        "usa_aggregation": "S+A; U excluded",
        "feature_blocks": {"S": n_genes, "U": n_genes, "A": n_genes},
        "selection": selection,
        "selected_cells": int(len(row_indices)),
        "selected_ntc_cells": int(sum(is_ntc[index] for index in row_indices)),
        "metadata_columns": metadata_columns + ["gene"],
        "ordering": "sorted original H5AD row order",
    }
    (output_dir / "adapter_manifest.json.partial").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "features.tsv.gz.partial").replace(output_dir / "features.tsv.gz")
    (output_dir / "metadata.tsv.gz.partial").replace(output_dir / "metadata.tsv.gz")
    (output_dir / "adapter_manifest.json.partial").replace(
        output_dir / "adapter_manifest.json"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
