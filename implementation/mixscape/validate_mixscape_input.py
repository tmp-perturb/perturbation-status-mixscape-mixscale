#!/usr/bin/env python3
"""Validate a Seurat-ready Mixscape input directory without loading counts."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def read_matrix_header(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        first = stream.readline().strip()
        if first != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Unexpected MatrixMarket header: {first}")
        line = stream.readline()
        while line.startswith("%"):
            line = stream.readline()
        rows, columns, nnz = map(int, line.split())
        return rows, columns, nnz


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    args = parser.parse_args()
    input_dir = args.input_dir
    manifest = json.loads((input_dir / "adapter_manifest.json").read_text())
    features = []
    with gzip.open(input_dir / "features.tsv.gz", "rt", encoding="utf-8") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if len(row) != 3:
                raise ValueError("Each feature row must have three columns")
            features.append(row)
    metadata_cells = 0
    metadata_names = []
    with gzip.open(input_dir / "metadata.tsv.gz", "rt", encoding="utf-8") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "cell_id":
            raise ValueError("Metadata index column must be cell_id")
        required = {"gene", "target_label", "is_ntc", "guide_id", "batch_id"}
        if not required.issubset(set(header[1:])):
            raise ValueError("Required Mixscape metadata columns are missing")
        for row in reader:
            metadata_cells += 1
            metadata_names.append(row[0])
    rows, columns, nnz = read_matrix_header(input_dir / "counts.mtx.gz")
    feature_names = [row[0] for row in features]
    if rows != len(features) or columns != metadata_cells:
        raise ValueError("Matrix dimensions do not match features/metadata")
    if metadata_names != list(dict.fromkeys(metadata_names)):
        raise ValueError("Metadata cell IDs are not unique")
    if not all(row[2] == "Gene Expression" for row in features):
        raise ValueError("Non-expression features found")
    if any(name.endswith(("_S", "_U", "_A")) for name in feature_names):
        raise ValueError("USA suffixes remain in gene-level features")
    expected_shape = tuple(manifest["output_shape"])
    if (rows, columns) != expected_shape or nnz != int(manifest["counts_nnz"]):
        raise ValueError("Matrix header does not match adapter manifest")
    print(json.dumps({
        "status": "valid",
        "shape": [rows, columns],
        "nnz": nnz,
        "file_bytes": (input_dir / "counts.mtx.gz").stat().st_size,
        "metadata_columns": header[1:],
        "gene_level_features": len(feature_names),
        "usa_aggregation": manifest.get("usa_aggregation"),
    }, indent=2))


if __name__ == "__main__":
    main()
