#!/usr/bin/env python3
"""Build the shared CalcPerturbSig checkpoint and frozen target shards."""

from __future__ import annotations

import argparse
import json
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

from common import dump_yaml, load_yaml, make_archive, r_env, rscript, run, write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--data.h5ad", required=True)
    parser.add_argument("--data.spec", required=True)
    args = parser.parse_args()
    h5ad = Path(getattr(args, "data.h5ad")).resolve()
    spec = json.loads(Path(getattr(args, "data.spec")).read_text(encoding="utf-8"))
    project = Path(spec["project_root"])
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()

    with tempfile.TemporaryDirectory(prefix=f"pse-shared-{args.name}-") as temp_text:
        temp = Path(temp_text)
        shared = temp / "shared"
        input_dir = shared / "input"
        initial_output = shared / "initial_mixscape"
        checkpoints = shared / "initial_checkpoints"
        logs = shared / "logs"

        config = load_yaml(spec["mixscape_config"])
        config.update({
            "dataset": args.name,
            "source_h5ad": str(h5ad),
            "input_dir": str(input_dir),
            "output_dir": str(initial_output),
            "logging": {"directory": str(logs), "stage_log": str(logs / "stage_resources.jsonl")},
            "checkpoints": {"enabled": True, "directory": str(checkpoints)},
        })
        config_path = temp / "mixscape.yaml"
        dump_yaml(config, config_path)
        run([sys.executable, project / "03_mixscape" / "adapt_mixscape_input.py", "--config", config_path])
        run([sys.executable, project / "03_mixscape" / "validate_mixscape_input.py", input_dir])
        run(
            [rscript(project, "mixscape"), project / "03_mixscape" / "run_mixscape.R", config_path],
            log=logs / "mixscape.log",
            env=r_env(project, "mixscape"),
        )

        pooled = load_yaml(spec["pooled_config"])
        pooled.update({
            "dataset_prefix": args.name,
            "source_h5ad": str(h5ad),
            "source_checkpoint_b": str(checkpoints / "checkpoint_B.rds"),
            "input_dir": str(input_dir),
            "shard_root": str(shared / "shards"),
        })
        pooled_path = temp / "pooled.yaml"
        dump_yaml(pooled, pooled_path)
        shard_started = time.monotonic()
        run(
            [rscript(project, "mixscape"), project / "03_mixscape" / "target_sharding" / "build_target_shards.R", pooled_path],
            log=logs / "shard_build.log",
            env=r_env(project, "mixscape"),
        )
        shard_seconds = time.monotonic() - shard_started
        (shared / "shard_build.time.txt").write_text(
            f"Elapsed (wall clock) time (h:mm:ss or m:ss): {shard_seconds:.3f}\n"
            f"Maximum resident set size (kbytes): {resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}\n",
            encoding="utf-8",
        )
        shutil.copy2(pooled_path, shared / "pooled_config.yaml")
        # The unsharded classification output and source checkpoint are only
        # construction intermediates. Frozen target shards contain the complete
        # downstream contract, so excluding duplicates keeps the Omni artifact
        # compact without changing any method input.
        shutil.rmtree(initial_output)
        shutil.rmtree(checkpoints)
        make_archive(shared, output / "shared_signature.tar.gz", "shared")

    elapsed = time.monotonic() - started
    write_json({
        "dataset": args.name,
        "stage": "shared_signature",
        "status": "PASS",
        "wall_seconds": elapsed,
        "shard_build_wall_seconds": shard_seconds,
        "peak_child_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "technical_role": "shared CalcPerturbSig signature and frozen target shards",
    }, output / "shared_manifest.json")


if __name__ == "__main__":
    main()
