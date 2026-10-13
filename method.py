#!/usr/bin/env python3
"""Run Mixscape or Mixscale from the shared CalcPerturbSig cohort."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

from common import dump_yaml, extract_archive, load_yaml, make_archive, prepare_r_library, rscript, run, write_json

HERE = Path(__file__).resolve().parent
MIXSCAPE = HERE / "implementation" / "mixscape"
MIXSCALE = HERE / "implementation" / "mixscale"
VENDOR = HERE / "vendor"


def _parallel(commands: list[tuple[list[object], Path, dict[str, str]]], workers: int) -> None:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run, command, log=log, env=env) for command, log, env in commands]
        for future in futures:
            future.result()


def run_mixscape(spec: dict, shared: Path, work: Path, workers: int, r_env: dict[str, str]) -> Path:
    pooled = load_yaml(MIXSCAPE / "sharding" / f'{spec["dataset"]}_pooled.yaml')
    root = work / "mixscape"
    commands = []
    for checkpoint in sorted((shared / "shards" / "checkpoints").glob("shard_*.rds")):
        shard = checkpoint.stem
        config = {
            "module": "M03", "method": "Mixscape", "dataset": f'{spec["dataset"]}_{shard}',
            "execution_mode": "checkpoint_only", "source_h5ad": spec["h5ad"],
            "input_dir": str(shared / "input"), "output_dir": str(root / "outputs" / shard),
            "resume_checkpoint": str(checkpoint), "selection": pooled["selection"],
            "seed": int(pooled["seed"]) + int(shard.rsplit("_", 1)[1]),
            "logging": {"directory": str(root / "logs" / shard), "stage_log": str(root / "logs" / shard / "stage_resources.jsonl")},
            "checkpoints": {"enabled": False, "directory": str(root / "unused" / shard)},
            "mixscape": pooled["mixscape"],
        }
        config_path = root / "configs" / f"{shard}.yaml"
        dump_yaml(config, config_path)
        commands.append(([rscript(), MIXSCAPE / "run_mixscape.R", config_path], root / "logs" / f"{shard}.log", r_env))
    _parallel(commands, workers)
    run([sys.executable, MIXSCAPE / "target_sharding" / "merge_pooled_shards.py", root, shared / "input" / "metadata.tsv.gz"])
    return root / "merged" / "mixscape_metadata.tsv.gz"


def run_mixscale(spec: dict, shared: Path, work: Path, workers: int, r_env: dict[str, str]) -> Path:
    config = load_yaml(MIXSCALE / "config" / f'{spec["dataset"]}_pooled.yaml')
    root = work / "mixscale"
    config.update({
        "dataset": spec["dataset"], "source_shard_root": str(shared / "shards"),
        "canonical_metadata": str(shared / "input" / "metadata.tsv.gz"),
        "output_root": str(root), "max_workers": workers,
    })
    config_path = work / "mixscale.yaml"
    dump_yaml(config, config_path)
    commands = []
    for checkpoint in sorted((shared / "shards" / "checkpoints").glob("shard_*.rds")):
        shard = checkpoint.stem
        commands.append(([rscript(), MIXSCALE / "sharding" / "run_shard.R", config_path, checkpoint, root / "shards" / shard], root / "logs" / f"{shard}.log", r_env))
    _parallel(commands, workers)
    run([rscript(), MIXSCALE / "sharding" / "finalize.R", config_path, root], env=r_env)
    return root / "merged" / "mixscale_metadata.tsv.gz"


RUNNERS = {"mixscape": run_mixscape, "mixscale": run_mixscale}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--data.h5ad", required=True)
    parser.add_argument("--data.spec", required=True)
    parser.add_argument("--shared_signature.bundle", required=True)
    parser.add_argument("--method", required=True, choices=sorted(RUNNERS))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(getattr(args, "data.spec")).read_text(encoding="utf-8"))
    spec["h5ad"] = str(Path(getattr(args, "data.h5ad")).resolve())
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=f"pse-{args.method}-{args.name}-") as temp_text:
        temp = Path(temp_text)
        archives = [VENDOR / "presto-3c97180d.tar.gz"]
        if args.method == "mixscape":
            archives.insert(0, VENDOR / "seurat-084d9e40.tar.gz")
        else:
            archives.append(VENDOR / "mixscale-3d3a2cd1.tar.gz")
        r_env = prepare_r_library(temp, archives)
        extract_archive(Path(getattr(args, "shared_signature.bundle")), temp)
        result = RUNNERS[args.method](spec, temp / "shared", temp / "work", args.workers, r_env)
        if not result.is_file():
            raise FileNotFoundError(result)
        result_out = output / "status_result.tsv.gz"
        if result.suffix == ".gz":
            shutil.copy2(result, result_out)
        else:
            with result.open("rb") as source, gzip.open(result_out, "wb") as destination:
                shutil.copyfileobj(source, destination)
        method_root = result.parent.parent
        make_archive(method_root, output / "method_aux.tar.gz", args.method)
    elapsed = time.monotonic() - started
    write_json({
        "dataset": spec["dataset"], "method": args.method, "status": "PASS",
        "wall_seconds": elapsed, "workers": args.workers,
        "peak_child_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "result_file": "status_result.tsv.gz", "aux_archive": "method_aux.tar.gz",
    }, output / "method_manifest.json")


if __name__ == "__main__":
    main()
