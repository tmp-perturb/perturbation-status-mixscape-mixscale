#!/usr/bin/env python3
"""Run one perturbation-status method from a common frozen cohort."""

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

from common import dump_yaml, extract_archive, load_yaml, make_archive, r_env, rscript, run, write_json


def _parallel(commands: list[tuple[list[object], Path, dict[str, str]]], workers: int) -> None:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run, command, log=log, env=env) for command, log, env in commands]
        for future in futures:
            future.result()


def run_mixscape(project: Path, spec: dict, shared: Path, work: Path, workers: int) -> Path:
    pooled = load_yaml(spec["pooled_config"])
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
        commands.append(([rscript(project, "mixscape"), project / "03_mixscape" / "run_mixscape.R", config_path], root / "logs" / f"{shard}.log", r_env(project, "mixscape")))
    _parallel(commands, workers)
    run([sys.executable, project / "03_mixscape" / "target_sharding" / "merge_pooled_shards.py", root, shared / "input" / "metadata.tsv.gz"])
    return root / "merged" / "mixscape_metadata.tsv.gz"


def run_mixscale(project: Path, spec: dict, shared: Path, work: Path, workers: int) -> Path:
    config = load_yaml(spec["mixscale_config"])
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
        commands.append(([rscript(project, "mixscale"), project / "04_mixscale" / "sharding" / "run_shard.R", config_path, checkpoint, root / "shards" / shard], root / "logs" / f"{shard}.log", r_env(project, "mixscale")))
    _parallel(commands, workers)
    run([rscript(project, "mixscale"), project / "04_mixscale" / "sharding" / "finalize.R", config_path, root], env=r_env(project, "mixscale"))
    return root / "merged" / "mixscale_metadata.tsv.gz"


def run_ps(project: Path, spec: dict, shared: Path, work: Path, workers: int) -> Path:
    input_config = load_yaml(spec["ps_input_config"])
    input_dir = work / "ps_input"
    input_config.update({
        "dataset": spec["dataset"], "source_h5ad": spec["h5ad"],
        "input_dir": str(input_dir), "shared_sparse_input_dir": str(shared / "input"),
        "seurat_rds": str(input_dir / "ps_input.rds"),
    })
    input_path = work / "ps_input.yaml"
    dump_yaml(input_config, input_path)
    run([sys.executable, project / "05_ps" / "adapt_ps_input.py", "--config", input_path])
    run([rscript(project, "ps"), project / "05_ps" / "build_ps_seurat.R", "--config", input_path], env=r_env(project, "ps"))
    config = load_yaml(spec["ps_run_config"])
    output = work / "ps"
    config.update({
        "dataset": spec["dataset"], "input_rds": str(input_dir / "ps_input.rds"),
        "barcode_table": str(input_dir / "barcode.tsv.gz"), "output_dir": str(output),
        "work_dir": str(work / "ps_work"), "max_workers": workers,
    })
    config.pop("reference_result", None)
    config_path = work / "ps_run.yaml"
    dump_yaml(config, config_path)
    run([sys.executable, project / "05_ps" / "sharding" / "run_scalable.py", "--config", config_path], env=r_env(project, "ps"))
    return output / "ps_metadata.tsv"


def run_baseline(project: Path, spec: dict, shared: Path, work: Path, workers: int) -> Path:
    config = load_yaml(spec["baseline_config"])
    output = work / "target_expression"
    config.update({
        "dataset": spec["dataset"], "checkpoint_dir": str(shared / "shards" / "checkpoints"),
        "target_manifest": str(shared / "shards" / "target_manifest.tsv"),
        "metadata": str(shared / "input" / "metadata.tsv.gz"), "output_dir": str(output),
        "max_workers": workers,
    })
    config_path = work / "baseline.yaml"
    dump_yaml(config, config_path)
    run([sys.executable, project / "06_target_expression_baseline" / "run_baseline.py", "--config", config_path], env=r_env(project, "mixscape"))
    return output / "target_expression_baseline.tsv"


RUNNERS = {"mixscape": run_mixscape, "mixscale": run_mixscale, "ps": run_ps, "target_expression": run_baseline}


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
    project = Path(spec["project_root"])
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=f"pse-{args.method}-{args.name}-") as temp_text:
        temp = Path(temp_text)
        extract_archive(Path(getattr(args, "shared_signature.bundle")), temp)
        result = RUNNERS[args.method](project, spec, temp / "shared", temp / "work", args.workers)
        if not result.is_file():
            raise FileNotFoundError(result)
        result_out = output / "status_result.tsv.gz"
        if result.suffix == ".gz":
            shutil.copy2(result, result_out)
        else:
            with result.open("rb") as source, gzip.open(result_out, "wb") as destination:
                shutil.copyfileobj(source, destination)
        method_root = result.parent.parent if args.method in {"mixscape", "mixscale"} else result.parent
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
