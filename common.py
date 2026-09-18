"""Shared execution helpers for perturbation-status Omnibenchmark modules."""

from __future__ import annotations

import json
import os
import subprocess
import tarfile
from pathlib import Path
from typing import Any

import yaml


THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected YAML mapping: {path}")
    return value


def dump_yaml(value: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def run(command: list[str | Path], *, log: Path | None = None, env: dict[str, str] | None = None) -> None:
    command = [str(item) for item in command]
    print("+ " + " ".join(command), flush=True)
    merged_env = {**os.environ, **THREAD_ENV, **(env or {})}
    if log is None:
        subprocess.run(command, check=True, env=merged_env)
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as stream:
        subprocess.run(command, check=True, env=merged_env, stdout=stream, stderr=subprocess.STDOUT)


def rscript(project: Path, environment: str) -> Path:
    path = project / "environments" / environment / "env" / "bin" / "Rscript"
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def r_env(project: Path, environment: str) -> dict[str, str]:
    prefix = project / "environments" / environment / "env"
    return {
        "R_LIBS_USER": str(prefix / "lib" / "R" / "library"),
        "PATH": f"{prefix / 'bin'}:{os.environ.get('PATH', '')}",
    }


def make_archive(source: Path, destination: Path, arcname: str) -> None:
    with tarfile.open(destination, "w:gz", compresslevel=1) as archive:
        archive.add(source, arcname=arcname, recursive=True)


def extract_archive(source: Path, destination: Path) -> None:
    with tarfile.open(source, "r:gz") as archive:
        root = destination.resolve()
        for member in archive.getmembers():
            target = (destination / member.name).resolve()
            if root not in target.parents and target != root:
                raise ValueError(f"Unsafe archive member: {member.name}")
        archive.extractall(destination)


def write_json(value: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
