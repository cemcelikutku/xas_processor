from __future__ import annotations
import os
import re
import fnmatch
from pathlib import Path
import numpy as np


def natural_key(s: str):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s)]


def collect_xasd_files(input_dir: str | Path) -> list[Path]:
    input_dir = Path(input_dir).expanduser().resolve()
    files: list[Path] = []
    for root_dir, _, names in os.walk(input_dir):
        for fname in fnmatch.filter(names, "*.xasd"):
            files.append(Path(root_dir) / fname)
    return sorted(files, key=lambda p: natural_key(p.name))


def load_xasd(path: str | Path) -> dict:
    """Read an ASTRA .xasd file.

    Expected numeric columns:
    E, Theta, dt, I0, I1, I2, IF, FDT, Ir
    """
    path = Path(path)
    data = np.genfromtxt(path, comments="#", delimiter=",")
    if data.ndim != 2 or data.shape[1] < 7:
        raise ValueError(f"Could not parse numeric data from {path}")
    return {
        "path": path,
        "filename": path.name,
        "energy": data[:, 0],
        "theta": data[:, 1],
        "dt": data[:, 2],
        "I0": data[:, 3],
        "I1": data[:, 4],
        "I2": data[:, 5],
        "IF": data[:, 6],
        "FDT": data[:, 7] if data.shape[1] > 7 else None,
        "Ir": data[:, 8] if data.shape[1] > 8 else None,
    }


def split_replicate_suffix(filename: str) -> tuple[str, int | None]:
    stem = Path(filename).stem
    m = re.match(r"^(.*)_(\d+)$", stem)
    if m:
        return m.group(1), int(m.group(2))
    return stem, None


def sanitize_name(name: str) -> str:
    name = name.strip()
    name = re.sub(r"[^\w\-.]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def build_default_output_directory(input_dir: str | Path) -> Path:
    input_dir = Path(input_dir).expanduser().resolve()
    return input_dir.parent / f"{input_dir.name}-processed"
