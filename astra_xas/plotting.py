from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _finite_xy(energy, mu):
    energy = np.asarray(energy, dtype=float)
    mu = np.asarray(mu, dtype=float)
    mask = np.isfinite(energy) & np.isfinite(mu)
    return energy[mask], mu[mask]


def plot_overview(records: Iterable[dict], output_path: str | Path, title: str, y_label: str, energy_range: tuple[float, float] | None = None):
    """Save one overview plot containing all spectra in *records*."""
    records = list(records)
    if not records:
        return None
    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    plotted = 0
    for rec in records:
        x, y = _finite_xy(rec["energy"], rec["mu"])
        if energy_range is not None:
            lo, hi = energy_range
            mask = (x >= lo) & (x <= hi)
            x, y = x[mask], y[mask]
        if len(x) == 0:
            continue
        ax.plot(x, y, linewidth=1.2, label=rec.get("label", "spectrum"))
        plotted += 1
    if plotted == 0:
        plt.close(fig)
        return None
    ax.set_title(title)
    ax.set_xlabel("Energy / eV")
    ax.set_ylabel(y_label)
    if plotted <= 18:
        ax.legend(fontsize=7, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def plot_replicate_qc(
    group_name: str,
    replicate_records: list[dict],
    average_record: dict,
    output_path: str | Path,
    energy_range: tuple[float, float] | None = None,
    y_label: str = "Normalized intensity",
):
    """Save a QC plot showing individual replicate spectra and their average."""
    if not replicate_records:
        return None
    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    plotted = 0
    for rec in replicate_records:
        x, y = _finite_xy(rec["energy"], rec["mu"])
        if energy_range is not None:
            lo, hi = energy_range
            mask = (x >= lo) & (x <= hi)
            x, y = x[mask], y[mask]
        if len(x) == 0:
            continue
        ax.plot(x, y, linewidth=0.8, alpha=0.55, label=rec.get("label", "replicate"))
        plotted += 1

    x, y = _finite_xy(average_record["energy"], average_record["mu"])
    if energy_range is not None:
        lo, hi = energy_range
        mask = (x >= lo) & (x <= hi)
        x, y = x[mask], y[mask]
    if len(x) > 0:
        ax.plot(x, y, linewidth=2.2, label=average_record.get("label", "average"))
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return None

    ax.set_title(f"Replicate QC: {group_name}")
    ax.set_xlabel("Energy / eV")
    ax.set_ylabel(y_label)
    if len(replicate_records) <= 12:
        ax.legend(fontsize=7, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path
