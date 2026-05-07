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
    title: str | None = None,
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

    ax.set_title(title or f"Replicate QC: {group_name}")
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


def plot_drift(
    shift_records: list[dict],
    output_path: str | Path,
    warn_threshold_eV: float = 2.0,
    quality_threshold: float = 0.7,
) -> None:
    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(10, 4))

    reliable = [rec for rec in shift_records if rec["alignment_quality"] >= quality_threshold]
    uncertain = [rec for rec in shift_records if rec["alignment_quality"] < quality_threshold]
    labelled = False

    if reliable:
        ax.scatter(
            [rec["scan_index"] for rec in reliable],
            [rec["shift_eV"] for rec in reliable],
            marker="o",
            color="steelblue",
            zorder=3,
            label=f"reliable (quality ≥ {quality_threshold:.2f})",
        )
        labelled = True
    if uncertain:
        ax.scatter(
            [rec["scan_index"] for rec in uncertain],
            [rec["shift_eV"] for rec in uncertain],
            marker="o",
            facecolors="none",
            edgecolors="tomato",
            zorder=3,
            label=f"uncertain (quality < {quality_threshold:.2f})",
        )
        labelled = True

    ax.axhline(0, color="grey", linewidth=0.8, linestyle="-")
    ax.axhline(warn_threshold_eV, color="red", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(-warn_threshold_eV, color="red", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_title("Energy shift per scan — drift tracker")
    ax.set_xlabel("Scan index")
    ax.set_ylabel("Energy shift (eV)")
    if labelled:
        ax.legend()
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_detector_health_overview(
    records: list[dict],
    output_path: str | Path,
    channels: list[tuple[str, str]],
    title: str = "Detector health overview",
    energy_range: tuple[float, float] | None = None,
) -> dict:
    """Save a stacked detector-channel QC plot using individual scan traces."""
    output_path = Path(output_path)
    included_channels = []
    skipped_channels = []
    channel_data = []

    for key, label in channels:
        traces = []
        for rec in records:
            values = rec.get(key)
            if values is None:
                continue
            x, y = _finite_xy(rec["energy"], values)
            if energy_range is not None:
                lo, hi = energy_range
                mask = (x >= lo) & (x <= hi)
                x, y = x[mask], y[mask]
            if len(x) == 0:
                continue
            traces.append((x, y, rec.get("label", rec.get("filename", "scan"))))
        if traces:
            included_channels.append(label)
            channel_data.append((label, traces))
        else:
            skipped_channels.append(label)

    if not channel_data:
        return {"path": None, "channels": included_channels, "skipped": skipped_channels}

    plt = _setup_matplotlib()
    n_panels = len(channel_data)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=(10, max(2.2 * n_panels, 4)),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    legend_count = max(len(traces) for _, traces in channel_data)
    for ax, (label, traces) in zip(axes, channel_data):
        for x, y, trace_label in traces:
            ax.plot(
                x,
                y,
                linewidth=0.75,
                alpha=0.55,
                label=trace_label if legend_count <= 12 else None,
            )
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.25)
        if legend_count <= 12:
            ax.legend(fontsize=6, loc="best")

    axes[-1].set_xlabel("Energy / eV")
    subtitle = "individual scan traces"
    if skipped_channels:
        subtitle += "; skipped: " + ", ".join(skipped_channels)
    fig.suptitle(f"{title} ({subtitle})")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return {"path": output_path, "channels": included_channels, "skipped": skipped_channels}


def plot_analysis_signal_qc(
    records: list[dict],
    output_path: str | Path,
    signal_label: str,
    title: str = "Analysis signal QC: individual scan traces before normalization",
    energy_range: tuple[float, float] | None = None,
) -> dict:
    """Save individual analysis-signal traces and an optional diagnostic average."""
    output_path = Path(output_path)
    traces = []
    skipped = []

    for rec in records:
        values = rec.get("signal")
        label = rec.get("label", rec.get("filename", "scan"))
        if values is None:
            skipped.append(f"{label}: missing signal")
            continue
        x, y = _finite_xy(rec["energy"], values)
        if energy_range is not None:
            lo, hi = energy_range
            mask = (x >= lo) & (x <= hi)
            x, y = x[mask], y[mask]
        if len(x) == 0:
            skipped.append(f"{label}: no finite signal in plot range")
            continue
        traces.append((x, y, label))

    if not traces:
        return {"path": None, "signal": signal_label, "n_traces": 0, "skipped": skipped}

    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(10, 5))

    for x, y, label in traces:
        ax.plot(
            x,
            y,
            linewidth=0.8,
            alpha=0.55,
            label=label if len(traces) <= 12 else None,
        )

    average_plotted = False
    if len(traces) > 1:
        master_x = traces[0][0]
        interp_values = []
        for x, y, _ in traces:
            if len(x) < 2:
                continue
            y_interp = np.interp(master_x, x, y, left=np.nan, right=np.nan)
            interp_values.append(y_interp)
        if interp_values:
            interp_array = np.asarray(interp_values, dtype=float)
            valid_counts = np.sum(np.isfinite(interp_array), axis=0)
            sums = np.nansum(interp_array, axis=0)
            avg = np.full_like(master_x, np.nan, dtype=float)
            valid_avg = valid_counts > 0
            avg[valid_avg] = sums[valid_avg] / valid_counts[valid_avg]
            valid = np.isfinite(master_x) & np.isfinite(avg)
            if np.any(valid):
                ax.plot(master_x[valid], avg[valid], color="black", linewidth=2.0, label="average", zorder=4)
                average_plotted = True

    ax.set_title(title)
    ax.set_xlabel("Energy / eV")
    ax.set_ylabel(signal_label)
    ax.grid(True, alpha=0.25)
    if len(traces) <= 12 or average_plotted:
        ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return {
        "path": output_path,
        "signal": signal_label,
        "n_traces": len(traces),
        "skipped": skipped,
        "average_plotted": average_plotted,
    }
