from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

from .config import AstraConfig


@dataclass
class SelfAbsorptionResult:
    sample: str
    status: str
    severity: str
    ratio_fluo_over_trans: float
    fluo_white_line_amplitude: float
    trans_white_line_amplitude: float
    threshold_used: float
    sensitivity: str
    white_line_window_eV: tuple[float, float]
    continuum_window_eV: tuple[float, float]
    n_replicates_used: int
    note: str
    fluo_norm: np.ndarray | None = field(default=None, repr=False)
    trans_norm: np.ndarray | None = field(default=None, repr=False)


def get_self_absorption_threshold(config: AstraConfig) -> float:
    sensitivity = str(getattr(config, "self_absorption_sensitivity", "normal")).strip().lower()
    thresholds = {
        "relaxed": 0.75,
        "normal": 0.85,
        "strict": 0.92,
    }
    if sensitivity == "custom":
        try:
            return float(getattr(config, "self_absorption_custom_threshold", 0.85))
        except (TypeError, ValueError):
            return 0.85
    return thresholds.get(sensitivity, thresholds["normal"])


def get_self_absorption_sensitivity(config: AstraConfig) -> str:
    sensitivity = str(getattr(config, "self_absorption_sensitivity", "normal")).strip().lower()
    if sensitivity in {"relaxed", "normal", "strict", "custom"}:
        return sensitivity
    return "normal"


def _format_float(value: float, fmt: str) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not np.isfinite(value):
        return "nan"
    return format(value, fmt)


def _smooth_spectrum(y: np.ndarray, window_length: int = 9, polyorder: int = 2) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(y)
    n_finite = int(np.count_nonzero(finite))
    if n_finite < 5:
        return y.copy()

    try:
        from scipy.signal import savgol_filter
    except Exception:
        return y.copy()

    window = min(int(window_length), n_finite)
    if window % 2 == 0:
        window -= 1
    if window < 3 or window <= polyorder:
        return y.copy()

    x = np.arange(len(y), dtype=float)
    filled = np.interp(x, x[finite], y[finite])
    try:
        smoothed = savgol_filter(filled, window_length=window, polyorder=polyorder)
    except Exception:
        return y.copy()
    smoothed[~finite] = np.nan
    return np.asarray(smoothed, dtype=float)


def _normalization_ok(pe: dict) -> bool:
    edge_step = pe.get("edge_step")
    norm = pe.get("norm")
    try:
        edge_step = float(edge_step)
    except (TypeError, ValueError):
        return False
    if not np.isfinite(edge_step) or edge_step <= 0:
        return False
    if norm is None:
        return False
    return np.isfinite(np.asarray(norm, dtype=float)).any()


def measure_white_line_amplitude(
    energy,
    norm_mu,
    e0,
    wl_min,
    wl_max,
    cont_min,
    cont_max,
    min_points: int = 5,
    use_continuum_fallback: bool = True,
) -> dict:
    energy = np.asarray(energy, dtype=float)
    norm_mu = _smooth_spectrum(np.asarray(norm_mu, dtype=float))
    wl_lo = float(e0) + float(wl_min)
    wl_hi = float(e0) + float(wl_max)
    cont_lo = float(e0) + float(cont_min)
    cont_hi = float(e0) + float(cont_max)

    wl_mask = (energy >= wl_lo) & (energy <= wl_hi) & np.isfinite(energy) & np.isfinite(norm_mu)
    n_wl = int(np.count_nonzero(wl_mask))
    if n_wl < int(min_points):
        return {
            "ok": False,
            "amplitude": float("nan"),
            "white_line_level": float("nan"),
            "continuum_level": float("nan"),
            "n_wl": n_wl,
            "n_cont": 0,
            "note": "skipped: insufficient finite points in white-line window",
        }

    cont_mask = (energy >= cont_lo) & (energy <= cont_hi) & np.isfinite(energy) & np.isfinite(norm_mu)
    n_cont = int(np.count_nonzero(cont_mask))
    note = ""
    if n_cont >= int(min_points):
        continuum_level = float(np.nanmedian(norm_mu[cont_mask]))
    elif use_continuum_fallback:
        continuum_level = 1.0
        note = "continuum window unavailable; used continuum_level = 1.0 fallback"
    else:
        return {
            "ok": False,
            "amplitude": float("nan"),
            "white_line_level": float("nan"),
            "continuum_level": float("nan"),
            "n_wl": n_wl,
            "n_cont": n_cont,
            "note": "skipped: insufficient finite points in continuum window",
        }

    white_line_level = float(np.nanpercentile(norm_mu[wl_mask], 95))
    amplitude = white_line_level - continuum_level
    return {
        "ok": bool(np.isfinite(amplitude)),
        "amplitude": float(amplitude),
        "white_line_level": white_line_level,
        "continuum_level": float(continuum_level),
        "n_wl": n_wl,
        "n_cont": n_cont,
        "note": note,
    }


def classify_self_absorption_ratio(ratio: float, threshold: float) -> tuple[str, str]:
    if not np.isfinite(ratio):
        return "skipped", "none"
    if ratio >= threshold:
        return "ok", "none"
    if ratio >= 0.75:
        return "flagged", "mild"
    if ratio >= 0.60:
        return "flagged", "moderate"
    return "flagged", "strong"


def _skipped_result(
    group_name: str,
    config: AstraConfig,
    n_replicates_used: int,
    note: str,
    fluo_amp: float = float("nan"),
    trans_amp: float = float("nan"),
    fluo_norm=None,
    trans_norm=None,
) -> SelfAbsorptionResult:
    return SelfAbsorptionResult(
        sample=group_name,
        status="skipped",
        severity="none",
        ratio_fluo_over_trans=float("nan"),
        fluo_white_line_amplitude=fluo_amp,
        trans_white_line_amplitude=trans_amp,
        threshold_used=get_self_absorption_threshold(config),
        sensitivity=get_self_absorption_sensitivity(config),
        white_line_window_eV=(
            float(config.e0) + float(getattr(config, "self_absorption_wl_min", 0.0)),
            float(config.e0) + float(getattr(config, "self_absorption_wl_max", 35.0)),
        ),
        continuum_window_eV=(
            float(config.e0) + float(getattr(config, "self_absorption_cont_min", 50.0)),
            float(config.e0) + float(getattr(config, "self_absorption_cont_max", 150.0)),
        ),
        n_replicates_used=int(n_replicates_used),
        note=note,
        fluo_norm=fluo_norm,
        trans_norm=trans_norm,
    )


def evaluate_self_absorption_for_group(
    group_name: str,
    energy,
    fluo_mu_average,
    trans_mu_average,
    config: AstraConfig,
    n_replicates_used: int,
    normalize_func: Callable,
) -> SelfAbsorptionResult:
    energy = np.asarray(energy, dtype=float)
    fluo_mu_average = np.asarray(fluo_mu_average, dtype=float)
    trans_mu_average = np.asarray(trans_mu_average, dtype=float)
    if not np.isfinite(trans_mu_average).any():
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            "skipped: usable transmission signal not available",
        )

    try:
        pe_fluo = normalize_func(energy, fluo_mu_average, config)
        pe_trans = normalize_func(energy, trans_mu_average, config)
    except Exception:
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            "skipped: diagnostic normalization failed",
        )

    fluo_norm = np.asarray(pe_fluo.get("norm", np.full_like(energy, np.nan)), dtype=float)
    trans_norm = np.asarray(pe_trans.get("norm", np.full_like(energy, np.nan)), dtype=float)
    if not _normalization_ok(pe_fluo) or not _normalization_ok(pe_trans):
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            "skipped: diagnostic normalization failed",
            fluo_norm=fluo_norm,
            trans_norm=trans_norm,
        )

    min_points = int(getattr(config, "self_absorption_min_points", 5))
    wl_min = float(getattr(config, "self_absorption_wl_min", 0.0))
    wl_max = float(getattr(config, "self_absorption_wl_max", 35.0))
    cont_min = float(getattr(config, "self_absorption_cont_min", 50.0))
    cont_max = float(getattr(config, "self_absorption_cont_max", 150.0))

    fluo = measure_white_line_amplitude(
        energy, fluo_norm, config.e0, wl_min, wl_max, cont_min, cont_max, min_points=min_points
    )
    if not fluo["ok"]:
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            fluo["note"],
            fluo_amp=fluo["amplitude"],
            fluo_norm=fluo_norm,
            trans_norm=trans_norm,
        )

    trans = measure_white_line_amplitude(
        energy, trans_norm, config.e0, wl_min, wl_max, cont_min, cont_max, min_points=min_points
    )
    if not trans["ok"]:
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            trans["note"],
            fluo_amp=fluo["amplitude"],
            trans_amp=trans["amplitude"],
            fluo_norm=fluo_norm,
            trans_norm=trans_norm,
        )

    trans_amp = float(trans["amplitude"])
    min_trans_amp = float(getattr(config, "self_absorption_min_trans_amp", 0.03))
    if not np.isfinite(trans_amp) or trans_amp <= 0:
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            "skipped: invalid transmission white-line amplitude",
            fluo_amp=fluo["amplitude"],
            trans_amp=trans_amp,
            fluo_norm=fluo_norm,
            trans_norm=trans_norm,
        )
    if trans_amp < min_trans_amp:
        return _skipped_result(
            group_name,
            config,
            n_replicates_used,
            "skipped: transmission white-line amplitude below minimum reliable threshold",
            fluo_amp=fluo["amplitude"],
            trans_amp=trans_amp,
            fluo_norm=fluo_norm,
            trans_norm=trans_norm,
        )

    threshold = get_self_absorption_threshold(config)
    ratio = float(fluo["amplitude"]) / trans_amp
    status, severity = classify_self_absorption_ratio(ratio, threshold)
    notes = [note for note in (fluo["note"], trans["note"]) if note]
    if status == "flagged":
        notes.insert(
            0,
            "possible self-absorption: fluorescence white-line amplitude suppressed relative to transmission",
        )
    else:
        notes.insert(0, "OK")

    return SelfAbsorptionResult(
        sample=group_name,
        status=status,
        severity=severity,
        ratio_fluo_over_trans=float(ratio),
        fluo_white_line_amplitude=float(fluo["amplitude"]),
        trans_white_line_amplitude=trans_amp,
        threshold_used=threshold,
        sensitivity=get_self_absorption_sensitivity(config),
        white_line_window_eV=(float(config.e0) + wl_min, float(config.e0) + wl_max),
        continuum_window_eV=(float(config.e0) + cont_min, float(config.e0) + cont_max),
        n_replicates_used=int(n_replicates_used),
        note="; ".join(notes),
        fluo_norm=fluo_norm,
        trans_norm=trans_norm,
    )


def write_self_absorption_flags(results: list[SelfAbsorptionResult], output_dir: str | Path, config: AstraConfig) -> Path:
    output_dir = Path(output_dir)
    path = output_dir / "ASTRA_self_absorption_flags.dat"
    threshold = get_self_absorption_threshold(config)
    wl_window = (
        float(config.e0) + float(getattr(config, "self_absorption_wl_min", 0.0)),
        float(config.e0) + float(getattr(config, "self_absorption_wl_max", 35.0)),
    )
    cont_window = (
        float(config.e0) + float(getattr(config, "self_absorption_cont_min", 50.0)),
        float(config.e0) + float(getattr(config, "self_absorption_cont_max", 150.0)),
    )
    with path.open("w", encoding="utf-8") as f:
        f.write("# ASTRA self-absorption diagnostic\n")
        f.write(f"# ASTRA version: {getattr(config, 'version', 'N/A')}\n")
        f.write("# Diagnostic name: possible fluorescence self-absorption / over-absorption flag\n")
        f.write(f"# Sensitivity: {get_self_absorption_sensitivity(config)}\n")
        f.write(f"# Threshold: {threshold:.4f}\n")
        f.write(f"# White-line window: {wl_window[0]:.6f}-{wl_window[1]:.6f} eV\n")
        f.write(f"# Continuum window: {cont_window[0]:.6f}-{cont_window[1]:.6f} eV\n")
        f.write(f"# Minimum transmission white-line amplitude: {float(getattr(config, 'self_absorption_min_trans_amp', 0.03)):.6g}\n")
        f.write("# This is a heuristic diagnostic and not a correction. It does not modify data.\n")
        f.write("# It should not be treated as proof that self-absorption is confirmed.\n")
        f.write(
            "# sample\tstatus\tseverity\tratio_fluo_over_trans\tfluo_white_line_amplitude\t"
            "trans_white_line_amplitude\tthreshold_used\tsensitivity\twhite_line_window_eV\t"
            "continuum_window_eV\tn_replicates_used\tnote\n"
        )
        for result in results:
            wl = f"{result.white_line_window_eV[0]:.6f}-{result.white_line_window_eV[1]:.6f}"
            cont = f"{result.continuum_window_eV[0]:.6f}-{result.continuum_window_eV[1]:.6f}"
            f.write(
                f"{result.sample}\t{result.status}\t{result.severity}\t"
                f"{_format_float(result.ratio_fluo_over_trans, '.4f')}\t"
                f"{_format_float(result.fluo_white_line_amplitude, '.6g')}\t"
                f"{_format_float(result.trans_white_line_amplitude, '.6g')}\t"
                f"{_format_float(result.threshold_used, '.4f')}\t{result.sensitivity}\t"
                f"{wl}\t{cont}\t{result.n_replicates_used}\t{result.note}\n"
            )
    return path


def plot_self_absorption_qc(
    result: SelfAbsorptionResult,
    energy,
    fluo_norm,
    trans_norm,
    output_path: str | Path,
) -> Path | None:
    if result.status not in {"ok", "flagged"}:
        return None
    energy = np.asarray(energy, dtype=float)
    fluo_norm = np.asarray(fluo_norm, dtype=float)
    trans_norm = np.asarray(trans_norm, dtype=float)
    mask_f = np.isfinite(energy) & np.isfinite(fluo_norm)
    mask_t = np.isfinite(energy) & np.isfinite(trans_norm)
    if not mask_f.any() or not mask_t.any():
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(energy[mask_f], fluo_norm[mask_f], color="darkorange", linewidth=1.4, label="normalized fluorescence diagnostic")
    ax.plot(energy[mask_t], trans_norm[mask_t], color="steelblue", linewidth=1.4, label="normalized sample transmission diagnostic")
    ax.axvspan(result.white_line_window_eV[0], result.white_line_window_eV[1], color="gold", alpha=0.18, label="white-line window")
    ax.axvspan(result.continuum_window_eV[0], result.continuum_window_eV[1], color="grey", alpha=0.14, label="continuum window")
    ax.set_title(f"Self-absorption QC: {result.sample}")
    ax.set_xlabel("Energy / eV")
    ax.set_ylabel("Normalized diagnostic μ(E)")
    ratio = _format_float(result.ratio_fluo_over_trans, ".4f")
    ax.text(
        0.02,
        0.98,
        f"ratio={ratio}, threshold={result.threshold_used:.4f}\nstatus={result.status}, severity={result.severity}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "lightgrey", "alpha": 0.9},
    )
    ax.legend(fontsize=7, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path
