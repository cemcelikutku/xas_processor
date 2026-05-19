"""Per-scan validation, QC, and signal computation.

This module owns the per-scan layer of the AstraXAS pipeline:

  - Build the canonical entry dict from a raw ``load_xasd`` scan
    (``_entry_from_scan``).
  - Compute analysis / alignment signal specs from config.
  - Run per-scan validation (channel coverage, alignment-signal structure,
    energy-range overlap) via ``_validate_single_scan``.
  - Detect per-channel detector jumps (``detect_detector_jumps`` and the
    ``_detect_jumps_for_entry`` aggregator).
  - Compose the above into a single ``process_single_scan`` entry point
    that returns a ``SingleScanResult`` consumable by callers that need
    per-scan QC without running the full batch pipeline.

Architectural rule (Phase 2.1): this module must NOT import from
``astra_xas.processor`` (absolute or relative). ``processor.py`` may import
from here. Downstream consumers (``beamtime/watcher.py``,
``beamtime/groups.py``, ``manifest.py``, etc.) continue to import the
per-scan helpers from ``astra_xas.processor`` via compatibility
re-exports added in the same phase.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import AstraConfig
from .io import split_replicate_suffix
from .signals import compute_signals


PRIMARY_DETECTOR_JUMP_CHANNELS = {"I0", "I1", "I2", "IF"}
FDT_DETECTOR_JUMP_CHANNELS = {"FDT"}
RAW_DETECTOR_JUMP_CHANNELS = PRIMARY_DETECTOR_JUMP_CHANNELS | FDT_DETECTOR_JUMP_CHANNELS
DERIVED_DETECTOR_JUMP_CHANNELS = {"IF_over_I0", "ln_I0_I1", "ln_I1_I2"}


# ---------------------------------------------------------------------------
# Tolerant config accessors
# ---------------------------------------------------------------------------

def _config_bool(config, name: str, default: bool = False) -> bool:
    value = getattr(config, name, default)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _config_float(config, name: str, default: float) -> float:
    try:
        return float(getattr(config, name, default))
    except (TypeError, ValueError):
        return float(default)


# ---------------------------------------------------------------------------
# Signal specs
# ---------------------------------------------------------------------------

def _analysis_signal_spec(config: AstraConfig) -> tuple[str, str]:
    mode = getattr(config, "analysis_mode", "fluo")
    if mode == "trans":
        return "mu_trans", "ln(I0/I1)"
    if mode == "ref":
        return "mu_ref", "ln(I1/I2)"
    return "mu_fluo", "IF/I0"


def _alignment_signal_spec(alignment_source: str, config: AstraConfig) -> tuple[str, str, str]:
    mode = "ref" if alignment_source == "inline_ref" else getattr(config, "foil_alignment_mode", "trans")
    mapping = {
        "trans": ("mu_trans", "ln(I0/I1)"),
        "ref": ("mu_ref", "ln(I1/I2)"),
        "fluo": ("mu_fluo", "IF/I0"),
    }
    signal_key, signal_label = mapping.get(mode, mapping["trans"])
    return mode, signal_key, signal_label


def _required_channels_for_signal(mode: str) -> list[tuple[str, str]]:
    if mode == "trans":
        return [("I0", "I0"), ("I1", "I1")]
    if mode == "ref":
        return [("I1", "I1"), ("I2", "I2")]
    return [("IF", "IF"), ("I0", "I0")]


# ---------------------------------------------------------------------------
# Entry assembly
# ---------------------------------------------------------------------------

def _entry_from_scan(scan: dict, config: AstraConfig, path: Path | None = None) -> dict:
    sigs = compute_signals(scan, config)
    source_path = Path(path) if path is not None else Path(scan.get("path", scan["filename"]))
    base_name, replicate_id = split_replicate_suffix(scan["filename"])
    return {
        "path": source_path,
        "filename": scan["filename"],
        "is_foil": config.foil_keyword.lower() in scan["filename"].lower() if config.foil_keyword else False,
        "energy": sigs["energy"],
        "mu_trans": sigs["mu_trans"],
        "mu_ref": sigs["mu_ref"],
        "mu_fluo": sigs["mu_fluo"],
        "I0": scan.get("I0"),
        "I1": scan.get("I1"),
        "I2": scan.get("I2"),
        "IF": scan.get("IF"),
        "FDT": scan.get("FDT"),
        "Ir": scan.get("Ir"),
        "base_name": base_name,
        "replicate_id": replicate_id,
    }


# ---------------------------------------------------------------------------
# Per-scan validation primitives
# ---------------------------------------------------------------------------

def _range_overlap_status(name: str, selected_range: tuple[float, float], data_range: tuple[float, float]) -> str | None:
    lo, hi = sorted((float(selected_range[0]), float(selected_range[1])))
    data_lo, data_hi = data_range
    if hi < data_lo or lo > data_hi:
        return (
            f"{name} {lo:.6g}-{hi:.6g} eV does not overlap data energy range "
            f"{data_lo:.6g}-{data_hi:.6g} eV."
        )
    if lo < data_lo or hi > data_hi:
        return (
            f"{name} {lo:.6g}-{hi:.6g} eV only partially overlaps data energy range "
            f"{data_lo:.6g}-{data_hi:.6g} eV."
        )
    return None


def _channel_validation_messages(entry: dict, key: str, label: str, require_positive: bool = False) -> tuple[list[str], list[str]]:
    fatal_errors = []
    warnings_out = []
    values = entry.get(key)
    filename = entry.get("filename", "unknown")
    if values is None:
        fatal_errors.append(f"{filename}: required channel {label} is missing.")
        return fatal_errors, warnings_out

    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        warnings_out.append(f"{filename}: required channel {label} contains no finite values.")
        return fatal_errors, warnings_out

    if np.all(finite == 0):
        warnings_out.append(f"{filename}: required channel {label} is all zeros.")
    if require_positive and np.any(finite <= 0):
        warnings_out.append(f"{filename}: required channel {label} contains non-positive values used in division/log calculations.")

    if finite.size >= 2:
        span = float(np.nanmax(finite) - np.nanmin(finite))
        scale = max(float(np.nanmax(np.abs(finite))), 1e-30)
        if span < scale * 1e-9:
            warnings_out.append(f"{filename}: required channel {label} is nearly flat.")
    return fatal_errors, warnings_out


def _alignment_structure_warning(entry: dict, signal_key: str, signal_label: str, config: AstraConfig) -> str | None:
    energy = np.asarray(entry.get("energy"), dtype=float)
    signal = entry.get(signal_key)
    if signal is None:
        return f"{entry.get('filename', 'unknown')}: alignment signal {signal_label} is missing."
    signal = np.asarray(signal, dtype=float)
    lo, hi = config.align_window
    mask = (energy >= lo) & (energy <= hi) & np.isfinite(energy) & np.isfinite(signal)
    if np.count_nonzero(mask) < 10:
        return (
            f"{entry.get('filename', 'unknown')}: alignment signal {signal_label} has fewer than "
            f"10 finite points in alignment window {lo:g}-{hi:g} eV."
        )
    ew = energy[mask]
    sw = signal[mask]
    order = np.argsort(ew)
    ew = ew[order]
    sw = sw[order]
    try:
        derivative = np.gradient(sw, ew)
    except Exception as exc:
        return f"{entry.get('filename', 'unknown')}: alignment signal {signal_label} derivative could not be evaluated: {exc}"
    signal_scale = np.nanmax(sw) - np.nanmin(sw)
    derivative_range = np.nanmax(derivative) - np.nanmin(derivative)
    floor = max(signal_scale * 1e-6, 1e-15)
    if (
        not np.isfinite(signal_scale)
        or not np.isfinite(derivative_range)
        or signal_scale < 1e-15
        or derivative_range < floor
    ):
        return (
            f"{entry.get('filename', 'unknown')}: selected alignment signal {signal_label} has weak "
            f"structure in alignment window {lo:g}-{hi:g} eV; alignment may be unreliable."
        )
    return None


# NOTE: Mirrors the per-scan portion of _validate_processing_inputs() in
# processor.py for the single-scan case. The batch validator stays in
# processor.py because it aggregates energies across multiple entries.
# A future phase may split _validate_processing_inputs into per-scan and
# per-batch primitives. Tracked in ROADMAP.md Phase 2+.
def _validate_single_scan(entry: dict, config: AstraConfig) -> tuple[list[str], list[str]]:
    """Per-scan validation. Returns (warnings, fatal_errors).

    Designed to match _validate_processing_inputs([entry], config) exactly
    for representative inputs; covered by an equivalence test.
    """
    warnings_out: list[str] = []
    fatal_errors: list[str] = []

    energy = entry.get("energy")
    finite_energy: np.ndarray | None = None
    if energy is not None:
        arr = np.asarray(energy, dtype=float)
        if np.isfinite(arr).any():
            finite_energy = arr[np.isfinite(arr)]

    if finite_energy is not None and finite_energy.size > 0:
        data_range = (float(np.nanmin(finite_energy)), float(np.nanmax(finite_energy)))
        ranges = [
            ("Plot range", (config.plot_energy_min, config.plot_energy_max)),
            ("Alignment window", config.align_window),
            ("Pre-edge range", (config.e0 + config.pre1, config.e0 + config.pre2)),
            ("Normalization range", (config.e0 + config.norm1, config.e0 + config.norm2)),
        ]
        for name, selected_range in ranges:
            message = _range_overlap_status(name, selected_range, data_range)
            if message is not None:
                warnings_out.append(message)
    else:
        warnings_out.append("No finite energy values found for validation.")

    mode = getattr(config, "analysis_mode", "fluo")
    required_channels = _required_channels_for_signal(mode)
    positive_channels_analysis = {"I0", "I1", "I2"} if mode in {"trans", "ref"} else {"I0"}
    if not entry.get("is_foil"):
        for key, label in required_channels:
            fatal, warn = _channel_validation_messages(
                entry, key, label, require_positive=key in positive_channels_analysis
            )
            fatal_errors.extend(fatal)
            warnings_out.extend(warn)

    alignment_source = getattr(config, "alignment_source", "separate_foil")
    if alignment_source == "inline_ref":
        # Mirrors the batch validator: in inline_ref mode the alignment-structure
        # check fires on every non-foil (sample) entry.
        if not entry.get("is_foil"):
            for key, label in _required_channels_for_signal("ref"):
                fatal, warn = _channel_validation_messages(
                    entry, key, label, require_positive=True
                )
                fatal_errors.extend(fatal)
                warnings_out.extend(warn)
            message = _alignment_structure_warning(entry, "mu_ref", "ln(I1/I2)", config)
            if message is not None:
                warnings_out.append(message)
    else:
        # separate_foil: alignment-structure check fires only on foil entries.
        if entry.get("is_foil"):
            alignment_mode = getattr(config, "foil_alignment_mode", "trans")
            signal_key, signal_label = {
                "trans": ("mu_trans", "ln(I0/I1)"),
                "ref": ("mu_ref", "ln(I1/I2)"),
                "fluo": ("mu_fluo", "IF/I0"),
            }.get(alignment_mode, ("mu_trans", "ln(I0/I1)"))
            positive_channels_align = (
                {"I0", "I1", "I2"} if alignment_mode in {"trans", "ref"} else {"I0"}
            )
            for key, label in _required_channels_for_signal(alignment_mode):
                fatal, warn = _channel_validation_messages(
                    entry, key, label, require_positive=key in positive_channels_align
                )
                fatal_errors.extend(fatal)
                warnings_out.extend(warn)
            message = _alignment_structure_warning(entry, signal_key, signal_label, config)
            if message is not None:
                warnings_out.append(message)

    warnings_out = list(dict.fromkeys(warnings_out))
    fatal_errors = list(dict.fromkeys(fatal_errors))
    return warnings_out, fatal_errors


# ---------------------------------------------------------------------------
# Detector jumps (single-channel + per-entry aggregator)
# ---------------------------------------------------------------------------

def detect_detector_jumps(
    energy: np.ndarray,
    channel: np.ndarray,
    channel_name: str,
    config: AstraConfig,
    filename: str,
) -> list[dict]:
    energy = np.asarray(energy, dtype=float)
    channel = np.asarray(channel, dtype=float)
    finite_mask = np.isfinite(energy) & np.isfinite(channel)
    energy_c = energy[finite_mask]
    channel_c = channel[finite_mask]
    if len(energy_c) < 20:
        return []
    if np.nanstd(channel_c) == 0:
        return []

    dch = np.abs(np.diff(channel_c))
    if len(dch) == 0:
        return []

    threshold = _config_float(config, "detector_jump_threshold", 10.0)
    min_relative = _config_float(config, "detector_jump_min_relative", 0.05)
    noise_mad = np.median(np.abs(dch - np.median(dch)))
    severity_noise_mad = max(float(noise_mad), 1e-12)
    if noise_mad < 1e-12:
        noise_mad = np.std(dch) + 1e-12
    if not np.isfinite(noise_mad) or noise_mad <= 0:
        return []

    records = []
    candidate_indices = np.where(dch > threshold * noise_mad)[0]
    for i in candidate_indices:
        i = int(i)
        look_ahead = min(5, len(dch) - i - 1)
        recovery_window = dch[i + 1 : i + 1 + look_ahead]

        is_spike = len(recovery_window) == 0
        ch_diff_at_i = channel_c[i + 1] - channel_c[i]
        if not is_spike:
            for j in range(len(recovery_window)):
                ch_diff_at_j = channel_c[i + 2 + j] - channel_c[i + 1 + j]
                if (
                    ch_diff_at_j * ch_diff_at_i < 0
                    and abs(ch_diff_at_j) > (threshold / 3.0) * noise_mad
                ):
                    is_spike = True
                    break
        if not is_spike:
            continue

        signal_at_i = max(abs(channel_c[i]), abs(channel_c[i + 1]), 1e-12)
        relative_jump = dch[i] / signal_at_i
        if relative_jump < min_relative:
            continue

        ratio = dch[i] / severity_noise_mad
        if ratio > 5.0 * threshold:
            severity = "high"
        elif ratio > 2.0 * threshold:
            severity = "medium"
        else:
            severity = "low"

        energy_at_jump = float(energy_c[i])
        e0 = config.e0
        inside_plot_window = config.plot_energy_min <= energy_at_jump <= config.plot_energy_max
        inside_alignment_window = config.align_window_min <= energy_at_jump <= config.align_window_max
        inside_preedge_window = (e0 + config.pre1) <= energy_at_jump <= (e0 + config.pre2)
        inside_norm_window = (e0 + config.norm1) <= energy_at_jump <= (e0 + config.norm2)

        notes = []
        if inside_alignment_window:
            notes.append("in alignment window")
        if inside_preedge_window:
            notes.append("in pre-edge window")
        if inside_norm_window:
            notes.append("in norm window")
        if channel_name in ("IF_over_I0", "ln_I0_I1", "ln_I1_I2"):
            notes.append("derived signal: edge region jumps may be real features")

        records.append({
            "_index": i,
            "filename": filename,
            "channel": channel_name,
            "energy_eV": energy_at_jump,
            "jump_size": float(dch[i]),
            "relative_jump": float(relative_jump),
            "severity": severity,
            "inside_plot_window": bool(inside_plot_window),
            "inside_alignment_window": bool(inside_alignment_window),
            "inside_preedge_window": bool(inside_preedge_window),
            "inside_norm_window": bool(inside_norm_window),
            "note": "; ".join(notes) if notes else "",
        })

    deduped = []
    for record in records:
        if not deduped or record["_index"] - deduped[-1]["_index"] > 3:
            deduped.append(record)
        elif record["jump_size"] > deduped[-1]["jump_size"]:
            deduped[-1] = record
    for record in deduped:
        record.pop("_index", None)
    return deduped


# NOTE: Similar logic exists in astra_xas.beamtime.watcher._add_derived_channels.
# Phase 2.2 will migrate watcher.py to use single_scan.py and remove the duplicate.
def _add_derived_channels(channels: dict) -> None:
    try:
        i0 = channels.get("I0")
        i1 = channels.get("I1")
        i2 = channels.get("I2")
        if_ = channels.get("IF")
        if i0 is not None and if_ is not None:
            i0 = np.asarray(i0, dtype=float)
            if_ = np.asarray(if_, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                channels["IF_over_I0"] = np.where(i0 > 0, if_ / i0, np.nan)
        if i0 is not None and i1 is not None:
            i0 = np.asarray(i0, dtype=float)
            i1 = np.asarray(i1, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                channels["ln_I0_I1"] = np.where(i1 > 0, np.log(i0 / i1), np.nan)
        if i1 is not None and i2 is not None:
            i1 = np.asarray(i1, dtype=float)
            i2 = np.asarray(i2, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                channels["ln_I1_I2"] = np.where(i2 > 0, np.log(i1 / i2), np.nan)
    except Exception:
        return


def _detect_jumps_for_entry(entry: dict, config: AstraConfig) -> list[dict]:
    """Run detect_detector_jumps over raw + derived channels for one entry.

    Returns a flat list of all jump records (raw channels first, then
    derived). Channels that are missing or contain no finite values are
    skipped silently — matching the existing per-channel behaviour in
    process_folder and watcher._count_detector_jumps.
    """
    raw_channels = ("I0", "I1", "I2", "IF", "FDT")
    channels: dict = {name: entry.get(name) for name in raw_channels}
    _add_derived_channels(channels)

    energy = entry.get("energy")
    filename = entry.get("filename", "unknown")
    all_records: list[dict] = []
    for name, values in channels.items():
        if values is None:
            continue
        try:
            arr = np.asarray(values, dtype=float)
        except Exception:
            continue
        if not np.isfinite(arr).any():
            continue
        records = detect_detector_jumps(energy, arr, name, config, filename)
        all_records.extend(records)
    return all_records


# ---------------------------------------------------------------------------
# Public per-scan API
# ---------------------------------------------------------------------------

@dataclass
class SingleScanResult:
    filename: str
    is_foil: bool
    entry: dict
    energy: np.ndarray
    analysis_signal: np.ndarray
    analysis_signal_label: str
    qc_status: str              # "ok" | "warn" | "reject"
    qc_warnings: list[str]
    qc_errors: list[str]
    detector_jumps: list[dict]
    metrics: dict


def process_single_scan(
    scan: dict,
    config: AstraConfig,
    filename: str | None = None,
) -> SingleScanResult:
    """Per-scan validation, QC, and signal computation.

    No file I/O. No alignment (relational; happens in batch). No grouping.
    No normalization. No mutation of caller-owned ``scan`` data: if a
    ``filename`` override is supplied, a shallow copy of ``scan`` is used.

    ``qc_status`` rule:
      - ``"reject"`` if any ``qc_errors``
      - else ``"warn"`` if any ``qc_warnings``
      - else ``"ok"``
    Detector jumps are reported in ``detector_jumps`` and
    ``metrics["n_detector_jumps"]`` but do NOT influence ``qc_status``.

    ``is_foil`` is derived from the filename via ``config.foil_keyword``.
    Manifest-supplied ``is_foil`` overrides are applied by callers after
    this function returns (existing Phase 1.2 behaviour, unchanged).
    """
    if filename is not None:
        scan = {**scan, "filename": filename}
    elif "filename" not in scan:
        raise ValueError(
            "scan dict missing 'filename' key; pass filename= explicitly"
        )

    entry = _entry_from_scan(scan, config)
    qc_warnings, qc_errors = _validate_single_scan(entry, config)
    detector_jumps = _detect_jumps_for_entry(entry, config)

    signal_key, signal_label = _analysis_signal_spec(config)
    energy_value = entry.get("energy")
    energy_arr = (
        np.asarray(energy_value, dtype=float)
        if energy_value is not None
        else np.array([], dtype=float)
    )
    signal_value = entry.get(signal_key)
    signal_arr = (
        np.asarray(signal_value, dtype=float)
        if signal_value is not None
        else np.array([], dtype=float)
    )

    if qc_errors:
        qc_status = "reject"
    elif qc_warnings:
        qc_status = "warn"
    else:
        qc_status = "ok"

    n_points = int(energy_arr.size)
    finite_energy_mask = (
        np.isfinite(energy_arr) if energy_arr.size > 0 else np.array([], dtype=bool)
    )
    n_finite_energy = int(np.count_nonzero(finite_energy_mask))
    if n_finite_energy > 0:
        finite_energy_values = energy_arr[finite_energy_mask]
        energy_min: float | None = float(finite_energy_values.min())
        energy_max: float | None = float(finite_energy_values.max())
        energy_range: list[float] | None = [energy_min, energy_max]
    else:
        energy_min = None
        energy_max = None
        energy_range = None

    channels_present: list[str] = []
    for ch in ("I0", "I1", "I2", "IF", "FDT", "Ir"):
        values = entry.get(ch)
        if values is None:
            continue
        try:
            arr = np.asarray(values, dtype=float)
        except (TypeError, ValueError):
            continue
        if np.isfinite(arr).any():
            channels_present.append(ch)

    if signal_arr.size > 0:
        analysis_finite_fraction = float(
            np.count_nonzero(np.isfinite(signal_arr)) / signal_arr.size
        )
    else:
        analysis_finite_fraction = 0.0

    metrics = {
        "n_points": n_points,
        "n_points_finite_energy": n_finite_energy,
        "energy_min": energy_min,
        "energy_max": energy_max,
        "energy_range_eV": energy_range,
        "channels_present": channels_present,
        "analysis_signal_finite_fraction": analysis_finite_fraction,
        "n_detector_jumps": int(len(detector_jumps)),
        "n_validation_warnings": int(len(qc_warnings)),
        "n_validation_errors": int(len(qc_errors)),
    }

    return SingleScanResult(
        filename=str(entry.get("filename", "")),
        is_foil=bool(entry.get("is_foil", False)),
        entry=entry,
        energy=energy_arr,
        analysis_signal=signal_arr,
        analysis_signal_label=signal_label,
        qc_status=qc_status,
        qc_warnings=qc_warnings,
        qc_errors=qc_errors,
        detector_jumps=detector_jumps,
        metrics=metrics,
    )
