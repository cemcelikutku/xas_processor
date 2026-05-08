from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from larch import Group
from larch.xafs import pre_edge

from .config import AstraConfig
from .io import collect_xasd_files, load_xasd, split_replicate_suffix, sanitize_name, build_default_output_directory
from .signals import compute_signals, get_signal
from .alignment import find_best_shift
from .grouping import group_samples
from .export import save_two_col
from .plotting import plot_overview, plot_replicate_qc, plot_drift, plot_detector_health_overview, plot_analysis_signal_qc


AUTO_DEGLITCH_WARNING = (
    "Automatic deglitching is intended for narrow point-like spikes. "
    "Use manual range deglitching for broad artifacts."
)

SHIFT_CONVENTION = (
    "positive shift_eV means +shift_eV is added to the scan energy before "
    "interpolation/averaging; shifts are relative to the alignment anchor."
)
PRIMARY_DETECTOR_JUMP_CHANNELS = {"I0", "I1", "I2", "IF"}
FDT_DETECTOR_JUMP_CHANNELS = {"FDT"}
RAW_DETECTOR_JUMP_CHANNELS = PRIMARY_DETECTOR_JUMP_CHANNELS | FDT_DETECTOR_JUMP_CHANNELS
DERIVED_DETECTOR_JUMP_CHANNELS = {"IF_over_I0", "ln_I0_I1", "ln_I1_I2"}


def interpolate_to_grid(E_source, mu_source, E_target, kind="linear"):
    f = interp1d(E_source, mu_source, kind=kind, bounds_error=False, fill_value=np.nan)
    return f(E_target)


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


def _relative_output_path(path: str | Path, output_dir: Path) -> str:
    path = Path(path)
    try:
        return path.relative_to(output_dir).as_posix()
    except ValueError:
        return path.as_posix()


def _write_alignment_metadata_header(
    f,
    alignment_source: str,
    alignment_signal_label: str,
    alignment_anchor_info: dict,
) -> None:
    f.write(f"# Alignment source: {alignment_source}\n")
    f.write(f"# Alignment signal: {alignment_signal_label}\n")
    f.write(f"# Alignment anchor mode: {alignment_anchor_info['mode']}\n")
    f.write(f"# Alignment anchor file: {alignment_anchor_info.get('path') or 'N/A'}\n")
    f.write(f"# Alignment anchor status: {alignment_anchor_info['status']}\n")
    f.write("# Absolute calibration: not guaranteed unless the alignment anchor was externally calibrated.\n")
    f.write(f"# Shift convention: {SHIFT_CONVENTION}\n")


def _as_clean_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        out = []
        for chunk in value.replace(",", "\n").splitlines():
            item = chunk.strip()
            if item:
                out.append(item)
        return out
    try:
        return [str(item).strip() for item in value if str(item).strip()]
    except TypeError:
        item = str(value).strip()
        return [item] if item else []


def _records_have_channel(records: list[dict], key: str) -> bool:
    for rec in records:
        values = rec.get(key)
        if values is None:
            continue
        values = np.asarray(values, dtype=float)
        if np.isfinite(values).any():
            return True
    return False


def _detector_health_channels(config: AstraConfig, records: list[dict]) -> list[tuple[str, str]]:
    mode = getattr(config, "analysis_mode", "fluo")
    if mode == "trans":
        channels = [("I0", "I0"), ("I1", "I1")]
    elif mode == "ref":
        channels = [("I1", "I1"), ("I2", "I2"), ("mu_ref", "ln(I1/I2)")]
    else:
        channels = [("I0", "I0"), ("IF", "IF"), ("FDT", "FDT")]

    if mode != "ref" and _records_have_channel(records, "mu_ref"):
        channels.append(("mu_ref", "ln(I1/I2)"))
    return channels


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


def _validate_processing_inputs(entries: list[dict], config: AstraConfig) -> tuple[list[str], list[str]]:
    warnings_out: list[str] = []
    fatal_errors: list[str] = []

    energies = [
        np.asarray(e.get("energy"), dtype=float)
        for e in entries
        if e.get("energy") is not None and np.isfinite(np.asarray(e.get("energy"), dtype=float)).any()
    ]
    if energies:
        finite_energy = np.concatenate([arr[np.isfinite(arr)] for arr in energies])
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
    positive_channels = {"I0", "I1", "I2"} if mode in {"trans", "ref"} else {"I0"}
    for entry in entries:
        if entry.get("is_foil"):
            continue
        for key, label in required_channels:
            fatal, warn = _channel_validation_messages(entry, key, label, require_positive=key in positive_channels)
            fatal_errors.extend(fatal)
            warnings_out.extend(warn)

    alignment_source = getattr(config, "alignment_source", "separate_foil")
    if alignment_source == "inline_ref":
        alignment_entries = [e for e in entries if not e.get("is_foil")]
        for entry in alignment_entries:
            for key, label in _required_channels_for_signal("ref"):
                fatal, warn = _channel_validation_messages(entry, key, label, require_positive=True)
                fatal_errors.extend(fatal)
                warnings_out.extend(warn)
            message = _alignment_structure_warning(entry, "mu_ref", "ln(I1/I2)", config)
            if message is not None:
                warnings_out.append(message)
    else:
        foil_entries = [e for e in entries if e.get("is_foil")]
        alignment_mode = getattr(config, "foil_alignment_mode", "trans")
        signal_key, signal_label = {
            "trans": ("mu_trans", "ln(I0/I1)"),
            "ref": ("mu_ref", "ln(I1/I2)"),
            "fluo": ("mu_fluo", "IF/I0"),
        }.get(alignment_mode, ("mu_trans", "ln(I0/I1)"))
        positive_channels = {"I0", "I1", "I2"} if alignment_mode in {"trans", "ref"} else {"I0"}
        for entry in foil_entries:
            for key, label in _required_channels_for_signal(alignment_mode):
                fatal, warn = _channel_validation_messages(entry, key, label, require_positive=key in positive_channels)
                fatal_errors.extend(fatal)
                warnings_out.extend(warn)
            message = _alignment_structure_warning(entry, signal_key, signal_label, config)
            if message is not None:
                warnings_out.append(message)

    # Preserve order while removing duplicates from overlapping analysis/alignment checks.
    warnings_out = list(dict.fromkeys(warnings_out))
    fatal_errors = list(dict.fromkeys(fatal_errors))
    return warnings_out, fatal_errors


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


def _detector_jump_sort_key(record: dict):
    severity_rank = {"high": 0, "medium": 1, "low": 2}
    return (
        severity_rank.get(record.get("severity"), 99),
        str(record.get("filename", "")),
        float(record.get("energy_eV", 0.0)),
    )


def _annotate_detector_jump_summary_inclusion(records: list[dict]) -> None:
    for record in records:
        channel = record.get("channel", "")
        severity = record.get("severity", "low")
        inside_alignment = bool(record.get("inside_alignment_window"))
        relative_jump = float(record.get("relative_jump", 0.0))
        include = True
        reasons = []

        if channel in DERIVED_DETECTOR_JUMP_CHANNELS:
            include = False
            reasons.append("excluded from main summary: derived signal")
        elif channel in FDT_DETECTOR_JUMP_CHANNELS:
            include = False
            reasons.append("reported separately as FDT diagnostic spike")
        elif channel not in RAW_DETECTOR_JUMP_CHANNELS:
            include = False
            reasons.append("excluded from main summary: non-standard channel")

        if include and inside_alignment:
            include = severity == "high" and relative_jump >= 0.5
            if not include:
                reasons.append("excluded from main summary: edge/alignment window")

        record["include_in_summary"] = bool(include)
        if reasons:
            note = record.get("note", "")
            reason_note = "; ".join(reasons)
            record["note"] = f"{note}; {reason_note}" if note else reason_note


def _write_detector_jumps(path: Path, records: list[dict], config: AstraConfig) -> None:
    threshold = _config_float(config, "detector_jump_threshold", 10.0)
    with path.open("w", encoding="utf-8") as f:
        f.write("# ASTRA detector jump diagnostics\n")
        f.write("# Generated by AstraXAS\n")
        f.write("# Jump detected using point-to-point MAD spike detector\n")
        f.write(
            f"# Severity: low={threshold:g}-{2 * threshold:g}x MAD, "
            f"medium={2 * threshold:g}-{5 * threshold:g}x MAD, "
            f"high=>{5 * threshold:g}x MAD\n"
        )
        f.write(f"# Threshold multiplier: {threshold}\n")
        f.write(f"# Min relative jump: {_config_float(config, 'detector_jump_min_relative', 0.05)}\n")
        f.write("# Spike vs step discrimination: recovery window = 5 points\n")
        f.write("# include_in_summary: True only for significant primary raw-channel jumps emphasized in ASTRA_processing_report.txt\n")
        f.write("# NOTE: This file is diagnostic only. No data was modified.\n")
        f.write("#\n")
        f.write(
            "# filename\tchannel\tenergy_eV\tjump_size\trelative_jump\tseverity\t"
            "inside_plot_window\tinside_alignment_window\tinside_preedge_window\t"
            "inside_norm_window\tinclude_in_summary\tnote\n"
        )
        for record in sorted(records, key=_detector_jump_sort_key):
            f.write(
                f"{record['filename']}\t{record['channel']}\t{record['energy_eV']:.6f}\t"
                f"{record['jump_size']:.8g}\t{record['relative_jump']:.8g}\t"
                f"{record['severity']}\t{record['inside_plot_window']}\t"
                f"{record['inside_alignment_window']}\t{record['inside_preedge_window']}\t"
                f"{record['inside_norm_window']}\t{record.get('include_in_summary', False)}\t"
                f"{record['note']}\n"
            )


def _detector_jump_energy_regions(records: list[dict]) -> str:
    if not records:
        return "None"
    energies = sorted(float(record["energy_eV"]) for record in records)
    regions = []
    start = energies[0]
    prev = energies[0]
    count = 1
    for energy in energies[1:]:
        if energy - prev <= 5.0:
            prev = energy
            count += 1
        else:
            if count == 1:
                regions.append(f"{start:.1f} eV")
            else:
                regions.append(f"{start:.1f}-{prev:.1f} eV ({count})")
            start = prev = energy
            count = 1
    if count == 1:
        regions.append(f"{start:.1f} eV")
    else:
        regions.append(f"{start:.1f}-{prev:.1f} eV ({count})")
    return ", ".join(regions)


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


def _same_file(path_a, path_b) -> bool:
    if path_a is None or path_b is None:
        return False
    try:
        return Path(path_a).expanduser().resolve() == Path(path_b).expanduser().resolve()
    except OSError:
        return False


def _load_selected_alignment_anchor(config: AstraConfig, alignment_source: str) -> tuple[dict | None, dict]:
    anchor_mode = getattr(config, "alignment_anchor_mode", "first_scan")
    _, signal_key, signal_label = _alignment_signal_spec(alignment_source, config)
    info = {
        "mode": anchor_mode,
        "path": None,
        "status": "not requested",
        "signal": signal_label,
    }
    if anchor_mode == "first_scan":
        info["status"] = "using first scan automatically"
        return None, info
    if anchor_mode != "selected_file":
        raise RuntimeError("alignment_anchor_mode must be 'first_scan' or 'selected_file'.")

    anchor_path_value = getattr(config, "alignment_anchor_path", None)
    if not anchor_path_value:
        raise RuntimeError("Selected alignment anchor mode requires alignment_anchor_path.")
    anchor_path = Path(anchor_path_value).expanduser()
    if not anchor_path.exists():
        raise RuntimeError(f"Selected alignment anchor file does not exist: {anchor_path}")
    if not anchor_path.is_file():
        raise RuntimeError(f"Selected alignment anchor path is not a file: {anchor_path}")

    try:
        anchor_path = anchor_path.resolve()
        scan = load_xasd(anchor_path)
        anchor_entry = _entry_from_scan(scan, config, path=anchor_path)
    except Exception as exc:
        raise RuntimeError(f"Could not load selected alignment anchor file {anchor_path}: {exc}") from exc

    signal_mode, signal_key, signal_label = _alignment_signal_spec(alignment_source, config)
    positive_channels = {"I0", "I1", "I2"} if signal_mode in {"trans", "ref"} else {"I0"}
    validation_errors = []
    for key, label in _required_channels_for_signal(signal_mode):
        fatal, warn = _channel_validation_messages(anchor_entry, key, label, require_positive=key in positive_channels)
        validation_errors.extend(fatal)
        validation_errors.extend([message for message in warn if " is nearly flat." not in message])
    structure_warning = _alignment_structure_warning(anchor_entry, signal_key, signal_label, config)
    if structure_warning is not None:
        validation_errors.append(structure_warning)
    if validation_errors:
        message = "Selected alignment anchor validation failed:\n" + "\n".join(f"- {err}" for err in validation_errors)
        raise RuntimeError(message)

    info = {
        "mode": "selected_file",
        "path": str(anchor_path),
        "status": "loaded and validated",
        "signal": signal_label,
    }
    return anchor_entry, info


def _safe_attr(group: Group, name: str, default=None):
    try:
        return getattr(group, name)
    except AttributeError:
        return default


def _run_pre_edge(energy, mu, config: AstraConfig) -> dict:
    """Run Larch pre_edge with the current config and return useful arrays/metadata.

    This function is called on the already-merged (averaged) μ(E) spectrum, so that
    normalization is performed once on the clean averaged signal. This matches the
    recommended Athena workflow: merge μ(E) first, then normalize the merged spectrum.
    """
    energy = np.asarray(energy, dtype=float)
    mu = np.asarray(mu, dtype=float)
    g = Group()
    pre_edge(
        energy,
        mu,
        group=g,
        e0=config.e0,
        step=config.step,
        pre1=config.pre1,
        pre2=config.pre2,
        norm1=config.norm1,
        norm2=config.norm2,
        nnorm=config.nnorm,
        nvict=config.nvict,
        make_flat=config.make_flat,
    )

    pre_edge_fit = _safe_attr(g, "pre_edge", None)
    if pre_edge_fit is None:
        pre_edge_fit = np.zeros_like(mu, dtype=float)

    norm = np.asarray(_safe_attr(g, "norm", np.full_like(mu, np.nan)), dtype=float)
    flat = np.asarray(_safe_attr(g, "flat", norm), dtype=float)
    bkgcorr = np.asarray(mu, dtype=float) - np.asarray(pre_edge_fit, dtype=float)

    edge_step = _safe_attr(g, "edge_step", None)
    if edge_step is None:
        edge_step = _safe_attr(g, "norm_c0", None)
    try:
        edge_step = float(edge_step)
    except (TypeError, ValueError):
        edge_step = float("nan")

    e0 = _safe_attr(g, "e0", config.e0)
    try:
        e0 = float(e0)
    except (TypeError, ValueError):
        e0 = float(config.e0)

    return {
        "group": g,
        "bkgcorr": bkgcorr,
        "norm": norm,
        "flat": flat,
        "edge_step": edge_step,
        "e0": e0,
    }


def auto_deglitch(energy, mu, config):
    energy = np.asarray(energy, dtype=float)
    mu = np.asarray(mu, dtype=float)
    if not getattr(config, "enable_auto_deglitch", False):
        return energy, mu, []

    window_radius = int(getattr(config, "deglitch_window", 5))
    if window_radius < 2:
        window_radius = 2
    threshold = float(getattr(config, "deglitch_threshold", 5.0))
    min_energy = getattr(config, "deglitch_min_energy", None)
    max_energy = getattr(config, "deglitch_max_energy", None)

    n = len(energy)
    if n < 5 or len(mu) != n:
        return energy, mu, []

    flags = np.zeros(n, dtype=bool)

    local_point_residuals = np.full(n, np.nan, dtype=float)
    for i in range(1, n - 1):
        if not np.all(np.isfinite([energy[i - 1], energy[i], energy[i + 1], mu[i - 1], mu[i], mu[i + 1]])):
            continue
        if energy[i - 1] == energy[i + 1]:
            continue
        expected = np.interp(energy[i], [energy[i - 1], energy[i + 1]], [mu[i - 1], mu[i + 1]])
        local_point_residuals[i] = mu[i] - expected

    finite_residuals = local_point_residuals[np.isfinite(local_point_residuals)]
    if finite_residuals.size == 0:
        return energy, mu, []

    finite_mu = mu[np.isfinite(mu)]
    amplitude_floor = np.nanmax(np.abs(finite_mu)) * 1e-12 if finite_mu.size else 0.0
    eps = max(np.finfo(float).eps, amplitude_floor)

    for i in range(1, n - 1):
        if min_energy is not None and energy[i] < min_energy:
            continue
        if max_energy is not None and energy[i] > max_energy:
            continue
        if not np.isfinite(local_point_residuals[i]):
            continue

        start = max(0, i - window_radius)
        end = min(n, i + window_radius + 1)
        window_residuals = local_point_residuals[start:end]
        window_residuals = np.delete(window_residuals, i - start)
        window_residuals = window_residuals[np.isfinite(window_residuals)]
        if window_residuals.size < 3:
            continue

        local_med = np.median(window_residuals)
        local_mad = np.median(np.abs(window_residuals - local_med))
        scale = 1.4826 * local_mad
        if not np.isfinite(scale) or scale <= 0:
            scale = np.median(np.abs(window_residuals - local_med))
        if not np.isfinite(scale) or scale <= 0:
            local_mu = mu[start:end]
            local_mu = np.delete(local_mu, i - start)
            local_mu = local_mu[np.isfinite(local_mu)]
            local_range = np.ptp(local_mu) if local_mu.size else 0.0
            scale = max(local_range * 0.05, eps)

        left_jump = mu[i] - mu[i - 1]
        right_jump = mu[i + 1] - mu[i]
        opposite_jumps = left_jump * right_jump < 0
        returned_to_baseline = abs(mu[i + 1] - mu[i - 1]) <= max(3.0 * scale, abs(local_point_residuals[i]) * 0.5)
        large_residual = abs(local_point_residuals[i] - local_med) > threshold * scale
        large_jumps = min(abs(left_jump), abs(right_jump)) > threshold * scale * 0.5

        if large_residual and large_jumps and opposite_jumps and returned_to_baseline:
            flags[i] = True

    neighbour_flags = np.zeros(n, dtype=bool)
    neighbour_flags[1:] |= flags[:-1]
    neighbour_flags[:-1] |= flags[1:]
    flags &= ~neighbour_flags
    flags[0] = False
    flags[-1] = False

    flagged_energies = energy[flags].tolist()
    if not np.any(flags):
        return energy, mu, []

    good = ~flags
    if np.count_nonzero(good) < 2:
        return energy, mu, flagged_energies

    cleaned_mu = np.array(mu, copy=True)
    cleaned_mu[flags] = np.interp(energy[flags], energy[good], mu[good])

    return energy, cleaned_mu, flagged_energies


def manual_deglitch_range(energy, mu, config):
    """Replace mu values in a specified energy range by interpolation from nearby neighbours."""
    energy = np.asarray(energy, dtype=float)
    mu = np.asarray(mu, dtype=float)
    
    if not getattr(config, "enable_manual_deglitch_range", False):
        return energy, mu, []

    min_energy = getattr(config, "manual_deglitch_min_energy", None)
    max_energy = getattr(config, "manual_deglitch_max_energy", None)
    margin_points = int(getattr(config, "manual_deglitch_margin_points", 5))
    
    if margin_points < 1:
        margin_points = 1
    if min_energy is None or max_energy is None or min_energy >= max_energy:
        return energy, mu, []

    n = len(energy)
    if n == 0:
        return energy, mu, []

    # Identify points in the deglitch range
    in_range = (energy >= min_energy) & (energy <= max_energy)
    range_indices = np.where(in_range)[0]
    
    if len(range_indices) == 0:
        return energy, mu, []

    # Get leftmost and rightmost indices in range
    left_idx = range_indices[0]
    right_idx = range_indices[-1]

    # Find margin points: points outside the range to use for interpolation
    left_margin_start = max(0, left_idx - margin_points)
    right_margin_end = min(n, right_idx + margin_points + 1)

    # Collect neighbours for interpolation
    left_neighbours_mask = (np.arange(n) >= left_margin_start) & (np.arange(n) < left_idx)
    right_neighbours_mask = (np.arange(n) > right_idx) & (np.arange(n) < right_margin_end)
    neighbours_mask = left_neighbours_mask | right_neighbours_mask

    if not np.any(neighbours_mask):
        # No neighbours available, cannot interpolate
        return energy, mu, []

    neighbour_indices = np.where(neighbours_mask)[0]
    if len(neighbour_indices) < 2:
        # Not enough neighbours for meaningful interpolation
        return energy, mu, []

    cleaned_mu = np.array(mu, copy=True)
    manually_replaced = energy[in_range].tolist()

    try:
        # Try cubic spline
        spline = CubicSpline(energy[neighbours_mask], mu[neighbours_mask], extrapolate=False)
        interpolated = spline(energy[in_range])
        if np.any(np.isnan(interpolated)):
            raise ValueError("spline produced NaN values")
        cleaned_mu[in_range] = interpolated
    except Exception:
        # Fallback to linear interpolation
        cleaned_mu[in_range] = np.interp(energy[in_range], energy[neighbours_mask], mu[neighbours_mask])

    return energy, cleaned_mu, manually_replaced


def save_detector_raw_entry(entry: dict, output_dir: Path) -> Path:
    detector_dir = output_dir / "detector_raw"
    detector_dir.mkdir(parents=True, exist_ok=True)

    stem = Path(entry["filename"]).stem
    out = detector_dir / f"{sanitize_name(stem)}_detector_raw.dat"

    columns = [
        ("energy_eV", entry.get("energy")),
        ("I0", entry.get("I0")),
        ("I1", entry.get("I1")),
        ("I2", entry.get("I2")),
        ("IF", entry.get("IF")),
        ("FDT", entry.get("FDT")),
        ("Ir", entry.get("Ir")),
        ("mu_trans_lnI0I1", entry.get("mu_trans")),
        ("mu_ref_lnI1I2", entry.get("mu_ref")),
        ("mu_fluo_IFI0", entry.get("mu_fluo")),
    ]
    columns = [(name, values) for name, values in columns if values is not None]
    names = [name for name, _ in columns]
    arrays = [np.asarray(values, dtype=float) for _, values in columns]
    data = np.column_stack(arrays)
    header = " ".join(names)
    comments = (
        f"# ASTRA detector raw export\n"
        f"# Source file: {entry['filename']}\n"
        f"# Columns: {header}\n"
    )
    np.savetxt(out, data, header=comments + header, comments="")
    return out


def process_folder(input_dir: str | Path, output_dir: str | Path | None = None, config: AstraConfig | None = None, log=print):
    config = config or AstraConfig()
    input_dir = Path(input_dir).expanduser().resolve()
    if output_dir is None:
        output_dir = build_default_output_directory(input_dir)
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    plots_dir = output_dir / "plots"
    overview_dir = plots_dir / "overview"
    replicate_qc_dir = plots_dir / "replicate_qc"
    plots_enabled = (
        getattr(config, "save_detector_health_overview_plot", True)
        or getattr(config, "save_analysis_signal_qc_plot", True)
        or getattr(config, "save_detector_raw_overview_plot", False)
        or getattr(config, "save_processed_overview_plot", getattr(config, "save_raw_overview_plot", True))
        or getattr(config, "save_bkgcorr_overview_plot", False)
        or getattr(config, "save_norm_overview_plot", True)
        or getattr(config, "save_processed_mu_replicate_qc_plot", True)
        or getattr(config, "save_replicate_qc_plots", True)
        or getattr(config, "save_drift_plot", False)
    )
    if plots_enabled:
        plots_dir.mkdir(parents=True, exist_ok=True)
        overview_dir.mkdir(parents=True, exist_ok=True)
        if (
            getattr(config, "save_processed_mu_replicate_qc_plot", True)
            or getattr(config, "save_replicate_qc_plots", True)
        ):
            replicate_qc_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    validation_warnings: list[str] = []
    log(f"Starting {config.version}")
    log(f"Input directory: {input_dir}")
    log(f"Output directory: {output_dir}")

    files = collect_xasd_files(input_dir)
    if not files:
        raise RuntimeError(f"No .xasd files found in: {input_dir}")

    entries = []
    for path in files:
        scan = load_xasd(path)
        entries.append(_entry_from_scan(scan, config, path=path))

    detector_raw_files = []
    for e in entries:
        try:
            detector_raw_files.append(save_detector_raw_entry(e, output_dir))
        except Exception as exc:
            warning = f"Could not save detector raw for {e['filename']}: {exc}"
            warnings.append(warning)
            log("WARNING: " + warning)
    if detector_raw_files:
        log(f"Saved detector raw exports: {len(detector_raw_files)} file(s) in {output_dir / 'detector_raw'}")

    # Manual exclusions before alignment.
    exclude_filenames = set(_as_clean_list(getattr(config, "exclude_filenames", ())))
    exclude_contains = _as_clean_list(getattr(config, "exclude_filename_contains", ()))
    excluded_entries = []
    auto_shift_rejected_entries = []
    auto_outlier_entries = []

    if exclude_filenames or exclude_contains:
        kept_entries = []
        for e in entries:
            fname = e["filename"]
            exact_match = fname in exclude_filenames
            contains_match = any(pattern in fname for pattern in exclude_contains)
            if exact_match or contains_match:
                reason = "exact filename" if exact_match else "filename contains"
                e["exclude_reason"] = reason
                excluded_entries.append(e)
                log(f"Excluding scan: {fname} ({reason})")
            else:
                kept_entries.append(e)
        entries = kept_entries

    if not entries:
        raise RuntimeError("No .xasd files remain after applying manual exclusions.")

    validation_warnings, validation_errors = _validate_processing_inputs(entries, config)
    for warning in validation_warnings:
        log("VALIDATION WARNING: " + warning)
    if validation_errors:
        message = "Input validation failed:\n" + "\n".join(f"- {err}" for err in validation_errors)
        log("ERROR: " + message)
        raise RuntimeError(message)

    alignment_source = getattr(config, "alignment_source", "separate_foil")
    selected_anchor_entry, alignment_anchor_info = _load_selected_alignment_anchor(config, alignment_source)
    _, _, alignment_signal_label = _alignment_signal_spec(alignment_source, config)
    alignment_anchor_path = alignment_anchor_info.get("path")
    shift_records = []

    def append_shift_record(filename: str, shift_value: float, quality_value: float):
        shift_records.append({
            "filename": filename,
            "shift_eV": shift_value,
            "alignment_quality": quality_value,
            "scan_index": len(shift_records),
        })

    def warn_alignment_quality(filename: str, err: float, quality: float):
        if not np.isfinite(err) and quality == 0.0:
            warnings.append(
                f"Alignment skipped (unusable reference or moving spectrum): "
                f"{filename} — shift set to 0.0"
            )
        elif np.isfinite(err) and quality < config.alignment_quality_warn_threshold:
            warnings.append(
                f"Low alignment quality: {filename} "
                f"quality={quality:.3f} "
                f"(threshold={config.alignment_quality_warn_threshold})"
            )

    if alignment_source == "separate_foil":
        foil_entries = [e for e in entries if e["is_foil"]]
        if not foil_entries:
            raise RuntimeError(
                f"No foil files found. Put '{config.foil_keyword}' in foil filenames "
                f"or change foil_keyword."
            )
        first_foil = foil_entries[0]
        if selected_anchor_entry is None:
            E_foil_ref = first_foil["energy"]
            foil_ref_signal = get_signal(first_foil, config.foil_alignment_mode)
            foil_shift_map = {
                first_foil["filename"]: {
                    "shift_eV": 0.0,
                    "fit_error": 0.0,
                    "alignment_quality": 1.0,
                }
            }
            log(f"Reference foil: {first_foil['filename']} shift = 0.0000 eV")
            foils_to_align = foil_entries[1:]
        else:
            E_foil_ref = selected_anchor_entry["energy"]
            foil_ref_signal = get_signal(selected_anchor_entry, config.foil_alignment_mode)
            foil_shift_map = {}
            foils_to_align = foil_entries
            log(f"Selected alignment anchor: {alignment_anchor_path} ({alignment_signal_label})")

        for foil in foils_to_align:
            if selected_anchor_entry is not None and _same_file(foil.get("path"), alignment_anchor_path):
                shift, err, quality = 0.0, 0.0, 1.0
            else:
                shift, err, quality = find_best_shift(
                    E_foil_ref,
                    foil_ref_signal,
                    foil["energy"],
                    get_signal(foil, config.foil_alignment_mode),
                    window=config.align_window,
                    bounds=config.shift_bounds,
                    grid_points=config.alignment_grid_points,
                )
            foil_shift_map[foil["filename"]] = {
                "shift_eV": shift,
                "fit_error": err,
                "alignment_quality": quality,
            }
            if not (err == 0.0 and quality == 1.0):
                warn_alignment_quality(foil["filename"], err, quality)
            if abs(shift) > config.warn_shift_abs_eV:
                warnings.append(f"Large foil shift: {foil['filename']} = {shift:.4f} eV")
            log(f"Foil {foil['filename']}: shift = {shift:.4f} eV, fit_error = {err:.6g}, quality = {quality:.3f}")

        current_foil_name = first_foil["filename"]
        current_shift = foil_shift_map[current_foil_name]["shift_eV"]
        current_quality = foil_shift_map[current_foil_name]["alignment_quality"]
        sample_before_first_foil = True
        for e in entries:
            if e["is_foil"]:
                current_foil_name = e["filename"]
                current_shift = foil_shift_map[current_foil_name]["shift_eV"]
                current_quality = foil_shift_map[current_foil_name]["alignment_quality"]
                sample_before_first_foil = False
            elif sample_before_first_foil:
                warnings.append(f"Sample appears before first foil and will use first foil shift: {e['filename']}")
            e["assigned_foil"] = current_foil_name
            e["energy_shift_eV"] = current_shift
            e["alignment_quality"] = current_quality
            append_shift_record(e["filename"], current_shift, current_quality)
        sample_entries = [e for e in entries if not e["is_foil"]]

    elif alignment_source == "inline_ref":
        sample_entries = [e for e in entries if not e["is_foil"]]
        if not sample_entries:
            raise RuntimeError("No sample files found.")
        ref_scan = sample_entries[0]
        if selected_anchor_entry is None:
            E_ref = ref_scan["energy"]
            mu_ref_signal = ref_scan["mu_ref"]
            ref_name = f"inline_ref:{ref_scan['filename']}"
            log(f"Inline reference scan: {ref_scan['filename']} shift = 0.0000 eV")
        else:
            E_ref = selected_anchor_entry["energy"]
            mu_ref_signal = selected_anchor_entry["mu_ref"]
            ref_name = f"inline_ref_anchor:{Path(alignment_anchor_path).name}"
            log(f"Selected alignment anchor: {alignment_anchor_path} ({alignment_signal_label})")
        foil_entries = []
        foil_shift_map = {}
        if selected_anchor_entry is None:
            foil_shift_map[ref_scan["filename"]] = {
                "shift_eV": 0.0,
                "fit_error": 0.0,
                "alignment_quality": 1.0,
            }

        for e in sample_entries:
            if selected_anchor_entry is None and e is ref_scan:
                shift, err, quality = 0.0, 0.0, 1.0
            elif selected_anchor_entry is not None and _same_file(e.get("path"), alignment_anchor_path):
                shift, err, quality = 0.0, 0.0, 1.0
            else:
                shift, err, quality = find_best_shift(
                    E_ref,
                    mu_ref_signal,
                    e["energy"],
                    e["mu_ref"],
                    window=config.align_window,
                    bounds=config.shift_bounds,
                    grid_points=config.alignment_grid_points,
                )
            e["assigned_foil"] = ref_name
            e["energy_shift_eV"] = shift
            e["alignment_quality"] = quality
            foil_shift_map[e["filename"]] = {
                "shift_eV": shift,
                "fit_error": err,
                "alignment_quality": quality,
            }
            append_shift_record(e["filename"], shift, quality)
            if not (err == 0.0 and quality == 1.0):
                warn_alignment_quality(e["filename"], err, quality)
            if abs(shift) > config.warn_shift_abs_eV:
                warnings.append(f"Large inline reference shift: {e['filename']} = {shift:.4f} eV")
            log(f"Inline ref {e['filename']}: shift = {shift:.4f} eV, fit_error = {err:.6g}, quality = {quality:.3f}")
    else:
        raise RuntimeError(f"Unknown alignment_source='{alignment_source}'. Use 'separate_foil' or 'inline_ref'.")

    if not sample_entries:
        raise RuntimeError("No sample files found.")

    all_jump_records: list[dict] = []
    detector_jump_diagnostic = {"status": "disabled", "error": ""}
    if _config_bool(config, "enable_detector_jump_warnings", True):
        try:
            for scan in entries:
                if scan.get("is_foil"):
                    continue
                shifted_energy = np.asarray(scan.get("energy"), dtype=float) + float(scan.get("energy_shift_eV", 0.0))
                raw_channels_to_check = {
                    "I0": scan.get("I0"),
                    "I1": scan.get("I1"),
                    "I2": scan.get("I2"),
                    "IF": scan.get("IF"),
                    "FDT": scan.get("FDT"),
                }
                try:
                    i0 = scan.get("I0")
                    i1 = scan.get("I1")
                    i2 = scan.get("I2")
                    if_ = scan.get("IF")
                    if i0 is not None and if_ is not None:
                        i0 = np.asarray(i0, dtype=float)
                        if_ = np.asarray(if_, dtype=float)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            raw_channels_to_check["IF_over_I0"] = np.where(i0 > 0, if_ / i0, np.nan)
                    if i0 is not None and i1 is not None:
                        i0 = np.asarray(i0, dtype=float)
                        i1 = np.asarray(i1, dtype=float)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            raw_channels_to_check["ln_I0_I1"] = np.where(i1 > 0, np.log(i0 / i1), np.nan)
                    if i1 is not None and i2 is not None:
                        i1 = np.asarray(i1, dtype=float)
                        i2 = np.asarray(i2, dtype=float)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            raw_channels_to_check["ln_I1_I2"] = np.where(i2 > 0, np.log(i1 / i2), np.nan)
                except Exception:
                    pass

                for ch_name, ch_data in raw_channels_to_check.items():
                    if ch_data is None:
                        continue
                    try:
                        records = detect_detector_jumps(
                            shifted_energy,
                            ch_data,
                            ch_name,
                            config,
                            scan.get("filename", "unknown"),
                        )
                        all_jump_records.extend(records)
                    except Exception as exc:
                        detector_jump_diagnostic = {"status": "error", "error": str(exc)}
                        log(f"WARNING: Detector jump diagnostic failed for {scan.get('filename', 'unknown')} {ch_name}: {exc}")

            detector_jump_path = output_dir / "ASTRA_detector_jumps.dat"
            if all_jump_records:
                _annotate_detector_jump_summary_inclusion(all_jump_records)
                _write_detector_jumps(detector_jump_path, all_jump_records, config)
                detector_jump_diagnostic = {
                    "status": "found" if detector_jump_diagnostic["status"] != "error" else "error",
                    "error": detector_jump_diagnostic["error"],
                }
            else:
                if detector_jump_path.exists():
                    detector_jump_path.unlink()
                if detector_jump_diagnostic["status"] != "error":
                    detector_jump_diagnostic = {"status": "none", "error": ""}
        except Exception as exc:
            detector_jump_diagnostic = {"status": "error", "error": str(exc)}
            warnings.append(f"Detector jump diagnostics failed: {exc}")
            log("WARNING: Detector jump diagnostics failed: " + str(exc))
    else:
        try:
            detector_jump_path = output_dir / "ASTRA_detector_jumps.dat"
            if detector_jump_path.exists():
                detector_jump_path.unlink()
        except Exception:
            pass

    enable_shift_rejection = _config_bool(config, "enable_shift_rejection", False)
    reject_shift_abs_eV = _config_float(config, "reject_shift_abs_eV", 3.0)
    if enable_shift_rejection:
        kept_samples = []
        for s in sample_entries:
            shift = float(s.get("energy_shift_eV", 0.0))
            if abs(shift) > reject_shift_abs_eV:
                s["exclude_reason"] = f"auto shift rejection |shift|>{reject_shift_abs_eV:g} eV"
                auto_shift_rejected_entries.append(s)
                warnings.append(f"Auto-excluded by shift safety: {s['filename']} shift={shift:.4f} eV")
                log(f"Auto-excluding by shift safety: {s['filename']} shift={shift:.4f} eV")
            else:
                kept_samples.append(s)
        sample_entries = kept_samples
        if not sample_entries:
            raise RuntimeError("No sample files remain after shift safety rejection.")

    groups = group_samples(sample_entries)
    plot_files = []
    plot_energy_range = (config.plot_energy_min, config.plot_energy_max)
    detector_health_plot_info = {"path": None, "channels": [], "skipped": []}
    analysis_signal_qc_info = {"path": None, "signal": "", "n_traces": 0, "skipped": []}
    drift_plot_info = {"path": None, "status": "disabled", "reason": ""}

    if getattr(config, "save_detector_health_overview_plot", True):
        detector_health_records = [
            {
                "energy": s["energy"],
                "I0": s.get("I0"),
                "I1": s.get("I1"),
                "I2": s.get("I2"),
                "IF": s.get("IF"),
                "FDT": s.get("FDT"),
                "mu_ref": s.get("mu_ref"),
                "label": s["filename"],
            }
            for s in sample_entries
        ]
        detector_health_plot_info = plot_detector_health_overview(
            detector_health_records,
            overview_dir / "detector_health_overview.png",
            _detector_health_channels(config, detector_health_records),
            title="Detector health overview",
            energy_range=plot_energy_range,
        )
        if detector_health_plot_info["path"] is not None:
            plot_files.append(detector_health_plot_info["path"])
            log(f"Saved plot: {detector_health_plot_info['path']}")

    if getattr(config, "save_analysis_signal_qc_plot", True):
        signal_key, signal_label = _analysis_signal_spec(config)
        analysis_signal_records = [
            {
                "energy": s["energy"],
                "signal": s.get(signal_key),
                "label": s["filename"],
            }
            for s in sample_entries
        ]
        analysis_signal_qc_info = plot_analysis_signal_qc(
            analysis_signal_records,
            overview_dir / "analysis_signal_qc.png",
            signal_label,
            energy_range=plot_energy_range,
        )
        if analysis_signal_qc_info["path"] is not None:
            plot_files.append(analysis_signal_qc_info["path"])
            log(f"Saved plot: {analysis_signal_qc_info['path']}")

    if config.save_drift_plot:
        overview_dir.mkdir(parents=True, exist_ok=True)
        drift_path = overview_dir / "drift_tracker.png"
        try:
            plot_drift(
                shift_records,
                drift_path,
                warn_threshold_eV=config.warn_shift_abs_eV,
                quality_threshold=config.alignment_quality_warn_threshold,
            )
            plot_files.append(drift_path)
            drift_plot_info = {"path": drift_path, "status": "created", "reason": ""}
            log(f"Saved plot: {drift_path}")
        except Exception as exc:
            reason = str(exc)
            drift_plot_info = {"path": None, "status": "failed", "reason": reason}
            warning = f"Could not save drift tracker plot: {reason}"
            warnings.append(warning)
            log("WARNING: " + warning)

    with open(output_dir / "ASTRA_excluded_scans.dat", "w", encoding="utf-8") as f:
        f.write("# filename\tbase_name\treplicate_id\treason\n")
        if excluded_entries:
            for e in excluded_entries:
                rep = "None" if e["replicate_id"] is None else str(e["replicate_id"])
                reason = e.get("exclude_reason", "manual")
                f.write(f"{e['filename']}\t{e['base_name']}\t{rep}\t{reason}\n")
        else:
            f.write("# None\n")

    with open(output_dir / "ASTRA_auto_rejected_scans.dat", "w", encoding="utf-8") as f:
        f.write("# filename\tbase_name\treplicate_id\treason\n")
        if auto_shift_rejected_entries:
            for e in auto_shift_rejected_entries:
                rep = "None" if e.get("replicate_id") is None else str(e.get("replicate_id"))
                reason = e.get("exclude_reason", "auto")
                f.write(f"{e['filename']}\t{e['base_name']}\t{rep}\t{reason}\n")
        else:
            f.write("# None\n")

    with open(output_dir / "ASTRA_energy_shifts.dat", "w", encoding="utf-8") as f:
        _write_alignment_metadata_header(f, alignment_source, alignment_signal_label, alignment_anchor_info)
        f.write("# filename\tbase_name\treplicate_id\tassigned_foil\tenergy_shift_eV\talignment_quality\n")
        for s in sample_entries:
            rep = "None" if s["replicate_id"] is None else str(s["replicate_id"])
            f.write(
                f"{s['filename']}\t{s['base_name']}\t{rep}\t{s['assigned_foil']}\t"
                f"{s['energy_shift_eV']:.6f}\t{s['alignment_quality']:.6f}\n"
            )

    with open(output_dir / "ASTRA_foil_alignment.dat", "w", encoding="utf-8") as f:
        _write_alignment_metadata_header(f, alignment_source, alignment_signal_label, alignment_anchor_info)
        f.write("# alignment_quality: Pearson r between z-scored derivatives at best shift.\n")
        f.write("#   Near 1.0 = reliable. Below 0.7 = suspect.\n")
        f.write("#   not-finite fit_error with quality=0.0 means alignment was skipped.\n")
        f.write("# foil_filename\tshift_eV\tfit_error\talignment_quality\n")
        for foil_name, info in foil_shift_map.items():
            f.write(
                f"{foil_name}\t{info['shift_eV']:.6f}\t{info['fit_error']:.8g}\t"
                f"{info['alignment_quality']:.6f}\n"
            )

    group_summary_rows = []
    group_norm_metadata_rows = []
    detector_raw_plot_records = []
    processed_plot_records = []
    bkgcorr_plot_records = []
    norm_plot_records = []
    total_auto_deglitched_points = 0
    total_manual_range_interpolated_points = 0
    processed_mu_replicate_qc_created = 0
    processed_mu_replicate_qc_skipped: list[str] = []
    normalized_replicate_qc_created = 0
    normalized_replicate_qc_skipped: list[str] = []

    for (base_name, assigned_foil), scans in groups.items():
        log(f"Processing group: {base_name} ({len(scans)} scan(s))")
        first_scan = scans[0]
        master_energy = first_scan["energy"] + first_scan["energy_shift_eV"]
        processed_spectra = []
        source_filenames = []
        detector_group_if = []
        deglitch_records = []  # Collect deglitch info per replicate

        for s in scans:
            shifted_energy = s["energy"] + s["energy_shift_eV"]
            mu = get_signal(s, config.analysis_mode)
            shifted_energy, mu, flagged_energies = auto_deglitch(shifted_energy, mu, config)
            if flagged_energies:
                log(f"Auto-interpolated {len(flagged_energies)} isolated spike point(s) in {s['filename']}")
            
            shifted_energy, mu, manually_replaced = manual_deglitch_range(shifted_energy, mu, config)
            if manually_replaced:
                log(f"Manual range-interpolated {len(manually_replaced)} point(s) in {s['filename']}")
            
            deglitch_records.append({
                "filename": s["filename"],
                "n_flagged_auto": len(flagged_energies),
                "flagged_energies_auto": flagged_energies,
                "n_replaced_manual": len(manually_replaced),
                "replaced_energies_manual": manually_replaced,
            })
            mu_interp = interpolate_to_grid(shifted_energy, mu, master_energy, kind=config.interp_kind)
            valid = np.isfinite(mu_interp)
            if np.count_nonzero(valid) < 0.9 * len(master_energy):
                warnings.append(f"Insufficient interpolation overlap; skipped {s['filename']}")
                continue
            if not np.all(valid):
                mu_interp = np.interp(master_energy, master_energy[valid], mu_interp[valid])
            processed_spectra.append(mu_interp)
            source_filenames.append(s["filename"])

            if s.get("IF") is not None:
                if_interp = interpolate_to_grid(shifted_energy, s["IF"], master_energy, kind=config.interp_kind)
                if np.isfinite(if_interp).any():
                    detector_group_if.append(if_interp)

        if not processed_spectra:
            warnings.append(f"No valid spectra remained for group {base_name}")
            continue

        processed_array = np.asarray(processed_spectra, dtype=float)

        # Optional automatic outlier detection on processed signal before normalization.
        enable_auto_outlier_detection = _config_bool(config, "enable_auto_outlier_detection", False)
        outlier_threshold = _config_float(config, "outlier_rms_threshold", 0.08)
        keep_mask = np.ones(len(processed_array), dtype=bool)
        if enable_auto_outlier_detection and len(processed_array) >= 3:
            median_spectrum = np.nanmedian(processed_array, axis=0)
            scale = float(np.nanmax(median_spectrum) - np.nanmin(median_spectrum))
            if not np.isfinite(scale) or scale <= 0:
                scale = 1.0
            for idx, spectrum in enumerate(processed_array):
                score = float(np.sqrt(np.nanmean(((spectrum - median_spectrum) / scale) ** 2)))
                if score > outlier_threshold:
                    keep_mask[idx] = False
                    auto_outlier_entries.append({
                        "filename": source_filenames[idx],
                        "base_name": base_name,
                        "replicate_id": scans[idx].get("replicate_id") if idx < len(scans) else None,
                        "score": score,
                        "threshold": outlier_threshold,
                    })
                    warnings.append(f"Auto-outlier flagged and skipped: {source_filenames[idx]} RMS={score:.5f}")
                    log(f"Auto-outlier skipped: {source_filenames[idx]} RMS={score:.5f}")
            if not keep_mask.any():
                warnings.append(f"All scans in group {base_name} were flagged as outliers; keeping original group.")
                keep_mask[:] = True

        processed_array = processed_array[keep_mask]
        source_filenames = [name for name, keep in zip(source_filenames, keep_mask) if keep]
        deglitch_records = [rec for rec, keep in zip(deglitch_records, keep_mask) if keep]
        total_auto_deglitched_points += sum(r["n_flagged_auto"] for r in deglitch_records)
        total_manual_range_interpolated_points += sum(r["n_replaced_manual"] for r in deglitch_records)

        # Recommended Athena workflow: merge μ(E) first, then normalize once on the
        # clean averaged spectrum. This gives a more reliable edge step and background
        # estimate than normalizing noisy individual replicates and averaging the results.
        processed_avg = np.nanmean(processed_array, axis=0)

        pe = _run_pre_edge(master_energy, processed_avg, config)
        bkgcorr_avg = pe["bkgcorr"]
        norm_avg = pe["norm"]
        flat_avg = pe["flat"]
        mean_edge_step = pe["edge_step"]
        std_edge_step = float("nan")   # single merged spectrum — no per-replicate spread
        mean_norm_e0 = pe["e0"]

        # Still compute per-replicate normalized spectra for the QC plot only,
        # so users can visually inspect replicate consistency.
        norm_reps = []
        for spectrum in processed_array:
            pe_rep = _run_pre_edge(master_energy, spectrum, config)
            norm_reps.append(pe_rep["norm"])
        norm_array = np.asarray(norm_reps, dtype=float)

        output_base = sanitize_name(base_name)

        if getattr(config, "save_processed_mu_replicate_qc_plot", True):
            if len(processed_array) > 1:
                try:
                    processed_replicate_records = [
                        {"energy": master_energy, "mu": processed_array[i], "label": source_filenames[i]}
                        for i in range(len(source_filenames))
                    ]
                    path = plot_replicate_qc(
                        group_name=output_base,
                        replicate_records=processed_replicate_records,
                        average_record={"energy": master_energy, "mu": processed_avg, "label": "average"},
                        output_path=replicate_qc_dir / f"{output_base}_processed_mu_replicate_qc.png",
                        energy_range=plot_energy_range,
                        y_label="Processed μ(E) before normalization",
                        title=f"Processed μ(E) replicate QC: {output_base}",
                    )
                    if path is not None:
                        plot_files.append(path)
                        processed_mu_replicate_qc_created += 1
                        log(f"Saving processed μ(E) replicate QC plot: {output_base}")
                    else:
                        processed_mu_replicate_qc_skipped.append(f"{output_base}: no finite processed μ(E) traces in plot range")
                except Exception as exc:
                    warning = f"Could not save processed μ(E) replicate QC plot for {output_base}: {exc}"
                    warnings.append(warning)
                    processed_mu_replicate_qc_skipped.append(f"{output_base}: {exc}")
                    log("WARNING: " + warning)
            else:
                processed_mu_replicate_qc_skipped.append(f"{output_base}: fewer than 2 valid processed scans")

        comments = (
            f"# ASTRA XAS Processor version: {config.version}\n"
            f"# Group base name: {base_name}\n"
            f"# Assigned foil: {assigned_foil}\n"
            f"# Analysis mode: {config.analysis_mode}\n"
            f"# Alignment source: {alignment_source}\n"
            f"# Normalization order nnorm: {config.nnorm}\n"
            f"# Normalize-before-average: False (merge-then-normalize)\n"
            f"# Auto-deglitch enabled: {_config_bool(config, 'enable_auto_deglitch', False)}\n"
            f"# Deglitch method: interpolate\n"
            f"# Deglitch threshold: {getattr(config, 'deglitch_threshold', 'N/A')}\n"
            f"# Deglitch window: {getattr(config, 'deglitch_window', 'N/A')}\n"
            f"# {AUTO_DEGLITCH_WARNING}\n"
            f"# Manual deglitch range enabled: {_config_bool(config, 'enable_manual_deglitch_range', False)}\n"
            f"# Manual deglitch range: {getattr(config, 'manual_deglitch_min_energy', 'N/A')} - {getattr(config, 'manual_deglitch_max_energy', 'N/A')} eV\n"
            f"# Manual deglitch margin points: {getattr(config, 'manual_deglitch_margin_points', 'N/A')}\n"
            f"# Mean edge_step: {mean_edge_step:.8g}\n"
            f"# Mean normalization E0: {mean_norm_e0:.8g}\n"
            f"# Shift rejection enabled: {enable_shift_rejection} threshold: {reject_shift_abs_eV} eV\n"
            f"# Auto outlier detection enabled: {enable_auto_outlier_detection} threshold: {outlier_threshold}\n"
            f"# Fluorescence multiplicative constant: {config.fluo_multiplicative_constant}\n"
            f"# e0: {config.e0}\n"
            f"# pre1/pre2/norm1/norm2: {config.pre1} {config.pre2} {config.norm1} {config.norm2}\n"
            f"# Number of scans averaged: {len(processed_array)}\n"
            "# Source scans:\n" + "".join(f"#   {name}\n" for name in source_filenames)
        )

        # Processed μ(E)-like signal. This replaces the old misleading *_raw.dat name.
        # True detector channels are exported separately under detector_raw/*.dat.
        save_two_col(output_dir / f"{output_base}_processed.dat", master_energy, processed_avg, "energy_eV   mu", comments)
        save_two_col(output_dir / f"{output_base}_bkgcorr.dat", master_energy, bkgcorr_avg, "energy_eV   bkgcorr", comments)
        save_two_col(output_dir / f"{output_base}_norm.dat", master_energy, norm_avg, "energy_eV   norm", comments)
        save_two_col(output_dir / f"{output_base}_flat.dat", master_energy, flat_avg, "energy_eV   flat", comments)

        processed_plot_records.append({"energy": master_energy, "mu": processed_avg, "label": output_base})
        bkgcorr_plot_records.append({"energy": master_energy, "mu": bkgcorr_avg, "label": output_base})
        norm_plot_records.append({"energy": master_energy, "mu": norm_avg, "label": output_base})
        if detector_group_if:
            detector_raw_plot_records.append({"energy": master_energy, "mu": np.nanmean(np.asarray(detector_group_if), axis=0), "label": output_base})

        if getattr(config, "save_replicate_qc_plots", True):
            if len(norm_array) > 1:
                try:
                    normalized_replicate_records = [
                        {"energy": master_energy, "mu": norm_array[i], "label": source_filenames[i]}
                        for i in range(len(source_filenames))
                    ]
                    path = plot_replicate_qc(
                        group_name=output_base,
                        replicate_records=normalized_replicate_records,
                        average_record={"energy": master_energy, "mu": norm_avg, "label": "average"},
                        output_path=replicate_qc_dir / f"{output_base}_normalized_replicate_qc.png",
                        energy_range=plot_energy_range,
                        y_label="Normalized intensity",
                    )
                    if path is not None:
                        plot_files.append(path)
                        normalized_replicate_qc_created += 1
                        log(f"Saving normalized replicate QC plot: {output_base}")
                    else:
                        normalized_replicate_qc_skipped.append(f"{output_base}: no finite normalized traces in plot range")
                except Exception as exc:
                    warning = f"Could not save normalized replicate QC plot for {output_base}: {exc}"
                    warnings.append(warning)
                    normalized_replicate_qc_skipped.append(f"{output_base}: {exc}")
                    log("WARNING: " + warning)
            else:
                normalized_replicate_qc_skipped.append(f"{output_base}: fewer than 2 valid normalized scans")

        group_filenames = set(source_filenames)
        detector_jump_count = sum(
            1
            for record in all_jump_records
            if record.get("filename") in group_filenames and record.get("include_in_summary")
        )
        group_summary_rows.append([
            output_base,
            base_name,
            assigned_foil,
            str(len(processed_array)),
            source_filenames[0],
            source_filenames[-1],
            str(detector_jump_count),
        ])
        std_str = "N/A" if not np.isfinite(std_edge_step) else f"{std_edge_step:.8g}"
        group_norm_metadata_rows.append([output_base, f"{mean_edge_step:.8g}", std_str, f"{mean_norm_e0:.8g}", str(config.nnorm)])

        # Write deglitch log if any deglitching occurred
        auto_enabled = _config_bool(config, "enable_auto_deglitch", False)
        manual_enabled = _config_bool(config, "enable_manual_deglitch_range", False)
        any_auto_flags = any(r["n_flagged_auto"] > 0 for r in deglitch_records)
        any_manual_replaced = any(r["n_replaced_manual"] > 0 for r in deglitch_records)
        
        if (auto_enabled and any_auto_flags) or (manual_enabled and any_manual_replaced):
            with open(output_dir / f"{output_base}_deglitch_log.dat", "w", encoding="utf-8") as f:
                f.write("# Deglitch summary for group: {}\n".format(base_name))
                f.write(f"# {AUTO_DEGLITCH_WARNING}\n")
                f.write("# filename\tmethod\tn_processed\tenergies\n")
                for rec in deglitch_records:
                    # Write auto deglitch records
                    if rec["n_flagged_auto"] > 0:
                        flagged_str = " ".join(f"{e:.6f}" for e in rec["flagged_energies_auto"])
                        f.write(f"{rec['filename']}\tauto_interpolate\t{rec['n_flagged_auto']}\t{flagged_str}\n")
                    # Write manual deglitch records
                    if rec["n_replaced_manual"] > 0:
                        replaced_str = " ".join(f"{e:.6f}" for e in rec["replaced_energies_manual"])
                        f.write(f"{rec['filename']}\tmanual_range_interpolate\t{rec['n_replaced_manual']}\t{replaced_str}\n")

    if plots_enabled:
        log(f"Generating plots in: {plots_dir}")

    if getattr(config, "save_detector_raw_overview_plot", False):
        path = plot_overview(
            detector_raw_plot_records,
            overview_dir / "aligned_averaged_IF_overview.png",
            "Aligned averaged IF detector signal",
            "IF detector signal",
            energy_range=plot_energy_range,
        )
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_processed_overview_plot", getattr(config, "save_raw_overview_plot", True)):
        path = plot_overview(processed_plot_records, overview_dir / "processed_mu_overview.png", "Processed μ(E) overview", "Processed μ(E)", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_bkgcorr_overview_plot", False):
        path = plot_overview(bkgcorr_plot_records, overview_dir / "background_corrected_overview.png", "Background-corrected overview", "Background-corrected μ(E)", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_norm_overview_plot", True):
        path = plot_overview(norm_plot_records, overview_dir / "normalized_overview.png", "Normalized overview", "Normalized intensity", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    with open(output_dir / "ASTRA_auto_outliers.dat", "w", encoding="utf-8") as f:
        f.write("# filename\tbase_name\treplicate_id\trms_score\tthreshold\n")
        if auto_outlier_entries:
            for e in auto_outlier_entries:
                rep = "None" if e.get("replicate_id") is None else str(e.get("replicate_id"))
                f.write(f"{e['filename']}\t{e['base_name']}\t{rep}\t{e['score']:.8g}\t{e['threshold']:.8g}\n")
        else:
            f.write("# None\n")

    with open(output_dir / "ASTRA_group_summary.dat", "w", encoding="utf-8") as f:
        f.write("# output_base\tbase_name\tassigned_foil\tn_scans_used\tfirst_scan\tlast_scan\tdetector_jumps\n")
        f.write("# detector_jumps: number of summary-level detector jump records among scans used in the group; data are not modified.\n")
        for row in group_summary_rows:
            f.write("\t".join(row) + "\n")

    with open(output_dir / "ASTRA_normalization_summary.dat", "w", encoding="utf-8") as f:
        f.write("# output_base\tmean_edge_step\tstd_edge_step\tmean_e0\tnnorm\n")
        for row in group_norm_metadata_rows:
            f.write("\t".join(row) + "\n")

    low_quality_count = sum(
        1
        for info in foil_shift_map.values()
        if not (info["fit_error"] == 0.0 and info["alignment_quality"] == 1.0)
        and info["alignment_quality"] < config.alignment_quality_warn_threshold
    )
    relative_plot_files = [_relative_output_path(path, output_dir) for path in plot_files]
    overview_plot_files = [path for path in relative_plot_files if path.startswith("plots/overview/")]
    replicate_qc_plot_files = [path for path in relative_plot_files if path.startswith("plots/replicate_qc/")]
    other_plot_files = [
        path
        for path in relative_plot_files
        if not path.startswith("plots/overview/") and not path.startswith("plots/replicate_qc/")
    ]

    with open(output_dir / "ASTRA_processing_report.txt", "w", encoding="utf-8") as f:
        f.write(f"{config.version}\n")
        f.write(f"Input directory: {input_dir}\nOutput directory: {output_dir}\n")
        f.write(f"Files found: {len(files)}\n")
        f.write("\nValidation warnings:\n")
        if validation_warnings:
            for warning in validation_warnings:
                f.write(f"- {warning}\n")
        else:
            f.write("Validation warnings: none\n")
        f.write("\nProcessing summary:\n")
        f.write(f"Files excluded manually: {len(excluded_entries)}\n")
        f.write(f"Files rejected by shift safety: {len(auto_shift_rejected_entries)}\n")
        f.write(f"Replicate outliers skipped: {len(auto_outlier_entries)}\n")
        f.write(f"Foils found: {len(foil_entries)}\n")
        f.write(f"Groups processed: {len(group_summary_rows)}\n")
        f.write(f"Alignment source: {alignment_source}\n")
        f.write(f"Alignment signal: {alignment_signal_label}\n")
        f.write(f"Alignment anchor mode: {alignment_anchor_info['mode']}\n")
        f.write(f"Alignment anchor file: {alignment_anchor_info.get('path') or 'N/A'}\n")
        f.write(f"Alignment anchor status: {alignment_anchor_info['status']}\n")
        f.write("Absolute calibration: not guaranteed unless the alignment anchor was externally calibrated.\n")
        f.write(f"Shift convention: {SHIFT_CONVENTION}\n")
        f.write(
            f"Low-quality alignments (quality < {config.alignment_quality_warn_threshold}): "
            f"{low_quality_count}\n"
        )
        f.write("\nDetector jump diagnostics:\n")
        if detector_jump_diagnostic["status"] == "disabled":
            f.write("Detector jump diagnostics: disabled\n")
        elif detector_jump_diagnostic["status"] == "error":
            f.write(
                "Detector jump diagnostics: ERROR - diagnostic failed "
                f"({detector_jump_diagnostic['error']})\n"
            )
            f.write("process_folder completed normally; no data was affected\n")
        elif detector_jump_diagnostic["status"] == "none":
            f.write("Detector jump diagnostics: none detected\n")
        else:
            summary_records = [record for record in all_jump_records if record.get("include_in_summary")]
            channel_counts = {}
            severity_counts = {"high": 0, "medium": 0, "low": 0}
            derived_count = 0
            edge_alignment_excluded = 0
            fdt_count = 0
            for record in all_jump_records:
                channel = record.get("channel", "unknown")
                if channel in DERIVED_DETECTOR_JUMP_CHANNELS:
                    derived_count += 1
                elif channel in FDT_DETECTOR_JUMP_CHANNELS:
                    fdt_count += 1
                elif record.get("inside_alignment_window") and not record.get("include_in_summary"):
                    edge_alignment_excluded += 1
            for record in summary_records:
                channel = record.get("channel", "unknown")
                channel_counts[channel] = channel_counts.get(channel, 0) + 1
                severity = record.get("severity", "low")
                severity_counts[severity] = severity_counts.get(severity, 0) + 1
            f.write(
                "  Significant primary raw-channel jumps outside edge/alignment window: "
                f"{len(summary_records)}\n"
            )
            f.write(f"  FDT diagnostic spikes: {fdt_count} entries, reported separately\n")
            f.write(f"  Full diagnostic entries: {len(all_jump_records)}\n")
            f.write("  Full diagnostic entries written to: ASTRA_detector_jumps.dat\n")
            f.write(
                "  Entries excluded from main summary: "
                f"derived-signal features={derived_count}; "
                f"edge/alignment-window features={edge_alignment_excluded}; "
                f"FDT diagnostic spikes={fdt_count}\n"
            )
            f.write("  Main affected channels:\n")
            if channel_counts:
                for channel in sorted(channel_counts):
                    f.write(f"    {channel}: {channel_counts[channel]} jump(s)\n")
            else:
                f.write("    None\n")
            f.write(f"  Main affected energy regions: {_detector_jump_energy_regions(summary_records)}\n")
            f.write("  By severity:\n")
            f.write(f"    high: {severity_counts.get('high', 0)}\n")
            f.write(f"    medium: {severity_counts.get('medium', 0)}\n")
            f.write(f"    low: {severity_counts.get('low', 0)}\n")
            f.write("  Jumps inside critical windows:\n")
            f.write(
                "    alignment window: "
                f"{sum(1 for record in summary_records if record.get('inside_alignment_window'))}\n"
            )
            f.write(
                "    pre-edge window:  "
                f"{sum(1 for record in summary_records if record.get('inside_preedge_window'))}\n"
            )
            f.write(
                "    norm window:      "
                f"{sum(1 for record in summary_records if record.get('inside_norm_window'))}\n"
            )
        f.write("\n")
        f.write(f"Auto-deglitched points: {total_auto_deglitched_points}\n")
        f.write(f"Manually range-interpolated points: {total_manual_range_interpolated_points}\n")
        f.write(f"Deglitch note: {AUTO_DEGLITCH_WARNING}\n")
        f.write(f"Normalization order nnorm: {config.nnorm}\n")
        f.write("Normalize-before-average: False (merge μ(E) first, normalize merged spectrum)\n")
        f.write(f"Plots folder: {plots_dir}\n")
        f.write(f"Overview plots folder: {overview_dir}\n")
        f.write(f"Replicate QC folder: {replicate_qc_dir}\n")
        f.write(f"Detector raw folder: {output_dir / 'detector_raw'}\n")
        f.write(f"Detector raw files created: {len(detector_raw_files)}\n")
        if getattr(config, "save_detector_health_overview_plot", True):
            if detector_health_plot_info["path"] is not None:
                channels = ", ".join(detector_health_plot_info["channels"])
                skipped = ", ".join(detector_health_plot_info["skipped"]) or "None"
                f.write(f"Detector health overview: created (channels: {channels}; skipped: {skipped})\n")
            else:
                skipped = ", ".join(detector_health_plot_info["skipped"]) or "None"
                f.write(f"Detector health overview: not created (no plottable channels; skipped: {skipped})\n")
        else:
            f.write("Detector health overview: disabled\n")
        if getattr(config, "save_analysis_signal_qc_plot", True):
            if analysis_signal_qc_info["path"] is not None:
                skipped = "; ".join(analysis_signal_qc_info["skipped"]) or "None"
                avg_status = "yes" if analysis_signal_qc_info.get("average_plotted") else "no"
                f.write(
                    f"Analysis signal QC: created "
                    f"(signal: {analysis_signal_qc_info['signal']}; "
                    f"individual traces: {analysis_signal_qc_info['n_traces']}; "
                    f"average trace: {avg_status}; skipped/warnings: {skipped})\n"
                )
            else:
                skipped = "; ".join(analysis_signal_qc_info["skipped"]) or "no plottable traces"
                signal_name = analysis_signal_qc_info.get("signal") or _analysis_signal_spec(config)[1]
                f.write(f"Analysis signal QC: not created (signal: {signal_name}; reason: {skipped})\n")
        else:
            f.write("Analysis signal QC: disabled\n")
        if getattr(config, "save_processed_mu_replicate_qc_plot", True):
            f.write(
                f"Processed μ(E) replicate QC plots: created {processed_mu_replicate_qc_created}; "
                f"skipped {len(processed_mu_replicate_qc_skipped)}\n"
            )
            for skipped in processed_mu_replicate_qc_skipped:
                f.write(f"  - {skipped}\n")
        else:
            f.write("Processed μ(E) replicate QC plots: disabled\n")
        if getattr(config, "save_replicate_qc_plots", True):
            f.write(
                f"Normalized replicate QC plots: created {normalized_replicate_qc_created}; "
                f"skipped {len(normalized_replicate_qc_skipped)}\n"
            )
            for skipped in normalized_replicate_qc_skipped:
                f.write(f"  - {skipped}\n")
        else:
            f.write("Normalized replicate QC plots: disabled\n")
        if getattr(config, "save_detector_raw_overview_plot", False):
            f.write("Aligned averaged IF overview: plots/overview/aligned_averaged_IF_overview.png\n")
        else:
            f.write("Aligned averaged IF overview: disabled\n")
        if drift_plot_info["status"] == "created":
            f.write(f"Drift tracker: created ({_relative_output_path(drift_plot_info['path'], output_dir)})\n")
        elif drift_plot_info["status"] == "failed":
            f.write(f"Drift tracker: failed ({drift_plot_info['reason']})\n")
        else:
            f.write("Drift tracker: disabled\n")
        f.write("Processed μ(E) files use *_processed.dat; true detector channels use detector_raw/*_detector_raw.dat.\n")
        f.write(f"Total plots created: {len(relative_plot_files)}\n")
        f.write(f"Overview plots created: {len(overview_plot_files)}\n")
        for plot_file in overview_plot_files:
            f.write(f"- {plot_file}\n")
        f.write(f"Replicate QC plots created: {len(replicate_qc_plot_files)}\n")
        for plot_file in replicate_qc_plot_files:
            f.write(f"- {plot_file}\n")
        if other_plot_files:
            f.write(f"Other plots created: {len(other_plot_files)}\n")
            for plot_file in other_plot_files:
                f.write(f"- {plot_file}\n")
        f.write("\nProcessing warnings:\n")
        if warnings:
            for w in warnings:
                f.write(f"- {w}\n")
        else:
            f.write("- None\n")

    log("Done.")
    return {
        "output_dir": output_dir,
        "plots_dir": plots_dir,
        "detector_raw_dir": output_dir / "detector_raw",
        "detector_raw_files": len(detector_raw_files),
        "warnings": warnings,
        "validation_warnings": validation_warnings,
        "groups_processed": len(group_summary_rows),
        "manual_excluded": len(excluded_entries),
        "shift_rejected": len(auto_shift_rejected_entries),
        "auto_outliers": len(auto_outlier_entries),
        "auto_deglitched_points": total_auto_deglitched_points,
        "manual_range_interpolated_points": total_manual_range_interpolated_points,
    }
