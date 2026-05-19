from __future__ import annotations

from datetime import datetime
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from larch import Group
from larch.xafs import pre_edge

from ._config_utils import load_config_json
from .config import AstraConfig
from .io import collect_xasd_files, load_xasd, split_replicate_suffix, sanitize_name, build_default_output_directory
from .manifest import Manifest, ScanEntry, load_manifest
from .signals import compute_signals, get_signal
from .alignment import find_best_shift
from .edge_presets import get_edge_preset
from .grouping import group_samples
from .export import save_two_col
from .plotting import plot_overview, plot_replicate_qc, plot_drift, plot_detector_health_overview, plot_analysis_signal_qc
from .self_absorption import (
    SelfAbsorptionResult,
    evaluate_self_absorption_for_group,
    get_self_absorption_sensitivity,
    get_self_absorption_threshold,
    plot_self_absorption_qc,
    write_self_absorption_flags,
)

# Phase 2.1 compatibility re-exports. Removed in Phase 2.2 when consumers
# (beamtime/watcher.py, beamtime/groups.py) migrate to importing directly
# from astra_xas.single_scan.
from .single_scan import (
    SingleScanResult,
    process_single_scan,
    _entry_from_scan,
    _analysis_signal_spec,
    _alignment_signal_spec,
    _required_channels_for_signal,
    _range_overlap_status,
    _channel_validation_messages,
    _alignment_structure_warning,
    detect_detector_jumps,
    _config_bool,
    _config_float,
    RAW_DETECTOR_JUMP_CHANNELS,
    PRIMARY_DETECTOR_JUMP_CHANNELS,
    FDT_DETECTOR_JUMP_CHANNELS,
    DERIVED_DETECTOR_JUMP_CHANNELS,
)


AUTO_DEGLITCH_WARNING = (
    "Automatic deglitching is intended for narrow point-like spikes. "
    "Use manual range deglitching for broad artifacts."
)

SHIFT_CONVENTION = (
    "positive shift_eV means +shift_eV is added to the scan energy before "
    "interpolation/averaging; shifts are relative to the alignment anchor."
)


def interpolate_to_grid(E_source, mu_source, E_target, kind="linear"):
    f = interp1d(E_source, mu_source, kind=kind, bounds_error=False, fill_value=np.nan)
    return f(E_target)


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


def _detector_jump_category(record: dict) -> str:
    channel = record.get("channel", "")
    if record.get("include_in_summary"):
        return "primary_summary"
    if channel in FDT_DETECTOR_JUMP_CHANNELS:
        return "fdt_diagnostic"
    if channel in DERIVED_DETECTOR_JUMP_CHANNELS:
        return "derived_signal"
    if record.get("inside_alignment_window"):
        return "edge_or_alignment_excluded"
    if record.get("severity") == "low":
        return "low_severity_excluded"
    return "other_excluded"


def _detector_jump_event_map(records: list[dict], energy_gap_eV: float = 5.0) -> tuple[dict[int, str], list[dict]]:
    summary_records = sorted(
        [record for record in records if record.get("include_in_summary")],
        key=lambda record: float(record.get("energy_eV", 0.0)),
    )
    event_ids: dict[int, str] = {}
    events: list[dict] = []
    current: list[dict] = []

    for record in summary_records:
        if not current:
            current = [record]
            continue
        previous_energy = float(current[-1].get("energy_eV", 0.0))
        energy = float(record.get("energy_eV", 0.0))
        if energy - previous_energy <= energy_gap_eV:
            current.append(record)
        else:
            events.append({"records": current})
            current = [record]
    if current:
        events.append({"records": current})

    for idx, event in enumerate(events, start=1):
        event_id = f"E{idx:03d}"
        event_records = event["records"]
        energies = [float(record.get("energy_eV", 0.0)) for record in event_records]
        channels = sorted({str(record.get("channel", "")) for record in event_records})
        scans = sorted({str(record.get("filename", "")) for record in event_records})
        if len(channels) >= 3 and len(scans) >= 2:
            interpretation = "common multi-channel feature"
        elif len(channels) >= 2:
            interpretation = "multi-channel feature"
        else:
            interpretation = "single-channel feature"
        event.update({
            "event_id": event_id,
            "energy_center_eV": float(np.mean(energies)),
            "energy_min_eV": float(np.min(energies)),
            "energy_max_eV": float(np.max(energies)),
            "channels": ",".join(channels),
            "scans_affected": len(scans),
            "entries": len(event_records),
            "interpretation": interpretation,
        })
        for record in event_records:
            event_ids[id(record)] = event_id
    return event_ids, events


def _write_detector_jump_rows(f, records: list[dict], event_ids: dict[int, str]) -> None:
    f.write(
        "# event_id\tcategory\tfilename\tchannel\tenergy_eV\tjump_size\trelative_jump\t"
        "severity\tinside_plot_window\tinside_alignment_window\tinside_preedge_window\t"
        "inside_norm_window\tinclude_in_summary\tnote\n"
    )
    for record in sorted(records, key=_detector_jump_sort_key):
        f.write(
            f"{event_ids.get(id(record), '')}\t{_detector_jump_category(record)}\t"
            f"{record['filename']}\t{record['channel']}\t{record['energy_eV']:.6f}\t"
            f"{record['jump_size']:.8g}\t{record['relative_jump']:.8g}\t"
            f"{record['severity']}\t{record['inside_plot_window']}\t"
            f"{record['inside_alignment_window']}\t{record['inside_preedge_window']}\t"
            f"{record['inside_norm_window']}\t{record.get('include_in_summary', False)}\t"
            f"{record['note']}\n"
        )


def _assign_event_ids(records: list[dict]) -> None:
    """
    Assign event_id to each record dict in-place.
    Only records with include_in_summary=True receive a non-empty event_id.
    Events are clusters of include_in_summary=True records within
    ENERGY_CLUSTER_TOLERANCE eV of each other, sorted by energy.
    """
    ENERGY_CLUSTER_TOLERANCE = 5.0
    summary = sorted(
        [r for r in records if r.get("include_in_summary")],
        key=lambda r: float(r.get("energy_eV", 0.0)),
    )
    for r in records:
        r["event_id"] = ""
    if not summary:
        return

    clusters: list[list[dict]] = []
    current_cluster: list[dict] = [summary[0]]
    for rec in summary[1:]:
        if float(rec["energy_eV"]) - float(current_cluster[-1]["energy_eV"]) <= ENERGY_CLUSTER_TOLERANCE:
            current_cluster.append(rec)
        else:
            clusters.append(current_cluster)
            current_cluster = [rec]
    clusters.append(current_cluster)

    for i, cluster in enumerate(clusters):
        event_id = f"E{i + 1:03d}"
        for rec in cluster:
            rec["event_id"] = event_id


def _write_detector_jumps(path: Path, records: list[dict], config: AstraConfig) -> None:
    _assign_event_ids(records)
    threshold = _config_float(config, "detector_jump_threshold", 10.0)
    min_relative = _config_float(config, "detector_jump_min_relative", 0.05)
    summary_records = [record for record in records if record.get("include_in_summary")]
    fdt_records = [record for record in records if record.get("channel") in FDT_DETECTOR_JUMP_CHANNELS]
    derived_records = [record for record in records if record.get("channel") in DERIVED_DETECTOR_JUMP_CHANNELS]
    other_records = [
        record
        for record in records
        if not record.get("include_in_summary")
        and record.get("channel") not in FDT_DETECTOR_JUMP_CHANNELS
        and record.get("channel") not in DERIVED_DETECTOR_JUMP_CHANNELS
    ]
    n_other = len(records) - len(summary_records) - len(fdt_records) - len(derived_records)
    total_scans = len({str(record.get("filename", "")) for record in records if record.get("filename")})

    severity_rank = {"low": 0, "medium": 1, "high": 2}
    table_header = (
        "filename\tchannel\tenergy_eV\tjump_size\trelative_jump\tseverity\t"
        "inside_plot_window\tinside_alignment_window\tinside_preedge_window\t"
        "inside_norm_window\tinclude_in_summary\tevent_id\tcategory\tnote\n"
    )

    def category(record: dict) -> str:
        channel = record.get("channel", "")
        if record.get("include_in_summary"):
            return "primary_summary"
        if channel in FDT_DETECTOR_JUMP_CHANNELS:
            return "fdt_diagnostic"
        if channel in DERIVED_DETECTOR_JUMP_CHANNELS:
            return "derived_signal"
        if channel in PRIMARY_DETECTOR_JUMP_CHANNELS:
            return "edge_or_alignment_excluded"
        return "other_excluded"

    def write_full_rows(f, rows: list[dict]) -> None:
        for record in rows:
            f.write(
                f"{record['filename']}\t{record['channel']}\t{record['energy_eV']:.6f}\t"
                f"{record['jump_size']:.8g}\t{record['relative_jump']:.8g}\t"
                f"{record['severity']}\t{record['inside_plot_window']}\t"
                f"{record['inside_alignment_window']}\t{record['inside_preedge_window']}\t"
                f"{record['inside_norm_window']}\t{record.get('include_in_summary', False)}\t"
                f"{record.get('event_id', '')}\t{category(record)}\t{record['note']}\n"
            )

    def filename_energy_key(record: dict):
        return (str(record.get("filename", "")), float(record.get("energy_eV", 0.0)))

    def section1_key(record: dict):
        return (
            str(record.get("event_id", "")),
            str(record.get("filename", "")),
            float(record.get("energy_eV", 0.0)),
        )

    def full_table_key(record: dict):
        return (not bool(record.get("include_in_summary")), *_detector_jump_sort_key(record))

    events = []
    for event_id in sorted({record.get("event_id", "") for record in summary_records if record.get("event_id")}):
        cluster = [record for record in summary_records if record.get("event_id") == event_id]
        energies = [float(record.get("energy_eV", 0.0)) for record in cluster]
        channels = ",".join(sorted({str(record.get("channel", "")) for record in cluster}))
        scans_affected = len({str(record.get("filename", "")) for record in cluster if record.get("filename")})
        max_severity = max((str(record.get("severity", "low")) for record in cluster), key=lambda s: severity_rank.get(s, -1))
        events.append({
            "event_id": event_id,
            "energy_center_eV": float(np.mean(energies)),
            "energy_min_eV": float(np.min(energies)),
            "energy_max_eV": float(np.max(energies)),
            "channels": channels,
            "scans_affected": scans_affected,
            "entries": len(cluster),
            "max_severity": max_severity,
        })

    with path.open("w", encoding="utf-8") as f:
        f.write("# ============================================================\n")
        f.write("# ASTRA DETECTOR JUMP DIAGNOSTICS\n")
        f.write("# Diagnostic only. No data was modified.\n")
        f.write("# Generated by AstraXAS\n")
        f.write("# ============================================================\n")
        f.write("#\n")
        f.write("# Detection method: point-to-point MAD spike detector\n")
        f.write(
            f"# Severity thresholds: low={threshold:g}-{2 * threshold:g}x MAD, "
            f"medium={2 * threshold:g}-{5 * threshold:g}x MAD, "
            f"high=>{5 * threshold:g}x MAD\n"
        )
        f.write(f"#   where T = detector_jump_threshold = {threshold}\n")
        f.write(f"# Min relative jump: {min_relative}\n")
        f.write("# Spike vs step discrimination: recovery window = 5 points\n")
        f.write("#\n")
        f.write("# Channel classification:\n")
        f.write("#   Primary raw channels (counted in main summary): I0, I1, I2, IF\n")
        f.write("#   Diagnostic channel (reported separately):       FDT\n")
        f.write("#   Derived signals (excluded from main summary):   IF_over_I0, ln_I0_I1, ln_I1_I2\n")
        f.write("#\n")
        f.write("# ============================================================\n")
        f.write("# SUMMARY\n")
        f.write("# ============================================================\n")
        f.write(f"# Full diagnostic entries:                   {len(records)}\n")
        f.write(f"# Significant primary raw-channel entries:   {len(summary_records)}   (include_in_summary=True)\n")
        f.write(f"# FDT diagnostic spikes:                     {len(fdt_records)}\n")
        f.write(f"# Derived-signal features excluded:          {len(derived_records)}\n")
        f.write(f"# Other excluded entries:                    {n_other}\n")
        f.write("#\n")
        f.write("# ============================================================\n")
        f.write("# GROUPED SIGNIFICANT EVENTS\n")
        f.write("# ============================================================\n")
        f.write("# Events are clusters of significant primary raw-channel entries\n")
        f.write("# within 5 eV of each other. One physical event may appear in\n")
        f.write("# multiple channels and/or multiple scans.\n")
        f.write("#\n")
        f.write("# event_id\tenergy_center_eV\tenergy_range_eV\tmax_severity\tchannels\tscans_affected\tentries\n")
        if events:
            for event in events:
                if event["energy_min_eV"] == event["energy_max_eV"]:
                    energy_range = f"{event['energy_min_eV']:.2f}"
                else:
                    energy_range = f"{event['energy_min_eV']:.2f}-{event['energy_max_eV']:.2f}"
                f.write(
                    f"# {event['event_id']}\t{event['energy_center_eV']:.2f}\t{energy_range}\t"
                    f"{event['max_severity']}\t{event['channels']}\t"
                    f"{event['scans_affected']}/{total_scans}\t{event['entries']}\n"
                )
        else:
            f.write("# No significant primary raw-channel events detected.\n")
        f.write("#\n")
        f.write(table_header)

        f.write("# ============================================================\n")
        f.write("# SECTION 1: Significant primary raw-channel jumps\n")
        f.write("# These entries have include_in_summary=True.\n")
        f.write("# Channels: I0, I1, I2, IF\n")
        f.write("# ============================================================\n")
        f.write("# Columns: filename channel energy_eV jump_size relative_jump severity event_id inside_alignment_window inside_preedge_window inside_norm_window note\n")
        if summary_records:
            write_full_rows(f, sorted(summary_records, key=section1_key))
        else:
            f.write("# None detected.\n")

        f.write("# ============================================================\n")
        f.write("# SECTION 2: FDT diagnostic spikes\n")
        f.write("# FDT is a diagnostic channel. These entries are not counted\n")
        f.write("# in the main detector-jump summary.\n")
        f.write("# ============================================================\n")
        f.write("# Columns: filename channel energy_eV jump_size relative_jump severity note\n")
        if fdt_records:
            write_full_rows(f, sorted(fdt_records, key=filename_energy_key))
        else:
            f.write("# None detected.\n")

        f.write("# ============================================================\n")
        f.write("# SECTION 3: Derived-signal sharp features\n")
        f.write("# Channels: IF_over_I0, ln_I0_I1, ln_I1_I2\n")
        f.write("# These are excluded from the main summary because sharp changes\n")
        f.write("# may reflect real XAS spectral features near absorption edges.\n")
        f.write("# ============================================================\n")
        f.write("# Columns: filename channel energy_eV jump_size relative_jump severity inside_alignment_window note\n")
        if derived_records:
            write_full_rows(f, sorted(derived_records, key=filename_energy_key))
        else:
            f.write("# None detected.\n")

        f.write("# ============================================================\n")
        f.write("# SECTION 4: Other excluded entries\n")
        f.write("# Primary raw-channel entries excluded from the main summary\n")
        f.write("# (e.g. inside edge/alignment window with low severity).\n")
        f.write("# ============================================================\n")
        f.write("# Columns: filename channel energy_eV jump_size relative_jump severity inside_alignment_window inside_preedge_window inside_norm_window note\n")
        if other_records:
            write_full_rows(f, sorted(other_records, key=filename_energy_key))
        else:
            f.write("# None.\n")

        f.write("# ============================================================\n")
        f.write("# SECTION 5: Full diagnostic table (machine-readable)\n")
        f.write("# All entries. Tab-separated. One row per detected feature.\n")
        f.write("# ============================================================\n")
        f.write("# Columns: filename channel energy_eV jump_size relative_jump severity inside_plot_window inside_alignment_window inside_preedge_window inside_norm_window include_in_summary event_id category note\n")
        if records:
            write_full_rows(f, sorted(records, key=full_table_key))
        else:
            f.write("# None detected.\n")


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


def _short_pdf_text(value, limit: int = 80) -> str:
    text = "" if value is None else str(value)
    text = text.replace("\n", " ").strip()
    if len(text) <= limit:
        return text
    return text[: max(limit - 3, 0)] + "..."


def _read_dat_table_preview(path: Path, max_rows: int = 12, max_cols: int = 8) -> tuple[list[str], list[list[str]]]:
    header: list[str] = []
    rows: list[list[str]] = []
    if not path.exists():
        return header, rows

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                clean = line.lstrip("#").strip()
                if not clean or clean.startswith("="):
                    continue
                parts = clean.split("\t") if "\t" in clean else clean.split()
                known_columns = {
                    "output_base",
                    "filename",
                    "foil_filename",
                    "energy_eV",
                    "mean_edge_step",
                    "shift_eV",
                    "event_id",
                    "sample",
                    "status",
                }
                if len(parts) > 1 and any(part in known_columns for part in parts):
                    header = parts
                continue

            parts = line.split("\t") if "\t" in line else line.split()
            if parts:
                rows.append(parts)
            if len(rows) >= max_rows:
                break

    if not header and rows:
        header = [f"col{i + 1}" for i in range(len(rows[0]))]

    if max_cols > 0:
        def trim(parts: list[str]) -> list[str]:
            if len(parts) <= max_cols:
                return parts
            return parts[:max_cols] + ["..."]

        header = trim(header)
        rows = [trim(row) for row in rows]

    return header, rows


def _read_detector_jump_pdf_summary(path: Path) -> list[str]:
    if not path.exists():
        return []
    lines: list[str] = []
    capture = False
    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("# GROUPED SIGNIFICANT EVENTS"):
                break
            if line.startswith("# SUMMARY"):
                capture = True
            if not capture:
                continue
            clean = line.lstrip("#").strip()
            if not clean or clean.startswith("="):
                continue
            lines.append(clean)
    return lines


def _write_pdf_qc_report(output_path: Path, context: dict) -> None:
    from xml.sax.saxutils import escape

    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        Image,
        KeepTogether,
        PageBreak,
        Paragraph,
        Preformatted,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    config = context["config"]
    output_dir = Path(context["output_dir"])
    overview_dir = Path(context["overview_dir"])
    replicate_qc_dir = Path(context["replicate_qc_dir"])
    alignment_anchor_info = context.get("alignment_anchor_info", {})
    edge_preset = get_edge_preset(getattr(config, "edge_preset_key", "custom"))
    edge_preset_ref_e0 = f"{edge_preset.e0_ref:.3f}" if edge_preset is not None else "N/A"
    edge_preset_note = getattr(config, "edge_preset_note", "") or (
        "Preset values are editable starting values, not absolute calibration."
        if getattr(config, "edge_preset_applied", False)
        else "N/A"
    )

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=36,
        rightMargin=36,
        topMargin=36,
        bottomMargin=36,
    )
    styles = getSampleStyleSheet()
    styles["BodyText"].fontSize = 9
    styles["BodyText"].leading = 11
    styles["Heading1"].spaceAfter = 8
    styles["Heading2"].spaceBefore = 6
    styles["Heading2"].spaceAfter = 6
    story = []

    def para(text, style_name: str = "BodyText"):
        safe = escape(str(text)).replace("\n", "<br/>")
        return Paragraph(safe, styles[style_name])

    def add_heading(text: str, level: int = 1) -> None:
        story.append(para(text, f"Heading{level}"))
        story.append(Spacer(1, 0.08 * inch))

    def add_list(title: str, items: list[str], empty_text: str = "- None") -> None:
        block = [para(title, "Heading3")]
        if items:
            for item in items[:25]:
                block.append(para(f"- {item}"))
            if len(items) > 25:
                block.append(para(f"- ... {len(items) - 25} more"))
        else:
            block.append(para(empty_text))
        block.append(Spacer(1, 0.08 * inch))
        if len(block) <= 10:
            story.append(KeepTogether(block))
        else:
            story.extend(block)

    def add_key_value_table(rows: list[tuple[str, object]]) -> None:
        data = [
            [para(_short_pdf_text(key, 45)), para(value)]
            for key, value in rows
        ]
        table = Table(data, colWidths=[2.05 * inch, doc.width - 2.05 * inch], hAlign="LEFT")
        table.setStyle(TableStyle([
            ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#f2f2f2")),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.append(table)
        story.append(Spacer(1, 0.14 * inch))

    def add_qc_status_box(status: str, findings: list[str]) -> None:
        color = {
            "Passed": colors.HexColor("#d8efe0"),
            "Check recommended": colors.HexColor("#fff0c7"),
            "Attention required": colors.HexColor("#f5d0cc"),
        }.get(status, colors.HexColor("#eeeeee"))
        findings_text = "<br/>".join(escape(f"- {finding}") for finding in findings) if findings else "- None"
        data = [
            [para("QC status", "Heading3"), para(status, "Heading3")],
            [para("Main findings"), Paragraph(findings_text, styles["BodyText"])],
        ]
        table = Table(data, colWidths=[1.7 * inch, doc.width - 1.7 * inch], hAlign="LEFT")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), color),
            ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor("#777777")),
            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#999999")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        story.append(KeepTogether([table, Spacer(1, 0.16 * inch)]))

    def add_image(path: Path, title: str, max_height: float = 4.6 * inch, caption: str | None = None) -> None:
        story.append(para(title, "Heading3"))
        if caption:
            story.append(para(caption))
        if not path.exists():
            story.append(para("not available"))
            story.append(Spacer(1, 0.08 * inch))
            return
        try:
            image = Image(str(path))
            scale = min(doc.width / image.imageWidth, max_height / image.imageHeight, 1.0)
            image.drawWidth = image.imageWidth * scale
            image.drawHeight = image.imageHeight * scale
            story.append(image)
        except Exception as exc:
            story.append(para(f"not available ({exc})"))
        story.append(Spacer(1, 0.12 * inch))

    def add_table_preview(path: Path, title: str, max_rows: int = 12, max_cols: int = 7) -> None:
        story.append(para(title, "Heading3"))
        if not path.exists():
            story.append(para("not available"))
            story.append(Spacer(1, 0.08 * inch))
            return
        header, rows = _read_dat_table_preview(path, max_rows=max_rows, max_cols=max_cols)
        if not rows:
            story.append(para("not available or no data rows"))
            story.append(Spacer(1, 0.08 * inch))
            return
        data = [[_short_pdf_text(cell, 28) for cell in header]]
        data.extend([[_short_pdf_text(cell, 28) for cell in row] for row in rows])
        col_width = doc.width / max(len(data[0]), 1)
        table = Table(data, colWidths=[col_width] * len(data[0]), repeatRows=1, hAlign="LEFT")
        table.setStyle(TableStyle([
            ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
            ("FONTSIZE", (0, 0), (-1, -1), 6.5),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e6edf5")),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(table)
        if len(rows) >= max_rows:
            story.append(para(f"Showing first {max_rows} rows. See {path.name} for full table."))
        else:
            story.append(para(f"Full table available in {path.name}."))
        story.append(Spacer(1, 0.12 * inch))

    def replicate_plot_title(path: Path) -> str:
        name = path.name
        if name.endswith("_processed_mu_replicate_qc.png"):
            sample = name.removesuffix("_processed_mu_replicate_qc.png")
            return f"Processed μ(E) replicate QC: {sample}"
        if name.endswith("_normalized_replicate_qc.png"):
            sample = name.removesuffix("_normalized_replicate_qc.png")
            return f"Normalized replicate QC: {sample}"
        return path.stem.replace("_", " ")

    jump_records = context.get("all_jump_records", [])
    detector_jump_diagnostic = context.get("detector_jump_diagnostic", {"status": "disabled"})
    summary_records = [record for record in jump_records if record.get("include_in_summary")]
    fdt_count = sum(1 for record in jump_records if record.get("channel") in FDT_DETECTOR_JUMP_CHANNELS)
    derived_count = sum(1 for record in jump_records if record.get("channel") in DERIVED_DETECTOR_JUMP_CHANNELS)
    self_absorption_summary = context.get("self_absorption_summary", {})
    self_absorption_results = context.get("self_absorption_results", [])
    validation_warnings_pdf = context.get("validation_warnings", [])
    processing_warnings_pdf = context.get("processing_warnings", [])
    low_quality_count_pdf = int(context.get("low_quality_count", 0) or 0)
    auto_deglitched_pdf = int(context.get("total_auto_deglitched_points", 0) or 0)
    manual_interpolated_pdf = int(context.get("total_manual_range_interpolated_points", 0) or 0)
    serious_terms = ("required channel", "no valid spectra", "failed", "error", "could not", "missing")
    serious_validation = any(
        "required channel" in str(warning).lower()
        and ("missing" in str(warning).lower() or "no finite" in str(warning).lower())
        for warning in validation_warnings_pdf
    )
    serious_processing = any(any(term in str(warning).lower() for term in serious_terms) for warning in processing_warnings_pdf)
    if detector_jump_diagnostic.get("status") == "error" or low_quality_count_pdf > 0 or serious_validation or serious_processing:
        qc_status = "Attention required"
    elif (
        validation_warnings_pdf
        or processing_warnings_pdf
        or summary_records
        or auto_deglitched_pdf > 0
        or manual_interpolated_pdf > 0
        or detector_jump_diagnostic.get("status") == "found"
    ):
        qc_status = "Check recommended"
    else:
        qc_status = "Passed"

    qc_findings = []
    if validation_warnings_pdf:
        first_validation = str(validation_warnings_pdf[0])
        qc_findings.append(f"Validation warning: {_short_pdf_text(first_validation, 120)}")
    else:
        qc_findings.append("Validation warnings: none")
    if summary_records:
        qc_findings.append("Significant primary raw-channel detector jumps detected; see detector jump summary.")
    else:
        qc_findings.append("Detector jumps: no significant primary raw-channel jumps.")
    if low_quality_count_pdf > 0:
        qc_findings.append(f"Alignment quality: {low_quality_count_pdf} low-quality alignment(s).")
    else:
        qc_findings.append("Alignment quality: no low-quality alignments.")
    if auto_deglitched_pdf or manual_interpolated_pdf:
        qc_findings.append(
            f"Deglitching/interpolation: {auto_deglitched_pdf} auto-deglitched point(s), "
            f"{manual_interpolated_pdf} manually interpolated point(s)."
        )
    else:
        qc_findings.append("Deglitching: no auto-deglitched or manually interpolated points.")

    add_heading("AstraXAS Processing and QC Report")
    story.append(para("Diagnostic report only. Spectra are not recomputed while building this PDF."))
    story.append(Spacer(1, 0.14 * inch))
    manifest_path_ctx = context.get("manifest_path")
    manifest_created_ctx = context.get("manifest_created_iso", "")
    manifest_row_value = (
        f"{manifest_path_ctx} (created {manifest_created_ctx or 'unknown'})"
        if manifest_path_ctx
        else "N/A (folder mode)"
    )
    add_key_value_table([
        ("ASTRA version", getattr(config, "version", "N/A")),
        ("Driven by manifest", manifest_row_value),
        ("Input directory", context.get("input_dir", "")),
        ("Output directory", output_dir),
        ("Processing date/time", context.get("processing_datetime", "")),
        ("Files found", context.get("files_found", "N/A")),
        ("Analysis mode", getattr(config, "analysis_mode", "N/A")),
        ("Alignment source", context.get("alignment_source", "N/A")),
        ("Alignment signal", context.get("alignment_signal_label", "N/A")),
        ("Alignment anchor mode", alignment_anchor_info.get("mode", "N/A")),
        ("Alignment anchor file", alignment_anchor_info.get("path") or "N/A"),
        ("Edge preset", getattr(config, "edge_preset_label", "Custom")),
        ("Preset applied", getattr(config, "edge_preset_applied", False)),
        ("Preset reference E0", edge_preset_ref_e0),
        ("Final E0 used", getattr(config, "e0", "N/A")),
        ("Preset note", edge_preset_note),
        ("Normalization order nnorm", getattr(config, "nnorm", "N/A")),
        ("Auto-deglitching", _config_bool(config, "enable_auto_deglitch", False)),
        ("Manual range deglitching", _config_bool(config, "enable_manual_deglitch_range", False)),
        ("Output folders", "plots/overview; plots/replicate_qc; detector_raw"),
    ])
    add_qc_status_box(qc_status, qc_findings[:5])

    add_heading("Warnings And Diagnostics", 2)
    add_list("Validation warnings", validation_warnings_pdf)
    add_list("Processing warnings", processing_warnings_pdf)
    detector_status = detector_jump_diagnostic.get("status", "disabled")
    if detector_status == "found":
        detector_summary = [
            f"Significant primary raw-channel jumps: {len(summary_records)}",
            f"Full diagnostic entries: {len(jump_records)}",
            f"FDT diagnostic spikes: {fdt_count}",
            f"Derived-signal features excluded: {derived_count}",
            "Full details: ASTRA_detector_jumps.dat",
        ]
    elif detector_status == "none":
        detector_summary = ["None detected"]
    elif detector_status == "error":
        detector_summary = [f"Failed ({detector_jump_diagnostic.get('error', '')})"]
    else:
        detector_summary = ["Disabled"]
    add_list("Detector jump diagnostics", detector_summary)
    if self_absorption_summary.get("enabled"):
        self_absorption_lines = [
            "Heuristic diagnostic only. It flags likely fluorescence self-absorption when the fluorescence white-line amplitude is suppressed relative to the simultaneously available sample transmission signal.",
            f"Sensitivity: {self_absorption_summary.get('sensitivity', 'normal')}",
            f"Threshold: {float(self_absorption_summary.get('threshold', 0.85)):.4f}",
            f"Groups checked: {self_absorption_summary.get('checked', 0)}",
            f"Groups flagged: {self_absorption_summary.get('flagged', 0)}",
            f"Groups skipped: {self_absorption_summary.get('skipped', 0)}",
        ]
        flagged = [result for result in self_absorption_results if result.status == "flagged"]
        for result in flagged[:8]:
            self_absorption_lines.append(
                f"{result.sample}: ratio={result.ratio_fluo_over_trans:.4f}, "
                f"threshold={result.threshold_used:.4f}, severity={result.severity}"
            )
        if len(flagged) > 8:
            self_absorption_lines.append(f"... {len(flagged) - 8} more flagged group(s)")
    else:
        self_absorption_lines = [
            self_absorption_summary.get("reason")
            or "disabled because analysis mode is not fluorescence"
        ]
    add_list("Self-absorption diagnostic", self_absorption_lines)
    add_key_value_table([
        ("Low-quality alignments", context.get("low_quality_count", 0)),
        ("Auto-deglitched points", context.get("total_auto_deglitched_points", 0)),
        ("Manually range-interpolated points", context.get("total_manual_range_interpolated_points", 0)),
    ])

    story.append(PageBreak())
    add_heading("Alignment And Drift", 2)
    add_key_value_table([
        ("Alignment source", context.get("alignment_source", "N/A")),
        ("Alignment signal", context.get("alignment_signal_label", "N/A")),
        ("Alignment anchor mode", alignment_anchor_info.get("mode", "N/A")),
        ("Alignment anchor status", alignment_anchor_info.get("status", "N/A")),
    ])
    story.append(para(f"Shift convention: {SHIFT_CONVENTION}"))
    story.append(Spacer(1, 0.08 * inch))
    add_image(overview_dir / "drift_tracker.png", "Drift tracker")

    story.append(PageBreak())
    add_heading("Detector QC", 2)
    add_image(overview_dir / "detector_health_overview.png", "Detector health overview")
    add_image(overview_dir / "analysis_signal_qc.png", "Analysis signal QC")
    jump_summary_path = output_dir / "ASTRA_detector_jumps.dat"
    jump_lines = _read_detector_jump_pdf_summary(jump_summary_path)
    story.append(para("Detector jump summary", "Heading3"))
    if jump_lines:
        story.append(Preformatted("\n".join(jump_lines[:18]), styles["Code"]))
        story.append(para("Full details: ASTRA_detector_jumps.dat"))
    else:
        story.append(para("not available"))

    story.append(PageBreak())
    add_heading("Spectrum Overview", 2)
    add_image(overview_dir / "processed_mu_overview.png", "Processed mu(E) overview")
    add_image(overview_dir / "background_corrected_overview.png", "Background-corrected overview")
    add_image(overview_dir / "normalized_overview.png", "Normalized overview")
    add_image(overview_dir / "aligned_averaged_IF_overview.png", "Aligned averaged IF detector signal")

    story.append(PageBreak())
    add_heading("Replicate QC", 2)
    replicate_paths = sorted(replicate_qc_dir.glob("*_processed_mu_replicate_qc.png"))
    replicate_paths.extend(sorted(replicate_qc_dir.glob("*_normalized_replicate_qc.png")))
    if replicate_paths:
        for idx, path in enumerate(replicate_paths):
            add_image(
                path,
                replicate_plot_title(path),
                max_height=3.8 * inch,
                caption=f"File: {_relative_output_path(path, output_dir)}",
            )
            if (idx + 1) % 2 == 0 and idx + 1 < len(replicate_paths):
                story.append(PageBreak())
    else:
        story.append(para("No replicate QC plots available."))

    story.append(PageBreak())
    add_heading("Main Output Files", 2)
    main_outputs = [
        "*_processed.dat",
        "*_norm.dat",
        "*_flat.dat, if created",
        "detector_raw/*_detector_raw.dat",
        "ASTRA_processing_report.txt",
        "ASTRA_energy_shifts.dat",
        "ASTRA_foil_alignment.dat, if created",
        "ASTRA_normalization_summary.dat",
        "plots/overview/",
        "plots/replicate_qc/",
    ]
    if (output_dir / "ASTRA_detector_jumps.dat").exists():
        main_outputs.insert(5, "ASTRA_detector_jumps.dat")
    if (output_dir / "ASTRA_self_absorption_flags.dat").exists():
        main_outputs.insert(6, "ASTRA_self_absorption_flags.dat")
    add_list("Main output files", main_outputs)

    add_heading("Summary Tables", 2)
    add_table_preview(output_dir / "ASTRA_normalization_summary.dat", "Normalization summary")
    add_table_preview(output_dir / "ASTRA_energy_shifts.dat", "Energy shifts", max_rows=12, max_cols=8)
    add_table_preview(output_dir / "ASTRA_self_absorption_flags.dat", "Self-absorption flags", max_rows=12, max_cols=8)
    story.append(para("Detector jump diagnostic header", "Heading3"))
    if jump_lines:
        story.append(Preformatted("\n".join(jump_lines[:18]), styles["Code"]))
    else:
        story.append(para("not available"))

    def draw_page_number(canvas, _doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.grey)
        canvas.drawRightString(_doc.pagesize[0] - 36, 18, f"Page {canvas.getPageNumber()}")
        canvas.restoreState()

    doc.build(story, onFirstPage=draw_page_number, onLaterPages=draw_page_number)


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


def _nanmean_without_warning(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return values
    counts = np.sum(np.isfinite(values), axis=0)
    sums = np.nansum(values, axis=0)
    out = np.full(values.shape[1], np.nan, dtype=float)
    valid = counts > 0
    out[valid] = sums[valid] / counts[valid]
    return out


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


def _build_entries_from_manifest(
    manifest: Manifest, config: AstraConfig, log=print
) -> tuple[list[Path], list[dict]]:
    """Triage manifest scans by status, verify files exist, and build entries.

    Returns (files, entries) parallel to the folder-mode block:
      - files: list of resolved .xasd paths actually loaded
      - entries: list of entry dicts in manifest array order, with
        ``base_name`` overridden to the manifest group id, ``is_foil``
        overridden by manifest.is_foil when non-null, and an injected
        ``_order_in_manifest`` so the grouping sort honours manifest order.
    """
    plan: list[tuple[Path, ScanEntry]] = []
    n_good = n_bad = n_excluded = n_unreviewed = 0
    missing: list[str] = []
    for scan in manifest.scans:
        if scan.status == "good":
            path = manifest.resolve_scan_path(scan)
            if not path.exists():
                missing.append(f"{scan.filename} -> {path}")
                continue
            plan.append((path, scan))
            n_good += 1
        elif scan.status == "bad":
            log(f"Manifest: skipping (bad data): {scan.filename}")
            n_bad += 1
        elif scan.status == "excluded":
            log(f"Manifest: skipping (excluded from analysis): {scan.filename}")
            n_excluded += 1
        elif scan.status == "unreviewed":
            n_unreviewed += 1

    if missing:
        raise RuntimeError(
            "Manifest references files that do not exist on disk:\n  "
            + "\n  ".join(missing)
        )
    if n_unreviewed:
        log(
            f"WARNING: {n_unreviewed} manifest scan(s) are 'unreviewed'; "
            f"run curation (or mark them 'good') before processing. "
            f"These scans were skipped for this run."
        )
    if n_good == 0:
        raise RuntimeError("Manifest has no scans with status='good'.")

    log(
        f"Manifest plan: {n_good} good, {n_bad} bad, "
        f"{n_excluded} excluded, {n_unreviewed} unreviewed."
    )

    files = [p for p, _ in plan]
    entries: list[dict] = []
    for i, (path, scan) in enumerate(plan):
        scan_data = load_xasd(path)
        entry = _entry_from_scan(scan_data, config, path=path)
        if scan.is_foil is not None:
            entry["is_foil"] = scan.is_foil
        if scan.group is not None:
            entry["base_name"] = scan.group
        entry["_order_in_manifest"] = i
        entries.append(entry)

    missing_groups = [
        entry["filename"]
        for entry, (_, scan) in zip(entries, plan)
        if not entry["is_foil"] and scan.group is None
    ]
    if missing_groups:
        raise RuntimeError(
            "Manifest mode requires non-foil scans with status='good' to have "
            "a non-null 'group' field. The following scans are missing groups:\n  "
            + "\n  ".join(missing_groups)
        )
    return files, entries


def process_folder(
    input_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    config: AstraConfig | None = None,
    *,
    session: str | Path | None = None,
    log=print,
):
    manifest: Manifest | None = None
    manifest_path: Path | None = None
    if session is not None:
        if input_dir is not None:
            raise RuntimeError(
                "process_folder: provide either input_dir or session, not both."
            )
        if config is not None:
            raise RuntimeError(
                "process_folder: provide either config or session, not both "
                "(the manifest specifies its own config)."
            )
        manifest_path = Path(session).expanduser().resolve()
        manifest = load_manifest(manifest_path)

        if manifest.config.source == "inline":
            raise RuntimeError(
                "config.source='inline' is reserved for v2; manifest schema v1 "
                "accepts but does not apply inline configs. Use config.source='path'."
            )

        config_path = manifest.resolve_config_path()
        if config_path is None or not config_path.exists():
            raise RuntimeError(
                f"Manifest references a config that does not exist: {config_path}"
            )
        try:
            config = load_config_json(config_path)
        except Exception as exc:
            raise RuntimeError(
                f"Could not load manifest config {config_path}: {exc}"
            ) from exc

        input_dir = manifest.resolve_base_dir()

        if _as_clean_list(getattr(config, "exclude_filenames", ())) or _as_clean_list(
            getattr(config, "exclude_filename_contains", ())
        ):
            log(
                "Manifest mode: bulk filename exclusions in config "
                "(exclude_filenames, exclude_filename_contains) are ignored; "
                "use manifest status='excluded' instead."
            )

        n_assigned_foil_overrides = sum(1 for s in manifest.scans if s.assigned_foil)
        if n_assigned_foil_overrides:
            log(
                f"Manifest field 'assigned_foil' is reserved for v2 and ignored "
                f"in this run (set on {n_assigned_foil_overrides} scan(s))."
            )

    config = config or AstraConfig()
    if input_dir is None:
        raise RuntimeError("process_folder: input_dir or session is required.")
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
        or (
            getattr(config, "analysis_mode", "fluo") == "fluo"
            and getattr(config, "enable_self_absorption_check", True)
            and getattr(config, "save_self_absorption_qc_plots", True)
        )
    )
    if plots_enabled:
        plots_dir.mkdir(parents=True, exist_ok=True)
        overview_dir.mkdir(parents=True, exist_ok=True)
        if (
            getattr(config, "save_processed_mu_replicate_qc_plot", True)
            or getattr(config, "save_replicate_qc_plots", True)
            or (
                getattr(config, "analysis_mode", "fluo") == "fluo"
                and getattr(config, "enable_self_absorption_check", True)
                and getattr(config, "save_self_absorption_qc_plots", True)
            )
        ):
            replicate_qc_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    validation_warnings: list[str] = []
    log(f"Starting {config.version}")
    if manifest is not None:
        log(f"Driven by session manifest: {manifest_path} (created {manifest.created_iso or 'unknown'})")
    log(f"Input directory: {input_dir}")
    log(f"Output directory: {output_dir}")

    if manifest is not None:
        files, entries = _build_entries_from_manifest(manifest, config, log=log)
    else:
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

    # Manual exclusions before alignment. Bypassed in manifest mode (status='excluded' replaces it).
    exclude_filenames: set[str] = set()
    exclude_contains: list[str] = []
    excluded_entries = []
    auto_shift_rejected_entries = []
    auto_outlier_entries = []

    if manifest is None:
        exclude_filenames = set(_as_clean_list(getattr(config, "exclude_filenames", ())))
        exclude_contains = _as_clean_list(getattr(config, "exclude_filename_contains", ()))
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
                "No foil files found. AstraXAS expected one or more files whose "
                f"names contain the foil_keyword (currently '{config.foil_keyword}') "
                "because alignment_source is 'separate_foil'.\n"
                "\n"
                "Possible fixes:\n"
                "  1. If your data includes a reference channel measured in every "
                "scan (typical for operando setups), set "
                "alignment_source = 'inline_ref' in your config or via the "
                "GUI/Python API.\n"
                "  2. If your foil scans use a different filename keyword "
                "(e.g., 'reference' or the element name), change foil_keyword "
                "in your config.\n"
                "  3. If you intended to include foil scans but forgot, add them "
                "to the input folder."
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
    self_absorption_results: list[SelfAbsorptionResult] = []
    self_absorption_should_run = (
        getattr(config, "analysis_mode", "fluo") == "fluo"
        and _config_bool(config, "enable_self_absorption_check", True)
    )
    self_absorption_summary = {
        "enabled": self_absorption_should_run,
        "reason": "",
        "threshold": get_self_absorption_threshold(config),
        "sensitivity": get_self_absorption_sensitivity(config),
        "checked": 0,
        "flagged": 0,
        "skipped": 0,
        "path": None,
    }
    if getattr(config, "analysis_mode", "fluo") != "fluo":
        self_absorption_summary["reason"] = "disabled because analysis mode is not fluorescence"
    elif not _config_bool(config, "enable_self_absorption_check", True):
        self_absorption_summary["reason"] = "disabled by user"

    for (base_name, assigned_foil), scans in groups.items():
        log(f"Processing group: {base_name} ({len(scans)} scan(s))")
        first_scan = scans[0]
        master_energy = first_scan["energy"] + first_scan["energy_shift_eV"]
        processed_spectra = []
        trans_spectra = []
        source_filenames = []
        detector_group_if = []
        deglitch_records = []  # Collect deglitch info per replicate

        for s in scans:
            raw_shifted_energy = s["energy"] + s["energy_shift_eV"]
            shifted_energy = raw_shifted_energy
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

            if self_absorption_should_run:
                trans_interp = np.full_like(master_energy, np.nan, dtype=float)
                try:
                    i0 = s.get("I0")
                    i1 = s.get("I1")
                    if i0 is not None and i1 is not None:
                        i0 = np.asarray(i0, dtype=float)
                        i1 = np.asarray(i1, dtype=float)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            mu_trans_diag = np.where(
                                np.isfinite(i0) & np.isfinite(i1) & (i0 > 0) & (i1 > 0),
                                np.log(i0 / i1),
                                np.nan,
                            )
                        candidate = interpolate_to_grid(raw_shifted_energy, mu_trans_diag, master_energy, kind=config.interp_kind)
                        valid_trans = np.isfinite(candidate)
                        if np.count_nonzero(valid_trans) >= 0.9 * len(master_energy):
                            if not np.all(valid_trans):
                                candidate = np.interp(master_energy, master_energy[valid_trans], candidate[valid_trans])
                            trans_interp = candidate
                except Exception:
                    pass
                trans_spectra.append(trans_interp)

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
        if self_absorption_should_run:
            trans_array = np.asarray(trans_spectra, dtype=float)[keep_mask]
        else:
            trans_array = None
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

        if self_absorption_should_run:
            try:
                trans_avg = _nanmean_without_warning(trans_array)
                result = evaluate_self_absorption_for_group(
                    output_base,
                    master_energy,
                    processed_avg,
                    trans_avg,
                    config,
                    n_replicates_used=len(processed_array),
                    normalize_func=_run_pre_edge,
                )
            except Exception as exc:
                result = SelfAbsorptionResult(
                    sample=output_base,
                    status="skipped",
                    severity="none",
                    ratio_fluo_over_trans=float("nan"),
                    fluo_white_line_amplitude=float("nan"),
                    trans_white_line_amplitude=float("nan"),
                    threshold_used=get_self_absorption_threshold(config),
                    sensitivity=get_self_absorption_sensitivity(config),
                    white_line_window_eV=(
                        config.e0 + getattr(config, "self_absorption_wl_min", 0.0),
                        config.e0 + getattr(config, "self_absorption_wl_max", 35.0),
                    ),
                    continuum_window_eV=(
                        config.e0 + getattr(config, "self_absorption_cont_min", 50.0),
                        config.e0 + getattr(config, "self_absorption_cont_max", 150.0),
                    ),
                    n_replicates_used=len(processed_array),
                    note=f"skipped: self-absorption diagnostic failed ({exc})",
                )
                log(f"WARNING: Self-absorption diagnostic failed for {output_base}: {exc}")
            self_absorption_results.append(result)
            if result.status == "flagged":
                warnings.append(
                    f"Possible self-absorption: {output_base} "
                    f"(ratio={result.ratio_fluo_over_trans:.3f}, "
                    f"threshold={result.threshold_used:.2f}, severity={result.severity})"
                )
            if (
                getattr(config, "save_self_absorption_qc_plots", True)
                and result.status in {"ok", "flagged"}
                and result.fluo_norm is not None
                and result.trans_norm is not None
            ):
                try:
                    path = plot_self_absorption_qc(
                        result,
                        master_energy,
                        result.fluo_norm,
                        result.trans_norm,
                        replicate_qc_dir / f"{output_base}_self_absorption_qc.png",
                    )
                    if path is not None:
                        plot_files.append(path)
                        log(f"Saving self-absorption QC plot: {output_base}")
                except Exception as exc:
                    warning = f"Could not save self-absorption QC plot for {output_base}: {exc}"
                    warnings.append(warning)
                    log("WARNING: " + warning)

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

    self_absorption_path = output_dir / "ASTRA_self_absorption_flags.dat"
    if self_absorption_should_run:
        try:
            path = write_self_absorption_flags(self_absorption_results, output_dir, config)
            self_absorption_summary.update({
                "path": path,
                "checked": sum(1 for result in self_absorption_results if result.status in {"ok", "flagged"}),
                "flagged": sum(1 for result in self_absorption_results if result.status == "flagged"),
                "skipped": sum(1 for result in self_absorption_results if result.status == "skipped"),
            })
            log(f"Saved self-absorption diagnostic: {path}")
        except Exception as exc:
            warning = f"Self-absorption diagnostic output failed: {exc}"
            warnings.append(warning)
            self_absorption_summary["reason"] = f"output failed: {exc}"
            log("WARNING: " + warning)
    else:
        try:
            if self_absorption_path.exists():
                self_absorption_path.unlink()
        except Exception:
            pass

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
    processing_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    edge_preset = get_edge_preset(getattr(config, "edge_preset_key", "custom"))
    edge_preset_ref_e0 = f"{edge_preset.e0_ref:.6g}" if edge_preset is not None else "N/A"
    edge_preset_note = getattr(config, "edge_preset_note", "") or (
        "Preset values are editable starting values, not absolute calibration."
        if getattr(config, "edge_preset_applied", False)
        else "N/A"
    )
    pdf_report_status = {"status": "disabled", "path": None, "error": ""}
    if getattr(config, "save_pdf_report", True):
        pdf_path = output_dir / "ASTRA_processing_and_QC_report.pdf"
        try:
            _write_pdf_qc_report(pdf_path, {
                "config": config,
                "input_dir": input_dir,
                "output_dir": output_dir,
                "processing_datetime": processing_datetime,
                "files_found": len(files),
                "alignment_source": alignment_source,
                "alignment_signal_label": alignment_signal_label,
                "alignment_anchor_info": alignment_anchor_info,
                "validation_warnings": validation_warnings,
                "processing_warnings": warnings,
                "detector_jump_diagnostic": detector_jump_diagnostic,
                "all_jump_records": all_jump_records,
                "self_absorption_summary": self_absorption_summary,
                "self_absorption_results": self_absorption_results,
                "low_quality_count": low_quality_count,
                "total_auto_deglitched_points": total_auto_deglitched_points,
                "total_manual_range_interpolated_points": total_manual_range_interpolated_points,
                "overview_dir": overview_dir,
                "replicate_qc_dir": replicate_qc_dir,
                "manifest_path": manifest_path,
                "manifest_created_iso": manifest.created_iso if manifest is not None else "",
            })
            pdf_report_status = {"status": "created", "path": pdf_path, "error": ""}
            log(f"Saved PDF QC report: {pdf_path}")
        except Exception as exc:
            pdf_report_status = {"status": "failed", "path": pdf_path, "error": str(exc)}
            log(f"WARNING: PDF QC report failed: {exc}")

    with open(output_dir / "ASTRA_processing_report.txt", "w", encoding="utf-8") as f:
        f.write(f"{config.version}\n")
        if manifest is not None:
            f.write(
                f"Driven by manifest: {manifest_path} "
                f"(created {manifest.created_iso or 'unknown'})\n"
            )
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
        f.write(f"Edge preset: {getattr(config, 'edge_preset_label', 'Custom')}\n")
        f.write(f"Preset applied: {getattr(config, 'edge_preset_applied', False)}\n")
        f.write(f"Preset reference E0: {edge_preset_ref_e0}\n")
        f.write(f"Final E0 used: {config.e0}\n")
        f.write(f"Preset note: {edge_preset_note}\n")
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
        f.write("\nSelf-absorption diagnostic:\n")
        f.write(
            "This is a heuristic diagnostic. It flags likely fluorescence self-absorption when "
            "the fluorescence white-line amplitude is suppressed relative to the simultaneously "
            "available sample transmission signal. It does not correct the spectrum and should "
            "not be treated as proof of self-absorption.\n"
        )
        if self_absorption_summary["enabled"]:
            f.write("Self-absorption diagnostic enabled: yes\n")
            f.write(f"Sensitivity: {self_absorption_summary['sensitivity']}\n")
            f.write(f"Threshold: {self_absorption_summary['threshold']:.4f}\n")
            f.write(f"Groups checked: {self_absorption_summary['checked']}\n")
            f.write(f"Groups flagged: {self_absorption_summary['flagged']}\n")
            f.write(f"Groups skipped: {self_absorption_summary['skipped']}\n")
            if self_absorption_summary.get("path") is not None:
                f.write(f"Output file: {_relative_output_path(self_absorption_summary['path'], output_dir)}\n")
            flagged = [result for result in self_absorption_results if result.status == "flagged"]
            if flagged:
                f.write("Flagged groups:\n")
                for result in flagged:
                    f.write(
                        f"- {result.sample}: ratio={result.ratio_fluo_over_trans:.4f}, "
                        f"threshold={result.threshold_used:.4f}, severity={result.severity}\n"
                    )
        elif getattr(config, "analysis_mode", "fluo") != "fluo":
            f.write("Self-absorption diagnostic disabled because analysis mode is not fluorescence.\n")
        else:
            f.write("Self-absorption diagnostic disabled by user.\n")
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
        if pdf_report_status["status"] == "created":
            f.write(f"PDF QC report: created ({_relative_output_path(pdf_report_status['path'], output_dir)})\n")
        elif pdf_report_status["status"] == "failed":
            f.write(f"PDF QC report: failed ({pdf_report_status['error']})\n")
        else:
            f.write("PDF QC report: disabled\n")
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
        "pdf_report": pdf_report_status["path"] if pdf_report_status["status"] == "created" else None,
        "self_absorption": self_absorption_summary,
    }
