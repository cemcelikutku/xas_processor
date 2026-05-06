from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
from larch import Group
from larch.xafs import pre_edge

from .config import AstraConfig
from .io import collect_xasd_files, load_xasd, split_replicate_suffix, sanitize_name, build_default_output_directory
from .signals import compute_signals, get_signal
from .alignment import find_best_shift
from .grouping import group_samples
from .export import save_two_col
from .plotting import plot_overview, plot_replicate_qc


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
    replicate_qc_dir = plots_dir / "replicate_qc"
    plots_enabled = (
        getattr(config, "save_detector_raw_overview_plot", False)
        or getattr(config, "save_processed_overview_plot", getattr(config, "save_raw_overview_plot", True))
        or getattr(config, "save_bkgcorr_overview_plot", False)
        or getattr(config, "save_norm_overview_plot", True)
        or getattr(config, "save_replicate_qc_plots", True)
    )
    if plots_enabled:
        plots_dir.mkdir(parents=True, exist_ok=True)
        if getattr(config, "save_replicate_qc_plots", True):
            replicate_qc_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    log(f"Starting {config.version}")
    log(f"Input directory: {input_dir}")
    log(f"Output directory: {output_dir}")

    files = collect_xasd_files(input_dir)
    if not files:
        raise RuntimeError(f"No .xasd files found in: {input_dir}")

    entries = []
    for path in files:
        scan = load_xasd(path)
        sigs = compute_signals(scan, config)
        is_foil = config.foil_keyword.lower() in scan["filename"].lower() if config.foil_keyword else False
        base_name, replicate_id = split_replicate_suffix(scan["filename"])
        entries.append({
            "path": path,
            "filename": scan["filename"],
            "is_foil": is_foil,
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
        })

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

    alignment_source = getattr(config, "alignment_source", "separate_foil")

    if alignment_source == "separate_foil":
        foil_entries = [e for e in entries if e["is_foil"]]
        if not foil_entries:
            raise RuntimeError(
                f"No foil files found. Put '{config.foil_keyword}' in foil filenames "
                f"or change foil_keyword."
            )
        first_foil = foil_entries[0]
        E_foil_ref = first_foil["energy"]
        foil_ref_signal = get_signal(first_foil, config.foil_alignment_mode)
        foil_shift_map = {first_foil["filename"]: {"shift_eV": 0.0, "fit_error": 0.0}}
        log(f"Reference foil: {first_foil['filename']} shift = 0.0000 eV")

        for foil in foil_entries[1:]:
            shift, err = find_best_shift(
                E_foil_ref,
                foil_ref_signal,
                foil["energy"],
                get_signal(foil, config.foil_alignment_mode),
                window=config.align_window,
                bounds=config.shift_bounds,
            )
            foil_shift_map[foil["filename"]] = {"shift_eV": shift, "fit_error": err}
            if abs(shift) > config.warn_shift_abs_eV:
                warnings.append(f"Large foil shift: {foil['filename']} = {shift:.4f} eV")
            log(f"Foil {foil['filename']}: shift = {shift:.4f} eV, fit_error = {err:.6g}")

        current_shift = 0.0
        current_foil_name = first_foil["filename"]
        sample_before_first_foil = True
        for e in entries:
            if e["is_foil"]:
                current_foil_name = e["filename"]
                current_shift = foil_shift_map[current_foil_name]["shift_eV"]
                sample_before_first_foil = False
            elif sample_before_first_foil:
                warnings.append(f"Sample appears before first foil and will use first foil shift: {e['filename']}")
            e["assigned_foil"] = current_foil_name
            e["energy_shift_eV"] = current_shift
        sample_entries = [e for e in entries if not e["is_foil"]]

    elif alignment_source == "inline_ref":
        sample_entries = [e for e in entries if not e["is_foil"]]
        if not sample_entries:
            raise RuntimeError("No sample files found.")
        ref_scan = sample_entries[0]
        E_ref = ref_scan["energy"]
        mu_ref_signal = ref_scan["mu_ref"]
        foil_entries = []
        foil_shift_map = {ref_scan["filename"]: {"shift_eV": 0.0, "fit_error": 0.0}}
        ref_name = f"inline_ref:{ref_scan['filename']}"
        log(f"Inline reference scan: {ref_scan['filename']} shift = 0.0000 eV")

        for e in sample_entries:
            if e is ref_scan:
                shift, err = 0.0, 0.0
            else:
                shift, err = find_best_shift(
                    E_ref,
                    mu_ref_signal,
                    e["energy"],
                    e["mu_ref"],
                    window=config.align_window,
                    bounds=config.shift_bounds,
                )
            e["assigned_foil"] = ref_name
            e["energy_shift_eV"] = shift
            foil_shift_map[e["filename"]] = {"shift_eV": shift, "fit_error": err}
            if abs(shift) > config.warn_shift_abs_eV:
                warnings.append(f"Large inline reference shift: {e['filename']} = {shift:.4f} eV")
            log(f"Inline ref {e['filename']}: shift = {shift:.4f} eV, fit_error = {err:.6g}")
    else:
        raise RuntimeError(f"Unknown alignment_source='{alignment_source}'. Use 'separate_foil' or 'inline_ref'.")

    if not sample_entries:
        raise RuntimeError("No sample files found.")

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
        f.write("# filename\tbase_name\treplicate_id\tassigned_foil\tenergy_shift_eV\n")
        for s in sample_entries:
            rep = "None" if s["replicate_id"] is None else str(s["replicate_id"])
            f.write(f"{s['filename']}\t{s['base_name']}\t{rep}\t{s['assigned_foil']}\t{s['energy_shift_eV']:.6f}\n")

    with open(output_dir / "ASTRA_foil_alignment.dat", "w", encoding="utf-8") as f:
        f.write("# foil_filename\tshift_eV\tfit_error\n")
        for foil_name, info in foil_shift_map.items():
            f.write(f"{foil_name}\t{info['shift_eV']:.6f}\t{info['fit_error']:.8g}\n")

    group_summary_rows = []
    group_norm_metadata_rows = []
    detector_raw_plot_records = []
    processed_plot_records = []
    bkgcorr_plot_records = []
    norm_plot_records = []
    plot_energy_range = (config.plot_energy_min, config.plot_energy_max)

    for (base_name, assigned_foil), scans in groups.items():
        log(f"Processing group: {base_name} ({len(scans)} scan(s))")
        first_scan = scans[0]
        master_energy = first_scan["energy"] + first_scan["energy_shift_eV"]
        processed_spectra = []
        source_filenames = []
        detector_group_if = []

        for s in scans:
            shifted_energy = s["energy"] + s["energy_shift_eV"]
            mu = get_signal(s, config.analysis_mode)
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

        comments = (
            f"# ASTRA XAS Processor version: {config.version}\n"
            f"# Group base name: {base_name}\n"
            f"# Assigned foil: {assigned_foil}\n"
            f"# Analysis mode: {config.analysis_mode}\n"
            f"# Alignment source: {alignment_source}\n"
            f"# Normalization order nnorm: {config.nnorm}\n"
            f"# Normalize-before-average: False (merge-then-normalize)\n"
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

        if getattr(config, "save_replicate_qc_plots", True) and len(norm_array) > 1:
            try:
                normalized_replicate_records = [
                    {"energy": master_energy, "mu": norm_array[i], "label": source_filenames[i]}
                    for i in range(len(source_filenames))
                ]
                plot_replicate_qc(
                    group_name=output_base,
                    replicate_records=normalized_replicate_records,
                    average_record={"energy": master_energy, "mu": norm_avg, "label": "average"},
                    output_path=replicate_qc_dir / f"{output_base}_replicate_qc.png",
                    energy_range=plot_energy_range,
                    y_label="Normalized intensity",
                )
                log(f"Saving replicate QC plot: {output_base}")
            except Exception as exc:
                warning = f"Could not save replicate QC plot for {output_base}: {exc}"
                warnings.append(warning)
                log("WARNING: " + warning)

        group_summary_rows.append([output_base, base_name, assigned_foil, str(len(processed_array)), source_filenames[0], source_filenames[-1]])
        std_str = "N/A" if not np.isfinite(std_edge_step) else f"{std_edge_step:.8g}"
        group_norm_metadata_rows.append([output_base, f"{mean_edge_step:.8g}", std_str, f"{mean_norm_e0:.8g}", str(config.nnorm)])

    plot_files = []
    if plots_enabled:
        log(f"Generating plots in: {plots_dir}")

    if getattr(config, "save_detector_raw_overview_plot", False):
        path = plot_overview(detector_raw_plot_records, plots_dir / "detector_raw_overview.png", "Detector raw overview", "IF detector signal", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_processed_overview_plot", getattr(config, "save_raw_overview_plot", True)):
        path = plot_overview(processed_plot_records, plots_dir / "processed_mu_overview.png", "Processed μ(E) overview", "Processed μ(E)", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_bkgcorr_overview_plot", False):
        path = plot_overview(bkgcorr_plot_records, plots_dir / "background_corrected_overview.png", "Background-corrected overview", "Background-corrected μ(E)", energy_range=plot_energy_range)
        if path is not None:
            plot_files.append(path); log(f"Saved plot: {path}")

    if getattr(config, "save_norm_overview_plot", True):
        path = plot_overview(norm_plot_records, plots_dir / "normalized_overview.png", "Normalized overview", "Normalized intensity", energy_range=plot_energy_range)
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
        f.write("# output_base\tbase_name\tassigned_foil\tn_scans_used\tfirst_scan\tlast_scan\n")
        for row in group_summary_rows:
            f.write("\t".join(row) + "\n")

    with open(output_dir / "ASTRA_normalization_summary.dat", "w", encoding="utf-8") as f:
        f.write("# output_base\tmean_edge_step\tstd_edge_step\tmean_e0\tnnorm\n")
        for row in group_norm_metadata_rows:
            f.write("\t".join(row) + "\n")

    with open(output_dir / "ASTRA_processing_report.txt", "w", encoding="utf-8") as f:
        f.write(f"{config.version}\n")
        f.write(f"Input directory: {input_dir}\nOutput directory: {output_dir}\n")
        f.write(f"Files found: {len(files)}\n")
        f.write(f"Files excluded manually: {len(excluded_entries)}\n")
        f.write(f"Files rejected by shift safety: {len(auto_shift_rejected_entries)}\n")
        f.write(f"Replicate outliers skipped: {len(auto_outlier_entries)}\n")
        f.write(f"Foils found: {len(foil_entries)}\n")
        f.write(f"Groups processed: {len(group_summary_rows)}\n")
        f.write(f"Normalization order nnorm: {config.nnorm}\n")
        f.write("Normalize-before-average: False (merge μ(E) first, normalize merged spectrum)\n")
        f.write(f"Plots folder: {plots_dir}\n")
        f.write(f"Detector raw folder: {output_dir / 'detector_raw'}\n")
        f.write(f"Detector raw files created: {len(detector_raw_files)}\n")
        f.write("Processed μ(E) files use *_processed.dat; true detector channels use detector_raw/*_detector_raw.dat.\n")
        f.write(f"Plots created: {len(plot_files)}\n")
        for plot_file in plot_files:
            f.write(f"- {plot_file}\n")
        f.write("\nWarnings:\n")
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
        "groups_processed": len(group_summary_rows),
        "manual_excluded": len(excluded_entries),
        "shift_rejected": len(auto_shift_rejected_entries),
        "auto_outliers": len(auto_outlier_entries),
    }
