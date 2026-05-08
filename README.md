# AstraXAS

**Open-source XAS preprocessing and visualization toolkit originally developed for ASTRA beamline `.xasd` data at the SOLARIS Synchrotron.**

AstraXAS provides automated workflows for X-ray absorption spectroscopy (XAS) preprocessing, including foil drift correction, alignment quality checks, optional deglitching, replicate alignment, scan merging, Athena-style normalization, automatic QC plots, and interactive spectrum visualization. Although originally developed around the ASTRA beamline at the SOLARIS Synchrotron, the processing workflow is adaptable to other XAS beamlines provided that compatible detector channels and energy-resolved scan formats are available.

---

## Features

- **Automatic foil drift correction** — aligns each scan to a reference foil (inline I₂ or separate foil files) using derivative-shape matching, a coarse global grid search, and local refinement
- **Alignment quality scoring** — reports a Pearson-r quality score for every alignment and warns when scans fall below a configurable quality threshold
- **Merge-then-normalize workflow** — averages raw μ(E) replicates first, then applies normalization once to the merged spectrum, following common Athena-style XAS preprocessing practice
- **Athena-compatible normalization** — uses Larch's `pre_edge` with full control over pre-edge range, normalization range, polynomial order, and E₀
- **Three analysis modes** — fluorescence (`IF/I0`), transmission (`ln(I0/I1)`), and reference (`ln(I1/I2)`)
- **Two alignment sources** — inline reference channel (I₂ measured in every scan) or separate foil files identified by a filename keyword
- **Pre-merge deglitching** — optional automatic interpolation of narrow detector spikes and manual range interpolation for broader inspected artifacts
- **Automatic outlier detection** — optionally flag and exclude replicates that deviate from the group mean by a configurable RMS threshold
- **Shift rejection** — optionally exclude replicates whose energy shift exceeds a threshold before merging
- **Validation warnings** — reports missing/weak channels, incompatible modes, and energy-window issues before processing decisions are made
- **Detector jump diagnostics** — optional diagnostic-only raw detector spike reporting with a structured `ASTRA_detector_jumps.dat` file
- **Self-absorption flag** — optional fluorescence-mode heuristic that flags likely self-absorption when the fluorescence white-line amplitude is suppressed relative to simultaneously available sample transmission; QC flag only, no correction
- **Detector raw export** — saves all raw detector channels (I0, I1, I2, IF, FDT, Ir) alongside processed outputs, plottable directly in the Spectrum Viewer
- **Automatic plots** — detector health overview, analysis signal QC, processed μ(E), background-corrected, and normalized overview plots; pre-normalization and normalized replicate QC plots; optional energy drift tracker
- **Edge presets** — editable starting templates for common ASTRA-accessible edges
- **PDF QC report** — optional single-file processing and quality-control report generated from existing outputs
- **Interactive Spectrum Viewer** — compare any `.dat` files side by side with Savitzky-Golay smoothing, raw/smoothed overlay, legend toggling, click-to-read energy values, and publication-ready figure export
- **JSON config system** — save and load processing parameters per edge or experiment type
- **CLI and Python API** — run headless from the command line or call `process_folder()` from a script

---

## Installation

AstraXAS requires Python 3.10+.

Clone the repository:

```bash
git clone https://github.com/cemcelikutku/AstraXAS.git
cd AstraXAS
```

Install the required Python dependencies manually:

```bash
pip install numpy scipy matplotlib xraylarch reportlab
```

On Ubuntu, `tkinter` may need to be installed separately:

```bash
sudo apt install python3-tk
```

**Main dependencies:**

- `numpy`, `scipy` — signal processing and alignment
- `xraylarch` — `pre_edge` normalization (Athena-compatible)
- `matplotlib` — automatic and interactive plots
- `reportlab` — optional PDF QC report generation
- `tkinter` — GUI toolkit

---

## Quick start

### GUI

```bash
python -m astra_xas.gui
```

### Command line

```bash
python -m astra_xas.cli /path/to/xasd/folder --mode fluo --e0 7121.030
```

Options:

| Flag | Default | Description |
|---|---|---|
| `input_dir` | required | Folder containing `.xasd` files |
| `-o`, `--output-dir` | `<input>-processed` | Output folder |
| `--mode` | `fluo` | Sample signal: `fluo`, `trans`, or `ref` |
| `--foil-mode` | `trans` | Foil alignment signal: `trans`, `ref`, or `fluo` |
| `--e0` | `7121.030` | Edge energy in eV |

### Python API

```python
from astra_xas import AstraConfig, process_folder

config = AstraConfig(
    analysis_mode="fluo",
    alignment_source="inline_ref",
    e0=7121.030,
    pre1=-229.74, pre2=-49.98,
    norm1=55.07,  norm2=227.22,
    nnorm=1,
    alignment_quality_warn_threshold=0.7,
    alignment_grid_points=50,
    save_drift_plot=True,
    plot_energy_min=7100.0,
    plot_energy_max=7160.0,
)

result = process_folder("/path/to/xasd/folder", config=config)
print(result["output_dir"])
```

---

## Workflow

```
.xasd files
    │
    ├─ Load detector channels (E, I0, I1, I2, IF, FDT, Ir)
    ├─ Run validation warnings on channels, modes, and energy windows
    ├─ Align each scan to foil reference via derivative-shape matching
    ├─ Score alignment quality and record energy drift
    ├─ [Optional] Reject scans with large alignment shifts
    ├─ [Optional] Run detector jump diagnostics on raw channels
    ├─ Compute μ(E) per scan (IF/I0, ln(I0/I1), or ln(I1/I2))
    ├─ [Optional] Deglitch aligned replicates before merging
    ├─ [Optional] Reject outlier replicates inside each group
    │
    ├─ Average aligned μ(E) replicates  ← merge first
    ├─ Run Larch pre_edge on merged spectrum  ← normalize once
    ├─ [Optional] Self-absorption flag (fluorescence mode only)
    │
    └─ Output
         ├─ <sample>_processed.dat   (merged processed μ(E))
         ├─ <sample>_bkgcorr.dat     (background-subtracted)
         ├─ <sample>_norm.dat        (normalized μ(E))
         ├─ <sample>_flat.dat        (flattened normalized μ(E))
         ├─ detector_raw/<scan>_detector_raw.dat  (all detector channels)
         ├─ plots/overview/          (dataset-level overview and QC plots)
         ├─ plots/replicate_qc/      (scan-to-scan replicate QC plots)
         ├─ ASTRA_energy_shifts.dat
         ├─ ASTRA_foil_alignment.dat
         ├─ ASTRA_processing_report.txt
         └─ [Optional] ASTRA_processing_and_QC_report.pdf
```

---

## Alignment and drift tracking

AstraXAS supports two alignment sources:

- `inline_ref` uses the reference channel (`ln(I1/I2)`) measured in each sample scan. The first sample scan is the zero-shift reference.
- `separate_foil` uses files whose names contain `foil_keyword` as reference foil scans. The first foil scan is the zero-shift reference, and each sample inherits the shift and quality of its most recent assigned foil scan.

By default, `alignment_anchor_mode="first_scan"` preserves this behavior. To compare multiple folders against the same internal reference, set `alignment_anchor_mode="selected_file"` and provide `alignment_anchor_path` pointing to a `.xasd` scan. The selected anchor file is loaded, validated, and used as the zero-shift alignment anchor for all scans in the folder. This separates the alignment source (`inline_ref` or `separate_foil`) from the alignment anchor (the file that defines zero shift).

Alignment is performed in the configured energy window (`align_window_min/max`) and bounded by `shift_bound_min/max`. The current engine sanitizes non-finite points, sorts spectra by energy, removes duplicate moving-spectrum energies, checks that the reference has usable derivative amplitude, then searches for the best shift using a coarse grid followed by local optimization. The moving spectrum is interpolated with a cubic spline before derivative comparison.

Each computed alignment writes:

- `shift_eV` — energy shift applied to the scan or foil.
- `fit_error` — residual mean-squared error between z-scored derivatives.
- `alignment_quality` — Pearson correlation between z-scored derivatives at the best shift. Values near `1.0` are reliable; values below `alignment_quality_warn_threshold` are suspect.

If alignment cannot be evaluated because the reference or moving spectrum is unusable, AstraXAS sets `shift_eV = 0.0`, `fit_error = NaN`, and `alignment_quality = 0.0`, then records an explicit warning.

Set `save_drift_plot=True` to write `plots/overview/drift_tracker.png`. Reliable scans are shown as filled blue circles, while low-quality scans are shown as open red circles. Red dashed horizontal lines mark `±warn_shift_abs_eV`.

The processing report and alignment shift tables record alignment source, alignment signal, alignment anchor mode, selected anchor file if any, whether the anchor loaded successfully, and the shift sign convention. A selected anchor provides a shared internal energy reference; it does not guarantee absolute energy calibration unless the anchor itself was externally calibrated.

---

## Detector health overview

When `save_detector_health_overview_plot=True`, AstraXAS writes `plots/overview/detector_health_overview.png`. This is a stacked diagnostic PNG that plots individual sample scan traces before normalization, so beam drops, detector jumps, spikes, or unstable channels are easier to spot than in averaged spectra.

The channel set is mode-aware:

- Fluorescence mode: `I0`, `IF`, and `FDT` when available.
- Transmission mode: `I0` and `I1`.
- Reference mode: `I1`, `I2`, and `ln(I1/I2)`.
- If reference-channel data are available in fluorescence or transmission mode, `ln(I1/I2)` is added as an extra panel.

Missing or non-plottable channels are skipped gracefully. `ASTRA_processing_report.txt` records whether the detector health overview was created and lists included and skipped channels.

---

## Analysis signal QC

When `save_analysis_signal_qc_plot=True`, AstraXAS writes `plots/overview/analysis_signal_qc.png`. This plot shows the actual per-scan analysis signal before final normalization, with individual scan traces kept visible and a black diagnostic average trace added when multiple traces can be overlaid.

The plotted signal follows `analysis_mode`:

- Fluorescence mode: `IF/I0`
- Transmission mode: `ln(I0/I1)`
- Reference mode: `ln(I1/I2)`

This plot is intended to answer whether detector-channel artifacts cancel out or survive in the signal that is actually processed. Missing or non-finite traces are skipped gracefully, and `ASTRA_processing_report.txt` records the plotted signal, number of individual traces, whether an average trace was added, and any skipped traces.

---

## Processed μ(E) replicate QC

When `save_processed_mu_replicate_qc_plot=True`, AstraXAS writes one pre-normalization replicate QC plot per sample group:

```text
plots/replicate_qc/<sample>_processed_mu_replicate_qc.png
```

This plot uses the aligned and interpolated processed μ(E) replicate spectra on the common group grid, after optional deglitching and outlier filtering, but before final Larch pre-edge normalization or flattening. Individual scans are shown with a thicker average trace so scan-to-scan consistency can be checked before normalization has a chance to hide differences.

This is different from `plots/overview/processed_mu_overview.png`: the overview shows only averaged group spectra, while processed μ(E) replicate QC shows individual scans plus the group average for each sample.

---

## Validation warnings

Before alignment and merging, AstraXAS runs a diagnostic validation pass on the scans that remain after manual exclusions. This pass does not modify data, deglitch, reject scans, or change processing results. It only records warnings when selected modes or energy ranges look inconsistent with the available data.

Validation checks include:

- Required channels for `analysis_mode`: `IF` and `I0` for fluorescence, `I0` and `I1` for transmission, `I1` and `I2` for reference.
- Required reference channels for `inline_ref` alignment.
- Required foil-alignment channels in `separate_foil` mode.
- All-zero, all-non-finite, nearly flat, or non-positive channels where division/log calculations need positive values.
- Plot range, alignment window, pre-edge range, and normalization range overlap with the data energy range.
- Weak structure in the selected alignment signal inside the alignment window.

Warnings are printed in the processing log and written to a dedicated `Validation warnings` section in `ASTRA_processing_report.txt`. If no validation warnings are found, the report writes `Validation warnings: none`.

---

## Detector jump diagnostics

When `enable_detector_jump_warnings=True`, AstraXAS runs a diagnostic-only detector jump check on raw detector channels after alignment shifts are known and before any deglitching or averaging. It can write `ASTRA_detector_jumps.dat` when spike-like jumps are detected, and it adds a conservative summary-level detector-jump count to `ASTRA_group_summary.dat`.

This check uses point-to-point MAD thresholding plus recovery-window spike-vs-step discrimination, so monotonic absorption-edge-like steps are not treated as detector jumps. The detailed table can include primary raw-channel, FDT, and derived-signal sharp features, but the main processing report emphasizes significant primary raw-channel jumps (`I0`, `I1`, `I2`, `IF`) outside the edge/alignment window. FDT diagnostic spikes and derived-signal edge features are reported separately and excluded from the main summary. It never modifies detector arrays, processed spectra, normalized spectra, alignment shifts, or plot data.

---

## Self-absorption diagnostic

When `analysis_mode="fluo"` and `enable_self_absorption_check=True`, AstraXAS can flag possible fluorescence self-absorption / over-absorption. The diagnostic compares the normalized white-line amplitude of the fluorescence signal against the normalized white-line amplitude of the simultaneously available sample transmission signal, `ln(I0/I1)`, and flags groups where the fluorescence white-line amplitude is suppressed.

This is a heuristic QC flag, not a correction. It does not modify processed spectra, normalized spectra, alignment, averaging, deglitching, or exported spectral values. The ratio can also be influenced by sample geometry, detector effects, noise, normalization quality, and chemical or sampling differences, so flagged groups should be inspected by the user.

The diagnostic runs only for fluorescence-mode sample groups and only uses accepted scans: scans that survive manual exclusions, optional shift rejection, and group-level outlier filtering. The diagnostic transmission average is built from the same accepted scans as the final fluorescence average. If usable `I0`/`I1` transmission data are not available, or if diagnostic normalization cannot be evaluated reliably, the group is marked `skipped`.

Sensitivity controls the fluorescence/transmission white-line amplitude ratio threshold:

| Sensitivity | Threshold |
|---|---:|
| `relaxed` | `0.75` |
| `normal` | `0.85` |
| `strict` | `0.92` |
| `custom` | `self_absorption_custom_threshold` |

Classification is based on `ratio_fluo_over_trans = fluo_white_line_amplitude / trans_white_line_amplitude`. Ratios greater than or equal to the selected threshold are `ok`; ratios below the threshold are `flagged`; invalid or insufficient diagnostics are `skipped`. Flag severity uses absolute ratio bands independent of the selected sensitivity: `mild` for ratios `>= 0.75`, `moderate` for ratios `>= 0.60`, and `strong` below `0.60`.

White-line and continuum windows are anchored to `config.e0`. Defaults are:

| Window | Default offset relative to E₀ |
|---|---|
| White-line | `0.0` to `35.0` eV |
| Continuum | `50.0` to `150.0` eV |

The normalized diagnostic spectra are measured using a smoothed 95th percentile in the white-line window minus the median continuum level. If the continuum window is unavailable but the white-line window is usable, AstraXAS falls back to a continuum level of `1.0` and records that note. The default minimum transmission white-line amplitude is `0.03`, and the default minimum number of finite points in a diagnostic window is `5`.

When the diagnostic runs, AstraXAS writes `ASTRA_self_absorption_flags.dat` even if every group is skipped. If the analysis mode is not fluorescence, or if `enable_self_absorption_check=False`, this file is not created. The table columns are:

```text
sample  status  severity  ratio_fluo_over_trans  fluo_white_line_amplitude
trans_white_line_amplitude  threshold_used  sensitivity  white_line_window_eV
continuum_window_eV  n_replicates_used  note
```

Rows with `status=ok` include measured amplitudes and an `OK` note. Rows with `status=flagged` include the ratio, threshold, severity, and a note that the fluorescence white-line amplitude is suppressed relative to transmission. Rows with `status=skipped` use `nan` where amplitudes or ratios are unavailable and include a clear skip reason.

When `save_self_absorption_qc_plots=True`, checked groups with `status=ok` or `status=flagged` also write:

```text
plots/replicate_qc/<sample>_self_absorption_qc.png
```

The QC plot overlays the normalized fluorescence diagnostic spectrum and normalized sample transmission diagnostic spectrum, with the white-line and continuum windows marked. Flagged groups are also added to the standard `Processing warnings` list and summarized in both `ASTRA_processing_report.txt` and `ASTRA_processing_and_QC_report.pdf`.

The GUI exposes this feature in the self-absorption diagnostic section, including enable/disable, sensitivity, custom threshold, white-line and continuum offsets, minimum transmission amplitude, minimum points, and QC plot saving.

---

## Deglitching

AstraXAS includes optional deglitching for scan-level artifacts. Deglitching operates on each aligned replicate before replicate averaging. The merge-then-normalize workflow is preserved: corrected μ(E) replicates are merged first, and Larch `pre_edge` normalization is applied once to the merged spectrum.

Two deglitching modes are available:

- **Automatic deglitching** interpolates isolated, narrow point-like detector spikes. It preserves the original energy grid and is intended for single-point excursions, not structured spectral features.
- **Manual range deglitching** interpolates a user-defined energy interval from neighboring points. It is intended for broader artifacts selected after visual inspection.

Automatic deglitching should be used conservatively. Broad artifacts, beamline disturbances, or multi-point distortions should be handled with manual range interpolation.

GUI usage:

1. Enable **Deglitching**.
2. Choose `automatic`, `manual`, or `both`.
3. For automatic deglitching, set the threshold, window half-width, and optional energy bounds.
4. For manual range deglitching, set the artifact energy range and the number of margin points used for interpolation.

Python examples:

```python
from astra_xas import AstraConfig, process_folder

# Automatic: narrow detector spikes only
config = AstraConfig(
    enable_auto_deglitch=True,
    deglitch_threshold=5.0,
    deglitch_window=5,
    deglitch_min_energy=7100.0,
    deglitch_max_energy=7160.0,
)

process_folder("/path/to/xasd/folder", config=config)
```

```python
from astra_xas import AstraConfig, process_folder

# Manual range: broader artifact selected after inspection
config = AstraConfig(
    enable_manual_deglitch_range=True,
    manual_deglitch_min_energy=7132.0,
    manual_deglitch_max_energy=7135.0,
    manual_deglitch_margin_points=5,
)

process_folder("/path/to/xasd/folder", config=config)
```

When deglitching changes points, AstraXAS writes a per-group `<sample>_deglitch_log.dat`. `ASTRA_processing_report.txt` summarizes automatic and manual deglitched point counts separately.

---

## File naming convention

AstraXAS groups `.xasd` files into samples by filename. Replicates are identified by a trailing numeric suffix:

```
Fe_sample_1.xasd   →  base name: Fe_sample,  replicate: 1
Fe_sample_2.xasd   →  base name: Fe_sample,  replicate: 2
Fe_foil_1.xasd     →  foil file (contains "foil" keyword)
```

The `foil_keyword` parameter (default: `"foil"`) determines which files are treated as reference foil scans for alignment.

---

## Output files

For each sample group, AstraXAS writes the following to the output directory:

| File | Contents |
|---|---|
| `<sample>_processed.dat` | Merged processed μ(E) before background removal |
| `<sample>_bkgcorr.dat` | Background-corrected μ(E) |
| `<sample>_norm.dat` | Normalized μ(E) |
| `<sample>_flat.dat` | Flattened normalized μ(E) |
| `detector_raw/<scan>_detector_raw.dat` | Raw detector channels for every individual scan |
| `plots/replicate_qc/<sample>_normalized_replicate_qc.png` | Normalized replicate overlay QC plot |
| `plots/replicate_qc/<sample>_processed_mu_replicate_qc.png` | Pre-normalization processed μ(E) replicate QC plot |
| `plots/replicate_qc/<sample>_self_absorption_qc.png` | Optional fluorescence self-absorption QC plot for checked groups when `save_self_absorption_qc_plots=True` |
| `plots/overview/detector_health_overview.png` | Mode-aware stacked detector QC plot using individual sample scan traces |
| `plots/overview/analysis_signal_qc.png` | Mode-aware per-scan analysis signal before final normalization, with optional average overlay |
| `plots/overview/aligned_averaged_IF_overview.png` | Optional aligned/interpolated/averaged IF detector signal by group |
| `plots/overview/processed_mu_overview.png` | All merged processed μ(E) spectra overlaid |
| `plots/overview/background_corrected_overview.png` | All background-corrected spectra overlaid |
| `plots/overview/normalized_overview.png` | All normalized spectra overlaid |
| `plots/overview/drift_tracker.png` | Optional scan-by-scan energy shift plot, written when `save_drift_plot=True` |
| `<sample>_deglitch_log.dat` | Deglitch point log, written only when deglitching modifies points |
| `ASTRA_detector_jumps.dat` | Diagnostic detector jump records, written only when jump-like spikes are detected |
| `ASTRA_self_absorption_flags.dat` | Fluorescence-mode possible self-absorption diagnostic table, written when `analysis_mode="fluo"` and `enable_self_absorption_check=True` |
| `ASTRA_excluded_scans.dat` | Manually excluded scans and exclusion reasons |
| `ASTRA_auto_rejected_scans.dat` | Scans rejected by optional shift rejection |
| `ASTRA_auto_outliers.dat` | Replicates rejected by optional group-level outlier detection |
| `ASTRA_processing_report.txt` | Full parameter log, validation warnings, processing warnings, plot file lists, replicate QC counts, per-group summary, low-quality alignment count, and deglitch point counts |
| `ASTRA_energy_shifts.dat` | Per-sample shift table with alignment-anchor metadata: filename, base name, replicate id, assigned foil/reference, shift, alignment quality |
| `ASTRA_foil_alignment.dat` | Per-foil or inline-reference alignment table with alignment-anchor metadata: filename, shift, fit error, alignment quality |
| `ASTRA_normalization_summary.dat` | Edge step, E₀, and normalization metadata per group |
| `ASTRA_group_summary.dat` | Sample names, foil assignments, replicate counts, and detector-jump summary counts |
| `ASTRA_processing_and_QC_report.pdf` | Optional single-file processing and QC report with metadata, QC status, warnings, summary tables, and generated QC plots |

All `.dat` files have a commented header listing parameters and column names, and spectral `.dat` files are directly loadable in the Spectrum Viewer. `ASTRA_energy_shifts.dat` and `ASTRA_foil_alignment.dat` include the alignment anchor context so shifts can be interpreted across folders.

`ASTRA_processing_report.txt` separates `Validation warnings` from `Processing warnings`, summarizes detector jump and self-absorption diagnostics, reports processed μ(E) replicate QC and normalized replicate QC counts separately, and lists created plots under `Overview plots created` and `Replicate QC plots created`.

When `save_pdf_report=True`, AstraXAS also writes `ASTRA_processing_and_QC_report.pdf`. The PDF collects run metadata, a compact QC status summary, validation and processing warnings, detector jump and self-absorption diagnostics, alignment/drift context, overview plots, replicate QC plots, main output-file locations, and compact previews of summary tables. It is diagnostic/reporting only: it uses already-created outputs and does not recompute or modify spectra. Missing plots or tables are skipped in the PDF as “not available”; processing still completes if PDF generation fails, and the failure is recorded in `ASTRA_processing_report.txt`.

Note: `plots/detector_raw_overview.png` was the former name for the aligned averaged IF overview. Current runs write only `plots/overview/aligned_averaged_IF_overview.png`. The normalized replicate QC plot was formerly named `<sample>_replicate_qc.png`; current runs write `<sample>_normalized_replicate_qc.png`.

---

## Edge Presets

AstraXAS includes an editable edge-preset library for common starting templates: Mg K, Al K, Si K, P K, S K, Cl K, K K, Ca K, Ti K, V K, Cr K, Mn K, Fe K, Co K, Ni K, Cu K, Zn K, Ga K, Ge K, As K, Se K, and Ce L3. Presets are not final scientific parameters and are not absolute calibration. They provide approximate reference/start E0 values plus relative pre-edge, normalization, alignment, and plot windows.

In the GUI, choosing a preset in the dropdown does not change any fields by itself. Click **Apply preset** to copy the template into the editable fields. After applying, E0, pre-edge range, normalization range, alignment window, and plot range remain normal editable inputs. The Spectrum Viewer uses the current processing plot range as its default energy range when opened, so applying a preset also updates the viewer's initial range.

When a preset is applied, AstraXAS records the preset label, whether it was applied, the preset reference E0, the final E0 used, and the preset note in `ASTRA_processing_report.txt` and `ASTRA_processing_and_QC_report.pdf`.

For scripts, presets are available from `astra_xas.edge_presets`:

```python
from astra_xas import AstraConfig
from astra_xas.edge_presets import apply_edge_preset_to_config

config = AstraConfig()
apply_edge_preset_to_config(config, "p_k")
# Edit config.e0, config.pre1/pre2, config.norm1/norm2, etc. as needed.
```

Old configurations and custom/manual workflows behave the same unless a preset is explicitly applied.

---

## Spectrum Viewer

The built-in Spectrum Viewer can open any `.dat` file produced by AstraXAS — including detector raw files — or any two-column text file.

**Key features:**

- Load multiple files and compare them in a single interactive plot
- Select plot channel from raw detector signals (`I0`, `I1`, `I2`, `IF`, `FDT`, `Ir`) or computed signals (`IF/I0`, `ln(I0/I1)`, `ln(I1/I2)`)
- Savitzky-Golay smoothing with configurable window and polynomial order
- Show raw and smoothed curves simultaneously
- Click legend entries to hide/show individual spectra
- Hover to see interpolated values; left-click to print energy and signal values for all visible curves
- Press `a` to show all hidden spectra; `n` to clear clicked markers
- Save publication-ready figures (PNG, PDF, SVG) at 300 dpi
- Fully configurable axis labels, title, legend position, figure size, and line width
- Save and reload viewer sessions (file list + all settings) as JSON

---

## Configuration

Processing parameters are stored in `AstraConfig` and saveable as JSON config files. The GUI exposes the common workflow controls, while some lower-level diagnostic thresholds are config-file/API options. Common parameters:

| Parameter | Description |
|---|---|
| `version` | Version string written to reports and output headers |
| `foil_keyword` | Filename substring used to identify separate foil scans; default `foil` |
| `analysis_mode` | Sample analysis signal: `fluo`, `trans`, or `ref` |
| `alignment_source` | Alignment source: `inline_ref` or `separate_foil` |
| `alignment_anchor_mode` | Zero-shift anchor selection: `first_scan` or `selected_file` |
| `alignment_anchor_path` | Path to selected `.xasd` anchor file when `alignment_anchor_mode="selected_file"` |
| `foil_alignment_mode` | Signal used for separate foil alignment: `trans`, `ref`, or `fluo` |
| `fluo_multiplicative_constant` | Scaling factor applied to IF before computing fluorescence μ(E) |
| `e0` | Edge energy in eV |
| `edge_preset_key` | Reporting metadata for the applied edge preset; default `custom` |
| `edge_preset_label` | Human-readable preset label recorded in reports |
| `edge_preset_applied` | Whether a preset was explicitly applied; default `False` |
| `edge_preset_note` | Preset note recorded in reports |
| `pre1`, `pre2` | Pre-edge fit range relative to E₀ (eV) |
| `norm1`, `norm2` | Post-edge normalization range relative to E₀ (eV) |
| `nnorm` | Normalization polynomial order (0, 1, or 2) |
| `step`, `nvict`, `make_flat` | Additional Larch `pre_edge` options passed through by the processor |
| `align_window_min/max` | Energy window used for foil alignment |
| `shift_bound_min/max` | Maximum allowed energy shift during alignment (eV) |
| `alignment_quality_warn_threshold` | Quality threshold below which alignment warnings are emitted; default `0.7` |
| `alignment_grid_points` | Number of coarse-search grid points before local alignment refinement; default `50` |
| `interp_kind` | Interpolation kind used for scan interpolation to the group grid |
| `warn_shift_abs_eV` | Shift magnitude that triggers a processing warning |
| `exclude_filenames` | Explicit filenames to exclude before processing |
| `exclude_filename_contains` | Filename substrings to exclude before processing |
| `enable_shift_rejection` | Exclude replicates with large energy shifts |
| `reject_shift_abs_eV` | Shift threshold for rejection (eV) |
| `enable_auto_outlier_detection` | Flag replicates that deviate from the group mean |
| `outlier_rms_threshold` | RMS deviation threshold for outlier detection |
| `enable_auto_deglitch` | Interpolate isolated narrow detector spikes before merging |
| `deglitch_threshold` | Robust local threshold for automatic spike detection |
| `deglitch_window` | Local half-window used by automatic deglitching |
| `deglitch_method` | Deglitch method name; current processing uses interpolation |
| `deglitch_min_energy/max_energy` | Optional energy bounds for automatic deglitching |
| `enable_manual_deglitch_range` | Interpolate a specified energy range before merging |
| `manual_deglitch_min_energy/max_energy` | Energy interval for manual range interpolation |
| `manual_deglitch_margin_points` | Neighboring points used to interpolate the manual range |
| `enable_detector_jump_warnings` | Run diagnostic-only raw detector jump warnings; default `True` |
| `detector_jump_threshold` | Point-to-point MAD multiplier for detector jump warnings; default `10.0` |
| `detector_jump_min_relative` | Minimum relative jump size for detector jump warnings; default `0.05` |
| `enable_self_absorption_check` | Enable fluorescence-mode possible self-absorption diagnostic; default `True` |
| `self_absorption_sensitivity` | Threshold preset: `relaxed`, `normal`, `strict`, or `custom`; default `normal` |
| `self_absorption_custom_threshold` | Custom ratio threshold used when sensitivity is `custom`; default `0.85` |
| `self_absorption_wl_min` | White-line window minimum offset from E₀; default `0.0` eV |
| `self_absorption_wl_max` | White-line window maximum offset from E₀; default `35.0` eV |
| `self_absorption_cont_min` | Continuum window minimum offset from E₀; default `50.0` eV |
| `self_absorption_cont_max` | Continuum window maximum offset from E₀; default `150.0` eV |
| `self_absorption_min_trans_amp` | Minimum reliable sample transmission white-line amplitude; default `0.03` |
| `self_absorption_min_points` | Minimum finite points required in diagnostic windows; default `5` |
| `save_detector_health_overview_plot` | Save `plots/overview/detector_health_overview.png`; default `True` |
| `save_analysis_signal_qc_plot` | Save `plots/overview/analysis_signal_qc.png`; default `True` |
| `save_detector_raw_overview_plot` | Save `plots/overview/aligned_averaged_IF_overview.png`; legacy config name retained for compatibility |
| `save_processed_overview_plot` | Save `plots/overview/processed_mu_overview.png` |
| `save_bkgcorr_overview_plot` | Save `plots/overview/background_corrected_overview.png` |
| `save_norm_overview_plot` | Save `plots/overview/normalized_overview.png` |
| `save_processed_mu_replicate_qc_plot` | Save pre-normalization processed μ(E) replicate QC plots; default `True` |
| `save_replicate_qc_plots` | Save per-group replicate QC plots |
| `save_raw_overview_plot` | Legacy alias retained for old configs |
| `save_drift_plot` | Save `plots/overview/drift_tracker.png` |
| `save_foil_alignment_plots` | Compatibility field; current processing does not generate a separate foil-alignment plot |
| `save_pdf_report` | Generate `ASTRA_processing_and_QC_report.pdf`; default `True` when ReportLab is available |
| `save_self_absorption_qc_plots` | Save `plots/replicate_qc/<sample>_self_absorption_qc.png` for checked groups; default `True` |
| `plot_energy_min/max` | Energy range used for automatic overview and QC plots |

In the GUI, the **Enable deglitching** checkbox and `automatic` / `manual` / `both` mode selector are translated into `enable_auto_deglitch` and `enable_manual_deglitch_range` in the saved configuration. The **Generate PDF QC report** checkbox controls `save_pdf_report`.

The GUI defaults to `inline_ref` alignment because many ASTRA datasets include the reference channel in every scan. The `AstraConfig` dataclass default remains `separate_foil`, so scripts should set `alignment_source` explicitly when a specific workflow is required.

Save a config file for each edge (Fe K, Cu K, etc.) and load it at the start of a session.

---

## Positioning relative to Athena and Fastosh

| Feature | Athena | Fastosh | AstraXAS |
|---|---|---|---|
| Automatic foil drift correction | Manual | Manual | ✅ Automatic |
| Merge-then-normalize workflow | Manual | Manual | ✅ Default workflow |
| Pre-merge deglitching | Manual | Manual | ✅ Automatic or manual range |
| Detector raw channel plotting | ✗ | ✗ | ✅ Built-in |
| Replicate QC plots | Manual | ✗ | ✅ Automatic |
| Outlier / shift rejection | ✗ | ✗ | ✅ Configurable |
| Self-absorption flag | ✗ | ✗ | ✅ Heuristic QC flag (no correction) |
| JSON config per edge | ✗ | ✗ | ✅ Save and load |
| Python scripting API | ✗ | ✗ | ✅ `process_folder()` |
| CLI | ✗ | ✗ | ✅ |
| Interactive spectrum viewer | ✓ | Limited | ✅ With session save |
| Publication figure export | ✓ | ✓ | ✅ PNG / PDF / SVG |
| EXAFS fitting | ✓ (Artemis) | Limited | Planned |

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

AstraXAS uses [xraylarch](https://xraypy.github.io/xraylarch/) for Athena-compatible normalization via `pre_edge`. Alignment follows the derivative-matching approach described in the [Athena User's Guide](https://bruceravel.github.io/demeter/documents/Athena/index.html), with added grid search, spline interpolation, and quality scoring in AstraXAS.
