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
- **Detector raw export** — saves all raw detector channels (I0, I1, I2, IF, FDT, Ir) alongside processed outputs, plottable directly in the Spectrum Viewer
- **Automatic plots** — overview plots for processed μ(E), background-corrected, and normalized spectra; per-group replicate QC plots; optional energy drift tracker
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
pip install numpy scipy matplotlib xraylarch
```

On Ubuntu, `tkinter` may need to be installed separately:

```bash
sudo apt install python3-tk
```

**Main dependencies:**

- `numpy`, `scipy` — signal processing and alignment
- `xraylarch` — `pre_edge` normalization (Athena-compatible)
- `matplotlib` — automatic and interactive plots
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
| `-o / --output-dir` | `<input>-processed` | Output folder |
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
    ├─ Compute μ(E) per scan (IF/I0, ln(I0/I1), or ln(I1/I2))
    ├─ Align each scan to foil reference via derivative-shape matching
    ├─ Score alignment quality and record energy drift
    ├─ [Optional] Deglitch aligned replicates before merging
    ├─ [Optional] Reject outlier replicates
    │
    ├─ Average aligned μ(E) replicates  ← merge first
    ├─ Run Larch pre_edge on merged spectrum  ← normalize once
    │
    └─ Output
         ├─ <sample>_processed.dat   (merged processed μ(E))
         ├─ <sample>_bkgcorr.dat     (background-subtracted)
         ├─ <sample>_norm.dat        (normalized μ(E))
         ├─ <sample>_flat.dat        (flattened normalized μ(E))
         ├─ detector_raw/<scan>.dat  (all detector channels)
         ├─ plots/                   (overview and QC plots)
         ├─ ASTRA_energy_shifts.dat
         ├─ ASTRA_foil_alignment.dat
         └─ ASTRA_processing_report.txt
```

---

## Alignment and drift tracking

AstraXAS supports two alignment sources:

- `inline_ref` uses the reference channel (`ln(I1/I2)`) measured in each sample scan. The first sample scan is the zero-shift reference.
- `separate_foil` uses files whose names contain `foil_keyword` as reference foil scans. The first foil scan is the zero-shift reference, and each sample inherits the shift and quality of its most recent assigned foil scan.

Alignment is performed in the configured energy window (`align_window_min/max`) and bounded by `shift_bound_min/max`. The current engine sanitizes non-finite points, sorts spectra by energy, removes duplicate moving-spectrum energies, checks that the reference has usable derivative amplitude, then searches for the best shift using a coarse grid followed by local optimization. The moving spectrum is interpolated with a cubic spline before derivative comparison.

Each computed alignment writes:

- `shift_eV` — energy shift applied to the scan or foil.
- `fit_error` — residual mean-squared error between z-scored derivatives.
- `alignment_quality` — Pearson correlation between z-scored derivatives at the best shift. Values near `1.0` are reliable; values below `alignment_quality_warn_threshold` are suspect.

If alignment cannot be evaluated because the reference or moving spectrum is unusable, AstraXAS sets `shift_eV = 0.0`, `fit_error = NaN`, and `alignment_quality = 0.0`, then records an explicit warning.

Set `save_drift_plot=True` to write `plots/drift_tracker.png`. Reliable scans are shown as filled blue circles, while low-quality scans are shown as open red circles. Red dashed horizontal lines mark `±warn_shift_abs_eV`.

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
| `detector_raw/<scan>.dat` | Raw detector channels for every individual scan |
| `plots/replicate_qc/<sample>_replicate_qc.png` | Replicate overlay QC plot |
| `plots/aligned_averaged_IF_overview.png` | Optional aligned/interpolated/averaged IF detector signal by group |
| `plots/processed_mu_overview.png` | All merged processed μ(E) spectra overlaid |
| `plots/background_corrected_overview.png` | All background-corrected spectra overlaid |
| `plots/normalized_overview.png` | All normalized spectra overlaid |
| `plots/drift_tracker.png` | Optional scan-by-scan energy shift plot, written when `save_drift_plot=True` |
| `<sample>_deglitch_log.dat` | Deglitch point log, written only when deglitching modifies points |
| `ASTRA_processing_report.txt` | Full parameter log, per-group summary, low-quality alignment count, warnings, and deglitch point counts |
| `ASTRA_energy_shifts.dat` | Per-sample shift table: filename, base name, replicate id, assigned foil/reference, shift, alignment quality |
| `ASTRA_foil_alignment.dat` | Per-foil or inline-reference alignment table: filename, shift, fit error, alignment quality |
| `ASTRA_normalization_summary.dat` | Edge step, E₀, and normalization metadata per group |
| `ASTRA_group_summary.dat` | Sample names, foil assignments, replicate counts |

All `.dat` files have a commented header listing parameters and column names, and spectral `.dat` files are directly loadable in the Spectrum Viewer.

Note: `plots/detector_raw_overview.png` was the former name for the aligned averaged IF overview. Current runs write only `plots/aligned_averaged_IF_overview.png`.

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

All processing parameters are exposed in the GUI and saveable as JSON config files. Common parameters:

| Parameter | Description |
|---|---|
| `e0` | Edge energy in eV |
| `pre1`, `pre2` | Pre-edge fit range relative to E₀ (eV) |
| `norm1`, `norm2` | Post-edge normalization range relative to E₀ (eV) |
| `nnorm` | Normalization polynomial order (0, 1, or 2) |
| `align_window_min/max` | Energy window used for foil alignment |
| `shift_bound_min/max` | Maximum allowed energy shift during alignment (eV) |
| `alignment_quality_warn_threshold` | Quality threshold below which alignment warnings are emitted; default `0.7` |
| `alignment_grid_points` | Number of coarse-search grid points before local alignment refinement; default `50` |
| `fluo_multiplicative_constant` | Scaling factor applied to IF before computing μ(E) |
| `enable_auto_deglitch` | Interpolate isolated narrow detector spikes before merging |
| `deglitch_threshold` | Robust local threshold for automatic spike detection |
| `deglitch_window` | Local half-window used by automatic deglitching |
| `deglitch_min_energy/max_energy` | Optional energy bounds for automatic deglitching |
| `enable_manual_deglitch_range` | Interpolate a specified energy range before merging |
| `manual_deglitch_min_energy/max_energy` | Energy interval for manual range interpolation |
| `manual_deglitch_margin_points` | Neighboring points used to interpolate the manual range |
| `enable_auto_outlier_detection` | Flag replicates that deviate from the group mean |
| `outlier_rms_threshold` | RMS deviation threshold for outlier detection |
| `enable_shift_rejection` | Exclude replicates with large energy shifts |
| `reject_shift_abs_eV` | Shift threshold for rejection (eV) |
| `save_detector_raw_overview_plot` | Save the aligned averaged IF overview plot; legacy config name retained for compatibility |
| `save_processed_overview_plot` | Save an overview plot of merged processed μ(E) spectra |
| `save_bkgcorr_overview_plot` | Save an overview plot of background-corrected spectra |
| `save_norm_overview_plot` | Save an overview plot of normalized spectra |
| `save_replicate_qc_plots` | Save per-group replicate QC plots |
| `save_drift_plot` | Save `plots/drift_tracker.png` |
| `plot_energy_min/max` | Energy range used for automatic overview and QC plots |

In the GUI, the **Enable deglitching** checkbox and `automatic` / `manual` / `both` mode selector are translated into `enable_auto_deglitch` and `enable_manual_deglitch_range` in the saved configuration.

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
