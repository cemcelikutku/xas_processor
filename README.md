# AstraXAS

**Open-source XAS preprocessing and visualization toolkit originally developed for ASTRA beamline `.xasd` data at the SOLARIS Synchrotron.**

AstraXAS provides automated workflows for X-ray absorption spectroscopy (XAS) preprocessing, including foil drift correction, replicate alignment, scan merging, Athena-style normalization, and interactive spectrum visualization. Although originally developed around the ASTRA beamline at the SOLARIS Synchrotron, the processing workflow is adaptable to other XAS beamlines provided that compatible detector channels and energy-resolved scan formats are available.

---

## Features

- **Automatic foil drift correction** — aligns each scan to a reference foil (inline I₂ or separate foil files) using derivative cross-correlation, with user-configurable energy windows and shift bounds
- **Merge-then-normalize workflow** — averages raw μ(E) replicates first, then applies normalization once to the merged spectrum, following common Athena-style XAS preprocessing practice
- **Athena-compatible normalization** — uses Larch's `pre_edge` with full control over pre-edge range, normalization range, polynomial order, and E₀
- **Three analysis modes** — fluorescence (`IF/I0`), transmission (`ln(I0/I1)`), and reference (`ln(I1/I2)`)
- **Two alignment sources** — inline reference channel (I₂ measured in every scan) or separate foil files identified by a filename keyword
- **Automatic outlier detection** — optionally flag and exclude replicates that deviate from the group mean by a configurable RMS threshold
- **Shift rejection** — optionally exclude replicates whose energy shift exceeds a threshold before merging
- **Detector raw export** — saves all raw detector channels (I0, I1, I2, IF, FDT, Ir) alongside processed outputs, plottable directly in the Spectrum Viewer
- **Automatic plots** — overview plots for raw μ(E), background-corrected, and normalized spectra; per-group replicate QC plots
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
    ├─ Align each scan to foil reference via derivative cross-correlation
    ├─ [Optional] Reject outlier replicates
    │
    ├─ Average aligned μ(E) replicates  ← merge first
    ├─ Run Larch pre_edge on merged spectrum  ← normalize once
    │
    └─ Output
         ├─ <sample>_raw.dat         (merged μ(E))
         ├─ <sample>_bkgcorr.dat     (background-subtracted)
         ├─ <sample>_norm.dat        (normalized μ(E))
         ├─ <sample>_flat.dat        (flattened normalized μ(E))
         ├─ detector_raw/<scan>.dat  (all detector channels)
         ├─ plots/                   (overview and QC plots)
         └─ ASTRA_processing_report.txt
```

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
| `<sample>_raw.dat` | Merged μ(E) before background removal |
| `<sample>_bkgcorr.dat` | Background-corrected μ(E) |
| `<sample>_norm.dat` | Normalized μ(E) |
| `<sample>_flat.dat` | Flattened normalized μ(E) |
| `detector_raw/<scan>.dat` | Raw detector channels for every individual scan |
| `plots/<sample>_qc.png` | Replicate overlay QC plot |
| `plots/overview_norm.png` | All normalized spectra overlaid |
| `plots/overview_processed.png` | All merged μ(E) spectra overlaid |
| `ASTRA_processing_report.txt` | Full parameter log and per-group summary |
| `ASTRA_normalization_summary.dat` | Edge step, E₀, and normalization metadata per group |
| `ASTRA_group_summary.dat` | Sample names, foil assignments, replicate counts |

All `.dat` files have a commented header listing parameters and column names, and are directly loadable in the Spectrum Viewer.

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
| `fluo_multiplicative_constant` | Scaling factor applied to IF before computing μ(E) |
| `enable_auto_outlier_detection` | Flag replicates that deviate from the group mean |
| `outlier_rms_threshold` | RMS deviation threshold for outlier detection |
| `enable_shift_rejection` | Exclude replicates with large energy shifts |
| `reject_shift_abs_eV` | Shift threshold for rejection (eV) |

Save a config file for each edge (Fe K, Cu K, etc.) and load it at the start of a session.

---

## Positioning relative to Athena and Fastosh

| Feature | Athena | Fastosh | AstraXAS |
|---|---|---|---|
| Automatic foil drift correction | Manual | Manual | ✅ Automatic |
| Merge-then-normalize workflow | Manual | Manual | ✅ Default workflow |
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

AstraXAS uses [xraylarch](https://xraypy.github.io/xraylarch/) for Athena-compatible normalization via `pre_edge`. Alignment is based on derivative cross-correlation following the approach described in the [Athena User's Guide](https://bruceravel.github.io/demeter/documents/Athena/index.html).
