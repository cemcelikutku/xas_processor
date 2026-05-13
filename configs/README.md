# Beamtime configs

Pre-prepared AstraConfig JSON files for the upcoming K K-edge operando electrochemistry experiment at ASTRA/SOLARIS. These are ready-to-use templates without paired example data. For reproducible end-to-end examples with paired .xasd data files, see [examples/](../examples/).

## Available configs

- **`k_k_operando_pd_foil.json`** — primary, separate_foil mode. Potassium K-edge (3608.4 eV) XANES
  in fluorescence mode at the ASTRA beamline (SOLARIS), with Pd foil
  scans in transmission used as the energy-drift reference. Intended
  for operando electrochemistry runs where applied potential varies
  between scans.
- **`k_k_operando_inline_ref.json`** — fallback, inline_ref mode. Alternate config for K K
  beamlines that measure an inline reference channel (I₁/I₂) during
  every sample scan instead of collecting separate foil scans. Use
  this when the beamline does NOT provide separate foil files and
  instead provides a reference signal measured simultaneously with
  each sample. Not applicable to the current ASTRA operando setup
  (which uses separate Pd foil scans), but kept as a template for
  other K K experiments.

## ⚠ Caution: verify which Pd L-edge will be measured

This config assumes the Pd foil is measured at the **L₁ edge
(3604 eV)**, which sits ~4 eV below K K and provides a clean drift
reference. Verify with beamline staff before the experiment which
Pd L-edge is actually being measured.

If they confirm **L₂ (3330 eV)** or **L₃ (3173 eV)** instead of L₁,
this config's alignment approach will not work correctly:

- The `align_window_min/max` (3588.4 to 3648.4 eV) is centered on
  L₁ and won't bracket the L₂ or L₃ derivative peaks.
- More fundamentally, the assumption that monochromator drift at
  the foil edge equals drift at the K K edge breaks down when the
  edges are hundreds of eV apart instead of ~4 eV apart.

In that case, the alignment_window in this config is wrong and
the drift correction approach is questionable: monochromator
drift over hundreds of eV is not necessarily uniform. Possible
responses: (a) ask beamline staff if Pd L₁ scans are also
feasible (most accurate fix); (b) widen the alignment window to
bracket the actual measured foil edge and accept that the K K
drift correction will be approximate; (c) post-process the data
manually using the foil scan information as a known reference.
The `k_k_operando_inline_ref.json` config is NOT a fix for this
scenario — it applies to a different beamline configuration
(inline reference channel) that is not available in the current
ASTRA spectroelectrochemical cell setup.

## K K + Pd L₁ alignment rationale

The Pd L₁ absorption edge sits at **3604 eV**, only ~4 eV below the
K K edge at **3608.4 eV**. For monochromator-drift correction, what
matters is that the foil edge falls within (or close to) the same
mono-angle range used for the sample edge: at this small separation,
the angular drift of the mono at K K is — to within experimental
uncertainty — the same drift seen at Pd L₁. That makes Pd foil an
excellent inline-energy standard for K K-edge work, even though it
isn't the K element itself.

The alignment window is set to ±30 eV around K E0
(`align_window_min: 3588.4`, `align_window_max: 3648.4`), which spans
both edge positions, and `shift_bound_min/max` is tightened to ±2 eV
because mono drift on a single beamtime should be well below that.

## Parameters that may need adjustment at the beamtime

- **`norm1` / `norm2`** (post-edge normalization range, currently
  E0+100 to E0+200 eV): real ASTRA scans may not extend to E0+200 eV,
  especially for short operando scans. If your scan range tops out
  earlier, lower `norm2` to match (and lower `norm1` proportionally
  to keep a reasonable post-edge window).
- **`foil_keyword`** (currently `"foil"`): if the Pd foil files at the
  beamtime are named with a different convention (e.g. `pd_ref_*` or
  the element name), change this to a substring that uniquely
  identifies them. Otherwise AstraXAS will not find the foil scans
  and will raise "No foil files found."

## Usage

From the CLI:

```
astra-xas /path/to/data --config configs/k_k_operando_pd_foil.json
```

From Python:

```python
from pathlib import Path
from astra_xas._config_utils import load_config_json
cfg = load_config_json(Path("configs/k_k_operando_pd_foil.json"))
```
