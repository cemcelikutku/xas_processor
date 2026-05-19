# Phase 2.1 — Dependency analysis for extracting `single_scan.py`

This document maps what `process_single_scan()` needs to call, traces transitive dependencies into a proposed file-ownership plan, and flags risks for the Session 2 implementer. **No files were modified.**

The target API is unchanged from the prompt:

```python
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

def process_single_scan(scan, config, filename=None) -> SingleScanResult: ...
```

`qc_status` derivation: `"reject"` if `qc_errors`; else `"warn"` if `qc_warnings`; else `"ok"`. Detector jumps are informational and do NOT influence `qc_status`. (This is a deliberate behavioural break from the watcher's current rule of `status="warn" if warnings or n_jumps`. The watcher is unchanged in Phase 2.1, so the divergence is bounded to the new path.)

---

## 1. Functions in `processor.py` that `process_single_scan()` will need to call

| Function | Lines | What it does | Shape |
|---|---|---|---|
| `_entry_from_scan(scan, config, path=None)` | 1259–1279 | Builds the canonical entry dict from a raw `load_xasd` scan: runs `compute_signals`, derives `base_name` / `replicate_id` from filename via `split_replicate_suffix`, sets `is_foil` from `config.foil_keyword`, copies raw channels (I0/I1/I2/IF/FDT/Ir). | **Per-scan** (pure). |
| `_analysis_signal_spec(config)` | 129–135 | Returns `(signal_key, signal_label)` for the analysis signal (`"mu_trans"/"ln(I0/I1)"`, etc.) based on `config.analysis_mode`. | **Per-scan** (config-only). |
| `_required_channels_for_signal(mode)` | 149–154 | Returns the list of `(channel_key, label)` pairs required for a given signal mode. | **Per-scan** (pure). |
| `_channel_validation_messages(entry, key, label, require_positive=False)` | 173–198 | Per-channel validation: missing/no-finite/all-zero/non-positive/nearly-flat. Returns `(fatal_errors, warnings)`. | **Per-scan**. |
| `_alignment_structure_warning(entry, signal_key, signal_label, config)` | 201–236 | Checks that the alignment signal has enough finite points and structure (derivative range) inside `config.align_window`. Returns a warning string or `None`. | **Per-scan**. |
| `_range_overlap_status(name, selected_range, data_range)` | 157–170 | Compares a config range (plot/align/pre-edge/norm) against an actual energy range; returns "does not overlap" / "partially overlaps" / `None`. | **Per-scan** when the data_range comes from one scan; reused by the batch validator on aggregated energies. |
| `_alignment_signal_spec(alignment_source, config)` | 138–146 | Returns `(mode, signal_key, signal_label)` for foil/inline alignment. Used by the validator's alignment-structure check. | **Per-scan** (config-only). |
| `_validate_processing_inputs(entries, config)` | 239–307 | Batch validator: aggregates energies → data_range → range-overlap checks; then per-entry channel + alignment-structure checks (with alignment-source-aware fanout). De-duplicates messages. | **Batch-shaped** (see §4). |
| `detect_detector_jumps(energy, channel, channel_name, config, filename)` | 310–414 | Single-channel jump detector. Returns a list of jump record dicts with severity, energy, jump_size, window flags. | **Per-scan, per-channel** (pure, no I/O). |
| `_config_bool(config, name, default=False)` | 50–54 | Tolerant truthy lookup. | **Per-scan** utility. |
| `_config_float(config, name, default)` | 57–61 | Tolerant float lookup. | **Per-scan** utility. |
| Constants: `RAW_DETECTOR_JUMP_CHANNELS`, `PRIMARY_DETECTOR_JUMP_CHANNELS`, `FDT_DETECTOR_JUMP_CHANNELS`, `DERIVED_DETECTOR_JUMP_CHANNELS` | 37–40 | Channel-class taxonomy used both by the per-scan jump scan and by the batch annotation function. | Per-scan taxonomy. |

`process_single_scan()` will also need a small loop over raw + derived channels that calls `detect_detector_jumps` for each. This loop is duplicated **twice today**: once in `process_folder` (lines ~1885–1935) and once in `beamtime/watcher.py::_count_detector_jumps` (with `_add_derived_channels`). Section 3 proposes consolidating it in `single_scan.py`.

---

## 2. Transitive dependencies

Stopping at module boundaries.

- `_entry_from_scan` → `compute_signals` (signals.py), `split_replicate_suffix` (io.py), `Path` (stdlib), `config.foil_keyword` / `config.fluo_multiplicative_constant`.
- `_analysis_signal_spec`, `_alignment_signal_spec`, `_required_channels_for_signal`, `_range_overlap_status` → no internal calls; pure functions of their args / `config`.
- `_channel_validation_messages` → numpy only; reads `entry[key]`, `entry["filename"]`.
- `_alignment_structure_warning` → numpy only; reads `entry["energy"]`, `entry[signal_key]`, `entry["filename"]`, `config.align_window`.
- `_validate_processing_inputs` → calls `_required_channels_for_signal`, `_channel_validation_messages`, `_alignment_structure_warning`, `_range_overlap_status`, `_alignment_signal_spec` (implicit via the inline mapping it duplicates inside its else branch). Reads `config.analysis_mode`, `config.alignment_source`, `config.foil_alignment_mode`, `config.plot_energy_*`, `config.align_window`, `config.e0`, `config.pre1/pre2/norm1/norm2`.
- `detect_detector_jumps` → numpy only; calls `_config_float` for two thresholds. Reads `config.e0`, `config.plot_energy_*`, `config.align_window_*`, `config.pre1/pre2/norm1/norm2`.
- `_config_bool`, `_config_float` → stdlib only.

**No transitive path from any of these reaches** `process_folder`, `alignment.find_best_shift`, `_run_pre_edge`, the grouping module, `save_two_col`, `plotting.*`, `self_absorption.*`, or `manifest.*`. The per-scan subgraph is cleanly bounded.

External deps used by the per-scan subgraph: `numpy`, `pathlib.Path`, `astra_xas.config.AstraConfig`, `astra_xas.signals.{compute_signals, get_signal}`, `astra_xas.io.split_replicate_suffix`.

---

## 3. Proposed file ownership after Phase 2.1

### Moves to `single_scan.py`

- `SingleScanResult` (new dataclass, exposed as public API)
- `process_single_scan(scan, config, filename=None)` (new public function)
- `_entry_from_scan` (rename to `_entry_from_scan` kept, OR expose as `entry_from_scan` — recommend **keep the leading underscore** to avoid widening the public surface; callers already import the private name)
- `_analysis_signal_spec`
- `_alignment_signal_spec`
- `_required_channels_for_signal`
- `_range_overlap_status`
- `_channel_validation_messages`
- `_alignment_structure_warning`
- `detect_detector_jumps` (public function, also moves)
- `_config_bool`, `_config_float` (small utilities; the alternative is a separate `_config_utils.py` but those already exist for JSON loading — keep these alongside the validators that use them)
- Constants: `RAW_DETECTOR_JUMP_CHANNELS`, `PRIMARY_DETECTOR_JUMP_CHANNELS`, `FDT_DETECTOR_JUMP_CHANNELS`, `DERIVED_DETECTOR_JUMP_CHANNELS`
- **New helper**: `detect_jumps_for_entry(entry, config) -> list[dict]` — applies `detect_detector_jumps` over the raw + derived channels for one entry. Consolidates logic currently duplicated in `process_folder` and `watcher._count_detector_jumps`. (Watcher is not modified in Phase 2.1, so its copy stays; future phase consolidates.)
- **New helper**: `_add_derived_channels(channels: dict) -> None` — moved from `watcher.py` only as a *copy* (watcher is not modified per the rule). The watcher's copy stays; `single_scan.py` owns its own copy. They are identical today; Phase 2.2 would unify.
- **New helper**: `validate_single_scan(entry, config) -> tuple[list[str], list[str]]` (warnings, fatal_errors). See §4.

### Stays in `processor.py`

- `process_folder` and all of its inline pipeline logic.
- `_validate_processing_inputs` (batch validator) — see §4.
- Batch detector-jump annotation/aggregation: `_detector_jump_sort_key`, `_annotate_detector_jump_summary_inclusion`, `_detector_jump_category`, `_detector_jump_event_map`, `_assign_event_ids`, `_write_detector_jump_rows`, `_write_detector_jumps`, `_detector_jump_energy_regions`.
- Alignment + anchor loading: `_load_selected_alignment_anchor`, `_same_file`.
- Group-merging utilities: `_run_pre_edge`, `_nanmean_without_warning`, `_safe_attr`.
- Spectrum-modifying batch operations: `auto_deglitch`, `manual_deglitch_range`.
- I/O & report writers: `save_detector_raw_entry`, `_write_alignment_metadata_header`, `_write_pdf_qc_report`, `_relative_output_path`, `_short_pdf_text`, `_read_dat_table_preview`, `_read_detector_jump_pdf_summary`.
- Constants tied to reporting: `AUTO_DEGLITCH_WARNING`, `SHIFT_CONVENTION`.
- Utility: `interpolate_to_grid` (used by group merging, not by `process_single_scan`).
- Helper: `_as_clean_list` (only used to parse the bulk-exclusion config fields, which are folder-mode concern).

### Compatibility re-exports from `processor.py` (zero-effort)

Add a re-export block near the top of `processor.py`:

```python
from .single_scan import (
    _entry_from_scan,
    detect_detector_jumps,
    _analysis_signal_spec,
    _alignment_signal_spec,
    _required_channels_for_signal,
    _range_overlap_status,
    _channel_validation_messages,
    _alignment_structure_warning,
    _config_bool,
    _config_float,
    PRIMARY_DETECTOR_JUMP_CHANNELS,
    FDT_DETECTOR_JUMP_CHANNELS,
    RAW_DETECTOR_JUMP_CHANNELS,
    DERIVED_DETECTOR_JUMP_CHANNELS,
)
```

Callers that import these names from `processor` (watcher.py uses `_entry_from_scan, _validate_processing_inputs, detect_detector_jumps`; beamtime/groups.py uses `_entry_from_scan, _run_pre_edge, interpolate_to_grid`) continue to work without modification. `_validate_processing_inputs` stays defined in `processor.py` directly — no re-export needed.

---

## 4. Validation strategy

**Recommendation: Option 3 — create a new per-scan validation helper in `single_scan.py` that preserves current behaviour exactly for the single-scan case.**

### What that looks like

```python
# in single_scan.py
def validate_single_scan(entry, config) -> tuple[list[str], list[str]]:
    """Per-scan validation. Returns (warnings, fatal_errors).
    Mirrors _validate_processing_inputs called with [entry] but skips
    the cross-entry aggregation and dedup that only matter for batch runs.
    """
```

It would:

1. Compute `data_range` from this one entry's finite energies and run `_range_overlap_status` against the four config ranges (plot, align, pre-edge, norm) — identical to what the batch path does with a single-entry energy list.
2. Skip foils (mirroring `if entry.get("is_foil"): continue`).
3. Run `_channel_validation_messages` per required channel for `analysis_mode`, with `require_positive` matching the batch logic.
4. Run the alignment-source-aware fanout (`inline_ref` vs `separate_foil`) — same as the batch path, but operating on this one entry only. Foil entries in `separate_foil` mode get the alignment checks; sample entries in `inline_ref` mode get them.
5. Dedup the resulting warnings/errors with `list(dict.fromkeys(...))` for parity.

### Why Option 3 over Option 1 or 2

- **Option 1 (move `_validate_processing_inputs` entirely)** is wrong because the batch function is genuinely batch-shaped — it aggregates energies across all entries to compute one `data_range`. Calling it with `[entry]` works today only because list-of-one is a trivial degenerate case. Moving the whole thing to `single_scan.py` puts a function with batch interface in a module whose job is per-scan; would force `process_folder` to call back into `single_scan` for batch logic, blurring ownership.
- **Option 2 (split the function)** sounds clean but the per-entry and batch parts are interleaved (range checks come first, then per-entry loops, with the alignment fanout repeating per-entry inside the foil branch). Cleanly splitting it would require restructuring the function meaningfully — the kind of change with a high risk of subtle behaviour drift in folder mode (which is the load-bearing path validated by the K K-edge dataset and existing tests).
- **Option 3** preserves the batch function untouched. The new per-scan helper is built from the same primitives (`_range_overlap_status`, `_channel_validation_messages`, `_alignment_structure_warning`) that the batch function uses, all of which move to `single_scan.py`. So both validators share the leaf checks; only the orchestration differs. This is the lowest-risk change.

### Risks of Option 3

- **Drift between the two validators.** If someone later adds a new warning class to one and forgets the other, behaviour diverges. Mitigation: Session 2 should add a regression test that asserts `validate_single_scan(entry, config) == _validate_processing_inputs([entry], config)` for representative cases (clean scan, missing channel, near-flat channel, foil scan in separate_foil mode, etc.). Cheap to write, catches drift.
- **The watcher's `_validate_single_scan` (in watcher.py line 33) still calls the batch function via re-export.** That's fine — watcher is unchanged, and the batch function behaves correctly on a singleton list. A future phase can switch watcher to call the new per-scan helper.

---

## 5. Circular import check

Proposed import graph after Phase 2.1:

```
single_scan.py
    imports: numpy, dataclasses, pathlib, .config, .signals, .io
    imported by: processor.py (direct), beamtime/watcher.py (transitively via processor re-exports)

processor.py
    imports: ..., .single_scan, .manifest, .grouping, .alignment, .signals,
             .edge_presets, .self_absorption, .plotting, .export, ._config_utils
    imported by: cli.py, beamtime/watcher.py, beamtime/groups.py

watcher.py (unchanged)
    imports: processor (gets re-exported per-scan helpers transparently)

beamtime/groups.py (unchanged)
    imports: processor (_entry_from_scan re-exported; _run_pre_edge and
             interpolate_to_grid stay in processor)
```

**Walk-through:**

- `single_scan.py` imports only `numpy`, stdlib, `.config`, `.signals`, `.io`. None of those import from `processor` (verified: `signals.py` imports `.config` only; `io.py` imports stdlib + numpy; `config.py` imports only `_version`).
- `processor.py` imports `single_scan` — fine, downstream of it.
- `watcher.py` and `beamtime/groups.py` import from `processor` — fine, they're consumers of the re-exports.
- `manifest.py` does not touch any of this graph.

**No cycles.** The single new edge added is `processor.py → single_scan.py`; everything else is unchanged.

---

## 6. Risks and unknowns for Session 2

1. **The `filename` parameter override** (`process_single_scan(scan, config, filename=None)`) interacts with `_entry_from_scan`, which reads `scan["filename"]` directly. The contract says "no mutation of caller-owned scan data," so when `filename` is provided the implementer must build a shallow-copy of `scan` with the override, not mutate the input dict in place. Same applies to the `is_foil` derivation inside `_entry_from_scan`.

2. **`_entry_from_scan` always reads `scan["filename"]`.** If `filename=None` and `scan["filename"]` is missing, today's code raises `KeyError`. `process_single_scan` should either require `scan["filename"]` to be present OR require `filename` to be passed; document this explicitly. (The watcher passes a scan dict from `load_xasd`, which guarantees the filename key.)

3. **`is_foil` is derived from filename inside `_entry_from_scan`** via `config.foil_keyword`. There is no overriding hook in the current function. The Phase 1.2 manifest integration overrides `is_foil` *after* `_entry_from_scan` returns, in `_build_entries_from_manifest`. `process_single_scan` will need to do likewise if it ever needs to honour a caller-supplied `is_foil` — but the documented signature has no such parameter, so v1 of `single_scan.py` can stick with filename-based detection. Flag this if the manifest-mode integration ever wants to call `process_single_scan`.

4. **Detector jump iteration** spans raw channels (`I0`, `I1`, `I2`, `IF`, `FDT`) plus three derived channels (`IF_over_I0`, `ln_I0_I1`, `ln_I1_I2`). The derived-channel computation lives in two places today: inlined in `process_folder` (lines 1899–1917) and as `_add_derived_channels` in `watcher.py`. The watcher version handles divide-by-zero via `np.where(I0 > 0, ..., np.nan)`. The `process_folder` inline copy is identical. Session 2 should reuse a single helper (in `single_scan.py`) and *not* try to also fix the duplicate in watcher.py (out of scope per the rule).

5. **`_validate_processing_inputs` reads many config fields** and has implicit assumptions about alignment-mode mappings (it duplicates the mapping that `_alignment_signal_spec` provides). The duplication looks intentional (the function evolved before `_alignment_signal_spec` existed). Resist the urge to clean it up while moving leaf helpers; folder-mode behaviour parity is the priority.

6. **`_alignment_structure_warning` in `separate_foil` mode is invoked on foil entries only.** In `inline_ref` mode it runs on every sample entry. This is not symmetric. The new `validate_single_scan` must preserve this asymmetry — i.e., a non-foil entry called with `alignment_source="separate_foil"` should NOT emit alignment-structure warnings, but a non-foil entry called with `alignment_source="inline_ref"` SHOULD. Mirror the batch logic exactly.

7. **`detect_detector_jumps` requires ≥20 finite points** and a non-zero stddev. Short scans silently return `[]`. `metrics["n_points_finite"]` is a useful field to expose in `SingleScanResult.metrics` so a caller can disambiguate "no jumps because clean" from "no jumps because the scan was too short to evaluate."

8. **The proposed `qc_status` rule excludes jumps from the status decision** ("Detector jumps do NOT automatically affect qc_status"). The watcher today sets `status="warn" if warnings or n_jumps` — this matches the roadmap Phase 2 plan ("detector jump diagnostics become informational only; status reflects data acquisition QC"). The new path adopts the future behaviour; the watcher remains on the legacy behaviour until Phase 2 watcher rework. Calling code in Session 2 should not silently expect the old behaviour.

9. **`metrics: dict`** is unspecified. Recommend Session 2 populate at minimum:
   - `n_points`: int
   - `n_points_finite_energy`: int
   - `energy_range_eV`: tuple[float, float] or None
   - `channels_present`: list[str] (channels with at least one finite value)
   - `analysis_signal_finite_fraction`: float in [0, 1]

   Keep it a flat dict (no nesting), keys snake_case, values JSON-serializable so the manifest's `qc_summary` can later mirror a subset of these.

10. **`scan` dict shape contract.** `load_xasd` produces keys `{path, filename, energy, theta, dt, I0, I1, I2, IF, FDT, Ir}`. `process_single_scan` must tolerate the FDT/Ir being `None` (they are when the .xasd has fewer than 9 columns) — `_entry_from_scan` already handles this via `.get`. The new metrics dict should not assume FDT is present.

11. **Removing `_entry_from_scan` from `processor.py`'s top-level definitions** while keeping the name importable (via re-export) means `from astra_xas.processor import _entry_from_scan` still works for `watcher.py` and `beamtime/groups.py`. Anyone using `astra_xas.processor._entry_from_scan` through a fully-qualified attribute lookup (e.g., in tests via `monkeypatch.setattr("astra_xas.processor._entry_from_scan", ...)`) will also still work — `from .single_scan import _entry_from_scan` binds the name on the processor module. Grep before merging to be safe; `tests/` does not currently appear to patch this attribute.

---

## Bottom line for Session 2

- The per-scan subgraph is cleanly extractable. No circular import risk.
- 13 functions / 4 constants move; 3 new helpers are added (`SingleScanResult`, `process_single_scan`, `validate_single_scan`, plus `detect_jumps_for_entry` and a private `_add_derived_channels`).
- `_validate_processing_inputs` stays put; `validate_single_scan` is the new per-scan path.
- A handful of compatibility re-exports in `processor.py` make `watcher.py` and `beamtime/groups.py` work without modification.
- The biggest hidden gotcha: the alignment-structure asymmetry between `separate_foil` and `inline_ref` modes in the batch validator. Mirror it exactly.
- Smallest follow-up debt: `_add_derived_channels` will be duplicated between `single_scan.py` and `watcher.py` until Phase 2 watcher rework. Acceptable.
