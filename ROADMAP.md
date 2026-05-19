# AstraXAS Roadmap

This document tracks the post-beamtime development roadmap for AstraXAS following the May 2026 beamtime at ASTRA/SOLARIS. It is a living document — updated as work progresses and assumptions are validated or revised.

## Background

AstraXAS reached v0.4.1 (Sprint 4) by May 2026 with:
- Working offline pipeline validated on Fe K, Si K, P K, and K K-edge data
- Preview Beamtime Mode for live monitoring during measurements
- Installable Python package with CLI, GUI, and Python API
- Documented operating procedure for beamtime use

A real beamtime at SOLARIS revealed both the strengths of the offline pipeline (publication-quality merged spectra on real K K-edge data) and limitations of Beamtime Mode in production conditions. Feedback from a beamline scientist after a 1.5-hour demo provided concrete user-requested features.

This roadmap responds to both the lessons learned and the user feedback, with a target of a SoftwareX publication describing the tool and its design.

## Strategic direction

Three architectural shifts inform this roadmap:

1. **Beamtime Mode as live logbook, not live processor.** Real-time analysis was the original vision; experience suggests "real-time logbook with progressive lightweight QC" is more useful and more achievable. Full scientific processing belongs offline.

2. **Session manifest as the unifying artifact.** A single curated record of which scans to process, how to group them, and what comments apply — produced by Beamtime Mode, consumed by the offline pipeline, editable by users between the two.

3. **Interoperability over reimplementation.** Don't try to replace Athena/Artemis. Make AstraXAS produce outputs that flow cleanly into those tools (XDI format, Larch-compatible data, etc.).

## Publication target

SoftwareX paper describing AstraXAS's design and validation, submitted approximately 6 months after the start of Sprint 5 work. The paper's contribution is the workflow design (live monitoring → curation → offline processing via session manifest) and the validation case studies on real synchrotron data, not the algorithms themselves.

## Phase 1: Session manifest foundation

**Goal:** Introduce the session manifest concept as the central abstraction. Make the offline pipeline able to read manifests as input alongside the current "folder" mode.

**Sub-tasks:**

- Design the manifest schema (JSON format). Fields include: session metadata, config path, scan list with per-scan group/replicate/status/comments. Schema lives in `astra_xas/manifest.py`.
- Implement `load_manifest()` and `save_manifest()` functions with schema validation.
- Modify `process_folder()` to optionally accept a manifest path. When provided, only process scans marked as `status="good"`, use manifest-supplied groupings, ignore filename-based grouping.
- Add CLI flag: `astra-xas --session session.json -o output/`.
- Maintain backward compatibility: existing folder-based mode still works.
- Add tests covering manifest parsing, validation, and offline pipeline integration.
- Document manifest format in a new `docs/manifest.md` file.

**Acceptance criteria:**

- An offline run can be driven by a manifest file
- Manifest schema is documented and validated
- All existing tests pass; new tests cover manifest functionality
- A real K K-edge dataset from the beamtime can be processed via a hand-edited manifest

**Why this phase first:** Defines the schema that Phase 2 will produce. Avoids designing Beamtime Mode in a way that's incompatible with offline consumption.

## Phase 2: Beamtime Mode rework

**Goal:** Refactor Beamtime Mode around the session manifest. Apply lessons from the SOLARIS beamtime: handle slow file writes properly, distinguish meaningful from noisy diagnostics, support comments and annotations.

**Sub-tasks:**

- Refactor the watcher to produce a session manifest as scans arrive (auto-populated defaults: status="auto-detected", inferred grouping from filename, no comments).
- Implement robust file-completion handling. Combines the work-in-progress branch `watcher-file-write-robustness` with the Option D refinement (longer stability window, minimum file size threshold). Test against synthetic slow-write fixtures.
- Reframe dashboard status meaning. "warn" should not fire on every scan due to XANES feature false positives in the detector jump diagnostic. Specifically: detector jump diagnostics become informational only; status reflects data acquisition QC (does it parse, are channels finite, is there an edge feature), not scientific result QC.
- Add comments support. Dashboard allows users to attach text comments to scans or to the session as a whole. Comments persist via the manifest, survive dashboard regeneration.
- Investigate subdirectory handling. Currently the watcher ignores subdirectories; some beamlines (including ASTRA) organize scans into subfolders. Decide whether to support recursive watching or improve the error message.
- Add structured output files to live mode. At minimum, write `ASTRA_detector_jumps.dat` incrementally so users can investigate warnings without running the offline pipeline.

**Acceptance criteria:**

- A real beamtime can be monitored without systematic rejections (regression test: the SOLARIS scenarios we know failed should now succeed)
- The dashboard's status field is informative — not every scan is `warn`
- Comments can be attached to scans and persist correctly
- The session manifest produced by Beamtime Mode is valid input for offline processing
- Tests cover both happy-path live processing and the edge cases we discovered (slow writes, aborted scans, sub-minimum-size files)

**Why this phase second:** The infrastructure from Phase 1 enables this work. The manifest is now the central data structure; Beamtime Mode produces it, offline consumes it.

## Phase 3: Session curation UI

**Goal:** Make manifest editing accessible to non-technical users. The user shouldn't have to edit JSON by hand to mark scans as good or bad.

**Sub-tasks:**

- Extend the Beamtime dashboard with curation actions: mark scan good/bad, edit comment, reassign group, mark for exclusion.
- Persist curation state to the session manifest as the user clicks.
- Add a "review session" view that summarizes the manifest state: N good scans, M excluded, N groups, etc.
- Add a "generate offline command" action that prints the appropriate `astra-xas --session ...` command for the curated manifest.

**Acceptance criteria:**

- A user can curate a session entirely through the dashboard without editing JSON
- The curated manifest validates and runs cleanly through the offline pipeline
- The workflow from "scan arrives at beamline" to "curated processed output" is end-to-end usable

**Why this phase third:** Phases 1 and 2 produce a working system; Phase 3 makes it pleasant. Order matters — better to have curation work via manifest editing first, then add UI, than to design UI without a clear data model.

## Phase 4: Channel selection and analysis flexibility

**Goal:** Implement the beamline scientist's most thoughtful feature request — let users explicitly choose numerator and denominator channels for the analysis signal, Athena-style.

**Sub-tasks:**

- Refactor `AstraConfig` to express analysis signal as `analysis_numerator`, `analysis_denominator`, and `analysis_log` fields rather than a single `analysis_mode` string.
- Maintain backward compatibility: parse old configs with `analysis_mode="trans"` etc. into the new representation.
- Add validation that the chosen signal is computable (e.g., division-by-zero protection, sign handling for `log()`).
- Update GUI to expose channel selection (matching Athena's pattern as closely as feasible).
- Update offline pipeline alignment code to use whatever analysis signal the user chose, including non-standard combinations.
- Document the new configuration format in README and migration guide.

**Acceptance criteria:**

- A user can choose any sensible channel combination (`I0/I1`, `IF/I0`, `I1/I2`, etc.) and processing works correctly
- The "I0 only" use case that surfaced at SOLARIS is naturally supported as `numerator=I0, denominator=1.0, log=false`
- Existing configs continue to work unchanged

**Why this phase here:** Foundational change to the data model. Belongs before Phase 5+ feature additions that depend on this flexibility.

## Phase 5: Interoperability outputs

**Goal:** Make AstraXAS outputs flow cleanly into the broader XAS tool ecosystem.

**Sub-tasks:**

- Add XDI (XAS Data Interchange format) output for processed spectra. Specification at https://github.com/XraySpectroscopy/XAS-Data-Interchange.
- Add a configurable option to also emit Larch-compatible HDF5 (using XDI internally) for direct loading in Athena/Artemis.
- Add an `astra-xas-export` CLI command that converts existing `.dat` outputs to XDI without re-running the pipeline.

**Acceptance criteria:**

- An AstraXAS-processed dataset can be loaded directly in Athena
- The XDI files round-trip cleanly (load, save, compare)
- Documentation explains the interoperability story

**Why this phase here:** Independent of the architectural work above. Could be done in parallel by a collaborator, or slotted in whenever convenient.

## Phase 6: Detector deadtime correction

**Goal:** Implement deadtime correction for fluorescence detectors (the beamline scientist's request #2).

**Sub-tasks:**

- Implement standard deadtime correction formulas (paralyzable and non-paralyzable). Use the existing FDT column in `.xasd` files.
- Add config options for deadtime mode and detector type.
- Validate against known datasets where deadtime is significant.
- Document the correction in README and add to the PDF report.

**Acceptance criteria:**

- Deadtime correction can be enabled/disabled via config
- Validation case shows correction produces expected behavior on high-count-rate data
- Inactive by default; existing configs unaffected

## Phase 7: GUI integration for Beamtime Mode

**Goal:** Make Beamtime Mode accessible through the existing GUI (the beamline scientist's request #6).

**Sub-tasks:**

- Add a "Beamtime Mode" tab/panel to the existing GUI.
- Embed or link to the live dashboard.
- Provide GUI controls for starting/stopping the watcher, choosing the watch folder, selecting config.
- Surface manifest curation directly in the GUI.

**Acceptance criteria:**

- A user who has only ever used the GUI can run Beamtime Mode without dropping to the command line.
- The GUI and CLI versions produce identical manifest output.

**Why this phase last:** Polish work, not capability work. Better done once everything else is stable.

## Cross-cutting work

These items should happen alongside the phases above, not after:

### Testing
- The offline pipeline has 0 pytest coverage. Add tests during Phase 1 work for `process_folder()` and its components.
- Each phase adds tests for its own functionality.
- By end of Phase 3, target ≥60% code coverage on the core processing modules.

### CI/CD
- Add GitHub Actions for pytest on every push and PR.
- Run on both Linux and Windows.
- Target: passing CI badge by end of Phase 2.

### PyPI publication
- Publish AstraXAS to PyPI after Phase 2 completes (when the major Beamtime Mode rework has settled).
- Pre-release/alpha versions during development, full release at v0.5.0 or v1.0.0 when Phase 4 completes.

### Documentation
- Each phase updates README.md and the appropriate docs/ files.
- BEAMTIME_USAGE.md gets a "live QC observability" section after Phase 2.
- A new `docs/manifest.md` after Phase 1.
- A new `docs/CONTRIBUTING.md` at some point to enable collaboration.

### Validation
- Process the K K-edge data from the May 2026 beamtime through the new pipeline once Phase 2 is complete. This becomes the third validation case study in the README (alongside Fe K and Si K).

## Known deferred items

These came up during Sprint 5 discussions but are not yet scheduled:

- **Analysis signal QC plot fix.** Currently averages across all scans regardless of group — misleading when multiple groups present. Decided post-beamtime fix. Probably Phase 4 cross-cutting.
- **Detector jump diagnostic recalibration.** Current MAD-based detector overcounts real XANES features on low-energy edges. Phase 2 handles the surfacing of this (informational not status), but the underlying algorithm may need rework — deferred until we understand whether the algorithm is fixable or should be replaced.
- **Multi-config sessions.** Currently a beamtime session uses one config; in practice, scientists may change between sample types and want different configs. Deferred; consider during Phase 2 design.
- **Session boundary detection.** Right now a "session" is defined by when the watcher was running. Better definitions (time-based, config-based) deferred.
- **Outlier detection improvements.** Current outlier rejection is conservative; some users may want more aggressive options. Deferred.

## Updating this document

This roadmap is updated:
- At the start of each phase (clarify what we're about to do)
- At the end of each phase (record what actually happened, lessons learned)
- Whenever a major decision changes scope (e.g., publication target shifts, new beamtime feedback arrives)

The git history of this file serves as the development log of the project. Commits to this file should briefly explain what changed and why.

## Current status

- **Date:** 2026-05-15
- **Last beamtime:** SOLARIS K K-edge, May 13-14, 2026
- **Current version:** v0.4.1
- **Active phase:** Pre-Phase-1 (planning complete, implementation not yet started)
- **Next milestone:** Phase 1 design document and initial implementation
