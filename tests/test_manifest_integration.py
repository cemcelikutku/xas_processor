"""End-to-end tests for the offline pipeline driven by a session manifest.

These tests use synthetic .xasd files (not physical) to exercise the
process_folder + manifest plumbing without depending on real beamtime
data. The synthetic helper produces enough signal for the pipeline's
alignment and validation steps to succeed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from astra_xas.beamtime._synthetic import write_synthetic_xasd
from astra_xas.config import AstraConfig
from astra_xas.grouping import group_samples
from astra_xas.manifest import (
    ConfigRef,
    GroupEntry,
    InputSpec,
    Manifest,
    ScanEntry,
    save_manifest,
)
from astra_xas.processor import process_folder


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_BASE_CONFIG_OVERRIDES = {
    "alignment_source": "inline_ref",
    "save_pdf_report": False,
    "save_drift_plot": False,
    "save_self_absorption_qc_plots": False,
    "save_detector_health_overview_plot": False,
    "save_analysis_signal_qc_plot": False,
    "save_detector_raw_overview_plot": False,
    "save_processed_overview_plot": False,
    "save_bkgcorr_overview_plot": False,
    "save_norm_overview_plot": False,
    "save_processed_mu_replicate_qc_plot": False,
    "save_replicate_qc_plots": False,
    "enable_detector_jump_warnings": False,
    "enable_self_absorption_check": False,
}


def _fast_config_dict(**overrides) -> dict:
    out = dict(_BASE_CONFIG_OVERRIDES)
    out.update(overrides)
    return out


def _fast_config(**overrides) -> AstraConfig:
    cfg = AstraConfig()
    for k, v in _fast_config_dict(**overrides).items():
        setattr(cfg, k, v)
    return cfg


def _write_config_json(path: Path, **overrides) -> Path:
    path.write_text(json.dumps(_fast_config_dict(**overrides)), encoding="utf-8")
    return path


def _write_synthetic_files(directory: Path, names: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for i, name in enumerate(names):
        write_synthetic_xasd(directory / name, seed=i, e_max=7400.0)


def _capture_log():
    """Return (log_fn, log_list) so tests can both run silently and assert on log lines."""
    lines: list[str] = []

    def log(message=""):
        lines.append(str(message))

    return log, lines


def _read_source_scan_order(dat_path: Path) -> list[str]:
    """Parse a per-group _processed.dat header and return the source filenames in order."""
    found: list[str] = []
    in_block = False
    for raw in dat_path.read_text(encoding="utf-8").splitlines():
        if raw.startswith("# Source scans"):
            in_block = True
            continue
        if in_block:
            if raw.startswith("#   "):
                found.append(raw[len("#   "):].strip())
            elif raw.startswith("#"):
                # Another comment header after Source scans block ended
                continue
            else:
                break
    return found


# ---------------------------------------------------------------------------
# Test 1: backward compatibility — folder mode still works unchanged
# ---------------------------------------------------------------------------

def test_folder_mode_still_works_without_session(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["sample_001.xasd", "sample_002.xasd"])
    out = tmp_path / "out"
    log, lines = _capture_log()
    process_folder(
        input_dir=data_dir, output_dir=out, config=_fast_config(), log=log
    )
    report = (out / "ASTRA_processing_report.txt").read_text(encoding="utf-8")
    assert "Driven by manifest" not in report
    assert (out / "sample_processed.dat").exists()
    # No manifest provenance line in log either.
    assert not any("Driven by session manifest" in line for line in lines)


# ---------------------------------------------------------------------------
# Test 2: basic manifest-driven processing — happy path
# ---------------------------------------------------------------------------

def test_basic_manifest_driven_processing(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["scan_001.xasd", "scan_002.xasd"])
    config_path = _write_config_json(tmp_path / "config.json")

    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="sampleA", label="Sample A")],
        scans=[
            ScanEntry(filename="scan_001.xasd", path="scan_001.xasd",
                      status="good", group="sampleA"),
            ScanEntry(filename="scan_002.xasd", path="scan_002.xasd",
                      status="good", group="sampleA"),
        ],
    )
    manifest_path = tmp_path / "session.json"
    save_manifest(manifest, manifest_path)

    out = tmp_path / "out"
    log, lines = _capture_log()
    process_folder(session=manifest_path, output_dir=out, log=log)

    # Output group is named from the manifest id, not from the filename.
    assert (out / "sampleA_processed.dat").exists()
    assert not (out / "scan_processed.dat").exists()
    # Provenance shows up in the log and in the text report.
    assert any("Driven by session manifest" in line for line in lines)
    report = (out / "ASTRA_processing_report.txt").read_text(encoding="utf-8")
    assert "Driven by manifest" in report
    assert str(manifest_path) in report


# ---------------------------------------------------------------------------
# Test 3: status filtering — bad/excluded scans skipped with distinct messages
# ---------------------------------------------------------------------------

def test_status_filtering_distinguishes_bad_and_excluded(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(
        data_dir,
        ["good_001.xasd", "good_002.xasd", "broken_001.xasd", "leftover_001.xasd"],
    )
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="good_002.xasd", path="good_002.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="broken_001.xasd", path="broken_001.xasd",
                      status="bad", group="grp"),
            ScanEntry(filename="leftover_001.xasd", path="leftover_001.xasd",
                      status="excluded", group="grp"),
        ],
    )
    manifest_path = tmp_path / "session.json"
    save_manifest(manifest, manifest_path)

    out = tmp_path / "out"
    log, lines = _capture_log()
    process_folder(session=manifest_path, output_dir=out, log=log)

    assert any("skipping (bad data): broken_001.xasd" in line for line in lines)
    assert any(
        "skipping (excluded from analysis): leftover_001.xasd" in line
        for line in lines
    )
    sources = _read_source_scan_order(out / "grp_processed.dat")
    assert "good_001.xasd" in sources
    assert "good_002.xasd" in sources
    assert "broken_001.xasd" not in sources
    assert "leftover_001.xasd" not in sources


# ---------------------------------------------------------------------------
# Test 4: unreviewed scans skipped with warning, processing continues
# ---------------------------------------------------------------------------

def test_unreviewed_scans_skipped_with_warning(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(
        data_dir, ["good_001.xasd", "good_002.xasd", "later_001.xasd"]
    )
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="good_002.xasd", path="good_002.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="later_001.xasd", path="later_001.xasd",
                      status="unreviewed", group="grp"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")

    out = tmp_path / "out"
    log, lines = _capture_log()
    process_folder(session=tmp_path / "session.json", output_dir=out, log=log)

    warning_lines = [line for line in lines if "unreviewed" in line.lower()]
    assert any("1 manifest scan(s)" in line for line in warning_lines)
    # Good scans still processed
    sources = _read_source_scan_order(out / "grp_processed.dat")
    assert "later_001.xasd" not in sources
    assert "good_001.xasd" in sources


# ---------------------------------------------------------------------------
# Test 5: missing file is a hard error
# ---------------------------------------------------------------------------

def test_missing_file_is_hard_error(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["good_001.xasd"])
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="ghost_001.xasd", path="ghost_001.xasd",
                      status="good", group="grp"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")
    with pytest.raises(RuntimeError, match="ghost_001.xasd"):
        process_folder(
            session=tmp_path / "session.json",
            output_dir=tmp_path / "out",
            log=lambda *_: None,
        )


# ---------------------------------------------------------------------------
# Test 6: inline config is rejected with a v2 message
# ---------------------------------------------------------------------------

def test_inline_config_rejected_with_v2_message(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["good_001.xasd"])
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="inline", inline={"alignment_source": "inline_ref"}),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")
    with pytest.raises(RuntimeError, match="reserved for v2"):
        process_folder(
            session=tmp_path / "session.json",
            output_dir=tmp_path / "out",
            log=lambda *_: None,
        )


# ---------------------------------------------------------------------------
# Test 7: bulk exclusions in config are ignored (with a warning) in manifest mode
# ---------------------------------------------------------------------------

def test_bulk_exclusions_ignored_in_manifest_mode(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["good_001.xasd", "good_002.xasd"])
    # Config tells folder mode to exclude "good_001.xasd" by exact filename.
    config_path = _write_config_json(
        tmp_path / "config.json",
        exclude_filenames=["good_001.xasd"],
    )
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp"),
            ScanEntry(filename="good_002.xasd", path="good_002.xasd",
                      status="good", group="grp"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")

    out = tmp_path / "out"
    log, lines = _capture_log()
    process_folder(session=tmp_path / "session.json", output_dir=out, log=log)

    assert any(
        "bulk filename exclusions" in line and "ignored" in line for line in lines
    )
    sources = _read_source_scan_order(out / "grp_processed.dat")
    assert "good_001.xasd" in sources
    assert "good_002.xasd" in sources


# ---------------------------------------------------------------------------
# Test 8: assigned_foil override is ignored in v1 with a warning
# ---------------------------------------------------------------------------

def test_assigned_foil_ignored_with_warning(tmp_path):
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["good_001.xasd", "good_002.xasd"])
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="grp")],
        scans=[
            ScanEntry(filename="good_001.xasd", path="good_001.xasd",
                      status="good", group="grp",
                      assigned_foil="some_other_foil.xasd"),
            ScanEntry(filename="good_002.xasd", path="good_002.xasd",
                      status="good", group="grp"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")

    log, lines = _capture_log()
    process_folder(
        session=tmp_path / "session.json",
        output_dir=tmp_path / "out",
        log=log,
    )
    assert any(
        "'assigned_foil' is reserved for v2" in line for line in lines
    )


# ---------------------------------------------------------------------------
# Test 9: manifest groups override filename-based grouping
# ---------------------------------------------------------------------------

def test_manifest_group_overrides_filename_grouping(tmp_path):
    """Two scans with different filenames stems are merged because the manifest
    puts them in the same group."""
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["alpha_001.xasd", "beta_001.xasd"])
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="combined")],
        scans=[
            ScanEntry(filename="alpha_001.xasd", path="alpha_001.xasd",
                      status="good", group="combined"),
            ScanEntry(filename="beta_001.xasd", path="beta_001.xasd",
                      status="good", group="combined"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")

    out = tmp_path / "out"
    process_folder(
        session=tmp_path / "session.json", output_dir=out, log=lambda *_: None
    )
    # One merged output file, both filenames inside it.
    assert (out / "combined_processed.dat").exists()
    assert not (out / "alpha_processed.dat").exists()
    assert not (out / "beta_processed.dat").exists()
    sources = _read_source_scan_order(out / "combined_processed.dat")
    assert set(sources) == {"alpha_001.xasd", "beta_001.xasd"}


# ---------------------------------------------------------------------------
# Test 10: both session and input_dir provided -> error
# ---------------------------------------------------------------------------

def test_session_plus_input_dir_rejected(tmp_path):
    with pytest.raises(RuntimeError, match="not both"):
        process_folder(
            input_dir=tmp_path / "data",
            session=tmp_path / "session.json",
            log=lambda *_: None,
        )


# ---------------------------------------------------------------------------
# Test 11: both session and config provided -> error
# ---------------------------------------------------------------------------

def test_session_plus_config_rejected(tmp_path):
    with pytest.raises(RuntimeError, match="not both"):
        process_folder(
            session=tmp_path / "session.json",
            config=AstraConfig(),
            log=lambda *_: None,
        )


# ---------------------------------------------------------------------------
# Test 12: explicit replicate ordering — manifest array order [c, b, a]
# ---------------------------------------------------------------------------

def test_manifest_order_overrides_alphabetical_filename_order(tmp_path):
    """Filenames a, b, c sort alphabetically; manifest declares [c, b, a].
    The merged output's source-scan list must reflect manifest order."""
    data_dir = tmp_path / "data"
    _write_synthetic_files(data_dir, ["a_001.xasd", "b_001.xasd", "c_001.xasd"])
    config_path = _write_config_json(tmp_path / "config.json")
    manifest = Manifest(
        input=InputSpec(base_dir=str(data_dir)),
        config=ConfigRef(source="path", path=str(config_path)),
        groups=[GroupEntry(id="ordered")],
        scans=[
            ScanEntry(filename="c_001.xasd", path="c_001.xasd",
                      status="good", group="ordered"),
            ScanEntry(filename="b_001.xasd", path="b_001.xasd",
                      status="good", group="ordered"),
            ScanEntry(filename="a_001.xasd", path="a_001.xasd",
                      status="good", group="ordered"),
        ],
    )
    save_manifest(manifest, tmp_path / "session.json")

    out = tmp_path / "out"
    process_folder(
        session=tmp_path / "session.json", output_dir=out, log=lambda *_: None
    )
    sources = _read_source_scan_order(out / "ordered_processed.dat")
    assert sources == ["c_001.xasd", "b_001.xasd", "a_001.xasd"]


# ---------------------------------------------------------------------------
# Unit-level test: group_samples honors _order_in_manifest
# ---------------------------------------------------------------------------

def test_group_samples_honors_order_in_manifest():
    """Direct unit test for the grouping sort key change."""
    entries = [
        {
            "filename": "c_001.xasd",
            "base_name": "grp",
            "assigned_foil": "foil",
            "replicate_id": 1,
            "_order_in_manifest": 0,
        },
        {
            "filename": "a_001.xasd",
            "base_name": "grp",
            "assigned_foil": "foil",
            "replicate_id": 1,
            "_order_in_manifest": 1,
        },
        {
            "filename": "b_001.xasd",
            "base_name": "grp",
            "assigned_foil": "foil",
            "replicate_id": 1,
            "_order_in_manifest": 2,
        },
    ]
    groups = group_samples(entries)
    ((_, _), scans), = groups.items()
    assert [s["filename"] for s in scans] == [
        "c_001.xasd", "a_001.xasd", "b_001.xasd",
    ]


def test_group_samples_folder_mode_unchanged():
    """Without _order_in_manifest, sort falls back to replicate-suffix + natural key."""
    entries = [
        {"filename": "x_003.xasd", "base_name": "x", "assigned_foil": "f", "replicate_id": 3},
        {"filename": "x_001.xasd", "base_name": "x", "assigned_foil": "f", "replicate_id": 1},
        {"filename": "x_002.xasd", "base_name": "x", "assigned_foil": "f", "replicate_id": 2},
    ]
    groups = group_samples(entries)
    ((_, _), scans), = groups.items()
    assert [s["filename"] for s in scans] == [
        "x_001.xasd", "x_002.xasd", "x_003.xasd",
    ]
