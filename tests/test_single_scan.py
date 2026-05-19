"""Phase 2.1 architectural and behavioral tests for astra_xas.single_scan."""
from __future__ import annotations

import ast
import inspect
import numpy as np
import pytest

# Tests can import private helpers to validate behavior.
import astra_xas.single_scan
import astra_xas.processor
from astra_xas.config import AstraConfig
from astra_xas.io import load_xasd
from astra_xas.beamtime._synthetic import write_synthetic_xasd
from astra_xas.single_scan import (
    SingleScanResult,
    process_single_scan,
    _validate_single_scan,
    _entry_from_scan,
)
from astra_xas.processor import _validate_processing_inputs


def test_single_scan_does_not_import_processor():
    """single_scan.py must not depend on processor.py (absolute OR relative)."""
    source = inspect.getsource(astra_xas.single_scan)
    tree = ast.parse(source)

    forbidden = {"processor", "astra_xas.processor"}

    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            assert module not in forbidden, (
                f"single_scan.py must not import from processor.py "
                f"(found module={module!r} at line {node.lineno})"
            )
            assert not module.startswith("astra_xas.processor."), (
                f"single_scan.py must not import from astra_xas.processor.* "
                f"(found module={module!r} at line {node.lineno})"
            )

        if isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.name
                assert name not in forbidden, (
                    f"single_scan.py must not import processor "
                    f"(found import {name!r} at line {node.lineno})"
                )
                assert not name.startswith("astra_xas.processor."), (
                    f"single_scan.py must not import astra_xas.processor.* "
                    f"(found import {name!r} at line {node.lineno})"
                )


def test_no_duplicate_function_definitions_in_processor():
    """Moved functions must not have definitions left in processor.py."""
    source = inspect.getsource(astra_xas.processor)
    moved_names = [
        "_entry_from_scan",
        "_analysis_signal_spec",
        "_alignment_signal_spec",
        "_required_channels_for_signal",
        "_range_overlap_status",
        "_channel_validation_messages",
        "_alignment_structure_warning",
        "detect_detector_jumps",
        "_config_bool",
        "_config_float",
    ]
    for name in moved_names:
        assert f"\ndef {name}(" not in source, (
            f"Function '{name}' moved to single_scan.py but a definition "
            f"still exists in processor.py (would shadow the re-import)"
        )


def test_no_duplicate_constant_definitions_in_processor():
    """Moved detector-jump constants must not be redefined in processor.py."""
    source = inspect.getsource(astra_xas.processor)
    moved_constants = [
        "RAW_DETECTOR_JUMP_CHANNELS",
        "PRIMARY_DETECTOR_JUMP_CHANNELS",
        "FDT_DETECTOR_JUMP_CHANNELS",
        "DERIVED_DETECTOR_JUMP_CHANNELS",
    ]
    for name in moved_constants:
        assert f"\n{name} =" not in source, (
            f"Constant '{name}' moved to single_scan.py but is still "
            f"defined in processor.py (would shadow the re-import)"
        )


def test_process_single_scan_smoke(tmp_path):
    """Direct call to process_single_scan() returns a properly-shaped result."""
    xasd_path = tmp_path / "smoke_001.xasd"
    write_synthetic_xasd(xasd_path, seed=0)
    scan = load_xasd(xasd_path)
    config = AstraConfig()

    result = process_single_scan(scan, config)

    # Type and shape assertions
    assert isinstance(result, SingleScanResult)
    assert isinstance(result.filename, str)
    assert isinstance(result.is_foil, bool)
    assert isinstance(result.entry, dict)
    assert isinstance(result.energy, np.ndarray)
    assert isinstance(result.analysis_signal, np.ndarray)
    assert isinstance(result.analysis_signal_label, str)
    assert result.qc_status in ("ok", "warn", "reject")
    assert isinstance(result.qc_warnings, list)
    assert isinstance(result.qc_errors, list)
    assert isinstance(result.detector_jumps, list)
    assert isinstance(result.metrics, dict)

    # Metrics shape — exactly the 10 required keys
    required_metrics = {
        "n_points",
        "n_points_finite_energy",
        "energy_min",
        "energy_max",
        "energy_range_eV",
        "channels_present",
        "analysis_signal_finite_fraction",
        "n_detector_jumps",
        "n_validation_warnings",
        "n_validation_errors",
    }
    assert set(result.metrics.keys()) == required_metrics

    # Spot-check some types
    assert isinstance(result.metrics["n_points"], int)
    assert isinstance(result.metrics["channels_present"], list)
    assert isinstance(result.metrics["analysis_signal_finite_fraction"], float)


def test_validate_single_scan_matches_batch_for_single_entry(tmp_path):
    """For a single entry, _validate_single_scan and _validate_processing_inputs
    should produce equivalent results (same set of warnings and errors)."""
    xasd_path = tmp_path / "validate_001.xasd"
    write_synthetic_xasd(xasd_path, seed=0)
    scan = load_xasd(xasd_path)
    config = AstraConfig()

    entry = _entry_from_scan(scan, config)

    single_warnings, single_errors = _validate_single_scan(entry, config)
    batch_warnings, batch_errors = _validate_processing_inputs([entry], config)

    # Set comparison: same warnings and errors, possibly different order
    assert set(single_warnings) == set(batch_warnings), (
        f"Warnings differ:\n"
        f"single: {sorted(single_warnings)}\n"
        f"batch:  {sorted(batch_warnings)}"
    )
    assert set(single_errors) == set(batch_errors), (
        f"Errors differ:\n"
        f"single: {sorted(single_errors)}\n"
        f"batch:  {sorted(batch_errors)}"
    )
