from __future__ import annotations

import json
import os
import shutil
import threading
import time
from pathlib import Path

import numpy as np

from astra_xas import AstraConfig
from astra_xas.beamtime import watcher as watcher_module
from astra_xas.beamtime._synthetic import write_synthetic_xasd
from astra_xas.beamtime.groups import restore_group_registry, update_group_with_entry
from astra_xas.beamtime.replay import replay
from astra_xas.beamtime.session import read_session_log
from astra_xas.beamtime.watcher import watch
from astra_xas.io import load_xasd
from astra_xas.single_scan import SingleScanResult, _entry_from_scan


def _write_scenario(path: Path, source: Path, incoming: Path, interval_s: float = 0.1) -> None:
    path.write_text(
        "\n".join(
            [
                f"source_dir: {source}",
                f"target_dir: {incoming}",
                f"interval_s: {interval_s}",
                "jitter_s: 0.0",
                "shuffle: false",
                "inject: []",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _start_watcher(incoming: Path, output: Path, max_files=None):
    stop_event = threading.Event()
    thread = threading.Thread(
        target=watch,
        kwargs={
            "incoming_dir": incoming,
            "output_dir": output,
            "config": AstraConfig(),
            "stop_event": stop_event,
            "max_files": max_files,
        },
        daemon=True,
    )
    thread.start()
    return stop_event, thread


def _stop_watcher(stop_event: threading.Event, thread: threading.Thread) -> None:
    stop_event.set()
    thread.join(timeout=5.0)
    if thread.is_alive():
        raise AssertionError("Beamtime watcher thread did not stop.")


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _first_data_row(path: Path) -> list[float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        return [float(value) for value in line.split()]
    raise AssertionError(f"No data row found in {path}")


def _atomic_copy(src: Path, dst: Path) -> None:
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    shutil.copy2(src, tmp)
    tmp.replace(dst)


def test_beamtime_live_group_merge_two_groups(tmp_path):
    source = tmp_path / "source"
    incoming = tmp_path / "incoming"
    output = tmp_path / "output"
    source.mkdir()
    for seed, name in enumerate(
        ["Cu_sample_1.xasd", "Cu_sample_2.xasd", "Fe_sample_1.xasd", "Fe_sample_2.xasd"],
        start=1,
    ):
        write_synthetic_xasd(source / name, seed=seed)
    scenario = tmp_path / "scenario.yaml"
    _write_scenario(scenario, source, incoming)

    stop_event, thread = _start_watcher(incoming, output, max_files=4)
    replay(scenario, log=lambda *_: None)

    deadline = time.monotonic() + 20.0
    missing = ""
    while time.monotonic() < deadline:
        rows = read_session_log(output / "ASTRA_beamtime_session.log")
        cu_norm = output / "groups" / "Cu_sample_norm.dat"
        fe_norm = output / "groups" / "Fe_sample_norm.dat"
        cu_summary = output / "groups" / "Cu_sample_group_summary.json"
        fe_summary = output / "groups" / "Fe_sample_group_summary.json"
        cu_plot = output / "plots" / "group_qc" / "Cu_sample_replicate_qc.png"
        fe_plot = output / "plots" / "group_qc" / "Fe_sample_replicate_qc.png"
        dashboard = output / "index.html"
        html = dashboard.read_text(encoding="utf-8") if dashboard.exists() else ""
        ready = (
            len(rows) == 4
            and cu_norm.exists()
            and fe_norm.exists()
            and cu_summary.exists()
            and fe_summary.exists()
            and cu_plot.exists()
            and fe_plot.exists()
            and "Live groups" in html
            and "Cu_sample" in html
            and "Fe_sample" in html
        )
        if ready:
            break
        missing = (
            f"rows={len(rows)}/4, cu_norm={cu_norm.exists()}, fe_norm={fe_norm.exists()}, "
            f"cu_summary={cu_summary.exists()}, fe_summary={fe_summary.exists()}, "
            f"cu_plot={cu_plot.exists()}, fe_plot={fe_plot.exists()}, dashboard={'Live groups' in html}"
        )
        time.sleep(0.1)
    else:
        raise AssertionError(f"Timed out waiting for live group outputs: {missing}")
    _stop_watcher(stop_event, thread)

    summaries = sorted((output / "groups").glob("*_group_summary.json"))
    assert len(summaries) == 2
    for summary_path in summaries:
        summary = _read_json(summary_path)
        assert summary["last_merge_status"] == "ready"
        assert summary["n_accepted"] == 2
        row = _first_data_row(output / summary["output_files"]["norm"])
        assert len(row) == 2
        assert all(np.isfinite(row))


def test_beamtime_rejects_do_not_enter_group_state_and_duplicate_is_noop(tmp_path):
    source = tmp_path / "source"
    incoming = tmp_path / "incoming"
    output = tmp_path / "output"
    source.mkdir()
    write_synthetic_xasd(source / "Sn_sample_1.xasd", seed=1)
    (source / "Sn_sample_2.xasd").write_bytes(b"# not real data\n")
    write_synthetic_xasd(source / "Sn_sample_3.xasd", seed=3)
    scenario = tmp_path / "scenario.yaml"
    _write_scenario(scenario, source, incoming)

    stop_event, thread = _start_watcher(incoming, output, max_files=3)
    replay(scenario, log=lambda *_: None)

    summary_path = output / "groups" / "Sn_sample_group_summary.json"
    deadline = time.monotonic() + 20.0
    while time.monotonic() < deadline:
        rows = read_session_log(output / "ASTRA_beamtime_session.log")
        if len(rows) == 3 and summary_path.exists():
            summary = _read_json(summary_path)
            if summary.get("n_accepted") == 2:
                break
        time.sleep(0.1)
    else:
        raise AssertionError("Timed out waiting for Sn_sample group summary with 2 accepted scans.")
    _stop_watcher(stop_event, thread)

    summary = _read_json(summary_path)
    assert summary["accepted_filenames"] == ["Sn_sample_1.xasd", "Sn_sample_3.xasd"]
    assert summary["n_accepted"] == 2
    assert summary["last_merge_status"] == "ready"
    assert (output / "groups" / "Sn_sample_norm.dat").exists()

    before_mtime_ns = os.stat(summary_path).st_mtime_ns
    before_bytes = summary_path.read_bytes()
    registry = {}
    registry_lock = threading.Lock()
    restore_group_registry(output, registry, registry_lock, log=lambda *_: None)
    duplicate_path = incoming / "Sn_sample_1.xasd"
    entry = _entry_from_scan(load_xasd(duplicate_path), AstraConfig(), path=duplicate_path)
    update_group_with_entry(
        entry,
        duplicate_path,
        "ok",
        output,
        AstraConfig(),
        registry,
        registry_lock,
        log=lambda *_: None,
    )
    assert os.stat(summary_path).st_mtime_ns == before_mtime_ns
    assert summary_path.read_bytes() == before_bytes


def test_beamtime_group_state_restores_and_updates_after_restart(tmp_path):
    source = tmp_path / "source"
    incoming = tmp_path / "incoming"
    output = tmp_path / "output"
    source.mkdir()
    write_synthetic_xasd(source / "Cu_sample_1.xasd", seed=1)
    write_synthetic_xasd(source / "Cu_sample_2.xasd", seed=2)
    write_synthetic_xasd(source / "Cu_sample_3.xasd", seed=3)
    first_source = tmp_path / "first_source"
    first_source.mkdir()
    shutil.copy2(source / "Cu_sample_1.xasd", first_source / "Cu_sample_1.xasd")
    shutil.copy2(source / "Cu_sample_2.xasd", first_source / "Cu_sample_2.xasd")
    scenario = tmp_path / "scenario.yaml"
    _write_scenario(scenario, first_source, incoming)

    stop_a, thread_a = _start_watcher(incoming, output, max_files=2)
    replay(scenario, log=lambda *_: None)

    summary_path = output / "groups" / "Cu_sample_group_summary.json"
    deadline = time.monotonic() + 20.0
    while time.monotonic() < deadline:
        if summary_path.exists():
            summary = _read_json(summary_path)
            outputs = summary.get("output_files", {})
            if (
                summary.get("n_accepted") == 2
                and summary.get("last_merge_status") == "ready"
                and (output / outputs.get("norm", "")).exists()
                and (output / outputs.get("qc_plot", "")).exists()
            ):
                break
        time.sleep(0.1)
    else:
        raise AssertionError("Timed out waiting for first live group merge.")
    _stop_watcher(stop_a, thread_a)

    stop_b, thread_b = _start_watcher(incoming, output, max_files=None)
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline:
        dashboard = output / "index.html"
        html = dashboard.read_text(encoding="utf-8") if dashboard.exists() else ""
        if "Live groups" in html and "Cu_sample" in html:
            break
        time.sleep(0.1)
    else:
        raise AssertionError("Timed out waiting for restored group dashboard.")

    _atomic_copy(source / "Cu_sample_3.xasd", incoming / "Cu_sample_3.xasd")
    deadline = time.monotonic() + 15.0
    while time.monotonic() < deadline:
        summary = _read_json(summary_path)
        if summary.get("n_accepted") == 3 and summary.get("last_merge_status") == "ready":
            break
        time.sleep(0.1)
    else:
        raise AssertionError("Timed out waiting for restored group to accept third replicate.")
    _stop_watcher(stop_b, thread_b)

    summary = _read_json(summary_path)
    assert summary["n_accepted"] == 3
    assert summary["last_merge_status"] == "ready"
    assert summary["accepted_filenames"] == [
        "Cu_sample_1.xasd",
        "Cu_sample_2.xasd",
        "Cu_sample_3.xasd",
    ]
    assert (output / "groups" / "Cu_sample_processed.dat").exists()
    assert (output / "groups" / "Cu_sample_norm.dat").exists()
    assert (output / "groups" / "Cu_sample_flat.dat").exists()
    assert (output / "plots" / "group_qc" / "Cu_sample_replicate_qc.png").exists()


def test_watcher_status_is_ok_when_jumps_present_but_no_warnings(monkeypatch, tmp_path):
    """Phase 2.2 C behavior change: detector jumps do NOT turn status to warn.

    Even if process_single_scan reports n_detector_jumps > 0, as long as
    qc_warnings and qc_errors are empty, the watcher logs status='ok'.
    Jumps are still counted in the n_jumps column.
    """
    incoming = tmp_path / "incoming"
    incoming.mkdir()
    xasd_path = incoming / "clean_with_jumps.xasd"
    write_synthetic_xasd(xasd_path, seed=0)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    config = AstraConfig()

    scan = load_xasd(xasd_path)
    real_entry = _entry_from_scan(scan, config, path=xasd_path)

    fake_result = SingleScanResult(
        filename=real_entry["filename"],
        is_foil=bool(real_entry.get("is_foil", False)),
        entry=real_entry,
        energy=real_entry["energy"],
        analysis_signal=real_entry["mu_fluo"],
        analysis_signal_label="IF/I0",
        qc_status="ok",
        qc_warnings=[],
        qc_errors=[],
        detector_jumps=[
            {"channel": "I0", "energy_eV": 7012.0, "severity": "low",
             "jump_size": 0.05, "relative_jump": 0.01,
             "filename": real_entry["filename"],
             "inside_plot_window": True, "inside_alignment_window": False,
             "inside_preedge_window": False, "inside_norm_window": False, "note": ""},
            {"channel": "I0", "energy_eV": 7016.0, "severity": "low",
             "jump_size": 0.04, "relative_jump": 0.008,
             "filename": real_entry["filename"],
             "inside_plot_window": True, "inside_alignment_window": False,
             "inside_preedge_window": False, "inside_norm_window": False, "note": ""},
            {"channel": "I0", "energy_eV": 7018.0, "severity": "low",
             "jump_size": 0.03, "relative_jump": 0.006,
             "filename": real_entry["filename"],
             "inside_plot_window": True, "inside_alignment_window": False,
             "inside_preedge_window": False, "inside_norm_window": False, "note": ""},
        ],
        metrics={
            "n_points": int(real_entry["energy"].size),
            "n_points_finite_energy": int(real_entry["energy"].size),
            "energy_min": float(real_entry["energy"].min()),
            "energy_max": float(real_entry["energy"].max()),
            "energy_range_eV": [float(real_entry["energy"].min()),
                                float(real_entry["energy"].max())],
            "channels_present": ["I0", "IF"],
            "analysis_signal_finite_fraction": 1.0,
            "n_detector_jumps": 3,
            "n_validation_warnings": 0,
            "n_validation_errors": 0,
        },
    )

    monkeypatch.setattr(
        watcher_module,
        "process_single_scan",
        lambda scan, config: fake_result,
    )

    session_log = output_dir / "ASTRA_beamtime_session.log"
    registry: dict = {}
    registry_lock = threading.Lock()

    watcher_module._process_scan(
        path=xasd_path,
        config=config,
        session_log=session_log,
        output_dir=output_dir,
        registry=registry,
        registry_lock=registry_lock,
        log=lambda *args, **kwargs: None,
    )

    rows = read_session_log(session_log)
    assert len(rows) >= 1, "Expected at least one session row to be written"
    last_row = rows[-1]

    assert last_row["status"] == "ok", (
        f"Expected status='ok' (jumps should not trigger warn in new policy), "
        f"got {last_row['status']}"
    )
    assert int(last_row["n_jumps"]) == 3, (
        f"Expected n_jumps=3 (jumps still counted), got {last_row['n_jumps']}"
    )
    assert int(last_row["n_warnings"]) == 0, (
        f"Expected n_warnings=0, got {last_row['n_warnings']}"
    )
