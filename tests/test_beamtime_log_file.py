from __future__ import annotations

import threading
import time

from astra_xas import AstraConfig
from astra_xas.beamtime.__main__ import make_tee_log
from astra_xas.beamtime._synthetic import write_synthetic_xasd
from astra_xas.beamtime.replay import replay
from astra_xas.beamtime.session import read_session_log
from astra_xas.beamtime.watcher import watch


def test_beamtime_watch_tee_log_file(tmp_path):
    source = tmp_path / "source"
    incoming = tmp_path / "incoming"
    output = tmp_path / "output"
    source.mkdir()
    write_synthetic_xasd(source / "scan_001.xasd", seed=1)
    write_synthetic_xasd(source / "scan_002.xasd", seed=2)

    scenario = tmp_path / "scenario.yaml"
    scenario.write_text(
        "\n".join(
            [
                f"source_dir: {source}",
                f"target_dir: {incoming}",
                "interval_s: 0.1",
                "jitter_s: 0.0",
                "shuffle: false",
                "inject: []",
                "",
            ]
        ),
        encoding="utf-8",
    )

    log_path = tmp_path / "beamtime.log"
    tee_fn, tee_fh = make_tee_log(log_path)
    stop_event = threading.Event()
    try:
        thread = threading.Thread(
            target=watch,
            kwargs={
                "incoming_dir": incoming,
                "output_dir": output,
                "config": AstraConfig(),
                "log": tee_fn,
                "stop_event": stop_event,
                "max_files": 2,
            },
            daemon=True,
        )
        thread.start()
        replay(scenario, log=print)

        session_log = output / "ASTRA_beamtime_session.log"
        deadline = time.monotonic() + 15.0
        content = ""
        while time.monotonic() < deadline:
            rows = read_session_log(session_log)
            has_log = log_path.exists() and log_path.stat().st_size > 0
            content = log_path.read_text(encoding="utf-8") if has_log else ""
            if (
                len(rows) == 2
                and has_log
                and "scan_001.xasd" in content
                and "scan_002.xasd" in content
                and "status=" in content
            ):
                break
            time.sleep(0.1)
        else:
            raise AssertionError("Timed out waiting for session rows and tee log status lines.")

        stop_event.set()
        thread.join(timeout=5.0)
        if thread.is_alive():
            raise AssertionError("Beamtime watcher thread did not stop.")

        assert "scan_001.xasd" in content
        assert "scan_002.xasd" in content
        assert "status=" in content
        assert any(line.startswith("[") for line in content.splitlines())
    finally:
        stop_event.set()
        if "thread" in locals():
            thread.join(timeout=5.0)
        if tee_fh is not None:
            tee_fh.close()
