from __future__ import annotations

import hashlib
import json
import os
import queue
import threading
import time
from datetime import datetime
from pathlib import Path

from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

from astra_xas.config import AstraConfig
from astra_xas.io import load_xasd, natural_key
from astra_xas.single_scan import process_single_scan

from .dashboard import render_dashboard
from .groups import restore_group_registry, update_group_with_entry
from .plots import render_per_scan_plot
from .session import append_session_row, write_session_ended_marker


def _timestamp() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_checkpoint(path: Path) -> dict:
    if not path.exists():
        return {"processed": {}}
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return {"processed": {}}
    if not isinstance(data, dict) or not isinstance(data.get("processed"), dict):
        return {"processed": {}}
    return {"processed": dict(data["processed"])}


def _write_checkpoint(path: Path, checkpoint: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(checkpoint, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def _wait_size_stable(path: Path, log=print, timeout_s: float = 10.0) -> bool:
    deadline = time.monotonic() + timeout_s
    last_size = -1
    stable_count = 0
    while time.monotonic() < deadline:
        try:
            size = path.stat().st_size
        except OSError:
            size = -1
        if size > 0 and size == last_size:
            stable_count += 1
            if stable_count >= 2:
                return True
        else:
            stable_count = 0
        last_size = size
        time.sleep(0.2)
    log(f"WARNING: file did not settle within 10 s; skipped {path.name}")
    return False


def _process_scan(
    path: Path,
    config: AstraConfig,
    session_log: Path,
    output_dir: Path,
    registry,
    registry_lock,
    log=print,
) -> None:
    timestamp = _timestamp()
    filename = path.name

    try:
        scan = load_xasd(path)
    except Exception as exc:
        status = "reject"
        notes = f"load_failed: {type(exc).__name__}: {exc}"
        append_session_row(
            session_log,
            {
                "timestamp_iso": timestamp,
                "filename": filename,
                "status": status,
                "n_warnings": 0,
                "n_jumps": 0,
                "notes": notes,
            },
        )
        log(f"{timestamp} {filename} status={status} warns=0 jumps=0")
        render_dashboard(output_dir, log=log)
        return

    try:
        result = process_single_scan(scan, config)
    except Exception as exc:
        status = "reject"
        notes = f"pipeline_failed: {type(exc).__name__}: {exc}"
        append_session_row(
            session_log,
            {
                "timestamp_iso": timestamp,
                "filename": filename,
                "status": status,
                "n_warnings": 0,
                "n_jumps": 0,
                "notes": notes,
            },
        )
        log(f"{timestamp} {filename} status={status} warns=0 jumps=0")
        render_dashboard(output_dir, log=log)
        return

    status = result.qc_status
    warnings = result.qc_warnings
    fatal_errors = result.qc_errors
    n_jumps = int(result.metrics["n_detector_jumps"])
    n_warnings = len(warnings)
    entry = result.entry

    notes_parts: list[str] = []
    if fatal_errors:
        notes_parts.append("fatal: " + " | ".join(fatal_errors))
    if warnings:
        notes_parts.append("warnings: " + " | ".join(warnings))
    notes = "; ".join(notes_parts)

    append_session_row(
        session_log,
        {
            "timestamp_iso": timestamp,
            "filename": filename,
            "status": status,
            "n_warnings": n_warnings,
            "n_jumps": n_jumps,
            "notes": notes,
        },
    )

    plot_path = output_dir / "plots" / "beamtime" / f"{Path(filename).stem}.png"
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    render_per_scan_plot(
        entry,
        config,
        status,
        n_warnings,
        n_jumps,
        timestamp,
        plot_path,
        log=log,
    )
    if status in {"ok", "warn"}:
        update_group_with_entry(
            entry,
            path,
            status,
            output_dir,
            config,
            registry,
            registry_lock,
            log=log,
        )
    render_dashboard(output_dir, log=log)
    log(f"{timestamp} {filename} status={status} warns={n_warnings} jumps={n_jumps}")


def watch(
    incoming_dir: Path,
    output_dir: Path | None = None,
    config: AstraConfig | None = None,
    log=print,
    stop_event: threading.Event | None = None,
    max_files: int | None = None,
) -> None:
    incoming_dir = Path(incoming_dir).expanduser().resolve()
    output_dir = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else incoming_dir.parent / f"{incoming_dir.name}-beamtime"
    )
    config = config or AstraConfig()
    stop_event = stop_event or threading.Event()

    output_dir.mkdir(parents=True, exist_ok=True)
    session_dir = output_dir / "_astra_session"
    session_dir.mkdir(parents=True, exist_ok=True)
    beamtime_plot_dir = output_dir / "plots" / "beamtime"
    beamtime_plot_dir.mkdir(parents=True, exist_ok=True)
    groups_dir = output_dir / "groups"
    groups_dir.mkdir(parents=True, exist_ok=True)
    group_qc_dir = output_dir / "plots" / "group_qc"
    group_qc_dir.mkdir(parents=True, exist_ok=True)
    session_log = output_dir / "ASTRA_beamtime_session.log"
    checkpoint_path = session_dir / "checkpoint.json"
    checkpoint = _load_checkpoint(checkpoint_path)

    file_queue: queue.Queue[Path] = queue.Queue()
    queued_or_processing: set[Path] = set()
    set_lock = threading.Lock()
    count_lock = threading.Lock()
    registry: dict = {}
    registry_lock = threading.Lock()
    rows_written = 0
    restore_group_registry(output_dir, registry, registry_lock, log=log)

    def should_skip_checkpoint(path: Path) -> bool:
        try:
            digest = _sha256_file(path)
        except OSError:
            return False
        return checkpoint["processed"].get(path.name) == digest

    def enqueue(path: Path) -> None:
        path = Path(path)
        if path.suffix.lower() != ".xasd" or path.name.endswith(".tmp"):
            return
        abs_path = path.expanduser().resolve()
        if should_skip_checkpoint(abs_path):
            return
        with set_lock:
            if abs_path in queued_or_processing:
                return
            queued_or_processing.add(abs_path)
        file_queue.put(abs_path)
        if file_queue.qsize() > 100:
            log(f"WARNING: beamtime watch queue has {file_queue.qsize()} pending files")

    def max_reached() -> bool:
        with count_lock:
            return max_files is not None and rows_written >= max_files

    def worker() -> None:
        nonlocal rows_written
        while True:
            if stop_event.is_set() and file_queue.empty():
                return
            if max_reached():
                stop_event.set()
                return
            try:
                path = file_queue.get(timeout=0.5)
            except queue.Empty:
                continue
            try:
                if not should_skip_checkpoint(path) and _wait_size_stable(path, log=log):
                    _process_scan(path, config, session_log, output_dir, registry, registry_lock, log=log)
                    try:
                        checkpoint["processed"][path.name] = _sha256_file(path)
                        _write_checkpoint(checkpoint_path, checkpoint)
                    except OSError as exc:
                        log(f"WARNING: could not update checkpoint for {path.name}: {exc}")
                    with count_lock:
                        rows_written += 1
                        if max_files is not None and rows_written >= max_files:
                            stop_event.set()
            finally:
                with set_lock:
                    queued_or_processing.discard(path)
                file_queue.task_done()

    class Handler(FileSystemEventHandler):
        def on_created(self, event):
            if not event.is_directory:
                enqueue(Path(event.src_path))

        def on_moved(self, event):
            if not event.is_directory:
                enqueue(Path(event.dest_path))

    worker_thread = threading.Thread(target=worker, daemon=True)
    worker_thread.start()

    incoming_dir.mkdir(parents=True, exist_ok=True)
    observer = Observer()
    observer.schedule(Handler(), str(incoming_dir), recursive=False)
    observer.start()

    for path in sorted(incoming_dir.glob("*.xasd"), key=lambda p: natural_key(p.name)):
        enqueue(path)

    while not file_queue.empty() and not stop_event.is_set():
        time.sleep(0.1)
    render_dashboard(output_dir, log=log)

    try:
        while not stop_event.is_set():
            if max_reached():
                stop_event.set()
                break
            time.sleep(0.2)
    except KeyboardInterrupt:
        stop_event.set()
    finally:
        observer.stop()
        observer.join(timeout=2.0)
        deadline = time.monotonic() + 5.0
        while not file_queue.empty() and time.monotonic() < deadline:
            time.sleep(0.1)
        stop_event.set()
        worker_thread.join(timeout=max(0.1, deadline - time.monotonic()))
        write_session_ended_marker(session_log)
