from __future__ import annotations

import argparse
import builtins
import signal
import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, IO

from astra_xas._config_utils import load_config_json
from astra_xas.config import AstraConfig

from .replay import replay
from .watcher import watch


def make_tee_log(
    log_file_path: Path | None,
) -> tuple[Callable, IO | None]:
    """Build a log function that mirrors stdout to a file.

    Returns (log_function, file_handle_or_None). The caller is
    responsible for closing the file handle (if not None) in a finally
    block.
    """
    if log_file_path is None:
        return builtins.print, None
    try:
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        file_handle = log_file_path.open("a", encoding="utf-8", buffering=1)
    except Exception as exc:
        builtins.print(
            f"WARNING: could not open beamtime log file {log_file_path}: "
            f"{type(exc).__name__}: {exc}"
        )
        return builtins.print, None

    file_lock = threading.Lock()

    def tee_log(*args, **kwargs):
        builtins.print(*args, **kwargs)
        message = " ".join(str(arg) for arg in args)
        ts = datetime.now().isoformat(timespec="seconds")
        line = f"[{ts}] {message}\n"
        with file_lock:
            try:
                file_handle.write(line)
                file_handle.flush()
            except Exception:
                pass

    return tee_log, file_handle


def main() -> None:
    parser = argparse.ArgumentParser(description="AstraXAS Beamtime Mode utilities.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    replay_parser = subparsers.add_parser("replay", help="Replay a scenario into a watch folder.")
    replay_parser.add_argument("scenario", help="Path to replay scenario YAML.")
    replay_parser.add_argument(
        "-l",
        "--log-file",
        type=Path,
        default=None,
        help="Mirror stdout to this file (append mode, line-buffered).",
    )

    watch_parser = subparsers.add_parser("watch", help="Watch a folder for new .xasd scans.")
    watch_parser.add_argument("incoming_dir", help="Folder receiving .xasd scans.")
    watch_parser.add_argument("-o", "--output-dir", default=None, help="Output folder. Default: <incoming>-beamtime")
    watch_parser.add_argument("-c", "--config", default=None, help="Optional AstraConfig JSON file.")
    watch_parser.add_argument(
        "-l",
        "--log-file",
        type=Path,
        default=None,
        help="Mirror stdout to this file (append mode, line-buffered).",
    )

    args = parser.parse_args()
    if args.command == "replay":
        log_path = args.log_file.expanduser().resolve() if args.log_file is not None else None
        log_fn, log_fh = make_tee_log(log_path)
        try:
            replay(Path(args.scenario), log=log_fn)
        finally:
            if log_fh is not None:
                log_fh.close()
        return

    if args.command == "watch":
        config = load_config_json(Path(args.config)) if args.config else AstraConfig()
        stop_event = threading.Event()
        previous_handler = signal.getsignal(signal.SIGINT)
        log_path = args.log_file.expanduser().resolve() if args.log_file is not None else None
        log_fn, log_fh = make_tee_log(log_path)
        if log_fh is not None:
            log_fn(f"AstraXAS Beamtime watcher starting; tee log -> {log_path}")

        def _handle_sigint(signum, frame):
            stop_event.set()

        signal.signal(signal.SIGINT, _handle_sigint)
        try:
            watch(
                incoming_dir=Path(args.incoming_dir),
                output_dir=Path(args.output_dir) if args.output_dir else None,
                config=config,
                log=log_fn,
                stop_event=stop_event,
            )
        except KeyboardInterrupt:
            stop_event.set()
        finally:
            if log_fh is not None:
                log_fh.close()
            signal.signal(signal.SIGINT, previous_handler)


if __name__ == "__main__":
    main()
