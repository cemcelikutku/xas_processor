from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

from ._config_utils import load_config_json
from .config import AstraConfig
from .processor import process_folder


def _resolve_config(args: argparse.Namespace) -> AstraConfig:
    """Build AstraConfig from --config and explicit CLI flags.

    Precedence (lowest to highest):
        AstraConfig defaults < --config JSON values < explicit CLI flags
    """
    if args.config is not None:
        config = load_config_json(args.config)
    else:
        config = AstraConfig()

    if args.mode is not None:
        config.analysis_mode = args.mode
    if args.foil_mode is not None:
        config.foil_alignment_mode = args.foil_mode
    if args.e0 is not None:
        config.e0 = args.e0

    return config


def _warn_session_overrides(args: argparse.Namespace) -> None:
    """Emit per-flag warnings when CLI overrides are passed alongside --session.

    In session mode the manifest is the source of truth for config; CLI flags
    are accepted by argparse but ignored. We warn rather than reject so that
    shell history / scripted invocations still work.
    """
    suffix = "manifest config is the source of truth in session mode"
    if args.config is not None:
        print(f"WARNING: ignoring --config {args.config} ({suffix})", file=sys.stderr)
    if args.mode is not None:
        print(f"WARNING: ignoring --mode {args.mode} ({suffix})", file=sys.stderr)
    if args.foil_mode is not None:
        print(f"WARNING: ignoring --foil-mode {args.foil_mode} ({suffix})", file=sys.stderr)
    if args.e0 is not None:
        print(f"WARNING: ignoring --e0 {args.e0} ({suffix})", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Process ASTRA fluorescence XAS .xasd files with foil drift correction.")
    parser.add_argument(
        "input_dir",
        nargs="?",
        default=None,
        help="Folder containing .xasd files. Omit when --session is provided.",
    )
    parser.add_argument("-o", "--output-dir", default=None, help="Output folder. Default: <input>-processed")
    parser.add_argument(
        "-s",
        "--session",
        type=Path,
        default=None,
        help=(
            "Path to a session manifest JSON. When provided, drives processing "
            "from the manifest: input_dir, --config, and per-field overrides "
            "(--mode, --foil-mode, --e0) are ignored with a warning."
        ),
    )
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        default=None,
        help=(
            "Path to a JSON config file. Values from the file override "
            "AstraConfig defaults; explicit CLI flags (--mode, --foil-mode, "
            "--e0) further override the config file values."
        ),
    )
    parser.add_argument(
        "--mode",
        default=None,
        choices=["fluo", "trans", "ref"],
        help="Sample analysis mode. Default: fluo (or value from --config).",
    )
    parser.add_argument(
        "--foil-mode",
        default=None,
        choices=["trans", "ref", "fluo"],
        help="Signal used for foil alignment. Default: trans (or value from --config).",
    )
    parser.add_argument(
        "--e0",
        type=float,
        default=None,
        help="Edge energy in eV. Default: 7121.030 (or value from --config).",
    )
    args = parser.parse_args()

    if args.session is not None:
        if args.input_dir is not None:
            print(
                "ERROR: positional input_dir is not allowed when --session is provided.",
                file=sys.stderr,
            )
            sys.exit(2)
        _warn_session_overrides(args)
        process_folder(session=args.session, output_dir=args.output_dir)
        return

    if args.input_dir is None:
        parser.error("input_dir is required when --session is not provided.")

    try:
        config = _resolve_config(args)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"ERROR: could not load config file {args.config}: {exc}", file=sys.stderr)
        sys.exit(1)
    process_folder(args.input_dir, args.output_dir, config=config)


if __name__ == "__main__":
    main()
