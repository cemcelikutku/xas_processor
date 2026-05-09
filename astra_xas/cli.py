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


def main():
    parser = argparse.ArgumentParser(description="Process ASTRA fluorescence XAS .xasd files with foil drift correction.")
    parser.add_argument("input_dir", help="Folder containing .xasd files")
    parser.add_argument("-o", "--output-dir", default=None, help="Output folder. Default: <input>-processed")
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
    try:
        config = _resolve_config(args)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"ERROR: could not load config file {args.config}: {exc}", file=sys.stderr)
        sys.exit(1)
    process_folder(args.input_dir, args.output_dir, config=config)


if __name__ == "__main__":
    main()
