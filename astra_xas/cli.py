from __future__ import annotations
import argparse
from .config import AstraConfig
from .processor import process_folder


def main():
    parser = argparse.ArgumentParser(description="Process ASTRA fluorescence XAS .xasd files with foil drift correction.")
    parser.add_argument("input_dir", help="Folder containing .xasd files")
    parser.add_argument("-o", "--output-dir", default=None, help="Output folder. Default: <input>-processed")
    parser.add_argument("--mode", default="fluo", choices=["fluo", "trans", "ref"], help="Sample analysis mode")
    parser.add_argument("--foil-mode", default="trans", choices=["trans", "ref", "fluo"], help="Signal used for foil alignment")
    parser.add_argument("--e0", type=float, default=7121.030)
    args = parser.parse_args()
    config = AstraConfig(analysis_mode=args.mode, foil_alignment_mode=args.foil_mode, e0=args.e0)
    process_folder(args.input_dir, args.output_dir, config=config)


if __name__ == "__main__":
    main()
