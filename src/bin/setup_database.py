#!/usr/bin/env python
"""
Setup GVClass database by downloading required files.
"""

import argparse
import sys
from pathlib import Path

# Add repository root to path to import src modules
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.bin.gvclass_cli import (
    CliOutput,
    check_and_setup_database,
    load_config,
    resolve_database_path,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download and setup the GVClass database."
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config/gvclass_config.yaml",
        help="Configuration file (default: config/gvclass_config.yaml)",
    )
    parser.add_argument(
        "-d",
        "--database",
        help="Database path (overrides GVCLASS_DB and config)",
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="Replace a complete database with the pinned version "
        "(without it, a complete database is never replaced)",
    )
    parser.add_argument(
        "--with-screen",
        action="store_true",
        help="Also install the contig-screen add-on (gvclass screen --install)",
    )
    return parser.parse_args()


def main():
    """Download and setup the GVClass database."""
    print("Setting up GVClass database...")
    args = parse_args()
    repo_dir = Path(__file__).resolve().parents[2]
    output = CliOutput()
    config = load_config(args.config, repo_dir, output)
    target_path = resolve_database_path(args, config, repo_dir)

    try:
        if not check_and_setup_database(
            target_path, config, output, update=args.update
        ):
            print("❌ Database setup was not completed")
            return 1
        print(f"✅ Database successfully set up at: {target_path}")
        if args.with_screen:
            # Imported here: the screen code loads heavy libraries.
            from src.bin.gvclass_screen import install_bundle  # noqa: PLC0415

            return install_bundle(target_path)
        return 0
    except Exception as exc:
        print(f"❌ Error setting up database: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
