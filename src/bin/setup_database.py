#!/usr/bin/env python
"""
Setup GVClass database by downloading required files.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent))

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
        if not check_and_setup_database(target_path, config, output):
            print("❌ Database setup was not completed")
            return 1
        print(f"✅ Database successfully set up at: {target_path}")
        return 0
    except Exception as exc:
        print(f"❌ Error setting up database: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
