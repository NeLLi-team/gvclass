#!/usr/bin/env python
"""
Setup GVClass database by downloading required files.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.bin.gvclass_cli import CliOutput, load_config, resolve_database_path
from src.utils.database_manager import DatabaseManager


def parse_args():
    parser = argparse.ArgumentParser(description="Download and setup the GVClass database.")
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


def resolve_setup_database_path(args, repo_dir: Path) -> Path:
    config = load_config(args.config, repo_dir, CliOutput(plain_output=True))
    return resolve_database_path(args, config, repo_dir)


def main():
    """Download and setup the GVClass database."""
    print("Setting up GVClass database...")
    args = parse_args()
    repo_dir = Path(__file__).resolve().parents[2]
    target_path = resolve_setup_database_path(args, repo_dir)

    try:
        db_path = DatabaseManager.setup_database(str(target_path))
        print(f"✅ Database successfully set up at: {db_path}")
        return 0
    except Exception as exc:
        print(f"❌ Error setting up database: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
