#!/usr/bin/env python
"""
Setup GVClass database by downloading required files.
"""

import sys
from pathlib import Path

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils.database_manager import DatabaseManager


def main():
    """Download and setup the GVClass database."""
    print("Setting up GVClass database...")

    try:
        db_path = DatabaseManager.setup_database()
        print(f"✅ Database successfully set up at: {db_path}")
        return 0
    except Exception as e:
        print(f"❌ Error setting up database: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
