"""
Database management utilities for GVClass.
Handles downloading and validating the reference database.
"""

import os
import subprocess
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class DatabaseManager:
    """Manages GVClass reference database."""

    REQUIRED_FILES = [
        "models/combined.hmm",
        "ncldvApril24_labels.txt",
        "models_APRIL24--databaseApril24.cutoffs",
        "order_completeness.tab",
    ]

    DATABASE_URL = "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources.tar.gz"

    @classmethod
    def setup_database(cls, database_path: Optional[str] = None) -> Path:
        """
        Setup the GVClass database.

        Args:
            database_path: Custom database location. If None, uses 'resources' in current directory.

        Returns:
            Path to the database directory
        """
        # Determine database path
        if database_path:
            db_path = Path(database_path).resolve()
        else:
            db_path = Path.cwd() / "resources"

        logger.info(f"Setting up database at: {db_path}")

        # Check if database exists and is complete
        if db_path.exists():
            missing_files = cls._check_missing_files(db_path)
            if not missing_files:
                logger.info("Database already exists and is complete")
                return db_path
            else:
                logger.warning(f"Database exists but missing files: {missing_files}")

        # Download and extract database
        cls._download_database(db_path)

        # Verify all files are present
        missing_files = cls._check_missing_files(db_path)
        if missing_files:
            raise ValueError(
                "Database download failed. Missing files:\n" + "\n".join(missing_files)
            )

        logger.info("Database setup complete")
        return db_path

    @classmethod
    def _check_missing_files(cls, db_path: Path) -> list:
        """Check for missing required files."""
        missing = []
        for file_path in cls.REQUIRED_FILES:
            if not (db_path / file_path).exists():
                missing.append(file_path)
        return missing

    @classmethod
    def _download_database(cls, db_path: Path):
        """Download and extract the database."""
        logger.info(f"Downloading database from {cls.DATABASE_URL}")

        # Create directory
        db_path.mkdir(parents=True, exist_ok=True)

        # Save current directory
        original_dir = Path.cwd()

        try:
            # Change to database directory
            os.chdir(str(db_path))

            # Download using wget (could also use requests)
            logger.info("Downloading resources.tar.gz...")
            result = subprocess.run(
                ["wget", "-O", "resources.tar.gz", cls.DATABASE_URL],
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                # Try with curl if wget fails
                logger.warning("wget failed, trying curl...")
                result = subprocess.run(
                    ["curl", "-L", "-o", "resources.tar.gz", cls.DATABASE_URL],
                    capture_output=True,
                    text=True,
                )

            if result.returncode != 0:
                raise RuntimeError(f"Failed to download database: {result.stderr}")

            # Extract
            logger.info("Extracting database...")
            subprocess.run(
                ["tar", "-xzvf", "resources.tar.gz", "--strip-components=1"], check=True
            )

            # Clean up
            Path("resources.tar.gz").unlink()

        finally:
            # Return to original directory
            os.chdir(str(original_dir))

    @classmethod
    def validate_database(cls, db_path: Path) -> bool:
        """
        Validate that all required database files exist.

        Args:
            db_path: Path to database directory

        Returns:
            True if valid, False otherwise
        """
        missing_files = cls._check_missing_files(db_path)
        if missing_files:
            logger.error(f"Missing database files: {missing_files}")
            return False
        return True
