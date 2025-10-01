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
        "models/combined.hmm",  # Combined HMM file with all models
        "gvclassSeptember25_labels.tsv",
        "order_completeness.tab",
    ]

    DATABASE_SOURCES = [
        {
            "version": "v1.1.1",
            "url": "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources_v1_1_1.tar.gz",
            "filename": "resources_v1_1_1.tar.gz",
        },
        {
            "version": "v1.1.0",
            "url": "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources_v1_1_0.tar.gz",
            "filename": "resources_v1_1_0.tar.gz",
        },
    ]
    DATABASE_VERSION = DATABASE_SOURCES[0]["version"]

    @classmethod
    def get_database_version(cls, database_path: Path) -> str:
        """
        Get the version of the installed database.

        Args:
            database_path: Path to the database directory

        Returns:
            Database version string or 'unknown' if not found
        """
        version_file = database_path / "DB_VERSION"
        if version_file.exists():
            return version_file.read_text().strip()
        return "unknown"

    @classmethod
    def write_database_version(cls, database_path: Path) -> None:
        """Write the database version to a file."""
        version_file = database_path / "DB_VERSION"
        version_file.write_text(cls.DATABASE_VERSION)

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
        # Create directory
        db_path.mkdir(parents=True, exist_ok=True)

        original_dir = Path.cwd()
        errors = []

        try:
            os.chdir(str(db_path))
            for source in cls.DATABASE_SOURCES:
                version = source["version"]
                url = source["url"]
                filename = source["filename"]

                logger.info(f"Attempting download of database {version} from {url}")
                result = subprocess.run(
                    ["wget", "-O", filename, url],
                    capture_output=True,
                    text=True,
                )

                if result.returncode != 0:
                    logger.warning("wget failed, trying curl...")
                    result = subprocess.run(
                        ["curl", "-L", "-o", filename, url],
                        capture_output=True,
                        text=True,
                    )

                if result.returncode != 0:
                    err_msg = result.stderr.strip() or result.stdout.strip() or "unknown error"
                    errors.append(f"{version}: download failed ({err_msg})")
                    Path(filename).unlink(missing_ok=True)
                    continue

                try:
                    logger.info(f"Extracting database archive {filename}")
                    subprocess.run(
                        [
                            "tar",
                            "-xzvf",
                            filename,
                            "--strip-components=1",
                            "--no-same-owner",
                        ],
                        check=True,
                    )
                except subprocess.CalledProcessError as exc:
                    errors.append(f"{version}: extraction failed ({exc})")
                    Path(filename).unlink(missing_ok=True)
                    continue

                cls.DATABASE_VERSION = version
                cls.write_database_version(db_path)
                Path(filename).unlink(missing_ok=True)
                logger.info(f"Database version {version} installed")
                return

        finally:
            os.chdir(str(original_dir))

        raise RuntimeError(
            "Failed to download database. Attempts: " + "; ".join(errors)
            if errors
            else "Failed to download database: no sources succeeded"
        )

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
