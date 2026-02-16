"""
Database management utilities for GVClass.
Handles downloading and validating the reference database.
"""

import hashlib
import logging
import shutil
import subprocess
import tarfile
from pathlib import Path, PurePosixPath
from typing import Optional

logger = logging.getLogger(__name__)


class DatabaseManager:
    """Manages GVClass reference database."""

    REQUIRED_FILES = [
        "models/combined.hmm",  # Combined HMM file with all models
        "gvclassFeb26_labels.tsv",
        "order_completeness.tab",
    ]

    DATABASE_SOURCES = [
        {
            "version": "v1.2.1",
            "url": "https://zenodo.org/records/18662446/files/resources_v1_2_1.tar.gz?download=1",
            "filename": "resources_v1_2_1.tar.gz",
            "sha256": "94ad59100a158ca7ac760e9a1d35f01ff5e2eb07c270c4ca1b45fbe9e3e637eb",
        },
        {
            "version": "v1.2.0",
            "url": "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources_v1_2_0.tar.gz",
            "filename": "resources_v1_2_0.tar.gz",
            "sha256": None,
        },
        {
            "version": "v1.1.1",
            "url": "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources_v1_1_1.tar.gz",
            "filename": "resources_v1_1_1.tar.gz",
            "sha256": None,
        },
        {
            "version": "v1.1.0",
            "url": "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources_v1_1_0.tar.gz",
            "filename": "resources_v1_1_0.tar.gz",
            "sha256": None,
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
        db_path.mkdir(parents=True, exist_ok=True)
        errors = []

        for source in cls.DATABASE_SOURCES:
            version = source["version"]
            archive_path = db_path / source["filename"]

            try:
                cls._install_database_source(db_path, source, archive_path)
                logger.info(f"Database version {version} installed")
                return
            except Exception as exc:
                errors.append(f"{version}: {exc}")
                logger.warning(f"Database source {version} failed: {exc}")
            finally:
                archive_path.unlink(missing_ok=True)

        raise RuntimeError(
            "Failed to download database. Attempts: " + "; ".join(errors)
            if errors
            else "Failed to download database: no sources succeeded"
        )

    @classmethod
    def _install_database_source(cls, db_path: Path, source: dict, archive_path: Path):
        """Download, verify, and extract a single database source."""
        version = source["version"]
        url = source["url"]

        logger.info(f"Attempting download of database {version} from {url}")
        cls._download_archive(url, archive_path)
        cls._verify_archive_checksum(source, archive_path)
        cls._extract_archive(archive_path, db_path)

        cls.DATABASE_VERSION = version
        cls.write_database_version(db_path)

    @classmethod
    def _download_archive(cls, url: str, archive_path: Path) -> None:
        """Download database archive with wget/curl fallback."""
        commands = [
            ("wget", ["wget", "-O", str(archive_path), url]),
            ("curl", ["curl", "-L", "-o", str(archive_path), url]),
        ]
        errors = []

        for idx, (tool, command) in enumerate(commands):
            if idx > 0:
                logger.warning(f"{commands[idx - 1][0]} failed, trying {tool}...")
            result = cls._run_download_command(command, tool)
            if result.returncode == 0 and archive_path.exists() and archive_path.stat().st_size > 0:
                return

            archive_path.unlink(missing_ok=True)
            error_text = cls._format_command_error(result)
            if result.returncode == 0:
                errors.append(f"{tool} produced an empty archive")
            else:
                errors.append(f"{tool} failed ({error_text})")

        raise RuntimeError("download failed (" + "; ".join(errors) + ")")

    @staticmethod
    def _run_download_command(command: list, tool_name: str) -> subprocess.CompletedProcess:
        """Run a download command and normalize missing-binary errors."""
        try:
            return subprocess.run(command, capture_output=True, text=True)
        except FileNotFoundError:
            return subprocess.CompletedProcess(
                args=command,
                returncode=127,
                stdout="",
                stderr=f"{tool_name} not found",
            )

    @staticmethod
    def _format_command_error(result: subprocess.CompletedProcess) -> str:
        """Extract a readable subprocess error string."""
        return result.stderr.strip() or result.stdout.strip() or "unknown error"

    @classmethod
    def _verify_archive_checksum(cls, source: dict, archive_path: Path) -> None:
        """Verify archive SHA256 when a checksum is available for the source."""
        expected_sha256 = source.get("sha256")
        if not expected_sha256:
            logger.warning(
                "No SHA256 configured for database source %s; skipping checksum verification",
                source["version"],
            )
            return

        observed_sha256 = cls._compute_sha256(archive_path)
        if observed_sha256.lower() != expected_sha256.lower():
            raise RuntimeError(
                "checksum mismatch "
                f"(expected {expected_sha256.lower()}, got {observed_sha256.lower()})"
            )
        logger.info(f"SHA256 verified for source {source['version']}")

    @staticmethod
    def _compute_sha256(file_path: Path) -> str:
        """Compute SHA256 checksum for a file."""
        digest = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @classmethod
    def _extract_archive(cls, archive_path: Path, db_path: Path) -> None:
        """Extract archive contents with path traversal protection."""
        logger.info(f"Extracting database archive {archive_path.name}")
        try:
            with tarfile.open(archive_path, "r:gz") as tar:
                cls._extract_tar_members(tar, db_path, strip_components=1)
        except (tarfile.TarError, OSError, RuntimeError) as exc:
            raise RuntimeError(f"extraction failed ({exc})") from exc

    @classmethod
    def _extract_tar_members(
        cls, tar: tarfile.TarFile, db_path: Path, strip_components: int
    ) -> None:
        """Extract a tar archive to db_path while preventing path traversal."""
        base_path = db_path.resolve()
        for member in tar.getmembers():
            relative_path = cls._normalize_member_path(member.name, strip_components)
            if relative_path is None:
                continue
            if member.issym() or member.islnk():
                raise RuntimeError(f"archive member uses links: {member.name}")

            target_path = cls._resolve_safe_target(base_path, relative_path, member.name)
            if member.isdir():
                target_path.mkdir(parents=True, exist_ok=True)
                continue
            if not member.isfile():
                continue

            target_path.parent.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"failed to read archive member: {member.name}")
            with extracted, open(target_path, "wb") as out_f:
                shutil.copyfileobj(extracted, out_f)

    @staticmethod
    def _normalize_member_path(member_name: str, strip_components: int) -> Optional[Path]:
        """Normalize member path and apply strip-components behavior."""
        raw_parts = [part for part in PurePosixPath(member_name).parts if part not in ("", ".")]
        if len(raw_parts) <= strip_components:
            return None

        stripped_parts = raw_parts[strip_components:]
        if any(part == ".." for part in stripped_parts):
            raise RuntimeError(f"path traversal component detected: {member_name}")
        return Path(*stripped_parts)

    @staticmethod
    def _resolve_safe_target(base_path: Path, relative_path: Path, member_name: str) -> Path:
        """Resolve extraction target and ensure it stays under base_path."""
        target_path = (base_path / relative_path).resolve()
        try:
            target_path.relative_to(base_path)
        except ValueError as exc:
            raise RuntimeError(f"path traversal detected in archive member: {member_name}") from exc
        return target_path

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
