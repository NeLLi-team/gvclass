"""
Database management utilities for GVClass.
Handles downloading and validating the reference database.
"""

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path, PurePosixPath
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)


class DatabaseManager:
    """Manages GVClass reference database."""

    REMOTE_LOOKUP_TIMEOUT_SECONDS = 5
    ZENODO_RECORD_RE = re.compile(
        r"https?://(?:www\.)?zenodo\.org/(?:api/)?records/(?P<record_id>\d+)"
    )

    REQUIRED_FILES = [
        "models/combined.hmm",  # Combined HMM file with all models
        "gvclassFeb26_labels.tsv",
        "order_completeness.tab",
        "novelty_completeness_config.json",
        "novelty_strategy2_tiers.tsv",
        "novelty_strategy2_baselines.tsv",
        "novelty_strategy3_model_metadata.tsv",
        "novelty_strategy3_models.joblib",
        "contamination_model.joblib",
    ]

    DEFAULT_DATABASE_SOURCE = {
        "version": "v1.4.0",
        "url": "https://zenodo.org/records/18926264/files/resources_v1_4_0.tar.gz?download=1",
        "filename": "resources_v1_4_0.tar.gz",
        "sha256": "95e2d75b5229a33f1910e16849fc067f3a2f55db5d12fbc25fb593aa61d9f3da",
    }
    LEGACY_DATABASE_SOURCES = [
        {
            "version": "v1.2.2",
            "url": "https://zenodo.org/records/18675742/files/resources_v1_2_2.tar.gz?download=1",
            "filename": "resources_v1_2_2.tar.gz",
            "sha256": "47bd2822b1da8a4c30db259607a0e7ae1747a9df6ce554a98cc015fd7b654ace",
        },
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
    DATABASE_SOURCES = [DEFAULT_DATABASE_SOURCE, *LEGACY_DATABASE_SOURCES]
    DATABASE_VERSION = DEFAULT_DATABASE_SOURCE["version"]

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
            version = version_file.read_text().strip()
            if version:
                return cls.normalize_version_label(version)
        return "unknown"

    @classmethod
    def write_database_version(
        cls, database_path: Path, version: Optional[str] = None
    ) -> None:
        """Write the database version to a file."""
        version_file = database_path / "DB_VERSION"
        version_file.write_text(
            cls.normalize_version_label(version or cls.DATABASE_VERSION)
        )

    @staticmethod
    def normalize_version_label(version: Optional[str]) -> str:
        """Normalize version labels to the repo's preferred v-prefixed format."""
        normalized = str(version or "").strip()
        if not normalized:
            return "unknown"
        if normalized.lower() in {"unknown", "not installed"}:
            return normalized.lower()
        if normalized.startswith("v"):
            return normalized
        if normalized[0].isdigit():
            return f"v{normalized}"
        return normalized

    @classmethod
    def is_newer_version(cls, candidate_version: str, current_version: str) -> bool:
        """Compare two database version labels, falling back gracefully for unknowns."""
        candidate_parts = cls._parse_version(candidate_version)
        current_parts = cls._parse_version(current_version)

        if candidate_parts is not None and current_parts is not None:
            return candidate_parts > current_parts
        if cls.normalize_version_label(current_version) in {"unknown", "not installed"}:
            return cls.normalize_version_label(candidate_version) not in {
                "unknown",
                "not installed",
            }
        return cls.normalize_version_label(
            candidate_version
        ) != cls.normalize_version_label(current_version)

    @staticmethod
    def _parse_version(version: Optional[str]) -> Optional[tuple[int, ...]]:
        """Parse dotted numeric versions such as v1.4.0."""
        normalized = str(version or "").strip()
        match = re.search(r"(\d+(?:\.\d+)*)", normalized)
        if not match:
            return None
        return tuple(int(part) for part in match.group(1).split("."))

    @classmethod
    def get_effective_database_source(
        cls, preferred_source: Optional[dict] = None
    ) -> dict:
        """Return the highest-priority download source after env/config resolution."""
        sources = cls._resolve_database_sources(preferred_source)
        if not sources:
            raise ValueError("No database download source is configured")
        return sources[0]

    @classmethod
    def get_latest_database_source(
        cls, preferred_source: Optional[dict] = None
    ) -> Optional[dict]:
        """Resolve the latest Zenodo-backed database source from the effective source."""
        try:
            effective_source = cls.get_effective_database_source(preferred_source)
        except ValueError:
            return None

        record_id = cls._extract_zenodo_record_id(effective_source["url"])
        if not record_id:
            return None

        try:
            record = cls._fetch_json(f"https://zenodo.org/api/records/{record_id}")
            latest_url = record.get("links", {}).get("latest")
            latest_record = cls._fetch_json(latest_url) if latest_url else record
            return cls._source_from_zenodo_record(latest_record)
        except (
            HTTPError,
            URLError,
            TimeoutError,
            ValueError,
            json.JSONDecodeError,
        ) as exc:
            logger.warning(
                "Could not resolve latest database version from Zenodo: %s", exc
            )
            return None

    @classmethod
    def setup_database(
        cls,
        database_path: Optional[str] = None,
        preferred_source: Optional[dict] = None,
        force: bool = False,
    ) -> Path:
        """
        Setup the GVClass database.

        Args:
            database_path: Custom database location. If None, uses 'resources' in current directory.
            preferred_source: Optional download source override from config/env.
            force: Reinstall even when the current database appears complete.

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
        if db_path.exists() and not force:
            missing_files = cls._check_missing_files(db_path)
            if not missing_files:
                logger.info("Database already exists and is complete")
                return db_path
            else:
                logger.warning(f"Database exists but missing files: {missing_files}")
        elif db_path.exists() and force:
            logger.info("Forcing database reinstall at: %s", db_path)

        # Download and extract database
        cls._download_database(db_path, preferred_source=preferred_source)

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
    def _download_database(cls, db_path: Path, preferred_source: Optional[dict] = None):
        """Download and extract the database."""
        db_path.mkdir(parents=True, exist_ok=True)
        errors = []

        for source in cls._resolve_database_sources(preferred_source):
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
    def _resolve_database_sources(
        cls, preferred_source: Optional[dict] = None
    ) -> list[dict]:
        """Resolve download sources in priority order with config/env overrides first."""
        sources = []
        for candidate in (
            cls._source_from_env(),
            preferred_source,
            *cls.DATABASE_SOURCES,
        ):
            normalized = cls._normalize_source(candidate)
            if normalized and all(
                existing["url"] != normalized["url"] for existing in sources
            ):
                sources.append(normalized)
        return sources

    @classmethod
    def _source_from_env(cls) -> Optional[dict]:
        """Read an optional database source override from environment variables."""
        url = os.environ.get("GVCLASS_DB_URL", "").strip()
        if not url:
            return None

        return {
            "version": os.environ.get(
                "GVCLASS_DB_VERSION",
                cls.DEFAULT_DATABASE_SOURCE["version"],
            ).strip(),
            "url": url,
            "filename": os.environ.get("GVCLASS_DB_FILENAME", "").strip() or None,
            "sha256": os.environ.get("GVCLASS_DB_SHA256", "").strip() or None,
        }

    @classmethod
    def _normalize_source(cls, source: Optional[dict]) -> Optional[dict]:
        """Fill in missing download-source metadata and discard empty entries."""
        if not source:
            return None

        url = str(source.get("url", "")).strip()
        if not url:
            return None

        filename = str(source.get("filename", "")).strip()
        if not filename:
            filename = PurePosixPath(urlparse(url).path).name
        if not filename:
            raise ValueError(f"Could not derive archive filename from URL: {url}")

        version = cls.normalize_version_label(
            str(source.get("version", "")).strip()
            or cls.DEFAULT_DATABASE_SOURCE["version"]
        )
        sha256 = source.get("sha256")
        if sha256 is not None:
            sha256 = str(sha256).strip() or None

        return {
            "version": version,
            "url": url,
            "filename": filename,
            "sha256": sha256,
        }

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
        cls.write_database_version(db_path, version=version)

    @classmethod
    def _extract_zenodo_record_id(cls, url: str) -> Optional[str]:
        """Extract a Zenodo record identifier from a standard record/files URL."""
        match = cls.ZENODO_RECORD_RE.search(url or "")
        if not match:
            return None
        return match.group("record_id")

    @staticmethod
    def _fetch_json(url: str) -> dict:
        """Fetch a JSON document from a remote API."""
        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "gvclass-database-manager",
            },
        )
        with urlopen(
            request, timeout=DatabaseManager.REMOTE_LOOKUP_TIMEOUT_SECONDS
        ) as response:
            return json.loads(response.read().decode("utf-8"))

    @staticmethod
    def _fetch_text(url: str) -> str:
        """Fetch a plain-text document from a remote URL."""
        request = Request(url, headers={"User-Agent": "gvclass-database-manager"})
        with urlopen(
            request, timeout=DatabaseManager.REMOTE_LOOKUP_TIMEOUT_SECONDS
        ) as response:
            return response.read().decode("utf-8").strip()

    @classmethod
    def _source_from_zenodo_record(cls, record: dict) -> dict:
        """Build a download source descriptor from a Zenodo record payload."""
        files = record.get("files", [])
        archive_file = next(
            (
                candidate
                for candidate in files
                if candidate.get("key", "").endswith(".tar.gz")
                and not candidate.get("key", "").endswith(".tar.gz.sha256")
            ),
            None,
        )
        if archive_file is None:
            raise ValueError("Zenodo record does not include a database archive")

        archive_name = archive_file.get("key", "")
        archive_url = archive_file.get("links", {}).get("self")
        if not archive_url:
            raise ValueError(
                f"Zenodo record archive is missing a download link: {archive_name}"
            )

        sha256_file = next(
            (
                candidate
                for candidate in files
                if candidate.get("key") == f"{archive_name}.sha256"
            ),
            None,
        )
        sha256 = None
        if sha256_file:
            sha256_url = sha256_file.get("links", {}).get("self")
            if sha256_url:
                try:
                    sha256_text = cls._fetch_text(sha256_url)
                except (HTTPError, URLError, TimeoutError, ValueError) as exc:
                    logger.warning(
                        "Could not retrieve SHA256 sidecar for Zenodo record %s: %s",
                        record.get("id") or record.get("recid") or "unknown",
                        exc,
                    )
                else:
                    sha256 = sha256_text.split()[0] if sha256_text else None

        return {
            "version": cls.normalize_version_label(
                record.get("metadata", {}).get("version")
            ),
            "url": archive_url,
            "filename": archive_name,
            "sha256": sha256,
            "record_id": str(record.get("id") or record.get("recid") or "").strip()
            or None,
            "conceptrecid": str(record.get("conceptrecid") or "").strip() or None,
        }

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
            if (
                result.returncode == 0
                and archive_path.exists()
                and archive_path.stat().st_size > 0
            ):
                return

            archive_path.unlink(missing_ok=True)
            error_text = cls._format_command_error(result)
            if result.returncode == 0:
                errors.append(f"{tool} produced an empty archive")
            else:
                errors.append(f"{tool} failed ({error_text})")

        raise RuntimeError("download failed (" + "; ".join(errors) + ")")

    @staticmethod
    def _run_download_command(
        command: list, tool_name: str
    ) -> subprocess.CompletedProcess:
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

            target_path = cls._resolve_safe_target(
                base_path, relative_path, member.name
            )
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
    def _normalize_member_path(
        member_name: str, strip_components: int
    ) -> Optional[Path]:
        """Normalize member path and apply strip-components behavior."""
        raw_parts = [
            part for part in PurePosixPath(member_name).parts if part not in ("", ".")
        ]
        if len(raw_parts) <= strip_components:
            return None

        stripped_parts = raw_parts[strip_components:]
        if any(part == ".." for part in stripped_parts):
            raise RuntimeError(f"path traversal component detected: {member_name}")
        return Path(*stripped_parts)

    @staticmethod
    def _resolve_safe_target(
        base_path: Path, relative_path: Path, member_name: str
    ) -> Path:
        """Resolve extraction target and ensure it stays under base_path."""
        target_path = (base_path / relative_path).resolve()
        try:
            target_path.relative_to(base_path)
        except ValueError as exc:
            raise RuntimeError(
                f"path traversal detected in archive member: {member_name}"
            ) from exc
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
