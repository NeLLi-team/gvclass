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
import sys
import tarfile
import time
from pathlib import Path, PurePosixPath
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

_TRANSFER_CHUNK_SIZE = 1 << 20


class _ProgressReporter:
    """Render lightweight progress updates to stderr."""

    BAR_WIDTH = 24
    MIN_TTY_UPDATE_SECONDS = 0.2
    MIN_PERCENT_STEP = 5.0
    MIN_BYTE_STEP = 64 * 1024 * 1024

    def __init__(
        self,
        label: str,
        total_bytes: Optional[int] = None,
        stream=None,
    ) -> None:
        self.label = label
        self.total_bytes = total_bytes if total_bytes and total_bytes > 0 else None
        self.stream = stream if stream is not None else sys.stderr
        self.current_bytes = 0
        self._last_percent = -1.0
        self._last_bytes = -1
        self._last_update = 0.0
        self._is_tty = bool(getattr(self.stream, "isatty", lambda: False)())

    def start(self) -> None:
        self._emit(force=True)

    def advance(self, size: int) -> None:
        self.current_bytes += max(size, 0)
        self._emit()

    def finish(self) -> None:
        if self.total_bytes is not None:
            self.current_bytes = max(self.current_bytes, self.total_bytes)
        self._emit(force=True)
        if self._is_tty:
            self.stream.write("\n")
            self.stream.flush()

    def _emit(self, force: bool = False) -> None:
        if not force and not self._should_emit():
            return

        message = self._format_message()
        if self._is_tty:
            self.stream.write(f"\r{message:<80}")
        else:
            self.stream.write(f"{message}\n")
        self.stream.flush()

        self._last_update = time.monotonic()
        self._last_bytes = self.current_bytes
        self._last_percent = self._percent_complete()

    def _should_emit(self) -> bool:
        if self._is_tty:
            return time.monotonic() - self._last_update >= self.MIN_TTY_UPDATE_SECONDS

        if self.total_bytes is not None:
            return (
                self._percent_complete() - self._last_percent >= self.MIN_PERCENT_STEP
            )

        return (self.current_bytes - self._last_bytes) >= self.MIN_BYTE_STEP

    def _percent_complete(self) -> float:
        if self.total_bytes is None:
            return 0.0
        return min(100.0, (self.current_bytes / self.total_bytes) * 100.0)

    def _format_message(self) -> str:
        if self.total_bytes is None:
            return f"{self.label}: {self._format_bytes(self.current_bytes)}"

        percent = self._percent_complete()
        transferred = self._format_bytes(self.current_bytes)
        total = self._format_bytes(self.total_bytes)
        if self._is_tty:
            filled = int((percent / 100.0) * self.BAR_WIDTH)
            bar = "#" * filled + "-" * (self.BAR_WIDTH - filled)
            return f"{self.label}: [{bar}] {percent:5.1f}% " f"({transferred}/{total})"
        return f"{self.label}: {percent:5.1f}% ({transferred}/{total})"

    @staticmethod
    def _format_bytes(size: int) -> str:
        units = ("B", "KiB", "MiB", "GiB", "TiB")
        value = float(size)
        unit_index = 0
        while value >= 1024.0 and unit_index < len(units) - 1:
            value /= 1024.0
            unit_index += 1
        if unit_index == 0:
            return f"{int(value)} {units[unit_index]}"
        return f"{value:.1f} {units[unit_index]}"


_BASE_PART_SHA256 = (
    "77ca92456cda7fa7305505da0c821481e61fdbad49514a21804a39407c23a274",
    "0772c0ebe7994dd2e71bcda4f65312a988d0936c01137921bb3d1e2d5b45753b",
    "867f81b2830efe5c0f0a2a241e0ff8970183577e80167b83ecf2718079557c6b",
    "cea7e47941e9fdcbc365359299f935930b45b3be65aefd28825231710ae4450a",
    "956168746bb5e7271031fb2ce772d2a357168f8ddb5eca5ad7da366745edbe7e",
    "491ae08cbdda457e86a979629177891e37d080effb9b049457a5248441fe92d5",
    "ea582c563190b5c78defdeb05f52ea287725536172040980a3abb1e6938676ac",
    "e563d0d82be563b9837e62ca5c5b9fe9d37e7e97b2bc07d9c4fdd5127cfd933d",
    "eb026bfb98fc6598da427f3d73e9e8e5d01c30c17c1ca6c59ccd23f37dffec5a",
    "6bb6c7152822278e1021c7b424ac13e933b0601e6ba5bf72b043e8f6e727f0da",
    "726af269668b3a7ea574dfc3c27160d119c0fc8625fd5f820614b1b7431518ce",
    "b0ca7600c0f4be8a62668adf11c305daaf4d252c570ed789438e4e57b6e1b8d2",
    "962274cfe37e442b43679f05bc45a19bf43602b066449a251cb26e7ede8cb81f",
    "1e8b5b9aa0d6727bd5e1390e8da43e6612b59c7c064688ab798cf63cb4c414d0",
    "7c6fa8a96723784f9712debb479fac0c32ba2a1b0296febda25ac9709c42d0e3",
    "f08d7d95ac9625900ecc5141f294f7e3ef3ef25ef677e1c50e6bcb1000d511fe",
    "42f7fdf4e0f2d24dd1374937d603f899bbbb82604f3b46e24ce323844602e6ca",
    "863f2525f18d40d45cd6af0ab6ec4c84cffa4aa0637cf39b16d7516957528d4e",
    "84ecef7f7823ed4113eb8a7be7504b1dc7d6f68e245a196e59fea409c9af296e",
    "b43d3dee217b34dcb5fcb6d0161604d7f88a2e9dd748d8035b173d605b3e4915",
    "e743bee11b8cf0c5d31eaa447f598375323c279a459f21fa2f7a3d682e3e49f4",
    "afd188e40c29f124c883abeda8049a4f09d83431a300ae85f2881df9c6fa79d8",
    "59898a67cf4a9e4c8b2e85549a2dbdf734749810a53eb56009c11a2e5746a94c",
    "f03f507f7f924420d156fcfa3aa72aebef914fc22220fa7743940af656132fa9",
)
_SCREEN_PART_SHA256 = (
    "6a80f773759b1c4221c3ed81c8d14c5506066441ed5885cb1eab3eecb602ce33",
    "c2f02bf0cf04ce401021deeef579d278975b61fc4bc9e5efc8c445a85c291898",
    "098dfba722a63262f217e2ae3284757f16216d7ad1ccb9c8921a6b17f5409934",
    "44145426fe9b82f53b612e6c98cfb7e2870846f5c2b7c02c11ee02050e1645ed",
    "6dda3572cafddc6c7ecab56497243ce398705679477a8eaf06588dcf11f28b38",
    "0c31d659ecbd4080270829435fbe5879df652b27e253a84ac662dcb1552941b7",
    "51a771b7d1ec32a8cffc3e6045bf8352a4f9a3855b3c7ebc8adeb364dd722a3e",
    "f247e969f493eae29d031adb358ad323dd3771d2330e71027489cfaf0916356c",
    "935b9c4e7bc9a9e3b0e06f4dbcd72016645e4406ea6b1e7cb621aeee9e0441d7",
    "a231c9bb74bba1354c30aa0fb90f388da3af337b101ea068e7558808d3a49475",
    "f31cdeb38ac8901b418e01f267d5a4595b832ca4625e44c01f5b1baa080766d5",
    "1938fd4c78c91443dc777f552f15b6ef593833a9a76658bd60479b23659e0ca4",
    "fcb944d9933b0ed46969f7b0fbd828f6662e0c53c48006be64090043699b92d2",
    "0e574a674b1dad842d4e615770300af38de7be4056a19252f89ff3eda912db32",
    "c61b7c530a4131f2dae93b0caa68b8d5366c4954ae5421286bbe2050c40c5e99",
    "553abba8e31c1d7dfa9a20deb487a49aaf3223ca11891129aed4d7b41abe6186",
    "9c9a60a326206e031e060c0f21bead79e16c0a52e7124387209752d6ee121b42",
    "6bd6caa35956694ae6e28f3fdbc40eea0bb18e60a2b66aaeab5d6674007099ac",
    "784649f2b5365424d7cabf0790cc5620464c767b5be935fb1f38a89f6ce705b4",
    "1391b2d987d900b8168e03f322fc808d70d366c301423653e4c8c23790593f86",
    "671c14c61cb03264555027f1cdaa0120cbc79304c4efe7ad85eb8127308f4c51",
    "01d33902ee87cf7b697f5d25f0b4dda40c93549864b305deb2b6a1cbc8d00489",
    "b86f0a92239ff2ade2d2c8cf833af7c8c088cce14e2ddfe7d19d20b80c4371ff",
    "572b92496e9acb1e22fbb2c7205d9144eb2dde70b002f4a751f7b888c37a1a41",
    "3ceffbdd67d6bb43cb8a92b743e77d78c75fe366fdceef33834934f307cadf56",
    "1495766c4ae6bfcd78512a782c24f31edee050d63fce5ae25e373b739a484aa4",
    "f6015890a63f55250b55e3131c5366723bc41e5b12f529433051c87e98c1b805",
    "11ed0d9424b02e76abf8902834ddff035b73cf179b0aed539b5179e6ba328cf2",
    "3a3e84091e6557502315781401b216b9959519214f477028e35923d1289ca7c2",
    "1e791c11bf66fd9e4fe8393365d1f9ef6544886fc6b5ba56ca8955d7121aed59",
    "6d0dd79358a17ed3f0568652ce66e4eaf33f042696041533c1f50209599dd053",
    "fc71bca7fc61e90462813cc7a1517f442af75b33c7fe98204bdc2b1fc9896c8b",
    "1905df2f971940340958efdc198a7f20802d02d69e48e97dfd308069dac07fc9",
    "7d2b8d8dd60b611871de7a1f40c0e58235d86d50b16e44afbaf6755342c2caa9",
    "e838ba474ac152e87be721518d8d4d7233a438eb37fe2e98af3054b2e40ce681",
    "2d06015f62b80073989e8b6e190485a96c356dc5d1c2489b2c18d41c55b4d056",
    "d5f12b4152d5475a3eb0e71b4158aba107e467251e8ee17331f423694ac7c614",
    "a4e68c9c9d44d904e7b41973a27a13806af58b76621d473406f345a1124b5161",
    "dc8381e84edcf295441ab3cd9c49c3646c213d7445ecaee613ad1c8a1f79c760",
    "53388fed6a3780d45fe6b3e72029159f64ae0b4a287ccf5e1a83a0f6ebd31439",
    "c0af5ef2d94748bad41d9fd42d4b06a137c5f70d6302bf1cf7e81f0a4c4d8403",
    "9f0e5f40cbcc029b52a11ce0c2910ce9857be8daf40f75cffcb6cc33670f2b76",
    "66256abb736a8a2ecb54e40c3218e211a9e58363a5fc26b0b7219f6bf31bc64a",
    "9082eafe02c0d2e65ccdd8b89fc9c72a7e695134b0ea98f51b8c4b8e5c491b64",
    "3cba2c6a4263394dbb891744f4b70d8aec590c5b0c08d5def7800711014eec87",
    "afa2871ca2313c653bcf67ab695d54a81797a515ec6892cfa731138c440580f2",
    "f2a255f309add762ac8d9a947039326cacf8f9b6cc825797715e082d1b870486",
    "f623825ff69666c43f7587fd557018e5ab98d1159030842e7ecb56eb6f2dee39",
)


def split_source(
    *, record: int, filename: str, sha256: str, part_sha256: tuple[str, ...]
) -> dict:
    """Return the download source of a Zenodo archive that is stored in parts."""
    base = f"https://zenodo.org/records/{record}/files/{filename}"
    parts = [
        {"url": f"{base}.part{number:03d}?download=1", "sha256": digest}
        for number, digest in enumerate(part_sha256)
    ]
    return {
        "version": "v2.1.0",
        "url": parts[0]["url"],
        "filename": filename,
        "sha256": sha256,
        "parts": parts,
    }


class DatabaseManager:
    """Manages GVClass reference database."""

    REMOTE_LOOKUP_TIMEOUT_SECONDS = 5
    ZENODO_RECORD_RE = re.compile(
        r"https?://(?:www\.)?zenodo\.org/(?:api/)?records/(?P<record_id>\d+)"
    )

    REQUIRED_FILES = (
        "hmm/combined.hmm",
        "labels.tsv",
        "parquet/faa.parquet",
        "markers/order_completeness.tab",
        "markers/stats.tsv",
        "markers/annotations.tsv",
        "completeness/config.json",
        "completeness/tiers.tsv",
        "completeness/baselines.tsv",
        "completeness/model_metadata.tsv",
        "completeness/model.joblib",
        "contamination/model.joblib",
        "contamination/model.yaml",
    )
    REQUIRED_FILE_ALTERNATIVES = {
        "labels.tsv": ("labels.tsv", "parquet/labels/labels.parquet"),
    }

    #: The archive is stored as parts (the Zenodo server cuts a long request); each
    #: part and the joined archive are checked against these hashes.
    DEFAULT_DATABASE_SOURCE = split_source(
        record=23381313,
        filename="resources_v2_1_0.tar.gz",
        sha256="823023a0d9af82d6b01d7beeb1a12af422432c7c1e3a11b57f04130bfa7f5ee6",
        part_sha256=_BASE_PART_SHA256,
    )
    DATABASE_VERSION = DEFAULT_DATABASE_SOURCE["version"]

    #: Source of the optional contig-screen resources (``gvclass screen --install``);
    #: pinned like the database source.
    CONTIG_SCREEN_SOURCE: dict | None = split_source(
        record=23464780,
        filename="contig_screen_v2_1_0.tar.gz",
        sha256="49183cfcd818940e17b4690f13b70caede7bab79a027f91d18183a699fe17135",
        part_sha256=_SCREEN_PART_SHA256,
    )

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
            if file_path == "parquet/faa.parquet":
                parquet_faa = db_path / file_path
                legacy_faa = db_path / "database" / "faa"
                has_parquet_faa = (
                    parquet_faa.is_file() and parquet_faa.stat().st_size > 0
                )
                has_legacy_faa = any(
                    marker_file.is_file() and marker_file.stat().st_size > 0
                    for marker_file in legacy_faa.glob("*.faa")
                )
                if not has_parquet_faa and not has_legacy_faa:
                    missing.append(file_path)
                continue
            alternatives = cls.REQUIRED_FILE_ALTERNATIVES.get(file_path)
            if alternatives:
                if not any(
                    (db_path / candidate).exists() for candidate in alternatives
                ):
                    missing.append(file_path)
                continue
            if not (db_path / file_path).exists():
                missing.append(file_path)
        return missing

    @classmethod
    def _download_database(cls, db_path: Path, preferred_source: Optional[dict] = None):
        """Download and extract the database via an atomic staging swap.

        Extraction writes into ``<db_path>.new`` and only swaps into place
        after the archive is verified, fully extracted, and the DB_VERSION
        marker is written. A mid-extraction crash leaves the prior
        ``db_path`` (if any) untouched and the half-built staging directory
        is removed on the next attempt. The downloaded archive is parked in
        ``db_path.parent`` as a sibling so it never contaminates the
        staging tree.
        """
        db_path.parent.mkdir(parents=True, exist_ok=True)
        # Clean up any leftover staging directory from a prior failed attempt
        # before we start, so the swap logic sees a predictable state.
        staging_path = db_path.parent / f"{db_path.name}.new"
        if staging_path.exists():
            logger.warning(f"Removing stale database staging directory: {staging_path}")
            shutil.rmtree(staging_path, ignore_errors=True)

        errors = []
        for source in cls._resolve_database_sources(preferred_source):
            version = source["version"]
            # Archive lands outside db_path and outside staging_path, so the
            # staging tree contains only extracted content.
            archive_path = db_path.parent / source["filename"]

            try:
                cls._install_database_source_atomic(
                    db_path, source, archive_path, staging_path
                )
                logger.info(f"Database version {version} installed")
                return
            except Exception as exc:
                errors.append(f"{version}: {exc}")
                logger.warning(f"Database source {version} failed: {exc}")
            finally:
                archive_path.unlink(missing_ok=True)
                if staging_path.exists():
                    shutil.rmtree(staging_path, ignore_errors=True)

        raise RuntimeError(
            "Failed to download database. Attempts: " + "; ".join(errors)
            if errors
            else "Failed to download database: no sources succeeded"
        )

    @classmethod
    def _resolve_database_sources(
        cls, preferred_source: Optional[dict] = None
    ) -> list[dict]:
        """Resolve the one exact source selected by env, config, or the default."""
        candidate = (
            cls._source_from_env() or preferred_source or cls.DEFAULT_DATABASE_SOURCE
        )
        normalized = cls._normalize_source(candidate)
        return [normalized] if normalized else []

    @classmethod
    def _source_from_env(cls) -> Optional[dict]:
        """Read an optional database source override from environment variables."""
        url = os.environ.get("GVCLASS_DB_URL", "").strip()
        if not url:
            return None

        sha256 = os.environ.get("GVCLASS_DB_SHA256", "").strip() or None
        allow_unverified = os.environ.get(
            "GVCLASS_DB_ALLOW_UNVERIFIED", ""
        ).strip().lower() in {"1", "true", "yes", "on"}
        # Fail closed here, at source resolution, BEFORE _download_database can
        # catch a per-source failure and silently fall back to the default
        # source: an attacker-influenced GVCLASS_DB_URL with no checksum would
        # otherwise feed an unverified, pickle-bearing archive into extraction.
        if sha256 is None and not allow_unverified:
            raise RuntimeError(
                "GVCLASS_DB_URL is set but GVCLASS_DB_SHA256 is not. Refusing to "
                "download an unverified database archive. Set GVCLASS_DB_SHA256 to "
                "the expected SHA-256, or set GVCLASS_DB_ALLOW_UNVERIFIED=1 to "
                "explicitly opt in to an unverified download."
            )

        return {
            "version": os.environ.get(
                "GVCLASS_DB_VERSION",
                cls.DEFAULT_DATABASE_SOURCE["version"],
            ).strip(),
            "url": url,
            "filename": os.environ.get("GVCLASS_DB_FILENAME", "").strip() or None,
            "sha256": sha256,
            "allow_unverified": allow_unverified,
        }

    @classmethod
    def _normalize_source(cls, source: Optional[dict]) -> Optional[dict]:
        """Fill in missing download-source metadata and discard empty entries."""
        if not source:
            return None

        url = str(source.get("url", "")).strip()
        if not url:
            return None

        scheme = urlparse(url).scheme.lower()
        if scheme != "https":
            raise ValueError(
                f"Refusing database source with non-https URL scheme "
                f"'{scheme or 'none'}': {url}. Database archives must be "
                "fetched over https."
            )

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
        allow_unverified = source.get("allow_unverified") is True
        if sha256 is None and not allow_unverified:
            raise ValueError(
                f"Database source {url} has no SHA-256 checksum. "
                "Provide the expected checksum or set "
                "GVCLASS_DB_ALLOW_UNVERIFIED=1 for an explicit environment source."
            )
        if sha256 is not None and re.fullmatch(r"[0-9a-fA-F]{64}", sha256) is None:
            raise ValueError(f"Database source {url} has an invalid SHA-256 checksum")

        normalized = {
            "version": version,
            "url": url,
            "filename": filename,
            "sha256": sha256,
        }
        if "parts" in source:
            normalized["parts"] = cls._normalize_parts(source["parts"], url)
            if sha256 is None:
                raise ValueError(
                    f"Database source {url} has parts but no SHA-256 of the "
                    "joined archive"
                )
        if allow_unverified:
            normalized["allow_unverified"] = True
        return normalized

    @staticmethod
    def _normalize_parts(parts: object, url: str) -> list[dict]:
        """Validate the parts of a split archive: https URL and SHA-256 each."""
        if not isinstance(parts, list) or not parts:
            raise ValueError(f"Database source {url} lists no parts")
        normalized = []
        for part in parts:
            part_url = str(part.get("url", "")).strip()
            if urlparse(part_url).scheme.lower() != "https":
                raise ValueError(
                    f"Refusing database part with non-https URL: {part_url or 'none'}"
                )
            part_sha = str(part.get("sha256") or "").strip()
            if re.fullmatch(r"[0-9a-fA-F]{64}", part_sha) is None:
                raise ValueError(f"Database part {part_url} has no valid SHA-256")
            normalized.append({"url": part_url, "sha256": part_sha})
        return normalized

    @classmethod
    def _fetch_verified_part(cls, part: dict, number: int, part_path: Path) -> None:
        """Download one part and check its SHA-256; the file is the part on success."""
        cls._download_archive(part["url"], part_path)
        observed = cls._compute_sha256(part_path)
        if observed.lower() != part["sha256"].lower():
            raise RuntimeError(
                f"checksum mismatch in part{number:03d} "
                f"(expected {part['sha256'].lower()}, got {observed.lower()})"
            )

    @classmethod
    def _download_parts(cls, source: dict, archive_path: Path) -> None:
        """Download the parts of a split archive, verify each, join them in order."""
        part_paths: list[Path] = []
        try:
            with archive_path.open("wb") as joined:
                for number, part in enumerate(source["parts"]):
                    part_path = archive_path.with_name(
                        f"{archive_path.name}.part{number:03d}"
                    )
                    part_paths.append(part_path)
                    cls._fetch_verified_part(part, number, part_path)
                    with part_path.open("rb") as piece:
                        shutil.copyfileobj(piece, joined, _TRANSFER_CHUNK_SIZE)
                    part_path.unlink()
        except BaseException:
            archive_path.unlink(missing_ok=True)
            raise
        finally:
            for part_path in part_paths:
                part_path.unlink(missing_ok=True)

    @classmethod
    def _install_database_source_atomic(
        cls,
        db_path: Path,
        source: dict,
        archive_path: Path,
        staging_path: Path,
    ) -> None:
        """Download, verify, extract to staging, then atomically swap into db_path.

        Replaces the previous in-place extraction so a crash or verification
        failure cannot leave a half-installed database tree. See
        :meth:`_download_database` for caller-level orchestration.
        """
        version = source["version"]
        url = source["url"]

        logger.info(f"Attempting download of database {version} from {url}")
        if source.get("parts"):
            cls._download_parts(source, archive_path)
        else:
            cls._download_archive(url, archive_path)
        cls._verify_archive_checksum(source, archive_path)

        staging_path.mkdir(parents=True, exist_ok=False)
        cls._extract_archive(archive_path, staging_path)
        cls.write_database_version(staging_path, version=version)

        cls._atomic_swap_directory(staging_path, db_path)

    @classmethod
    def install_addon(cls, db_path: Path, name: str, source: dict) -> Path:
        """Install an optional archive as ``<db_path>/<name>``, replacing it whole.

        The archive holds one top-level directory. It is downloaded (in parts if the
        source lists them), checked against its SHA-256, extracted to a staging
        directory and swapped into place; a failure keeps the installed add-on.
        """
        if re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name) is None:
            raise ValueError(f"Invalid add-on name: {name!r}")
        normalized = cls._normalize_source(source)
        if normalized is None:
            raise ValueError("The add-on source has no URL")
        db_path = Path(db_path)
        archive_path = db_path / f".{name}.download"
        staging_path = db_path / f".{name}.new"
        shutil.rmtree(staging_path, ignore_errors=True)
        try:
            if normalized.get("parts"):
                cls._download_parts(normalized, archive_path)
            else:
                cls._download_archive(normalized["url"], archive_path)
            cls._verify_archive_checksum(normalized, archive_path)
            staging_path.mkdir(parents=True)
            cls._extract_archive(archive_path, staging_path)
            target = db_path / name
            cls._atomic_swap_directory(staging_path, target)
        finally:
            archive_path.unlink(missing_ok=True)
            shutil.rmtree(staging_path, ignore_errors=True)
        return target

    @classmethod
    def _atomic_swap_directory(cls, staging_path: Path, db_path: Path) -> None:
        """Swap ``staging_path`` into ``db_path`` with rollback on failure.

        If ``db_path`` already exists it is renamed aside, the staging tree
        is moved into place, and the old tree is removed. On any failure
        the previous database is restored. A parent-directory fsync on
        POSIX anchors the rename in the filesystem journal so a crash
        immediately after the swap still resolves to a durable state.
        """
        backup_path = db_path.parent / f"{db_path.name}.old"
        if backup_path.exists():
            logger.warning(f"Removing stale database backup: {backup_path}")
            shutil.rmtree(backup_path, ignore_errors=True)

        had_previous = db_path.exists()
        remove_backup = False
        try:
            if had_previous:
                os.replace(db_path, backup_path)
            os.replace(staging_path, db_path)
        except Exception:
            # Restore the previous database if we already moved it aside.
            if had_previous and backup_path.exists() and not db_path.exists():
                try:
                    os.replace(backup_path, db_path)
                except OSError as restore_exc:
                    logger.error(
                        f"Failed to restore previous database from backup: {restore_exc}"
                    )
                else:
                    remove_backup = True
            raise
        else:
            # Best-effort fsync of the parent directory so the rename is
            # durable on POSIX (no-op on platforms where O_RDONLY on dirs
            # is not supported).
            try:
                dir_fd = os.open(str(db_path.parent), os.O_RDONLY)
            except OSError:
                dir_fd = None
            if dir_fd is not None:
                try:
                    os.fsync(dir_fd)
                finally:
                    os.close(dir_fd)
            remove_backup = True
        finally:
            if remove_backup and backup_path.exists():
                shutil.rmtree(backup_path, ignore_errors=True)

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
        if any(
            re.search(r"\.tar\.gz\.part\d{3}$", candidate.get("key", ""))
            for candidate in files
        ):
            return cls._source_from_zenodo_parts(record)
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
        if sha256_file is None:
            raise ValueError(f"Zenodo archive {archive_name} has no SHA-256 sidecar")

        sha256_url = sha256_file.get("links", {}).get("self")
        if not sha256_url:
            raise ValueError(
                f"Zenodo archive {archive_name} has no SHA-256 sidecar download link"
            )
        sha256_text = cls._fetch_text(sha256_url)
        sha256 = sha256_text.split()[0] if sha256_text else ""
        if re.fullmatch(r"[0-9a-fA-F]{64}", sha256) is None:
            raise ValueError(
                f"Zenodo archive {archive_name} has an invalid SHA-256 sidecar"
            )

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
    def _source_from_zenodo_parts(cls, record: dict) -> dict:
        """Build a split-archive source from ``<name>.tar.gz.partNNN`` record files."""
        files = {f.get("key", ""): f for f in record.get("files", [])}
        pattern = re.compile(r"^(?P<name>.+\.tar\.gz)\.part(?P<number>\d{3})$")
        found = {}
        for key in files:
            match = pattern.match(key)
            if match:
                found.setdefault(match["name"], {})[int(match["number"])] = key
        if len(found) != 1:
            raise ValueError("Zenodo record holds parts of more than one archive")
        archive_name, numbered = next(iter(found.items()))
        if sorted(numbered) != list(range(len(numbered))):
            raise ValueError(
                f"Zenodo parts of {archive_name} do not have consecutive numbers"
            )
        whole = cls._sidecar_hashes(files, f"{archive_name}.sha256")
        part_hashes = cls._sidecar_hashes(files, f"{archive_name}.parts.sha256")
        parts = []
        for number in sorted(numbered):
            key = numbered[number]
            if key not in part_hashes:
                raise ValueError(f"Zenodo part {key} has no SHA-256 in the sidecar")
            parts.append(
                {"url": files[key]["links"]["self"], "sha256": part_hashes[key]}
            )
        if archive_name not in whole:
            raise ValueError(f"Zenodo archive {archive_name} has no SHA-256 sidecar")
        return {
            "version": cls.normalize_version_label(
                record.get("metadata", {}).get("version")
            ),
            "url": parts[0]["url"],
            "filename": archive_name,
            "sha256": whole[archive_name],
            "parts": parts,
            "record_id": str(record.get("id") or record.get("recid") or "").strip()
            or None,
            "conceptrecid": str(record.get("conceptrecid") or "").strip() or None,
        }

    @classmethod
    def _sidecar_hashes(cls, files: dict, key: str) -> dict[str, str]:
        """Return file name to SHA-256 from a ``sha256sum`` sidecar of a record."""
        sidecar = files.get(key)
        link = (sidecar or {}).get("links", {}).get("self")
        if not link:
            raise ValueError(f"Zenodo record has no SHA-256 sidecar {key}")
        hashes = {}
        for line in cls._fetch_text(link).splitlines():
            fields = line.split()
            if len(fields) == 2 and re.fullmatch(r"[0-9a-fA-F]{64}", fields[0]):
                hashes[fields[1].lstrip("*")] = fields[0]
        return hashes

    @classmethod
    def _download_archive(cls, url: str, archive_path: Path) -> None:
        """Download database archive with streamed progress reporting."""
        request = Request(url, headers={"User-Agent": "gvclass-database-manager"})
        try:
            with urlopen(request) as response, open(archive_path, "wb") as out_f:
                progress = _ProgressReporter(
                    "Downloading database archive",
                    total_bytes=cls._response_content_length(response),
                )
                progress.start()
                for chunk in iter(lambda: response.read(_TRANSFER_CHUNK_SIZE), b""):
                    out_f.write(chunk)
                    progress.advance(len(chunk))
                progress.finish()
        except (HTTPError, URLError, OSError, TimeoutError) as exc:
            archive_path.unlink(missing_ok=True)
            raise RuntimeError(f"download failed ({exc})") from exc

        if not archive_path.exists() or archive_path.stat().st_size == 0:
            archive_path.unlink(missing_ok=True)
            raise RuntimeError("download failed (empty archive)")

    @staticmethod
    def _response_content_length(response) -> Optional[int]:
        """Extract an integer content length from a urllib response."""
        content_length = response.headers.get("Content-Length")
        if not content_length:
            return None
        try:
            return int(content_length)
        except ValueError:
            return None

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

        observed_sha256 = cls._compute_sha256(
            archive_path, progress_label="Verifying archive checksum"
        )
        if observed_sha256.lower() != expected_sha256.lower():
            raise RuntimeError(
                "checksum mismatch "
                f"(expected {expected_sha256.lower()}, got {observed_sha256.lower()})"
            )
        logger.info(f"SHA256 verified for source {source['version']}")

    @staticmethod
    def _compute_sha256(file_path: Path, progress_label: Optional[str] = None) -> str:
        """Compute SHA256 checksum for a file."""
        digest = hashlib.sha256()
        progress = None
        if progress_label:
            progress = _ProgressReporter(
                progress_label,
                total_bytes=file_path.stat().st_size,
            )
            progress.start()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(_TRANSFER_CHUNK_SIZE), b""):
                digest.update(chunk)
                if progress is not None:
                    progress.advance(len(chunk))
        if progress is not None:
            progress.finish()
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
        prepared_members = []
        total_bytes = 0

        for member in tar.getmembers():
            relative_path = cls._normalize_member_path(member.name, strip_components)
            if relative_path is None:
                continue
            if member.issym() or member.islnk():
                raise RuntimeError(f"archive member uses links: {member.name}")

            target_path = cls._resolve_safe_target(
                base_path, relative_path, member.name
            )
            prepared_members.append((member, target_path))
            if member.isfile():
                total_bytes += member.size

        progress = _ProgressReporter(
            "Extracting database archive",
            total_bytes=total_bytes,
        )
        progress.start()

        for member, target_path in prepared_members:
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
                for chunk in iter(lambda: extracted.read(_TRANSFER_CHUNK_SIZE), b""):
                    out_f.write(chunk)
                    progress.advance(len(chunk))

        progress.finish()

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
