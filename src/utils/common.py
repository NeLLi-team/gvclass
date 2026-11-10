"""Shared helpers: logger setup, path validation and file hashing."""

import hashlib
import logging
from pathlib import Path

from .error_handling import FileError, ValidationError


def setup_logging(name: str, level: str = "INFO") -> logging.Logger:
    """Return the logger ``name`` at ``level`` with one stream handler."""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(handler)

    return logger


def validate_file_path(file_path: str | Path, *, must_exist: bool = True) -> Path:
    """Return ``file_path`` resolved.

    Raises:
        ValidationError: If the path is empty or contains ``..``.
        FileError: If ``must_exist`` is true and the file is missing.
    """
    if not file_path:
        raise ValidationError("File path cannot be empty", field="file_path")

    path = Path(file_path).resolve()

    if ".." in Path(file_path).parts:
        raise ValidationError(
            f"Path traversal detected in: {file_path}",
            field="file_path",
            value=file_path,
        )

    if must_exist and not path.exists():
        raise FileError(
            f"File not found: {file_path}",
            file_path=str(file_path),
            operation="validate_path",
        )

    return path


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 hex digest of ``path``, read in chunks.

    Shared by the resume log, the bundled-model integrity gate and the
    contig-screen resource check so all integrity checks hash the same way.
    """
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()
