"""Shared helpers: path validation and file hashing."""

import hashlib
from collections.abc import Callable
from pathlib import Path

from .error_handling import FileError, ValidationError

_HASH_CHUNK_BYTES = 1 << 20


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


def sha256_file(
    path: str | Path, *, on_chunk: Callable[[int], None] | None = None
) -> str:
    """Return the SHA-256 hex digest of ``path``, read in 1 MiB chunks.

    Shared by the resume log, the bundled-model integrity gate, the contig-screen
    resource check and the database download so all integrity checks hash the same
    way.

    Args:
        path: File to hash.
        on_chunk: Called with the size of each chunk read, for progress reports.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(_HASH_CHUNK_BYTES):
            digest.update(chunk)
            if on_chunk is not None:
                on_chunk(len(chunk))
    return digest.hexdigest()
