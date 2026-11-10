"""Keep archive version markers intact through database installation."""

import hashlib
import io
import json
import tarfile
from pathlib import Path

import pytest

from src.utils.database_manager import DatabaseManager


@pytest.mark.parametrize("marker", [b"v2.1.0\n", b"v2.1.0", b"2.1.0\n"])
def test_install_preserves_archived_version_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, marker: bytes
) -> None:
    """Installation must preserve every byte covered by the resource manifest."""
    manifest = {
        "DB_VERSION": {
            "sha256": hashlib.sha256(marker).hexdigest(),
            "size": len(marker),
        }
    }

    database = _install_archive(
        tmp_path,
        monkeypatch,
        {"DB_VERSION": marker, "RESOURCE_MANIFEST.json": json.dumps(manifest).encode()},
    )

    assert (database / "DB_VERSION").read_bytes() == marker
    assert json.loads((database / "RESOURCE_MANIFEST.json").read_bytes()) == manifest


@pytest.mark.parametrize("marker", [b"v2.0.0\n", b"\n"])
def test_install_rejects_version_mismatch_before_replacing_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, marker: bytes
) -> None:
    """A verified archive with the wrong version must not replace existing data."""
    database = tmp_path / "resources"
    database.mkdir()
    (database / "DB_VERSION").write_bytes(b"v1.0.0\n")
    (database / "existing.txt").write_bytes(b"keep me\n")

    with pytest.raises(RuntimeError, match="Database version mismatch"):
        _install_archive(tmp_path, monkeypatch, {"DB_VERSION": marker})

    assert (database / "DB_VERSION").read_bytes() == b"v1.0.0\n"
    assert (database / "existing.txt").read_bytes() == b"keep me\n"


def test_install_writes_missing_version_with_newline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Archives without a version marker receive the declared source version."""
    database = _install_archive(tmp_path, monkeypatch, {"content.txt": b"resource\n"})

    assert (database / "DB_VERSION").read_bytes() == b"v2.1.0\n"


def _install_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, bytes]
) -> Path:
    """Install a small local archive through the real verification and swap path."""
    archive = tmp_path / "resources.tar.gz"
    required = dict.fromkeys(DatabaseManager.REQUIRED_FILES, b"x\n")
    with tarfile.open(archive, "w:gz") as handle:
        for name, content in {**required, **files}.items():
            member = tarfile.TarInfo(f"resources/{name}")
            member.size = len(content)
            handle.addfile(member, io.BytesIO(content))
    source = {
        "url": "https://example.test/resources.tar.gz",
        "version": "v2.1.0",
        "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    }
    monkeypatch.setattr(
        DatabaseManager,
        "_download_archive",
        staticmethod(lambda url, path: path.write_bytes(archive.read_bytes())),
    )
    database = tmp_path / "resources"

    DatabaseManager._install_database(database, source)

    return database
