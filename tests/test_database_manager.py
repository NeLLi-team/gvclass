import io
import tarfile
from pathlib import Path

from src.utils import database_manager as database_manager_module
from src.utils.database_manager import DatabaseManager


def _create_complete_database(db_path: Path) -> None:
    (db_path / "database").mkdir(parents=True, exist_ok=True)
    for relative_path in DatabaseManager.REQUIRED_FILES:
        target = db_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test\n")
    (db_path / "DB_VERSION").write_text("v1.2.2\n")


def test_setup_database_skips_download_when_database_is_complete(tmp_path, monkeypatch):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    called = {"count": 0}

    def fake_download(cls, target_path, preferred_source=None):
        called["count"] += 1

    monkeypatch.setattr(
        DatabaseManager, "_download_database", classmethod(fake_download)
    )

    result = DatabaseManager.setup_database(str(db_path))

    assert result == db_path.resolve()
    assert called["count"] == 0


def test_setup_database_force_reinstalls_complete_database(tmp_path, monkeypatch):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    called = {"count": 0}

    def fake_download(cls, target_path, preferred_source=None):
        called["count"] += 1

    monkeypatch.setattr(
        DatabaseManager, "_download_database", classmethod(fake_download)
    )

    result = DatabaseManager.setup_database(str(db_path), force=True)

    assert result == db_path.resolve()
    assert called["count"] == 1


def test_get_latest_database_source_follows_zenodo_latest_record(monkeypatch):
    preferred_source = {
        "version": "v1.2.2",
        "url": "https://zenodo.org/records/18675742/files/resources_v1_2_2.tar.gz?download=1",
        "filename": "resources_v1_2_2.tar.gz",
    }
    latest_archive_url = (
        "https://zenodo.org/api/records/18926264/files/resources_v1_4_0.tar.gz/content"
    )
    latest_sha_url = "https://zenodo.org/api/records/18926264/files/resources_v1_4_0.tar.gz.sha256/content"
    payloads = {
        "https://zenodo.org/api/records/18675742": {
            "links": {"latest": "https://zenodo.org/api/records/18926264"}
        },
        "https://zenodo.org/api/records/18926264": {
            "id": 18926264,
            "conceptrecid": "18662445",
            "metadata": {"version": "1.4.0"},
            "files": [
                {
                    "key": "resources_v1_4_0.tar.gz",
                    "links": {"self": latest_archive_url},
                },
                {
                    "key": "resources_v1_4_0.tar.gz.sha256",
                    "links": {"self": latest_sha_url},
                },
            ],
        },
    }

    monkeypatch.setattr(
        DatabaseManager, "_fetch_json", staticmethod(lambda url: payloads[url])
    )
    monkeypatch.setattr(
        DatabaseManager,
        "_fetch_text",
        staticmethod(lambda url: "abc123 resources_v1_4_0.tar.gz\n"),
    )

    source = DatabaseManager.get_latest_database_source(preferred_source)

    assert source == {
        "version": "v1.4.0",
        "url": latest_archive_url,
        "filename": "resources_v1_4_0.tar.gz",
        "sha256": "abc123",
        "record_id": "18926264",
        "conceptrecid": "18662445",
    }


def test_download_archive_streams_progress_to_stderr(tmp_path, monkeypatch):
    archive_path = tmp_path / "resources.tar.gz"
    payload = b"abcdefghijklmnopqrstuvwxyz" * 1024

    class FakeResponse:
        def __init__(self, data: bytes):
            self._data = data
            self._offset = 0
            self.headers = {"Content-Length": str(len(data))}

        def read(self, size: int = -1) -> bytes:
            if self._offset >= len(self._data):
                return b""
            if size < 0:
                size = len(self._data) - self._offset
            start = self._offset
            end = min(len(self._data), start + size)
            self._offset = end
            return self._data[start:end]

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class CaptureStream(io.StringIO):
        def isatty(self) -> bool:
            return False

    stderr = CaptureStream()
    monkeypatch.setattr(
        database_manager_module,
        "urlopen",
        lambda request, timeout=None: FakeResponse(payload),
    )
    monkeypatch.setattr(database_manager_module.sys, "stderr", stderr)

    DatabaseManager._download_archive(
        "https://example.test/resources.tar.gz", archive_path
    )

    assert archive_path.read_bytes() == payload
    output = stderr.getvalue()
    assert "Downloading database archive" in output
    assert "100.0%" in output


def test_extract_archive_reports_progress(tmp_path, monkeypatch):
    archive_path = tmp_path / "resources.tar.gz"
    db_path = tmp_path / "resources"
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "labels.tsv").write_text("label\n")
    (source_dir / "nested").mkdir()
    (source_dir / "nested" / "stats.tsv").write_text("stats\n")

    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(source_dir, arcname="bundle")

    class CaptureStream(io.StringIO):
        def isatty(self) -> bool:
            return False

    stderr = CaptureStream()
    monkeypatch.setattr(database_manager_module.sys, "stderr", stderr)

    DatabaseManager._extract_archive(archive_path, db_path)

    assert (db_path / "labels.tsv").read_text() == "label\n"
    assert (db_path / "nested" / "stats.tsv").read_text() == "stats\n"
    output = stderr.getvalue()
    assert "Extracting database archive" in output
    assert "100.0%" in output
