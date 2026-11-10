import hashlib
import io
import os
import tarfile
import threading
import time
from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

from src.utils import database_manager as database_manager_module
from src.utils.database_manager import DatabaseManager


def _create_complete_database(db_path: Path) -> None:
    (db_path / "database").mkdir(parents=True, exist_ok=True)
    for relative_path in DatabaseManager.REQUIRED_FILES:
        target = db_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test\n")
    (db_path / "DB_VERSION").write_text("v2.0.0\n")


@pytest.fixture(autouse=True)
def _isolated_source_settings(monkeypatch):
    """No pause between retries, and no source override from the environment."""
    monkeypatch.setattr(database_manager_module, "DOWNLOAD_RETRY_PAUSES_S", (0, 0, 0))
    for name in ("GVCLASS_DB_URL", "GVCLASS_DB_SHA256", "GVCLASS_DB_VERSION"):
        monkeypatch.delenv(name, raising=False)


def make_database_archive(
    path: Path,
    *,
    omit: tuple[str, ...] = (),
    extra: dict[str, bytes] | None = None,
    version: str = "v9.9.9",
) -> bytes:
    """Write a small database archive with the required files; return its bytes."""
    members = dict.fromkeys(DatabaseManager.REQUIRED_FILES, b"x\n")
    members["DB_VERSION"] = f"{version}\n".encode()
    members.update(extra or {})
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            if name in omit:
                continue
            info = tarfile.TarInfo(f"resources/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    data = path.read_bytes()
    path.unlink()
    return data


def archive_source(data: bytes, version: str = "v9.9.9") -> dict:
    return {
        "version": version,
        "url": "https://example.test/resources.tar.gz",
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def serve_archive(monkeypatch, data: bytes) -> list[Path]:
    """Answer every download with ``data``; return the paths downloaded to."""
    asked: list[Path] = []

    def fake(url: str, path: Path) -> None:
        asked.append(path)
        path.write_bytes(data)

    monkeypatch.setattr(DatabaseManager, "_download_archive", staticmethod(fake))
    return asked


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
        self._offset = min(len(self._data), start + size)
        return self._data[start : self._offset]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_setup_database_skips_download_when_database_is_complete(tmp_path, monkeypatch):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    called = {"count": 0}

    def fake_install(cls, target_path, source):
        called["count"] += 1

    monkeypatch.setattr(DatabaseManager, "_install_database", classmethod(fake_install))

    result = DatabaseManager.setup_database(str(db_path))

    assert result == db_path.resolve()
    assert called["count"] == 0


def test_setup_database_force_reinstalls_complete_database(tmp_path, monkeypatch):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    called = {"count": 0}

    def fake_install(cls, target_path, source):
        called["count"] += 1

    monkeypatch.setattr(DatabaseManager, "_install_database", classmethod(fake_install))

    result = DatabaseManager.setup_database(str(db_path), force=True)

    assert result == db_path.resolve()
    assert called["count"] == 1


def test_database_completeness_accepts_parquet_labels(tmp_path):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    (db_path / "labels.tsv").unlink()
    parquet_labels = db_path / "parquet" / "labels" / "labels.parquet"
    parquet_labels.parent.mkdir(parents=True)
    parquet_labels.write_text("placeholder\n")

    assert DatabaseManager._check_missing_files(db_path) == []


def test_database_completeness_accepts_legacy_marker_faa(tmp_path):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    (db_path / "parquet" / "faa.parquet").unlink()
    legacy_faa = db_path / "database" / "faa" / "marker.faa"
    legacy_faa.parent.mkdir(parents=True, exist_ok=True)
    legacy_faa.write_text(">reference\nMA\n")

    assert DatabaseManager._check_missing_files(db_path) == []


def test_database_completeness_rejects_empty_marker_resources(tmp_path):
    db_path = tmp_path / "resources"
    _create_complete_database(db_path)
    (db_path / "parquet" / "faa.parquet").write_text("")
    (db_path / "database" / "faa").mkdir(parents=True, exist_ok=True)

    assert "parquet/faa.parquet" in DatabaseManager._check_missing_files(db_path)


def test_get_latest_database_source_follows_zenodo_latest_record(monkeypatch):
    preferred_source = {
        "version": "v1.2.2",
        "url": "https://zenodo.org/records/18675742/files/resources_v1_2_2.tar.gz?download=1",
        "sha256": "b" * 64,
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
        staticmethod(lambda url: f"{'a' * 64} resources_v1_4_0.tar.gz\n"),
    )

    source = DatabaseManager.get_latest_database_source(preferred_source)

    assert source == {
        "version": "v1.4.0",
        "url": latest_archive_url,
        "sha256": "a" * 64,
        "record_id": "18926264",
        "conceptrecid": "18662445",
    }


def test_setup_database_does_not_fallback_from_selected_source(tmp_path, monkeypatch):
    preferred_source = {
        "version": "v8.0.0",
        "url": "https://example.com/config.tar.gz",
        "sha256": "b" * 64,
    }
    attempted_urls = []

    def fail_install(cls, db_path, source):
        attempted_urls.append(source["url"])
        raise RuntimeError("selected source unavailable")

    monkeypatch.setattr(DatabaseManager, "_install_database", classmethod(fail_install))

    with pytest.raises(RuntimeError, match="selected source unavailable"):
        DatabaseManager.setup_database(
            str(tmp_path / "resources"), preferred_source=preferred_source
        )

    assert attempted_urls == ["https://example.com/config.tar.gz"]


def test_download_archive_streams_progress_to_stderr(tmp_path, monkeypatch):
    archive_path = tmp_path / "resources.tar.gz"
    payload = b"abcdefghijklmnopqrstuvwxyz" * 1024

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


def test_download_uses_a_sixty_second_timeout(tmp_path, monkeypatch):
    timeouts = []

    def fake_urlopen(request, timeout=None):
        timeouts.append(timeout)
        return FakeResponse(b"data")

    monkeypatch.setattr(database_manager_module, "urlopen", fake_urlopen)

    DatabaseManager._download_archive("https://example.test/a", tmp_path / "a")

    assert timeouts == [60]


def _flaky_urlopen(monkeypatch, failures: list[Exception], payload: bytes) -> list[int]:
    """Raise each failure in turn, then answer with ``payload``."""
    calls: list[int] = []

    def fake_urlopen(request, timeout=None):
        calls.append(1)
        if len(calls) <= len(failures):
            raise failures[len(calls) - 1]
        return FakeResponse(payload)

    monkeypatch.setattr(database_manager_module, "urlopen", fake_urlopen)
    return calls


def _http_error(code: int) -> HTTPError:
    return HTTPError("https://example.test/a", code, "boom", {}, None)


@pytest.mark.parametrize(
    "failure", [URLError("reset"), TimeoutError("slow"), _http_error(503)]
)
def test_single_archive_download_retries_transient_errors(
    tmp_path, monkeypatch, failure
):
    data = b"archive"
    calls = _flaky_urlopen(monkeypatch, [failure, failure], data)
    archive = tmp_path / "a.tar.gz"

    DatabaseManager._fetch_archive(archive_source(data), archive)

    assert archive.read_bytes() == data
    assert len(calls) == 3


def test_download_gives_up_after_three_retries(tmp_path, monkeypatch):
    calls = _flaky_urlopen(monkeypatch, [URLError("down")] * 10, b"")

    with pytest.raises(RuntimeError, match="download failed"):
        DatabaseManager._fetch_archive(archive_source(b"x"), tmp_path / "a")

    assert len(calls) == 4


def test_download_pauses_two_eight_and_thirty_seconds(tmp_path, monkeypatch):
    pauses = []
    monkeypatch.setattr(database_manager_module, "DOWNLOAD_RETRY_PAUSES_S", (2, 8, 30))
    monkeypatch.setattr(database_manager_module.time, "sleep", pauses.append)
    _flaky_urlopen(monkeypatch, [URLError("down")] * 10, b"")

    with pytest.raises(RuntimeError, match="download failed"):
        DatabaseManager._fetch_archive(archive_source(b"x"), tmp_path / "a")

    assert pauses == [2, 8, 30]


def test_download_does_not_retry_a_missing_file(tmp_path, monkeypatch):
    calls = _flaky_urlopen(monkeypatch, [_http_error(404)] * 5, b"")

    with pytest.raises(RuntimeError, match="404"):
        DatabaseManager._fetch_archive(archive_source(b"x"), tmp_path / "a")

    assert len(calls) == 1


def test_single_archive_download_retries_a_checksum_mismatch(tmp_path, monkeypatch):
    good = b"good archive"
    answers = [b"torn", good]
    calls: list[int] = []

    def fake(url: str, path: Path) -> None:
        calls.append(1)
        path.write_bytes(answers[len(calls) - 1])

    monkeypatch.setattr(DatabaseManager, "_download_archive", staticmethod(fake))
    archive = tmp_path / "a.tar.gz"

    DatabaseManager._fetch_archive(archive_source(good), archive)

    assert archive.read_bytes() == good
    assert len(calls) == 2


def test_source_dictionaries_carry_no_filename():
    default = DatabaseManager.get_effective_database_source()

    assert "filename" not in default
    assert "filename" not in DatabaseManager.CONTIG_SCREEN_SOURCE


def test_environment_source_carries_no_filename(monkeypatch):
    monkeypatch.setenv("GVCLASS_DB_URL", "https://example.com/db.tar.gz")
    monkeypatch.setenv("GVCLASS_DB_SHA256", "a" * 64)

    source = DatabaseManager.get_effective_database_source()

    assert "filename" not in source


def test_install_stages_in_a_private_directory_and_removes_it(tmp_path, monkeypatch):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    asked = serve_archive(monkeypatch, data)
    db = tmp_path / "work" / "resources"

    DatabaseManager.setup_database(str(db), preferred_source=archive_source(data))

    assert [path.parent.parent for path in asked] == [db.parent]
    assert asked[0].parent.name.startswith(".resources.install.")
    assert sorted(path.name for path in db.parent.iterdir()) == [
        ".resources.lock",
        "resources",
    ]
    assert DatabaseManager.get_database_version(db) == "v9.9.9"


def test_install_leaves_a_file_named_like_the_archive_untouched(tmp_path, monkeypatch):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    mine = tmp_path / "resources.tar.gz"
    mine.write_text("my own download")

    DatabaseManager.setup_database(str(db), preferred_source=archive_source(data))

    assert mine.read_text() == "my own download"


def test_a_missing_required_file_keeps_the_installed_database(tmp_path, monkeypatch):
    data = make_database_archive(
        tmp_path / "fixture.tar.gz", omit=("markers/stats.tsv",)
    )
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "sentinel.txt").write_text("keep me")

    with pytest.raises(ValueError, match=r"markers/stats\.tsv"):
        DatabaseManager.setup_database(
            str(db), preferred_source=archive_source(data), force=True
        )

    assert (db / "sentinel.txt").read_text() == "keep me"
    assert (db / "DB_VERSION").read_text() == "v2.0.0\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        ".resources.lock",
        "resources",
    ]


def test_a_base_update_keeps_the_contig_screen_add_on(tmp_path, monkeypatch):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "contig_screen").mkdir()
    (db / "contig_screen" / "manifest.json").write_text("{}")

    DatabaseManager.setup_database(
        str(db), preferred_source=archive_source(data), force=True
    )

    assert (db / "contig_screen" / "manifest.json").read_text() == "{}"
    assert DatabaseManager.get_database_version(db) == "v9.9.9"


def test_an_add_on_in_the_archive_wins_over_the_installed_one(tmp_path, monkeypatch):
    data = make_database_archive(
        tmp_path / "fixture.tar.gz", extra={"contig_screen/new.txt": b"new"}
    )
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "contig_screen").mkdir()
    (db / "contig_screen" / "old.txt").write_text("old")

    DatabaseManager.setup_database(
        str(db), preferred_source=archive_source(data), force=True
    )

    assert sorted(p.name for p in (db / "contig_screen").iterdir()) == ["new.txt"]


def test_a_failed_swap_returns_the_add_on_to_the_installed_database(
    tmp_path, monkeypatch
):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "contig_screen").mkdir()
    (db / "contig_screen" / "manifest.json").write_text("{}")

    def fail(cls, staged, target, **_kwargs):
        raise OSError("simulated swap failure")

    monkeypatch.setattr(DatabaseManager, "_atomic_swap_directory", classmethod(fail))

    with pytest.raises(OSError, match="simulated swap failure"):
        DatabaseManager.setup_database(
            str(db), preferred_source=archive_source(data), force=True
        )

    assert (db / "contig_screen" / "manifest.json").read_text() == "{}"


def _wait_for(condition, what: str, timeout_s: float = 10.0) -> None:
    deadline = time.monotonic() + timeout_s
    while not condition():
        assert time.monotonic() < deadline, f"timed out waiting for {what}"
        time.sleep(0.01)


def test_a_second_installer_waits_and_skips_the_download(tmp_path, monkeypatch, caplog):
    caplog.set_level("INFO", logger=database_manager_module.logger.name)
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    source = archive_source(data)
    db = tmp_path / "resources"
    downloads: list[Path] = []
    first_is_downloading = threading.Event()
    release_first = threading.Event()

    def slow_download(url: str, path: Path) -> None:
        downloads.append(path)
        first_is_downloading.set()
        assert release_first.wait(10)
        path.write_bytes(data)

    monkeypatch.setattr(
        DatabaseManager, "_download_archive", staticmethod(slow_download)
    )
    errors: list[BaseException] = []

    def install() -> None:
        try:
            DatabaseManager.setup_database(str(db), preferred_source=source)
        except BaseException as exc:  # noqa: BLE001  # the test reads it below
            errors.append(exc)

    first = threading.Thread(target=install)
    first.start()
    assert first_is_downloading.wait(10)
    second = threading.Thread(target=install)
    second.start()
    _wait_for(lambda: "Waiting for another" in caplog.text, "the second installer")

    staging_while_waiting = list(tmp_path.glob(".resources.install.*"))
    release_first.set()
    first.join(10)
    second.join(10)

    assert errors == []
    assert len(staging_while_waiting) == 1
    assert len(downloads) == 1
    assert DatabaseManager.get_database_version(db) == "v9.9.9"
    assert not list(tmp_path.glob(".resources.install.*"))


def test_a_complete_database_is_used_without_taking_the_lock(tmp_path):
    db = tmp_path / "resources"
    _create_complete_database(db)

    DatabaseManager.setup_database(str(db))

    assert not (tmp_path / ".resources.lock").exists()


def test_an_add_on_install_waits_for_the_add_on_lock_inside_the_database(
    tmp_path, monkeypatch, caplog
):
    caplog.set_level("INFO", logger=database_manager_module.logger.name)
    db = tmp_path / "resources"
    db.mkdir()
    holder_has_lock = threading.Event()
    release = threading.Event()

    def hold() -> None:
        with DatabaseManager._addon_lock(db, "contig_screen"):
            holder_has_lock.set()
            release.wait(10)

    holder = threading.Thread(target=hold)
    holder.start()
    assert holder_has_lock.wait(10)
    fetched: list[Path] = []
    monkeypatch.setattr(
        DatabaseManager,
        "_fetch_archive",
        classmethod(lambda cls, source, path: fetched.append(path)),
    )
    installer = threading.Thread(
        target=lambda: DatabaseManager.install_addon(
            db, "contig_screen", archive_source(b"x")
        )
    )
    installer.start()
    _wait_for(lambda: "Waiting for another" in caplog.text, "the add-on installer")

    still_waiting = list(fetched)
    release.set()
    holder.join(10)
    installer.join(10)

    assert still_waiting == []
    assert (db / ".contig_screen.lock").is_file()
    assert not (tmp_path / ".resources.lock").exists()


def test_an_add_on_install_needs_no_write_access_to_the_database_parent(
    tmp_path, monkeypatch
):
    """``gvclass screen --install`` in the image binds only the database directory."""
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    parent = tmp_path / "mount"
    db = parent / "screen-resources"
    db.mkdir(parents=True)
    parent.chmod(0o555)
    try:
        DatabaseManager.install_addon(db, "contig_screen", archive_source(data))
    finally:
        parent.chmod(0o755)

    assert (db / "contig_screen" / "DB_VERSION").is_file()
    assert [path.name for path in parent.iterdir()] == ["screen-resources"]


def test_a_base_update_waits_for_the_add_on_lock_before_it_moves_the_add_on(
    tmp_path, monkeypatch, caplog
):
    caplog.set_level("INFO", logger=database_manager_module.logger.name)
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "contig_screen").mkdir()
    (db / "contig_screen" / "manifest.json").write_text("{}")
    holder_has_lock = threading.Event()
    release = threading.Event()

    def hold() -> None:
        with DatabaseManager._addon_lock(db, "contig_screen"):
            holder_has_lock.set()
            release.wait(10)

    holder = threading.Thread(target=hold)
    holder.start()
    assert holder_has_lock.wait(10)
    errors: list[BaseException] = []

    def update() -> None:
        try:
            DatabaseManager.setup_database(
                str(db), preferred_source=archive_source(data), force=True
            )
        except BaseException as exc:  # noqa: BLE001  # the test reads it below
            errors.append(exc)

    updater = threading.Thread(target=update)
    updater.start()
    _wait_for(lambda: "Waiting for another" in caplog.text, "the base update")

    add_on_in_place = (db / "contig_screen" / "manifest.json").is_file()
    version_while_waiting = DatabaseManager.get_database_version(db)
    release.set()
    holder.join(10)
    updater.join(10)

    assert errors == []
    assert add_on_in_place
    assert version_while_waiting == "v2.0.0"
    assert (db / "contig_screen" / "manifest.json").read_text() == "{}"
    assert DatabaseManager.get_database_version(db) == "v9.9.9"


def _interrupt_the_staged_rename(monkeypatch, db: Path, *, break_restore: bool):
    """Raise ``KeyboardInterrupt`` when the new tree is renamed to ``db``.

    With ``break_restore`` the rename of the backup back to ``db`` fails too.
    """
    real_replace = os.replace

    def replace(source, target):
        if Path(target).resolve() == db.resolve():
            if Path(source).name == "tree":
                raise KeyboardInterrupt
            if break_restore and ".old." in Path(source).name:
                raise OSError("simulated restore failure")
        real_replace(source, target)

    monkeypatch.setattr(os, "replace", replace)


def _database_with_add_on(tmp_path: Path) -> Path:
    db = tmp_path / "resources"
    _create_complete_database(db)
    (db / "contig_screen").mkdir()
    (db / "contig_screen" / "manifest.json").write_text("{}")
    return db


def test_an_interrupt_between_the_two_renames_restores_database_and_add_on(
    tmp_path, monkeypatch
):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = _database_with_add_on(tmp_path)
    _interrupt_the_staged_rename(monkeypatch, db, break_restore=False)

    with pytest.raises(KeyboardInterrupt):
        DatabaseManager.setup_database(
            str(db), preferred_source=archive_source(data), force=True
        )

    assert (db / "contig_screen" / "manifest.json").read_text() == "{}"
    assert DatabaseManager.get_database_version(db) == "v2.0.0"
    assert not list(tmp_path.glob(".resources.old.*"))
    assert not list(tmp_path.glob(".resources.install.*"))


def test_an_interrupt_with_a_failed_restore_keeps_the_add_on_in_the_backup(
    tmp_path, monkeypatch
):
    data = make_database_archive(tmp_path / "fixture.tar.gz")
    serve_archive(monkeypatch, data)
    db = _database_with_add_on(tmp_path)
    _interrupt_the_staged_rename(monkeypatch, db, break_restore=True)

    with pytest.raises(KeyboardInterrupt):
        DatabaseManager.setup_database(
            str(db), preferred_source=archive_source(data), force=True
        )

    (backup,) = tmp_path.glob(".resources.old.*")
    assert (backup / "contig_screen" / "manifest.json").read_text() == "{}"
    assert (backup / "DB_VERSION").read_text() == "v2.0.0\n"
    assert not list(tmp_path.glob(".resources.install.*"))
