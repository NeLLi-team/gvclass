"""Split database archives: parts are verified one by one and joined in order."""

from __future__ import annotations

import hashlib
import io
import tarfile

import pytest

from src.utils.database_manager import DatabaseManager


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parts_source(chunks: list[bytes]) -> dict:
    whole = b"".join(chunks)
    return {
        "version": "v9.9.9",
        "url": "https://example.com/db.tar.gz.part000",
        "filename": "db.tar.gz",
        "sha256": sha(whole),
        "parts": [
            {
                "url": f"https://example.com/db.tar.gz.part{i:03d}",
                "sha256": sha(chunk),
            }
            for i, chunk in enumerate(chunks)
        ],
    }


def serve(monkeypatch: pytest.MonkeyPatch, payload: dict[str, bytes]) -> list[str]:
    asked: list[str] = []

    def fake(url: str, path) -> None:
        asked.append(url)
        path.write_bytes(payload[url])

    monkeypatch.setattr(DatabaseManager, "_download_archive", staticmethod(fake))
    return asked


def test_parts_are_joined_in_order(monkeypatch, tmp_path):
    chunks = [b"alpha", b"beta", b"gamma"]
    source = parts_source(chunks)
    asked = serve(
        monkeypatch, {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    )
    archive = tmp_path / "db.tar.gz"

    DatabaseManager._download_parts(source, archive)

    assert archive.read_bytes() == b"alphabetagamma"
    assert asked == [p["url"] for p in source["parts"]]
    assert not list(tmp_path.glob("*.part*"))


def test_a_corrupt_part_names_itself_and_leaves_no_archive(monkeypatch, tmp_path):
    chunks = [b"alpha", b"beta"]
    source = parts_source(chunks)
    payload = {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    payload[source["parts"][1]["url"]] = b"BETA"
    serve(monkeypatch, payload)
    archive = tmp_path / "db.tar.gz"

    with pytest.raises(RuntimeError, match="part001"):
        DatabaseManager._download_parts(source, archive)

    assert not archive.exists()
    assert not list(tmp_path.glob("*.part*"))


def test_install_of_a_split_source_checks_the_whole_archive(monkeypatch, tmp_path):
    chunks = [b"alpha", b"beta"]
    source = parts_source(chunks)
    source["sha256"] = sha(b"something else")
    serve(
        monkeypatch, {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    )
    archive = tmp_path / "db.tar.gz"

    with pytest.raises(RuntimeError, match="checksum mismatch"):
        DatabaseManager._install_database_source_atomic(
            tmp_path / "db", source, archive, tmp_path / "db.new"
        )


def test_normalize_keeps_the_parts():
    source = parts_source([b"a", b"b"])

    normalized = DatabaseManager._normalize_source(source)

    assert normalized["parts"] == source["parts"]


@pytest.mark.parametrize(
    "bad",
    [
        {"url": "http://example.com/p0", "sha256": "a" * 64},
        {"url": "https://example.com/p0", "sha256": "xyz"},
        {"url": "https://example.com/p0"},
    ],
)
def test_normalize_refuses_an_unsafe_part(bad):
    source = parts_source([b"a", b"b"])
    source["parts"][1] = bad

    with pytest.raises(ValueError, match=r"https|SHA-256|checksum"):
        DatabaseManager._normalize_source(source)


def test_normalize_refuses_empty_parts():
    source = parts_source([b"a"])
    source["parts"] = []

    with pytest.raises(ValueError, match="parts"):
        DatabaseManager._normalize_source(source)


def zenodo_record(names: list[str], sidecars: dict[str, str]) -> dict:
    files = [
        {"key": name, "links": {"self": f"https://zenodo.org/f/{name}"}}
        for name in names
    ]
    return {"id": 77, "metadata": {"version": "v9.9.9"}, "files": files}


def test_zenodo_record_with_parts_builds_a_split_source(monkeypatch):
    whole = "c" * 64
    part_hashes = ["a" * 64, "b" * 64]
    names = [
        "db.tar.gz.part000",
        "db.tar.gz.part001",
        "db.tar.gz.sha256",
        "db.tar.gz.parts.sha256",
    ]
    record = zenodo_record(names, {})
    texts = {
        "https://zenodo.org/f/db.tar.gz.sha256": f"{whole}  db.tar.gz\n",
        "https://zenodo.org/f/db.tar.gz.parts.sha256": (
            f"{part_hashes[0]}  db.tar.gz.part000\n"
            f"{part_hashes[1]}  db.tar.gz.part001\n"
        ),
    }
    monkeypatch.setattr(
        DatabaseManager, "_fetch_text", staticmethod(lambda url: texts[url])
    )

    source = DatabaseManager._source_from_zenodo_record(record)

    assert source["filename"] == "db.tar.gz"
    assert source["sha256"] == whole
    assert [p["sha256"] for p in source["parts"]] == part_hashes
    assert [p["url"] for p in source["parts"]] == [
        "https://zenodo.org/f/db.tar.gz.part000",
        "https://zenodo.org/f/db.tar.gz.part001",
    ]


def test_zenodo_record_with_a_missing_part_hash_is_refused(monkeypatch):
    names = [
        "db.tar.gz.part000",
        "db.tar.gz.part001",
        "db.tar.gz.sha256",
        "db.tar.gz.parts.sha256",
    ]
    record = zenodo_record(names, {})
    texts = {
        "https://zenodo.org/f/db.tar.gz.sha256": f"{'c' * 64}  db.tar.gz\n",
        "https://zenodo.org/f/db.tar.gz.parts.sha256": (
            f"{'a' * 64}  db.tar.gz.part000\n"
        ),
    }
    monkeypatch.setattr(
        DatabaseManager, "_fetch_text", staticmethod(lambda url: texts[url])
    )

    with pytest.raises(ValueError, match="part001"):
        DatabaseManager._source_from_zenodo_record(record)


def test_zenodo_record_with_a_gap_in_the_part_numbers_is_refused(monkeypatch):
    names = [
        "db.tar.gz.part000",
        "db.tar.gz.part002",
        "db.tar.gz.sha256",
        "db.tar.gz.parts.sha256",
    ]
    record = zenodo_record(names, {})
    texts = {
        "https://zenodo.org/f/db.tar.gz.sha256": f"{'c' * 64}  db.tar.gz\n",
        "https://zenodo.org/f/db.tar.gz.parts.sha256": (
            f"{'a' * 64}  db.tar.gz.part000\n{'b' * 64}  db.tar.gz.part002\n"
        ),
    }
    monkeypatch.setattr(
        DatabaseManager, "_fetch_text", staticmethod(lambda url: texts[url])
    )

    with pytest.raises(ValueError, match="consecutive"):
        DatabaseManager._source_from_zenodo_record(record)


def make_addon_archive(tmp_path, files):
    path = tmp_path / "addon.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(f"contig_screen/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


def addon_source(archive, chunk=200):
    data = archive.read_bytes()
    chunks = [data[i : i + chunk] for i in range(0, len(data), chunk)]
    source = parts_source(chunks)
    source["filename"] = archive.name
    return source, chunks


def test_an_addon_is_installed_beside_the_database(monkeypatch, tmp_path):
    archive = make_addon_archive(tmp_path, {"thresholds.json": b"{}", "m/x.bin": b"xx"})
    source, chunks = addon_source(archive)
    serve(
        monkeypatch, {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    )
    db = tmp_path / "db"
    db.mkdir()

    target = DatabaseManager.install_addon(db, "contig_screen", source)

    assert target == db / "contig_screen"
    assert (target / "thresholds.json").read_bytes() == b"{}"
    assert (target / "m" / "x.bin").read_bytes() == b"xx"
    assert not list(db.glob(".contig_screen*"))


def test_a_new_addon_replaces_the_old_one_whole(monkeypatch, tmp_path):
    db = tmp_path / "db"
    (db / "contig_screen").mkdir(parents=True)
    (db / "contig_screen" / "old.txt").write_text("old")
    archive = make_addon_archive(tmp_path, {"new.txt": b"new"})
    source, chunks = addon_source(archive)
    serve(
        monkeypatch, {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    )

    DatabaseManager.install_addon(db, "contig_screen", source)

    assert sorted(p.name for p in (db / "contig_screen").iterdir()) == ["new.txt"]


def test_a_failed_addon_download_keeps_the_installed_one(monkeypatch, tmp_path):
    db = tmp_path / "db"
    (db / "contig_screen").mkdir(parents=True)
    (db / "contig_screen" / "old.txt").write_text("old")
    archive = make_addon_archive(tmp_path, {"new.txt": b"new"})
    source, chunks = addon_source(archive)
    payload = {p["url"]: c for p, c in zip(source["parts"], chunks, strict=True)}
    payload[source["parts"][0]["url"]] = b"corrupt" + chunks[0][7:]
    serve(monkeypatch, payload)

    with pytest.raises(RuntimeError, match="part000"):
        DatabaseManager.install_addon(db, "contig_screen", source)

    assert (db / "contig_screen" / "old.txt").read_text() == "old"
    assert not list(db.glob(".contig_screen*"))


def test_an_addon_name_with_a_path_is_refused(tmp_path):
    with pytest.raises(ValueError, match="name"):
        DatabaseManager.install_addon(tmp_path, "../x", parts_source([b"a"]))


def test_the_built_in_database_source_is_the_split_archive_with_pinned_hashes():
    source = DatabaseManager._normalize_source(DatabaseManager.DEFAULT_DATABASE_SOURCE)

    assert source["filename"] == "resources_v2_1_0.tar.gz"
    assert len(source["sha256"]) == 64
    assert len(source["parts"]) == 24
    assert all(
        p["url"].startswith("https://zenodo.org/records/") for p in source["parts"]
    )
    assert [p["url"].split("/files/")[1].split("?")[0] for p in source["parts"]] == [
        f"resources_v2_1_0.tar.gz.part{i:03d}" for i in range(24)
    ]


def test_the_contig_screen_source_is_the_split_add_on_archive():
    source = DatabaseManager._normalize_source(DatabaseManager.CONTIG_SCREEN_SOURCE)

    assert source["filename"] == "contig_screen_v2_1_0.tar.gz"
    assert len(source["parts"]) == 48
    assert all(len(p["sha256"]) == 64 for p in source["parts"])
    records = {p["url"].split("/")[4] for p in source["parts"]}
    assert records == {"23464780"}
