"""The archive splitter writes parts and the SHA-256 sidecars the installer reads."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from scripts.split_archive import main, split_archive


def make_archive(tmp_path: Path, size: int) -> Path:
    archive = tmp_path / "db.tar.gz"
    archive.write_bytes(bytes(range(256)) * (size // 256) + b"x" * (size % 256))
    return archive


def test_parts_joined_equal_the_archive(tmp_path):
    archive = make_archive(tmp_path, 1000)

    parts = split_archive(archive, tmp_path / "out", part_bytes=300)

    assert [p.name for p in parts] == [f"db.tar.gz.part{i:03d}" for i in range(4)]
    assert b"".join(p.read_bytes() for p in parts) == archive.read_bytes()
    assert [p.stat().st_size for p in parts] == [300, 300, 300, 100]


def test_sidecars_list_the_whole_archive_and_every_part(tmp_path):
    archive = make_archive(tmp_path, 1000)
    out = tmp_path / "out"

    parts = split_archive(archive, out, part_bytes=400)

    whole = (out / "db.tar.gz.sha256").read_text().split()
    assert whole == [hashlib.sha256(archive.read_bytes()).hexdigest(), "db.tar.gz"]
    lines = (out / "db.tar.gz.parts.sha256").read_text().splitlines()
    assert [line.split()[1] for line in lines] == [p.name for p in parts]
    for line, part in zip(lines, parts, strict=True):
        assert line.split()[0] == hashlib.sha256(part.read_bytes()).hexdigest()


def test_an_archive_that_fits_one_part_gives_one_part(tmp_path):
    archive = make_archive(tmp_path, 100)

    parts = split_archive(archive, tmp_path / "out", part_bytes=1000)

    assert [p.name for p in parts] == ["db.tar.gz.part000"]


def test_a_part_size_of_zero_is_refused(tmp_path):
    archive = make_archive(tmp_path, 100)

    with pytest.raises(ValueError, match="part size"):
        split_archive(archive, tmp_path / "out", part_bytes=0)


def test_existing_output_is_not_overwritten(tmp_path):
    archive = make_archive(tmp_path, 100)
    out = tmp_path / "out"
    split_archive(archive, out, part_bytes=50)

    with pytest.raises(FileExistsError):
        split_archive(archive, out, part_bytes=50)


def test_main_takes_the_part_size_in_mebibytes(tmp_path, capsys):
    archive = make_archive(tmp_path, 3 * 1024 * 1024 + 5)

    main([str(archive), str(tmp_path / "out"), "--part-mib", "1"])

    assert len(list((tmp_path / "out").glob("db.tar.gz.part[0-9]*"))) == 4
    assert "4 parts" in capsys.readouterr().out
