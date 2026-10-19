"""Tests of the manifest and the integrity checks of the screen resources."""

from __future__ import annotations

import io
import json
import pickle
from pathlib import Path

import joblib
import numpy as np
import pytest

from src.contig_screen import layout
from tests.contig_screen.conftest import BundleFactory, write_bundle

MODEL = f"{layout.SCREEN_DIR}/models/boosting.joblib"
DATA = f"{layout.SCREEN_DIR}/ref/ref_genomes.tsv"


def test_open_reads_the_manifest_when_its_digest_matches(tmp_path: Path) -> None:
    digest = write_bundle(tmp_path, {DATA: b"a\tb\n"})

    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    assert resources.has(DATA)
    assert resources.path(DATA).read_bytes() == b"a\tb\n"


def test_open_refuses_a_manifest_with_another_digest(tmp_path: Path) -> None:
    write_bundle(tmp_path, {DATA: b"x"})

    with pytest.raises(layout.IntegrityError, match="manifest SHA-256 mismatch"):
        layout.ScreenResources.open(tmp_path, manifest_sha256="0" * 64)


def test_open_fails_closed_when_no_digest_is_pinned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_bundle(tmp_path, {DATA: b"x"})
    monkeypatch.setattr(layout, "MANIFEST_SHA256", None)

    with pytest.raises(layout.IntegrityError, match="no manifest digest"):
        layout.ScreenResources.open(tmp_path)


def test_open_uses_the_pinned_constant_without_an_argument(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    digest = write_bundle(tmp_path, {DATA: b"x"})
    monkeypatch.setattr(layout, "MANIFEST_SHA256", digest)

    assert layout.ScreenResources.open(tmp_path).has(DATA)


def test_open_refuses_a_manifest_that_changed_after_pinning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    digest = write_bundle(tmp_path, {DATA: b"x"})
    monkeypatch.setattr(layout, "MANIFEST_SHA256", digest)
    manifest = tmp_path / layout.MANIFEST_FILE
    manifest.write_bytes(manifest.read_bytes() + b" ")

    with pytest.raises(layout.IntegrityError, match="mismatch"):
        layout.ScreenResources.open(tmp_path)


def test_open_refuses_a_missing_manifest(tmp_path: Path) -> None:
    with pytest.raises(layout.IntegrityError, match="manifest not found"):
        layout.ScreenResources.open(tmp_path, manifest_sha256="0" * 64)


def test_parse_manifest_refuses_paths_that_leave_the_bundle() -> None:
    raw = json.dumps(
        {"version": 1, "files": {"../x": {"size": 1, "sha256": "0" * 64}}}
    ).encode()

    with pytest.raises(layout.IntegrityError, match="unsafe path"):
        layout.parse_manifest(raw)


def test_parse_manifest_refuses_another_version() -> None:
    with pytest.raises(layout.IntegrityError, match="version"):
        layout.parse_manifest(b'{"version": 2, "files": {}}')


def test_path_refuses_a_file_the_manifest_does_not_list(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(tmp_path, {DATA: b"x"})
    (tmp_path / "extra.txt").write_text("y")

    with pytest.raises(layout.IntegrityError, match="not listed"):
        resources.path("extra.txt")


def test_path_refuses_a_listed_file_that_is_missing(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(tmp_path, {DATA: b"x"})
    (tmp_path / DATA).unlink()

    with pytest.raises(layout.IntegrityError, match="missing"):
        resources.path(DATA)


def test_path_refuses_a_file_of_another_size(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(tmp_path, {DATA: b"abc"})
    (tmp_path / DATA).write_bytes(b"abcd")

    with pytest.raises(layout.IntegrityError, match="manifest lists 3"):
        resources.path(DATA)


def test_path_refuses_a_file_with_the_same_size_and_other_content(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(tmp_path, {MODEL: b"abc"})
    (tmp_path / MODEL).write_bytes(b"abd")

    with pytest.raises(layout.IntegrityError, match="SHA-256 mismatch"):
        resources.path(MODEL)


def test_a_file_over_the_size_limit_is_checked_by_size_unless_verify_all(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dmnd = f"{layout.SCREEN_DIR}/ref/ref7.dmnd"
    digest = write_bundle(tmp_path, {dmnd: b"abcd"})
    (tmp_path / dmnd).write_bytes(b"abce")
    monkeypatch.setattr(layout, "HASH_SIZE_LIMIT", 2)

    sized = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)
    full = layout.ScreenResources.open(
        tmp_path, manifest_sha256=digest, verify_all=True
    )

    assert sized.path(dmnd) == tmp_path / dmnd
    with pytest.raises(layout.IntegrityError, match="SHA-256 mismatch"):
        full.path(dmnd)


def test_a_pickle_file_over_the_size_limit_is_still_hashed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    digest = write_bundle(tmp_path, {MODEL: b"abcd"})
    (tmp_path / MODEL).write_bytes(b"abce")
    monkeypatch.setattr(layout, "HASH_SIZE_LIMIT", 2)
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    with pytest.raises(layout.IntegrityError, match="SHA-256 mismatch"):
        resources.path(MODEL)


def test_load_joblib_returns_the_object_of_a_verified_file(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    buffer = io.BytesIO()
    joblib.dump({"columns": ["a"]}, buffer)
    resources = make_resources(tmp_path, {MODEL: buffer.getvalue()})

    assert layout.load_joblib(resources, MODEL) == {"columns": ["a"]}


def test_load_joblib_never_runs_a_swapped_pickle(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    marker = tmp_path / "ran"
    resources = make_resources(tmp_path, {MODEL: b"original"})
    (tmp_path / MODEL).write_bytes(pickle.dumps(_Exploit(marker)))

    with pytest.raises(layout.IntegrityError):
        layout.load_joblib(resources, MODEL)
    assert not marker.exists()


class _Exploit:
    """A pickle whose loading writes a marker file."""

    def __init__(self, marker: Path) -> None:
        self.marker = marker

    def __reduce__(self) -> tuple[object, tuple[str, str]]:
        return (Path.write_text, (str(self.marker), "ran"))


def test_load_npz_returns_arrays_and_refuses_pickled_objects(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    plain, pickled = io.BytesIO(), io.BytesIO()
    np.savez(plain, mu=np.array([1.0, 2.0]), comp=np.array(["c_a", "c_b"]))
    np.savez(pickled, bad=np.array([object()], dtype=object))
    resources = make_resources(
        tmp_path, {"plain.npz": plain.getvalue(), "pickled.npz": pickled.getvalue()}
    )

    arrays = layout.load_npz(resources, "plain.npz")

    assert arrays["mu"].tolist() == [1.0, 2.0]
    assert arrays["comp"].tolist() == ["c_a", "c_b"]
    with pytest.raises(ValueError, match="pickle"):
        layout.load_npz(resources, "pickled.npz")


def test_load_thresholds_reads_both_stages(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    content = b'{"stage1": 0.5, "stage2": 0.25}'
    resources = make_resources(tmp_path, {layout.THRESHOLDS_FILE: content})

    assert layout.load_thresholds(resources) == {"stage1": 0.5, "stage2": 0.25}


def test_directory_checks_every_listed_file_and_refuses_an_empty_directory(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(
        tmp_path,
        {f"{layout.ESM_DIR}/config.json": b"{}", f"{layout.ESM_DIR}/vocab.txt": b"a"},
    )
    (tmp_path / layout.ESM_DIR / "vocab.txt").write_bytes(b"b")

    with pytest.raises(layout.IntegrityError, match="SHA-256 mismatch"):
        resources.directory(layout.ESM_DIR)
    with pytest.raises(layout.IntegrityError, match="lists no file"):
        resources.directory(layout.SET_MODEL_DIR)


def test_names_lists_the_files_directly_inside_a_directory(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    resources = make_resources(
        tmp_path, {"d/b.txt": b"1", "d/a.txt": b"2", "d/sub/c.txt": b"3", "e.txt": b"4"}
    )

    assert resources.names("d") == ["d/a.txt", "d/b.txt"]


def test_directory_refuses_a_file_that_the_manifest_does_not_list(
    tmp_path: Path,
) -> None:
    listed = f"{layout.ESM_DIR}/model.safetensors"
    digest = write_bundle(tmp_path, {listed: b"weights"})
    (tmp_path / layout.ESM_DIR / "pytorch_model.bin").write_bytes(b"other")
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    with pytest.raises(layout.IntegrityError, match=r"pytorch_model\.bin"):
        resources.directory(layout.ESM_DIR)


def test_directory_accepts_a_directory_with_exactly_the_listed_files(
    tmp_path: Path,
) -> None:
    files = {
        f"{layout.ESM_DIR}/model.safetensors": b"weights",
        f"{layout.ESM_DIR}/config.json": b"{}",
    }
    digest = write_bundle(tmp_path, files)
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    assert resources.directory(layout.ESM_DIR) == tmp_path / layout.ESM_DIR
