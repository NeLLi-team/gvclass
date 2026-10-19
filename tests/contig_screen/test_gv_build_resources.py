"""Tests of the giant-virus files in the contig_screen resource tree."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

from scripts import build_contig_screen_resources as build
from src.contig_screen import layout
from tests.contig_screen.test_build_resources import make_inputs


def with_giant_virus(tmp_path: Path) -> build.Inputs:
    gv_dir = tmp_path / "src" / "gv"
    gv_dir.mkdir(parents=True)
    (gv_dir / "model.joblib").write_bytes(b"gv model")
    (gv_dir / "thresholds.json").write_text('{"0.01": 0.74299}')
    base = make_inputs(tmp_path, with_shared=True)
    return dataclasses.replace(base, gv_dir=gv_dir)


def test_the_giant_virus_files_are_listed_with_the_shared_catalog_files(
    tmp_path: Path,
) -> None:
    out = tmp_path / "bundle"

    build.build_resources(with_giant_virus(tmp_path), out)

    listed = set(json.loads((out / layout.MANIFEST_FILE).read_text())["files"])
    assert set(layout.GV_FILES) <= listed
    assert (out / layout.GV_MODEL_FILE).read_bytes() == b"gv model"


def test_a_build_without_a_giant_virus_folder_lists_none_of_its_files(
    tmp_path: Path,
) -> None:
    out = tmp_path / "bundle"

    build.build_resources(make_inputs(tmp_path, with_shared=True), out)

    listed = set(json.loads((out / layout.MANIFEST_FILE).read_text())["files"])
    assert layout.GV_MODEL_FILE not in listed
    assert layout.GV_THRESHOLDS_FILE not in listed


def test_a_giant_virus_folder_without_a_model_is_an_error(tmp_path: Path) -> None:
    inputs = with_giant_virus(tmp_path)
    assert inputs.gv_dir is not None
    (inputs.gv_dir / "model.joblib").unlink()

    with pytest.raises(FileNotFoundError, match=r"model\.joblib"):
        build.build_resources(inputs, tmp_path / "bundle")


def test_main_takes_the_giant_virus_folder(tmp_path: Path) -> None:
    inputs = with_giant_virus(tmp_path)
    argv = [
        "--out-root", str(tmp_path / "bundle"),
        "--ref-dir", str(inputs.ref_dir),
        "--boosting", str(inputs.boosting),
        "--set-dir", str(inputs.set_dir),
        "--stage2-dir", str(inputs.stage2_dir),
        "--threshold1", "0.4",
        "--threshold2", "0.3",
        "--esm-dir", str(inputs.esm_dir),
        "--label-domain", str(inputs.label_domain),
        "--refs-dir", str(inputs.refs_dir),
        "--hmm", str(inputs.hmm),
        "--labels", str(inputs.labels),
        "--gv-dir", str(inputs.gv_dir),
    ]  # fmt: skip

    code = build.main(argv)

    assert code == 0
    assert (tmp_path / "bundle" / layout.GV_THRESHOLDS_FILE).is_file()
