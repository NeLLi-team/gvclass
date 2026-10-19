"""Tests of the script that assembles the contig_screen resource tree."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts import build_contig_screen_resources as build
from src.contig_screen import layout


def make_inputs(tmp_path: Path, *, with_shared: bool) -> build.Inputs:
    src = tmp_path / "src"
    (src / "ref").mkdir(parents=True)
    for name in build.REF_FILES:
        (src / "ref" / name).write_bytes(name.encode())
    (src / "set").mkdir()
    for name in (
        "model_hits_struct_k0.pt",
        "model_hits_struct_k1.pt",
        "model_hits_k0.pt",
    ):
        (src / "set" / name).write_bytes(name.encode())
    (src / "set" / "scaler_hits_struct.json").write_text("{}")
    (src / "stage2").mkdir()
    (src / "stage2" / "stage2.joblib").write_bytes(b"s2")
    (src / "stage2" / "norm.npz").write_bytes(b"norm")
    (src / "esm").mkdir()
    for name in ("config.json", "model.safetensors", "vocab.txt"):
        (src / "esm" / name).write_bytes(name.encode())
    (src / "refs").mkdir()
    (src / "refs" / "model_to_group.json").write_text("{}")
    (src / "refs" / "A.faa").write_text(">a\nMKT\n")
    (src / "refs" / "markers.tsv").write_text("ignored\n")
    (src / "boosting.joblib").write_bytes(b"gbm")
    (src / "label_domain.parquet").write_bytes(b"ld")
    (src / "combined.hmm").write_bytes(b"hmm")
    (src / "labels.parquet").write_bytes(b"labels")
    return build.Inputs(
        ref_dir=src / "ref",
        boosting=src / "boosting.joblib",
        set_dir=src / "set",
        stage2_dir=src / "stage2",
        esm_dir=src / "esm",
        label_domain=src / "label_domain.parquet",
        refs_dir=src / "refs",
        threshold1=0.4,
        threshold2=0.3,
        hmm=src / "combined.hmm" if with_shared else None,
        labels=src / "labels.parquet" if with_shared else None,
    )


def test_the_tree_follows_the_layout_and_the_manifest_lists_every_file(
    tmp_path: Path,
) -> None:
    out = tmp_path / "bundle"

    digest = build.build_resources(make_inputs(tmp_path, with_shared=True), out)

    manifest = json.loads((out / layout.MANIFEST_FILE).read_text())
    listed = set(manifest["files"])
    assert set(layout.FIXED_FILES) <= listed
    assert f"{layout.MARKER_REFS_DIR}/A.faa" in listed
    assert f"{layout.MARKER_REFS_DIR}/markers.tsv" not in listed
    assert f"{layout.SET_MODEL_DIR}/model_hits_struct_k1.pt" in listed
    assert f"{layout.SET_MODEL_DIR}/model_hits_k0.pt" not in listed
    assert f"{layout.ESM_DIR}/model.safetensors" in listed
    assert layout.MANIFEST_FILE not in listed
    assert (
        digest == hashlib.sha256((out / layout.MANIFEST_FILE).read_bytes()).hexdigest()
    )


def test_the_manifest_records_the_size_and_digest_of_each_file(tmp_path: Path) -> None:
    out = tmp_path / "bundle"
    build.build_resources(make_inputs(tmp_path, with_shared=True), out)

    entry = json.loads((out / layout.MANIFEST_FILE).read_text())["files"][
        layout.BOOSTING_FILE
    ]

    assert entry == {"size": 3, "sha256": hashlib.sha256(b"gbm").hexdigest()}


def test_the_built_tree_opens_with_the_returned_digest_and_has_the_thresholds(
    tmp_path: Path,
) -> None:
    out = tmp_path / "bundle"
    digest = build.build_resources(make_inputs(tmp_path, with_shared=True), out)

    resources = layout.ScreenResources.open(out, manifest_sha256=digest)

    assert layout.load_thresholds(resources) == {"stage1": 0.4, "stage2": 0.3}
    assert resources.path(layout.HMM_FILE).read_bytes() == b"hmm"


def test_linking_leaves_the_inputs_in_place(tmp_path: Path) -> None:
    out = tmp_path / "bundle"
    inputs = make_inputs(tmp_path, with_shared=True)

    build.build_resources(inputs, out, link=True)

    assert (out / layout.BOOSTING_FILE).is_symlink()
    assert (out / layout.BOOSTING_FILE).resolve() == inputs.boosting.resolve()


def test_the_bundle_must_have_the_hmm_and_labels_when_they_are_not_given(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError, match="--hmm and --labels"):
        build.build_resources(
            make_inputs(tmp_path, with_shared=False), tmp_path / "bundle"
        )


def test_a_bundle_that_already_has_the_shared_files_needs_no_option(
    tmp_path: Path,
) -> None:
    out = tmp_path / "bundle"
    for relative in (layout.HMM_FILE, layout.LABELS_FILE):
        (out / relative).parent.mkdir(parents=True, exist_ok=True)
        (out / relative).write_bytes(b"shared")

    digest = build.build_resources(make_inputs(tmp_path, with_shared=False), out)

    assert layout.ScreenResources.open(out, manifest_sha256=digest).has(layout.HMM_FILE)


def test_a_missing_input_is_an_error(tmp_path: Path) -> None:
    inputs = make_inputs(tmp_path, with_shared=True)
    inputs.boosting.unlink()

    with pytest.raises(FileNotFoundError, match=r"boosting\.joblib"):
        build.build_resources(inputs, tmp_path / "bundle")


def test_rebuilding_replaces_the_files_of_an_earlier_build(tmp_path: Path) -> None:
    out = tmp_path / "bundle"
    inputs = make_inputs(tmp_path, with_shared=True)
    first = build.build_resources(inputs, out)
    inputs.boosting.write_bytes(b"new model")

    second = build.build_resources(inputs, out)

    assert first != second
    assert (out / layout.BOOSTING_FILE).read_bytes() == b"new model"


def test_main_prints_the_digest_that_belongs_in_the_layout(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    inputs = make_inputs(tmp_path, with_shared=True)
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
    ]  # fmt: skip

    code = build.main(argv)

    manifest = (tmp_path / "bundle" / layout.MANIFEST_FILE).read_bytes()
    assert code == 0
    assert hashlib.sha256(manifest).hexdigest() in capsys.readouterr().out


def test_place_refuses_a_destination_that_is_its_own_source(tmp_path: Path) -> None:
    source = tmp_path / "hmm" / "combined.hmm"
    source.parent.mkdir()
    source.write_bytes(b"hmm")

    with pytest.raises(ValueError, match="same file"):
        build.place(source, source, link=False)

    assert source.read_bytes() == b"hmm"


def test_place_refuses_a_link_that_resolves_to_its_source(tmp_path: Path) -> None:
    source = tmp_path / "a.bin"
    source.write_bytes(b"x")
    alias = tmp_path / "alias.bin"
    alias.symlink_to(source)

    with pytest.raises(ValueError, match="same file"):
        build.place(source, alias, link=True)

    assert source.read_bytes() == b"x"
