from pathlib import Path

from src.bin import gvclass_apptainer


def test_apptainer_wrapper_uses_release_sylabs_image_by_default() -> None:
    assert (
        gvclass_apptainer.DEFAULT_IMAGE
        == "library://nelligroup-jgi/gvclass/gvclass:2.0.1"
    )


def test_apptainer_wrapper_pulls_https_image(
    tmp_path: Path, monkeypatch
) -> None:
    calls = []
    monkeypatch.setattr(gvclass_apptainer, "IMAGE_CACHE_DIR", tmp_path / "images")
    monkeypatch.setattr(
        gvclass_apptainer.subprocess,
        "call",
        lambda cmd: calls.append(cmd) or 0,
    )

    image = "https://example.org/gvclass_test.sif"
    local_image = gvclass_apptainer.ensure_local_image(image)

    assert calls == [
        [
            "apptainer",
            "pull",
            "-F",
            local_image,
            image,
        ]
    ]


def test_apptainer_wrapper_pulls_library_image_from_sylabs(
    tmp_path: Path, monkeypatch
) -> None:
    calls = []
    monkeypatch.setattr(gvclass_apptainer, "IMAGE_CACHE_DIR", tmp_path / "images")
    monkeypatch.setattr(
        gvclass_apptainer.subprocess,
        "call",
        lambda cmd: calls.append(cmd) or 0,
    )

    local_image = gvclass_apptainer.ensure_local_image(gvclass_apptainer.SYLABS_IMAGE)

    assert calls == [
        [
            "apptainer",
            "pull",
            "-F",
            local_image,
            "--library",
            gvclass_apptainer.PUBLIC_LIBRARY_URL,
            gvclass_apptainer.SYLABS_IMAGE,
        ]
    ]


def test_apptainer_wrapper_binds_resource_cache(
    tmp_path: Path, monkeypatch
) -> None:
    query = tmp_path / "query"
    output = tmp_path / "output"
    cache = tmp_path / "cache"
    query.mkdir()
    (query / "genome.faa").write_text(">p1\nMA\n")

    calls = []
    monkeypatch.setattr(gvclass_apptainer, "ensure_local_image", lambda image: image)
    monkeypatch.setattr(
        gvclass_apptainer.subprocess,
        "call",
        lambda cmd: calls.append(cmd) or 0,
    )

    rc = gvclass_apptainer.run_container(
        query,
        output,
        8,
        gvclass_apptainer.DEFAULT_IMAGE,
        cache,
        contigs_mode=True,
        min_length=30000,
        contigs_min_length=50000,
    )

    assert rc == 0
    assert cache.is_dir()
    assert calls
    command = calls[0]
    assert f"{query.resolve()}:/input:ro" in command
    assert f"{cache.resolve()}:/resource-cache" in command
    assert "GVCLASS_RESOURCE_CACHE=/resource-cache" in command
    assert gvclass_apptainer.DEFAULT_IMAGE in command
    assert "--contigs" in command
    assert command[command.index("--min-length") + 1] == "30000"
    assert command[command.index("--contigs-min-length") + 1] == "50000"
