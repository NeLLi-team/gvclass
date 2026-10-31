import runpy
import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest

from src.bin import gvclass_apptainer

ROOT = Path(__file__).resolve().parents[1]
SIF_LIST = (
    "ID |GROUP |LINK |SIF POSITION (start-end) |TYPE\n"
    "1 |1 |NONE |32768-33000 |Def.FILE\n"
    "2 |1 |NONE |49152-50000 |FS (Squashfs/*System/amd64)\n"
)
IMAGE = "https://example.org/gvclass_test.sif"
GHCR_IMAGE = "oras://ghcr.io/nelli-team/gvclass:2.1.1"
CUSTOM_LIBRARY_IMAGE = "library://sylabs.example.org/team/collection/gvclass:2.1.1"


@pytest.fixture(params=["gvclass", "src/bin/gvclass_apptainer.py"])
def image_launcher(request, tmp_path: Path, monkeypatch) -> Callable[[str], str]:
    """Load each independently distributed launcher with an isolated cache."""
    namespace = runpy.run_path(str(ROOT / request.param))
    ensure = namespace["ensure_local_image"]
    monkeypatch.setitem(ensure.__globals__, "IMAGE_CACHE_DIR", tmp_path / "images")
    return ensure


@pytest.fixture
def sif_listing(monkeypatch) -> list[list[str]]:
    """Stub only the external SIF table reader; real file bounds still apply."""
    calls = []

    def inspect(command, **kwargs):
        calls.append(command)
        assert command[:3] == ["apptainer", "sif", "list"]
        return subprocess.CompletedProcess(command, 0, SIF_LIST, "")

    monkeypatch.setattr(subprocess, "run", inspect)
    return calls


def cache_path(launcher: Callable[[str], str], image: str = IMAGE) -> Path:
    """Find the documented cache entry without running the launcher."""
    state = launcher.__globals__
    return state["IMAGE_CACHE_DIR"] / state["cache_filename_for_image"](image)


def test_apptainer_wrapper_uses_the_release_ghcr_image_by_default() -> None:
    assert gvclass_apptainer.DEFAULT_IMAGE == GHCR_IMAGE


def test_the_launcher_has_no_sylabs_constants() -> None:
    for name in ("SYLABS_IMAGE", "PUBLIC_LIBRARY_URL", "is_library_uri"):
        assert not hasattr(gvclass_apptainer, name), name


@pytest.mark.parametrize(
    "image",
    [
        IMAGE,
        GHCR_IMAGE,
        "docker://ghcr.io/nelli-team/gvclass:2.1.1",
        CUSTOM_LIBRARY_IMAGE,
    ],
    ids=["https", "oras", "docker", "library"],
)
def test_image_pull_is_staged_and_checked(
    image_launcher, sif_listing, monkeypatch, image: str
) -> None:
    calls = []
    cached = cache_path(image_launcher, image)

    def pull(command):
        calls.append(command)
        staged = Path(command[3])
        assert staged != cached
        assert staged.parent.parent == cached.parent
        assert not cached.exists()
        staged.write_bytes(b"\0" * 50000)
        return 0

    monkeypatch.setattr(subprocess, "call", pull)
    assert image_launcher(image) == str(cached)
    assert cached.stat().st_size == 50000
    assert len(calls) == 1
    assert calls[0] == ["apptainer", "pull", "-F", calls[0][3], image]
    assert sif_listing == [["apptainer", "sif", "list", calls[0][3]]]
    assert list(cached.parent.iterdir()) == [cached]


def test_valid_cached_image_is_checked_without_pull(
    image_launcher, sif_listing, monkeypatch
) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    cached.write_bytes(b"\0" * 50000)
    monkeypatch.setattr(subprocess, "call", lambda cmd: pytest.fail("Unexpected pull"))

    assert image_launcher(IMAGE) == str(cached)
    assert sif_listing == [["apptainer", "sif", "list", str(cached)]]


def test_invalid_cached_bytes_are_preserved(
    image_launcher, sif_listing, monkeypatch
) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    cached.write_bytes(b"not a sif\n")
    monkeypatch.setattr(subprocess, "call", lambda cmd: pytest.fail("Unexpected pull"))

    with pytest.raises(SystemExit, match=r"Invalid.*image"):
        image_launcher(IMAGE)
    assert cached.read_bytes() == b"not a sif\n"


@pytest.mark.parametrize(
    ("target_name", "target_exists"), [("missing.sif", False), ("user-image.sif", True)]
)
def test_cached_symlink_is_preserved_and_rejected(
    image_launcher,
    sif_listing,
    monkeypatch,
    tmp_path: Path,
    target_name: str,
    target_exists: bool,
) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    target = tmp_path / target_name
    (tmp_path / "user-image.sif").write_bytes(b"\0" * 50000)
    cached.symlink_to(target)
    monkeypatch.setattr(subprocess, "call", lambda cmd: pytest.fail("Unexpected pull"))

    with pytest.raises(SystemExit, match=r"Invalid.*image"):
        image_launcher(IMAGE)
    assert cached.is_symlink()
    assert target.exists() == target_exists
    assert not sif_listing


def test_cached_directory_is_preserved_and_rejected(
    image_launcher, sif_listing
) -> None:
    cached = cache_path(image_launcher)
    cached.mkdir(parents=True)
    (cached / "keep.txt").write_text("keep", encoding="utf-8")

    with pytest.raises(SystemExit, match=r"Invalid.*image"):
        image_launcher(IMAGE)
    assert (cached / "keep.txt").read_text(encoding="utf-8") == "keep"
    assert not sif_listing


@pytest.mark.parametrize("pull_status", [0, 7])
def test_failed_or_partial_pull_never_enters_cache(
    image_launcher, sif_listing, monkeypatch, pull_status: int
) -> None:
    cached = cache_path(image_launcher)

    def pull(command):
        Path(command[3]).write_bytes(b"partial")
        return pull_status

    monkeypatch.setattr(subprocess, "call", pull)
    with pytest.raises(SystemExit):
        image_launcher(IMAGE)
    assert not cached.exists()
    assert list(cached.parent.iterdir()) == []


def test_concurrent_valid_cache_entry_is_not_overwritten(
    image_launcher, sif_listing, monkeypatch
) -> None:
    cached = cache_path(image_launcher)
    existing = b"winner" + b"\0" * 49994

    def pull(command):
        Path(command[3]).write_bytes(b"\0" * 50000)
        cached.write_bytes(existing)
        return 0

    monkeypatch.setattr(subprocess, "call", pull)
    assert image_launcher(IMAGE) == str(cached)
    assert cached.read_bytes() == existing
    assert list(cached.parent.iterdir()) == [cached]


def test_concurrent_invalid_cache_entry_is_preserved(
    image_launcher, sif_listing, monkeypatch
) -> None:
    cached = cache_path(image_launcher)

    def pull(command):
        Path(command[3]).write_bytes(b"\0" * 50000)
        cached.write_bytes(b"user data")
        return 0

    monkeypatch.setattr(subprocess, "call", pull)
    with pytest.raises(SystemExit, match=r"Invalid.*image"):
        image_launcher(IMAGE)
    assert cached.read_bytes() == b"user data"
    assert list(cached.parent.iterdir()) == [cached]


@pytest.mark.parametrize(
    ("listing", "returncode"),
    [
        ("", 0),
        (SIF_LIST.replace("*System", "Data"), 0),
        (SIF_LIST, 1),
        (SIF_LIST.replace("49152-50000", "49152-49152"), 0),
    ],
)
def test_invalid_sif_listing_is_rejected(
    image_launcher, monkeypatch, listing: str, returncode: int
) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    cached.write_bytes(b"\0" * 50000)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda command, **kwargs: subprocess.CompletedProcess(
            command, returncode, listing, ""
        ),
    )

    with pytest.raises(SystemExit, match=r"Invalid.*image"):
        image_launcher(IMAGE)


def test_empty_metadata_object_is_allowed(image_launcher, monkeypatch) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    cached.write_bytes(b"\0" * 50000)
    listing = SIF_LIST.replace("32768-33000", "32768-32768")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, listing, ""),
    )

    assert image_launcher(IMAGE) == str(cached)


@pytest.mark.parametrize(
    "error",
    [FileNotFoundError("apptainer"), subprocess.TimeoutExpired("apptainer", 30)],
)
def test_reader_errors_have_a_clear_exit(image_launcher, monkeypatch, error) -> None:
    cached = cache_path(image_launcher)
    cached.parent.mkdir()
    cached.write_bytes(b"\0" * 50000)

    def inspect(command, **kwargs):
        raise error

    monkeypatch.setattr(subprocess, "run", inspect)
    with pytest.raises(SystemExit, match="Could not inspect Apptainer image"):
        image_launcher(IMAGE)
    assert cached.stat().st_size == 50000


def test_unsupported_cache_filesystem_has_a_clear_exit(
    image_launcher, sif_listing, monkeypatch
) -> None:
    cached = cache_path(image_launcher)

    def pull(command):
        Path(command[3]).write_bytes(b"\0" * 50000)
        return 0

    def install(*args, **kwargs):
        raise OSError("hard links unsupported")

    monkeypatch.setattr(subprocess, "call", pull)
    monkeypatch.setattr(Path, "hardlink_to", install)
    with pytest.raises(SystemExit, match="Could not install cached Apptainer image"):
        image_launcher(IMAGE)
    assert list(cached.parent.iterdir()) == []


def test_explicit_local_image_keeps_existing_behavior(
    image_launcher, tmp_path: Path
) -> None:
    local = str(tmp_path / "image with spaces.sif")
    assert image_launcher(local) == local


def test_apptainer_wrapper_binds_resource_cache(tmp_path: Path, monkeypatch) -> None:
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
        species_tree_combined=True,
        resume=True,
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
    assert "--species-tree-combined" in command
    assert "--resume" in command
    assert command[command.index("--min-length") + 1] == "30000"
    assert command[command.index("--contigs-min-length") + 1] == "50000"
