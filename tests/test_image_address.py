"""The image address of the Apptainer launchers: ``oras://`` pulls and their cache keys.

Both launchers (the ``gvclass`` script that users download with wget, and
``src/bin/gvclass_apptainer.py``) are copied into an empty folder and run as programs
with a fake ``apptainer`` on the PATH. It records every call, writes a file of 60 kB
for ``pull`` and prints a valid SIF table for ``sif list``, so the real staging, the
checks and the cache of the launcher run, with no repository and no network.
"""

from __future__ import annotations

import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
COPIES = {
    "script": REPO_ROOT / "gvclass",
    "module": REPO_ROOT / "src" / "bin" / "gvclass_apptainer.py",
}
FAKE_APPTAINER = """#!/bin/sh
printf '%s\\n' "$*" >> "$FAKE_LOG"
case "$1" in
  pull) head -c 60000 /dev/zero > "$3" ;;
  sif)
    printf 'ID |GROUP |LINK |SIF POSITION (start-end) |TYPE\\n'
    printf '1 |1 |NONE |32768-33000 |Def.FILE\\n'
    printf '2 |1 |NONE |49152-50000 |FS (Squashfs/*System/amd64)\\n'
    ;;
esac
"""
GHCR = "oras://ghcr.io/nelli-team/gvclass:2.1.1"
DOCKER = "docker://ghcr.io/nelli-team/gvclass:2.1.1"
CUSTOM_LIBRARY = "library://sylabs.example.org/team/collection/gvclass:2.1.1"


class Launcher:
    """A copy of one launcher in an empty folder, run with a fake apptainer."""

    def __init__(self, tmp_path: Path, source: Path) -> None:
        self.program = tmp_path / "standalone" / source.name
        self.program.parent.mkdir()
        shutil.copy(source, self.program)
        self.work = tmp_path / "work"
        self.work.mkdir()
        self.home = tmp_path / "home"
        self.home.mkdir()
        self.log = tmp_path / "apptainer.log"
        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        fake = fake_bin / "apptainer"
        fake.write_text(FAKE_APPTAINER, encoding="utf-8")
        fake.chmod(0o755)
        self.env = {
            "PATH": f"{fake_bin}:/usr/bin:/bin",
            "HOME": str(self.home),
            "FAKE_LOG": str(self.log),
        }
        self.query = tmp_path / "bins"
        self.query.mkdir()
        (self.query / "a.fna").write_text(">c\nACGT\n", encoding="utf-8")
        self.cache_dir = self.home / ".cache" / "gvclass" / "images"
        self.name_of = runpy.run_path(str(self.program))["cache_filename_for_image"]

    def classify(self, *image: str) -> list[str]:
        """Run a classification; ``image`` is ``--image ADDRESS`` or empty."""
        subprocess.run(
            [
                sys.executable,
                str(self.program),
                str(self.query),
                "-o",
                str(self.query.parent / "results"),
                "--resource-cache-dir",
                str(self.home / "resource-cache"),
                *(["--image", image[0]] if image else []),
            ],
            capture_output=True,
            text=True,
            check=True,
            cwd=self.work,
            env=self.env,
        )
        return self.log.read_text(encoding="utf-8").splitlines()

    def pulls(self) -> list[list[str]]:
        calls = self.log.read_text(encoding="utf-8").splitlines()
        return [call.split() for call in calls if call.startswith("pull ")]


@pytest.fixture(params=sorted(COPIES))
def launcher(request: pytest.FixtureRequest, tmp_path: Path) -> Launcher:
    return Launcher(tmp_path, COPIES[request.param])


@pytest.mark.parametrize("address", [GHCR, DOCKER, CUSTOM_LIBRARY], ids=lambda a: a[:5])
def test_a_remote_address_is_pulled_plainly_and_cached_under_its_key(
    launcher: Launcher, address: str
) -> None:
    calls = launcher.classify(address)

    cached = launcher.cache_dir / launcher.name_of(address)
    [pull] = launcher.pulls()
    assert pull[:3] == ["pull", "-F", pull[2]]
    assert pull[3:] == [address]
    assert Path(pull[2]).parent.parent == launcher.cache_dir
    assert cached.stat().st_size == 60000
    assert sorted(path.name for path in launcher.cache_dir.iterdir()) == [cached.name]
    assert calls[1] == f"sif list {pull[2]}"
    assert calls[-1].split()[calls[-1].split().index(str(cached)) + 1] == "/input"


def test_the_default_image_is_the_ghcr_address_of_this_version(
    launcher: Launcher,
) -> None:
    launcher.classify()

    [pull] = launcher.pulls()
    assert pull[-1] == GHCR
    assert "--library" not in pull


def test_a_cached_image_is_checked_and_not_pulled_again(launcher: Launcher) -> None:
    launcher.classify(GHCR)
    calls = launcher.classify(GHCR)

    cached = launcher.cache_dir / launcher.name_of(GHCR)
    assert len(launcher.pulls()) == 1
    assert f"sif list {cached}" in calls


def test_addresses_with_one_readable_name_get_separate_cache_files(
    launcher: Launcher,
) -> None:
    launcher.classify(GHCR)
    launcher.classify(DOCKER)

    files = sorted(path.name for path in launcher.cache_dir.iterdir())
    assert len(launcher.pulls()) == 2
    assert files == sorted([launcher.name_of(GHCR), launcher.name_of(DOCKER)])
    assert launcher.name_of(GHCR) != launcher.name_of(DOCKER)


def test_the_library_pull_has_no_default_endpoint_and_the_image_help_says_so(
    launcher: Launcher,
) -> None:
    launcher.classify(CUSTOM_LIBRARY)

    helped = subprocess.run(
        [sys.executable, str(launcher.program), "--help"],
        capture_output=True,
        text=True,
        check=True,
        cwd=launcher.work,
        env=launcher.env,
    ).stdout
    assert "--library" not in launcher.pulls()[0]
    assert "sylabs.io" not in helped
    assert "oras://" in helped
