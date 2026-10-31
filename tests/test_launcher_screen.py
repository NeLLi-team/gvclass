"""A classification run through the Apptainer launchers: ``--no-screen`` and the add-on.

The launcher exists twice: as the ``gvclass`` script, which users download with wget
and run in a folder without the repository, and as ``src/bin/gvclass_apptainer.py``.
Each test copies the file into an empty folder and runs it as a program with a fake
``apptainer`` on the PATH that records its arguments, so the path that a download
takes is the path that is tested. The ``screen`` command has its own tests in
``tests/test_container_screen.py``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FAKE_APPTAINER = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_APPTAINER_LOG"
"""
COPIES = {
    "script": REPO_ROOT / "gvclass",
    "module": REPO_ROOT / "src" / "bin" / "gvclass_apptainer.py",
}
ADDON_MOUNT = "/opt/gvclass/resources/contig_screen:ro"


class Launcher:
    """A copy of one launcher in an empty folder, run with a fake apptainer."""

    def __init__(self, tmp_path: Path, source: Path) -> None:
        self.work = tmp_path / "work"
        self.work.mkdir()
        self.log = tmp_path / "apptainer.args"
        self.program = tmp_path / "standalone" / source.name
        self.program.parent.mkdir()
        shutil.copy(source, self.program)
        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        fake = fake_bin / "apptainer"
        fake.write_text(FAKE_APPTAINER, encoding="utf-8")
        fake.chmod(0o755)
        home = tmp_path / "home"
        home.mkdir()
        self.env = {
            "PATH": f"{fake_bin}:/usr/bin:/bin",
            "HOME": str(home),
            "FAKE_APPTAINER_LOG": str(self.log),
        }

    def screen(self, *args: str, **env: str) -> list[str]:
        """Run the ``screen`` command and return the arguments that apptainer got."""
        subprocess.run(
            [sys.executable, str(self.program), "screen", *args, "--image", "x.sif"],
            capture_output=True,
            text=True,
            check=True,
            cwd=self.work,
            env={**self.env, **env},
        )
        return self.log.read_text(encoding="utf-8").splitlines()

    def classify(self, query: Path, cache: Path, *flags: str, **env: str) -> list[str]:
        """Run a classification and return the arguments that apptainer got."""
        subprocess.run(
            [
                sys.executable,
                str(self.program),
                str(query),
                "-o",
                str(query.parent / "results"),
                "--image",
                "x.sif",
                "--resource-cache-dir",
                str(cache),
                *flags,
            ],
            capture_output=True,
            text=True,
            check=True,
            cwd=self.work,
            env={**self.env, **env},
        )
        return self.log.read_text(encoding="utf-8").splitlines()


@pytest.fixture(params=sorted(COPIES))
def launcher(request: pytest.FixtureRequest, tmp_path: Path) -> Launcher:
    return Launcher(tmp_path, COPIES[request.param])


@pytest.fixture
def query(tmp_path: Path) -> Path:
    folder = tmp_path / "bins"
    folder.mkdir()
    (folder / "a.fna").write_text(">c\nACGT\n", encoding="utf-8")
    return folder


def install_addon(root: Path) -> Path:
    """Put an add-on manifest in ``root``; return the folder that is mounted."""
    folder = root / "contig_screen"
    folder.mkdir(parents=True)
    (folder / "manifest.json").write_text("{}", encoding="utf-8")
    return folder


def test_a_container_run_forwards_no_screen(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"

    without = launcher.classify(query, cache)
    with_flag = launcher.classify(query, cache, "--no-screen")

    assert "--no-screen" not in without
    assert "--no-screen" in with_flag
    assert with_flag[with_flag.index("x.sif") + 1] == "/input"


def test_a_container_run_mounts_the_addon_that_is_in_the_resource_cache(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    addon = install_addon(cache)

    recorded = launcher.classify(query, cache)

    assert f"{addon}:{ADDON_MOUNT}" in recorded
    assert recorded.index(f"{addon}:{ADDON_MOUNT}") < recorded.index("x.sif")


def test_a_container_run_without_an_addon_mounts_nothing_for_it(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()

    recorded = launcher.classify(query, cache)

    assert not any(item.endswith(ADDON_MOUNT) for item in recorded)


def test_a_folder_that_is_not_a_complete_addon_is_not_mounted(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    (cache / "contig_screen").mkdir(parents=True)

    recorded = launcher.classify(query, cache)

    assert not any(item.endswith(ADDON_MOUNT) for item in recorded)


def test_gvclass_db_selects_the_addon_like_the_screen_command_does(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    install_addon(cache)
    database = tmp_path / "db"
    addon = install_addon(database)

    recorded = launcher.classify(query, cache, GVCLASS_DB=str(database))

    assert f"{addon}:{ADDON_MOUNT}" in recorded
    assert f"{cache / 'contig_screen'}:{ADDON_MOUNT}" not in recorded


def test_gvclass_db_without_an_addon_hides_the_cache_addon(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    install_addon(cache)
    database = tmp_path / "db"
    database.mkdir()

    recorded = launcher.classify(query, cache, GVCLASS_DB=str(database))

    assert not any(item.endswith(ADDON_MOUNT) for item in recorded)


def test_the_two_launchers_choose_the_addon_mount_the_same_way() -> None:
    bodies = []
    for source in COPIES.values():
        text = source.read_text(encoding="utf-8")
        start = text.index("def _addon_mount(")
        bodies.append(text[start : text.index("\ndef ", start + 1)])

    assert bodies[0] == bodies[1]


IMAGE_DATABASE = "/opt/gvclass/resources"


def env_settings(recorded: list[str]) -> dict[str, str]:
    """Return the ``--env NAME=VALUE`` pairs of an apptainer command."""
    pairs = [recorded[i + 1] for i, item in enumerate(recorded) if item == "--env"]
    return dict(pair.split("=", 1) for pair in pairs)


@pytest.mark.parametrize("host_setting", ["set", "unset"])
def test_a_container_run_uses_the_database_of_the_image_whatever_the_host_has(
    launcher: Launcher, query: Path, tmp_path: Path, host_setting: str
) -> None:
    host_database = tmp_path / "host_db"
    install_addon(host_database)
    env = {"GVCLASS_DB": str(host_database)} if host_setting == "set" else {}

    recorded = launcher.classify(query, tmp_path / "cache", **env)

    assert env_settings(recorded)["GVCLASS_DB"] == IMAGE_DATABASE
    assert recorded.index("GVCLASS_DB=" + IMAGE_DATABASE) < recorded.index("x.sif")
    assert str(host_database) not in env_settings(recorded).values()


def test_the_host_database_still_decides_where_the_addon_comes_from(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    host_database = tmp_path / "host_db"
    addon = install_addon(host_database)

    recorded = launcher.classify(
        query, tmp_path / "cache", GVCLASS_DB=str(host_database)
    )

    assert f"{addon}:{ADDON_MOUNT}" in recorded
    assert env_settings(recorded)["GVCLASS_DB"] == IMAGE_DATABASE


def test_the_resource_cache_setting_of_a_container_run_is_unchanged(
    launcher: Launcher, query: Path, tmp_path: Path
) -> None:
    recorded = launcher.classify(query, tmp_path / "cache")

    assert env_settings(recorded) == {
        "GVCLASS_RESOURCE_CACHE": "/resource-cache",
        "GVCLASS_DB": IMAGE_DATABASE,
    }


@pytest.mark.parametrize("command", ["screen", "install"])
def test_the_screen_command_names_its_database_so_the_host_value_cannot_reach_it(
    launcher: Launcher, query: Path, tmp_path: Path, command: str
) -> None:
    host_database = tmp_path / "host_db"
    install_addon(host_database)
    cache = tmp_path / "cache"
    install_addon(cache)
    args = (
        [str(query), "-o", str(tmp_path / "out")]
        if command == "screen"
        else ["--install"]
    )

    recorded = launcher.screen(*args, "--resource-cache-dir", str(cache))
    with_host = launcher.screen(
        *args, "--resource-cache-dir", str(cache), GVCLASS_DB=str(host_database)
    )

    for call in (recorded, with_host):
        inner = call[call.index("screen") :]
        assert inner[inner.index("-d") + 1].startswith(("/opt/gvclass/", "/screen-"))
        assert str(host_database) not in inner
