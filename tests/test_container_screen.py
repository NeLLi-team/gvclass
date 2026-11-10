"""Behavioral checks for the launcher's ``screen`` command and recipe source staging."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock

import pytest

from tests.launcher_helpers import load_launcher

ROOT = Path(__file__).resolve().parents[1]
IMAGE_DATABASE = "/opt/gvclass/resources"
ADDON_BIND = f"{IMAGE_DATABASE}/contig_screen:ro"


@pytest.fixture
def launch(monkeypatch: pytest.MonkeyPatch) -> tuple[ModuleType, Mock]:
    """Load the launcher without running its entry point; capture the container call."""
    module = load_launcher()
    call = Mock(return_value=7)
    monkeypatch.setattr(subprocess, "call", call)
    monkeypatch.setattr(module, "ensure_local_image", lambda image: image)
    monkeypatch.delenv("GVCLASS_DB", raising=False)
    return module, call


@pytest.fixture(params=["file", "directory"])
def screen_input(
    request: pytest.FixtureRequest, tmp_path: Path
) -> tuple[Path, str, str]:
    """Supply file and directory inputs with their expected container paths."""
    source = tmp_path / "bins.fna"
    if request.param == "file":
        source.write_text(">x\nACGT\n", encoding="utf-8")
        return source, f"{tmp_path}:/input:ro", "/input/bins.fna"
    source.mkdir()
    return source, f"{source}:/input:ro", "/input"


@pytest.mark.parametrize("position", ["before", "between", "after"])
def test_screen_binds_defaults(
    launch: tuple, tmp_path: Path, screen_input: tuple[Path, str, str], position: str
) -> None:
    module, call = launch
    source, expected_bind, expected_input = screen_input
    cache = tmp_path / "cache"
    (cache / "contig_screen").mkdir(parents=True)
    (cache / "contig_screen/manifest.json").write_text("{}", encoding="utf-8")
    flags = ["--apptainer", "--image", "image.sif", "--resource-cache-dir", str(cache)]
    argv = {
        "before": [*flags, "screen", str(source)],
        "between": ["screen", *flags, str(source)],
        "after": ["screen", str(source), *flags],
    }[position]

    assert module.container_main(argv) == 7

    cmd = call.call_args.args[0]
    assert f"{cache}/contig_screen:{ADDON_BIND}" in cmd
    assert f"{source}_screen:/output" in cmd
    assert cmd[-4:] == ["-o", "/output", "-d", IMAGE_DATABASE]
    assert expected_bind in cmd
    assert cmd[cmd.index("screen") + 1] == expected_input


@pytest.mark.parametrize("via_env", [True, False])
def test_screen_with_an_addon_only_database_keeps_the_database_of_the_image(
    launch: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, via_env: bool
) -> None:
    module, call = launch
    database = tmp_path / "db"
    (database / "contig_screen").mkdir(parents=True)
    (database / "contig_screen/manifest.json").touch()
    source = tmp_path / "input"
    source.mkdir()
    options = ["-d", str(database)]
    if via_env:
        monkeypatch.setenv("GVCLASS_DB", str(database))
        options = []

    code = module.container_main(
        [
            "screen",
            str(source),
            "--image",
            "image.sif",
            "-t",
            "9",
            "-o",
            str(tmp_path / "out"),
            *options,
        ]
    )

    cmd = call.call_args.args[0]
    assert code == 7
    assert f"{database}/contig_screen:{ADDON_BIND}" in cmd
    assert not any(item.startswith(f"{database}:") for item in cmd)
    assert f"{tmp_path}/out:/output" in cmd
    assert cmd[-4:] == ["-d", IMAGE_DATABASE, "-t", "9"]


@pytest.mark.parametrize("setting", ["default", "cache", "env", "explicit"])
def test_install_writable_without_input(
    launch: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, setting: str
) -> None:
    module, call = launch
    root = tmp_path / setting
    monkeypatch.setattr(module, "DEFAULT_RESOURCE_CACHE_DIR", root)
    options = []
    if setting == "cache":
        options = ["--resource-cache-dir", str(root)]
    if setting == "env":
        monkeypatch.setenv("GVCLASS_DB", str(root))
    if setting == "explicit":
        monkeypatch.setenv("GVCLASS_DB", str(tmp_path / "ignored"))
        options = ["-d", str(root)]
    assert (
        module.container_main(["screen", "--install", "--image", "image.sif", *options])
        == 7
    )
    assert root.is_dir()
    assert call.call_args.args[0] == [
        "apptainer",
        "run",
        "-B",
        f"{root}:/screen-resources",
        "image.sif",
        "screen",
        "--install",
        "-d",
        "/screen-resources",
    ]


def test_screen_missing_bundle_never_executes(launch: tuple, tmp_path: Path) -> None:
    module, call = launch
    with pytest.raises(SystemExit, match="screen --install"):
        module.container_main(
            ["screen", str(tmp_path), "--resource-cache-dir", str(tmp_path)]
        )
    call.assert_not_called()


def test_screen_help_without_resources(
    launch: tuple, capsys: pytest.CaptureFixture
) -> None:
    module, call = launch
    with pytest.raises(SystemExit) as error:
        module.container_main(["--image", "missing.sif", "screen", "--help"])
    assert error.value.code == 0
    assert "--install" in capsys.readouterr().out
    call.assert_not_called()


def test_image_value_screen_stays_classification(launch: tuple, tmp_path: Path) -> None:
    module, call = launch
    code = module.container_main(
        [
            str(tmp_path),
            "--image",
            "screen",
            "--resource-cache-dir",
            str(tmp_path / "cache"),
        ]
    )
    assert code == 7
    cmd = call.call_args.args[0]
    assert cmd[cmd.index("screen") + 1] == "/input"
    assert cmd[-2:] == ["-o", "/output"]


def _stage_recipe(recipe: str, stage_root: Path) -> Path:
    """Copy the source allowlist exactly as the selected recipe declares it."""
    in_files = False
    for line in (ROOT / recipe).read_text(encoding="utf-8").splitlines():
        if line.startswith("%"):
            in_files = line == "%files"
            continue
        if line.startswith("COPY "):
            sources_and_dest = shlex.split(line)[1:]
        elif in_files and line.strip():
            sources_and_dest = shlex.split(line)
        else:
            continue
        *sources, destination = sources_and_dest
        target = stage_root / destination.lstrip("/")
        for source in sources:
            origin = ROOT / source
            resolved = target / origin.name if destination.endswith("/") else target
            resolved.parent.mkdir(parents=True, exist_ok=True)
            if origin.is_dir():
                shutil.copytree(origin, resolved)
            else:
                shutil.copy2(origin, resolved)
    return stage_root / "opt/gvclass"


@pytest.mark.parametrize(
    "recipe", ["containers/docker/Dockerfile", "containers/apptainer/gvclass.def"]
)
def test_recipe_staged_screen_help(recipe: str, tmp_path: Path) -> None:
    """Run only recipe-staged sources, without imports from the checkout."""
    staged = _stage_recipe(recipe, tmp_path)
    execution = subprocess.run(
        [sys.executable, str(staged / "gvclass"), "screen", "--help"],
        cwd=staged,
        env={
            **os.environ,
            "PYTHONPATH": str(staged),
            "GVCLASS_REPO": str(staged),
            # The staged tree has no environment of its own, and the test must not
            # start a pixi installation.
            "GVCLASS_RUNTIME_SELECTED": "1",
        },
        text=True,
        capture_output=True,
        check=False,
    )
    assert execution.returncode == 0, execution.stderr
    assert "--install" in execution.stdout
