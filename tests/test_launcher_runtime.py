"""The launcher: which runtime it picks, and what it passes to the container."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import NoReturn
from unittest.mock import Mock

import pytest

from src.bin.gvclass_cli import build_parser
from tests.launcher_helpers import LAUNCHER_PATH, load_launcher

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNING_PYTHON = Path("/usr/bin/python3")
PIXI = "/usr/local/bin/pixi"


@pytest.fixture
def launcher() -> ModuleType:
    return load_launcher()


@pytest.fixture
def fake_repo(tmp_path: Path) -> Path:
    """A folder that the launcher takes for a GVClass checkout, with its own Python."""
    repo = tmp_path / "repo"
    (repo / "src" / "bin").mkdir(parents=True)
    (repo / "src" / "bin" / "gvclass_cli.py").write_text("", encoding="utf-8")
    python = repo / ".pixi" / "envs" / "default" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("", encoding="utf-8")
    return repo


def decide(
    launcher: ModuleType,
    repo: Path,
    *,
    running: Path = RUNNING_PYTHON,
    pixi: str | None = PIXI,
    selected: bool = False,
) -> list[str] | None:
    return launcher.runtime_argv(
        repo,
        repo / "gvclass",
        ["bins", "-t", "4"],
        running_python=running,
        pixi=pixi,
        selected=selected,
    )


def test_the_launcher_runs_again_under_the_python_of_the_project(
    launcher: ModuleType, fake_repo: Path
) -> None:
    env_python = fake_repo / ".pixi" / "envs" / "default" / "bin" / "python"

    argv = decide(launcher, fake_repo)

    assert argv == [str(env_python), str(fake_repo / "gvclass"), "bins", "-t", "4"]


def test_the_launcher_stays_in_the_python_of_the_project(
    launcher: ModuleType, fake_repo: Path
) -> None:
    env_python = fake_repo / ".pixi" / "envs" / "default" / "bin" / "python"

    assert decide(launcher, fake_repo, running=env_python) is None


def test_the_launcher_uses_pixi_when_the_project_has_no_environment(
    launcher: ModuleType, fake_repo: Path
) -> None:
    (fake_repo / ".pixi" / "envs" / "default" / "bin" / "python").unlink()

    argv = decide(launcher, fake_repo)

    assert argv == [
        PIXI,
        "run",
        "--manifest-path",
        str(fake_repo / "pixi.toml"),
        "python",
        str(fake_repo / "gvclass"),
        "bins",
        "-t",
        "4",
    ]


def test_the_launcher_goes_on_without_environment_and_without_pixi(
    launcher: ModuleType, fake_repo: Path
) -> None:
    (fake_repo / ".pixi" / "envs" / "default" / "bin" / "python").unlink()

    assert decide(launcher, fake_repo, pixi=None) is None


@pytest.mark.parametrize("pixi", [PIXI, None])
def test_the_launcher_never_runs_again_after_it_chose_a_runtime(
    launcher: ModuleType, fake_repo: Path, pixi: str | None
) -> None:
    assert decide(launcher, fake_repo, pixi=pixi, selected=True) is None


class _Exec(Exception):  # noqa: N818  # stands for a process image that never returns
    def __init__(self, argv: list[str]) -> None:
        super().__init__(argv)
        self.argv = argv


def test_main_replaces_the_process_and_marks_the_runtime_as_chosen(
    launcher: ModuleType, fake_repo: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def replace_process(program: str, argv: list[str]) -> NoReturn:
        raise _Exec(argv)

    monkeypatch.setattr(launcher, "resolve_repo_dir", lambda script: fake_repo)
    monkeypatch.setattr(launcher.shutil, "which", lambda name: None)
    monkeypatch.setattr(os, "execvp", replace_process)
    # "0" lets the test restore the variable that the launcher sets.
    monkeypatch.setenv(launcher.RUNTIME_SELECTED_ENV, "0")

    with pytest.raises(_Exec) as raised:
        launcher.main(["bins"])

    assert raised.value.argv[0].endswith(".pixi/envs/default/bin/python")
    assert os.environ[launcher.RUNTIME_SELECTED_ENV] == "1"


def test_main_runs_the_cli_of_the_repository_once_the_runtime_is_chosen(
    launcher: ModuleType, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    monkeypatch.setenv(launcher.RUNTIME_SELECTED_ENV, "1")
    monkeypatch.setenv("GVCLASS_REPO", str(REPO_ROOT))

    with pytest.raises(SystemExit) as raised:
        launcher.main(["--help"])

    help_text = capsys.readouterr().out
    assert raised.value.code == 0
    assert "--resource-cache-dir" in help_text
    assert "--cluster-type" not in help_text


def test_the_script_prints_help_without_touching_a_database() -> None:
    result = subprocess.run(
        [sys.executable, str(LAUNCHER_PATH), "--help"],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "GVCLASS_RUNTIME_SELECTED": "1"},
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


@pytest.mark.parametrize(
    ("argv", "has_repo", "expected"),
    [
        ([], False, True),
        (["bins"], True, False),
        (["bins", "--apptainer"], True, True),
        (["bins", "--image", "x.sif"], True, True),
        (["bins", "--image=x.sif"], True, True),
        (["screen", "--install"], True, False),
    ],
)
def test_container_mode_needs_a_container_flag_or_no_repository(
    launcher: ModuleType, argv: list[str], has_repo: bool, expected: bool
) -> None:
    assert launcher._should_use_apptainer(argv, has_local_repo=has_repo) is expected


@pytest.fixture
def container_call(launcher: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Mock:
    call = Mock(return_value=0)
    monkeypatch.setattr(subprocess, "call", call)
    monkeypatch.setattr(launcher, "ensure_local_image", lambda image: image)
    return call


def test_unknown_options_pass_through_to_the_container_unchanged(
    launcher: ModuleType, container_call: Mock, tmp_path: Path
) -> None:
    query = tmp_path / "bins"
    query.mkdir()

    code = launcher.container_main(
        [
            str(query),
            "--image",
            "x.sif",
            "--resource-cache-dir",
            str(tmp_path / "cache"),
            "--future-flag",
            "value",
            "-t",
            "12",
            "--tree-method",
            "iqtree",
            "--min-length",
            "0",
            "--future-setting=3",
            "--plain-output",
        ]
    )

    command = container_call.call_args.args[0]
    tail = command[command.index("/input") + 1 :]
    assert code == 0
    assert tail[:2] == ["-o", "/output"]
    assert command[command.index("--future-flag") + 1] == "value"
    assert command[command.index("-t") + 1] == "12"
    assert command[command.index("--tree-method") + 1] == "iqtree"
    assert command[command.index("--min-length") + 1] == "0"
    assert sorted(tail[2:]) == sorted(
        [
            "--future-flag",
            "value",
            "-t",
            "12",
            "--tree-method",
            "iqtree",
            "--min-length",
            "0",
            "--future-setting=3",
            "--plain-output",
        ]
    )


def test_the_container_command_sets_no_thread_default(
    launcher: ModuleType, container_call: Mock, tmp_path: Path
) -> None:
    query = tmp_path / "bins"
    query.mkdir()

    launcher.container_main([str(query), "--image", "x.sif"])

    command = container_call.call_args.args[0]
    assert "-t" not in command
    assert "--threads" not in command


@pytest.mark.parametrize(
    "flags",
    [["-c", "x.yaml"], ["--config", "x.yaml"], ["-d", "/db"], ["--database=/db"]],
)
def test_options_of_the_native_run_are_rejected_in_container_mode(
    launcher: ModuleType, container_call: Mock, tmp_path: Path, flags: list[str]
) -> None:
    query = tmp_path / "bins"
    query.mkdir()

    with pytest.raises(SystemExit) as raised:
        launcher.container_main([str(query), "--image", "x.sif", *flags])

    assert raised.value.code == 2
    container_call.assert_not_called()


def test_a_single_file_needs_contigs_in_container_mode(
    launcher: ModuleType,
    container_call: Mock,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    source = tmp_path / "assembly.fna"
    source.write_text(">c\nACGT\n", encoding="utf-8")

    with pytest.raises(SystemExit) as raised:
        launcher.container_main([str(source), "--image", "x.sif"])

    assert raised.value.code == 2
    assert "Use --contigs" in capsys.readouterr().err
    container_call.assert_not_called()


def test_the_launcher_handles_every_cli_option_that_takes_a_value(
    launcher: ModuleType,
) -> None:
    """A new CLI option with a value must also be known to the launcher, or its value
    could be taken for the input path."""
    cli_value_options = {
        option
        for action in build_parser()._actions
        if action.nargs != 0
        for option in action.option_strings
    }
    handled = (
        set(launcher.FORWARDED_VALUE_OPTIONS)
        | set(launcher.OWN_VALUE_OPTIONS)
        | set(launcher.NATIVE_ONLY_VALUE_OPTIONS)
    )

    assert cli_value_options <= handled
    assert handled - cli_value_options == {"--image"}
