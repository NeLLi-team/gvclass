"""The container launcher passes mode flags on as typed; the CLI in the image decides.

The precedence of ``--extended`` over ``--mode-fast``, and of ``--no-sensitive``
over ``--sensitive``, is tested for the command line in ``tests/test_gvclass_cli.py``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import Mock

import pytest

from tests.launcher_helpers import load_launcher


@pytest.mark.parametrize(
    "flags",
    [
        [],
        ["--mode-fast"],
        ["--extended"],
        ["-e"],
        ["--extended", "--mode-fast"],
        ["--mode-fast", "--extended"],
        ["-e", "--mode-fast"],
        ["--mode-fast", "-e"],
        ["--sensitive", "--no-sensitive"],
    ],
)
def test_container_mode_flags_reach_the_image_as_typed(
    flags: list[str], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    launcher = load_launcher()
    call = Mock(return_value=0)
    monkeypatch.setattr(subprocess, "call", call)
    query = tmp_path / "input"
    query.mkdir()

    code = launcher.container_main(
        [
            str(query),
            "--image",
            "local.sif",
            "--resource-cache-dir",
            str(tmp_path / "cache"),
            *flags,
        ]
    )

    command = call.call_args.args[0]
    assert code == 0
    assert command[command.index("/input") + 1 :] == ["-o", "/output", *flags]
