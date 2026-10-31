"""Check native mode-flag precedence through both container launchers."""

from __future__ import annotations

import runpy
import subprocess
from pathlib import Path
from unittest.mock import Mock

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("launcher_path", ["gvclass", "src/bin/gvclass_apptainer.py"])
@pytest.mark.parametrize(
    ("flags", "expected"),
    [
        ([], []),
        (["--mode-fast"], ["--mode-fast"]),
        (["--extended"], ["--extended"]),
        (["-e"], ["--extended"]),
        (["--extended", "--mode-fast"], ["--extended"]),
        (["--mode-fast", "--extended"], ["--extended"]),
        (["-e", "--mode-fast"], ["--extended"]),
        (["--mode-fast", "-e"], ["--extended"]),
    ],
)
def test_container_mode_forwarding(
    launcher_path: str,
    flags: list[str],
    expected: list[str],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Extended wins in either order; an absent choice leaves native defaults intact."""
    namespace = runpy.run_path(str(ROOT / launcher_path))
    main = namespace.get("_standalone_apptainer_main", namespace["main"])
    call = Mock(return_value=0)
    monkeypatch.setattr(subprocess, "call", call)
    query = tmp_path / "input"
    query.mkdir()

    assert (
        main(
            [
                str(query),
                "--image",
                "local.sif",
                "--resource-cache-dir",
                str(tmp_path / "cache"),
                *flags,
            ]
        )
        == 0
    )

    command = call.call_args.args[0]
    assert [arg for arg in command if arg in {"--mode-fast", "--extended"}] == expected
