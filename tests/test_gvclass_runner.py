"""The runner entry point: read ``run_config.json``, log to a file, run the flow."""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

import pytest

from src.bin import gvclass_runner
from src.pipeline.run_config import RunConfig
from tests.pipeline_helpers import make_config

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def config_file(tmp_path: Path) -> Path:
    path = tmp_path / "run_config.json"
    path.write_text(make_config(tmp_path).to_json(), encoding="utf-8")
    return path


def test_main_runs_the_flow_on_the_configuration_in_the_file(
    config_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    isolated_logging: None,
) -> None:
    seen: list[RunConfig] = []
    monkeypatch.setattr(
        gvclass_runner,
        "gvclass_flow",
        lambda config: seen.append(config) or tmp_path / "summary.tsv",
    )

    code = gvclass_runner.main([str(config_file)])

    assert code == 0
    assert seen == [make_config(tmp_path)]


def test_main_logs_to_the_output_directory_and_to_stderr(
    config_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    isolated_logging: None,
) -> None:
    def flow(config: RunConfig) -> Path:
        logging.getLogger("gvclass_runner").info("flow started")
        return tmp_path / "summary.tsv"

    monkeypatch.setattr(gvclass_runner, "gvclass_flow", flow)

    gvclass_runner.main([str(config_file)])

    log = (tmp_path / "out" / "gvclass.log").read_text(encoding="utf-8")
    assert "flow started" in log
    assert "flow started" in capsys.readouterr().err


def test_main_logs_the_failure_and_returns_1(
    config_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    isolated_logging: None,
) -> None:
    def flow(config: RunConfig) -> Path:
        raise RuntimeError("tree failed")

    monkeypatch.setattr(gvclass_runner, "gvclass_flow", flow)

    code = gvclass_runner.main([str(config_file)])

    assert code == 1
    assert "tree failed" in capsys.readouterr().err
    assert "tree failed" in (tmp_path / "out" / "gvclass.log").read_text(
        encoding="utf-8"
    )


@pytest.mark.parametrize(
    ("flags", "level"), [([], logging.INFO), (["--verbose"], logging.DEBUG)]
)
def test_main_sets_the_log_level_from_the_verbose_flag(
    config_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    isolated_logging: None,
    flags: list[str],
    level: int,
) -> None:
    monkeypatch.setattr(gvclass_runner, "gvclass_flow", lambda config: tmp_path)

    gvclass_runner.main([str(config_file), *flags])

    assert logging.getLogger("gvclass_runner").level == level


def test_the_runner_does_not_load_click() -> None:
    code = (
        "import sys; sys.path.insert(0, sys.argv[1]); import src.bin.gvclass_runner; "
        "print('click' in sys.modules)"
    )

    result = subprocess.run(
        [sys.executable, "-c", code, str(REPO_ROOT)],
        capture_output=True,
        text=True,
        check=True,
    )

    assert result.stdout.strip() == "False"
