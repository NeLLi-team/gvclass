"""Run ``gvclass_cli.main`` with the database setup and the runner process replaced.

The fakes record what the command would have started, so a test sees the run
configuration and the runner command without any pipeline work.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

from src.bin import gvclass_cli
from src.bin.gvclass_cli import RUN_CONFIG_FILENAME, RunnerOutcome
from src.pipeline.run_config import RunConfig

SUCCESS = RunnerOutcome(return_code=0, tail=(), interrupted=False, peak_memory_mb=1.0)
VALID = {"success": True, "issues": [], "warnings": []}


@dataclass
class CliRun:
    """What one call of ``main`` did."""

    exit_code: int
    output_dir: Path
    runner_calls: list[dict[str, Any]] = field(default_factory=list)

    @property
    def config(self) -> RunConfig:
        """The run configuration that the command wrote for the runner."""
        text = (self.output_dir / RUN_CONFIG_FILENAME).read_text(encoding="utf-8")
        return RunConfig.from_json(text)


def run_cli(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *args: str,
    n_queries: int = 2,
    suffix: str = ".faa",
    body: str = ">p\nM\n",
    outcome: RunnerOutcome = SUCCESS,
    runner_error: Exception | None = None,
    validation: dict[str, Any] | None = None,
    database_check: Callable[..., bool] = lambda *a, **k: True,
) -> CliRun:
    """Call ``main`` on a folder of ``n_queries`` small files, with fakes."""
    queries = tmp_path / "queries"
    queries.mkdir(parents=True, exist_ok=True)
    for index in range(n_queries):
        (queries / f"q{index}{suffix}").write_text(body, encoding="utf-8")
    database = tmp_path / "db"
    database.mkdir(parents=True, exist_ok=True)
    output_dir = tmp_path / "out"
    run = CliRun(exit_code=-1, output_dir=output_dir)

    def fake_run_runner(command: list[str], **kwargs: Any) -> RunnerOutcome:
        run.runner_calls.append({"command": command, **kwargs})
        if runner_error is not None:
            raise runner_error
        return outcome

    monkeypatch.setattr(gvclass_cli, "check_and_setup_database", database_check)
    monkeypatch.setattr(gvclass_cli, "run_runner", fake_run_runner)
    monkeypatch.setattr(
        "src.utils.output_validator.validate_pipeline_outputs",
        lambda **kwargs: validation or VALID,
    )
    monkeypatch.setenv("GVCLASS_RESOURCE_CACHE", str(tmp_path / "cache"))
    monkeypatch.setenv("GVCLASS_PLAIN_OUTPUT", "1")
    argv = [str(queries), "-o", str(output_dir), "-d", str(database), *args]
    run.exit_code = gvclass_cli.main(argv)
    return run
