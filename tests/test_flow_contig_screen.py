"""The run flow screens every query after the queries finish, before the summary."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from src.contig_screen import layout
from src.pipeline import parallel_runner
from src.pipeline.contig_screen_step import ScreenOptions
from src.pipeline.run_config import RunConfig
from tests.fake_screen import install_fake_screen
from tests.pipeline_helpers import fake_flow, make_config, write_queries


class StepSpy:
    """Stands in for the screen step; records its call and sets ``screen_status``."""

    def __init__(self, output: Path) -> None:
        self.output = output
        self.queries: list[str] = []
        self.query_files: dict[str, Path] = {}
        self.options: ScreenOptions | None = None
        self.database: Path | None = None
        self.summary_existed = False

    def __call__(
        self,
        results: list[dict[str, Any]],
        query_files: dict[str, Path],
        output_dir: Path,
        database: Path,
        options: ScreenOptions,
        logger: object,
    ) -> None:
        self.queries = [result["query"] for result in results]
        self.query_files = dict(query_files)
        self.options = options
        self.database = database
        self.summary_existed = (self.output / "gvclass_summary.tsv").exists()
        for result in results:
            result["summary_data"]["screen_status"] = "ok"
            result["summary_data"]["screen_flagged_pct"] = 12.345


@pytest.fixture
def run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> RunConfig:
    fake_flow(monkeypatch)
    return make_config(tmp_path, threads=8, max_workers=None, threads_per_worker=None)


def read_rows(summary_file: Path) -> list[dict[str, str]]:
    with summary_file.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def refuse_to_process(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: object, **kwargs: object) -> None:
        msg = "resume must not process a completed query"
        raise AssertionError(msg)

    monkeypatch.setattr(parallel_runner, "process_single_query", refuse)


def flow(config: RunConfig, **options: Any) -> list[dict[str, str]]:
    return read_rows(parallel_runner.gvclass_flow(replace(config, **options)))


def install_addon(config: RunConfig) -> Path:
    """Put an add-on manifest where the flow looks for it; return the screen folder."""
    manifest = config.database_path / "contig_screen" / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("{}", encoding="utf-8")
    return config.output_dir / "contig_screen"


def test_every_query_is_screened_before_the_final_summary_is_written(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1", "q2"], ".fna")
    spy = StepSpy(run.output_dir)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    rows = flow(run)

    assert sorted(spy.queries) == ["q1", "q2"]
    assert not spy.summary_existed
    assert spy.database == run.database_path
    assert [row["screen_status"] for row in rows] == ["ok", "ok"]
    assert [row["screen_flagged_pct"] for row in rows] == ["12.35", "12.35"]


def test_queries_skipped_by_resume_are_screened_too(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1"], ".fna")
    flow(run)
    write_queries(run, ["q2", "q3"], ".fna")
    spy = StepSpy(run.output_dir)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    rows = flow(run, resume=True)

    assert sorted(spy.queries) == ["q1", "q2", "q3"]
    assert sorted(spy.query_files) == ["q1", "q2", "q3"]
    assert spy.options is not None
    assert spy.options.resume is True
    assert [row["screen_status"] for row in rows] == ["ok"] * 3


def test_a_run_in_which_resume_skips_every_query_is_screened(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1", "q2"], ".fna")
    flow(run)
    refuse_to_process(monkeypatch)
    spy = StepSpy(run.output_dir)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    rows = flow(run, resume=True)

    assert sorted(spy.queries) == ["q1", "q2"]
    assert [row["screen_status"] for row in rows] == ["ok", "ok"]


def test_the_step_gets_the_screen_options_and_the_worker_plan_of_the_run(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1", "q2"], ".fna")
    spy = StepSpy(run.output_dir)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    flow(run, contig_screen=False, contigs_mode=True)

    options = spy.options
    assert options is not None
    assert (options.enabled, options.contigs_mode, options.resume) == (
        False,
        True,
        False,
    )
    assert options.plan_workers(2) == (2, 4)
    assert options.plan_workers(1) == (1, 8)
    assert options.plan_workers(16) == (2, 4)


def test_without_the_addon_the_summary_says_not_installed_or_not_applicable(
    run: RunConfig,
) -> None:
    write_queries(run, ["bin_a"], ".fna")
    write_queries(run, ["bin_b"], ".faa")

    rows = {row["query"]: row for row in flow(run)}

    assert rows["bin_a"]["screen_status"] == "not_installed"
    assert rows["bin_b"]["screen_status"] == "not_applicable"
    assert rows["bin_a"]["screen_cellular_pct"] == "NaN"
    assert rows["bin_a"]["taxonomy_majority"] == "processed_bin_a"
    assert rows["bin_a"]["estimated_completeness"] == "91.25"


def test_a_screen_that_breaks_never_fails_the_run(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1"], ".fna")
    install_addon(run)

    def broken_run(*args: object, **kwargs: object) -> None:
        msg = "cannot start"
        raise OSError(msg)

    monkeypatch.setattr("src.pipeline.contig_screen_step.subprocess.Popen", broken_run)

    row = flow(run)[0]

    assert row["screen_status"] == "failed"
    assert row["taxonomy_majority"] == "processed_q1"


def test_resume_screens_skipped_queries_once_the_addon_is_installed(
    run: RunConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    install_fake_screen(tmp_path, monkeypatch)
    write_queries(run, ["ok_a", "ok_b"], ".fna")
    first = flow(run)
    screen_dir = install_addon(run)
    refuse_to_process(monkeypatch)

    second = flow(run, resume=True)

    assert [row["screen_status"] for row in first] == ["not_installed"] * 2
    assert [row["screen_status"] for row in second] == ["ok", "ok"]
    assert [row["screen_flagged_pct"] for row in second] == ["25.00", "25.00"]
    assert [row["taxonomy_majority"] for row in second] == [
        "processed_ok_a",
        "processed_ok_b",
    ]
    assert sorted(path.name for path in screen_dir.iterdir()) == ["ok_a", "ok_b"]


def test_resume_screens_a_skipped_query_and_a_processed_one_after_the_install(
    run: RunConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    install_fake_screen(tmp_path, monkeypatch)
    write_queries(run, ["ok_a"], ".fna")
    flow(run)
    write_queries(run, ["ok_b"], ".fna")
    install_addon(run)

    rows = flow(run, resume=True)

    assert [row["taxonomy_majority"] for row in rows] == [
        "processed_ok_a",
        "processed_ok_b",
    ]
    assert [row["screen_status"] for row in rows] == ["ok", "ok"]


def test_resume_without_the_screen_keeps_saved_results_and_reports_disabled(
    run: RunConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    install_fake_screen(tmp_path, monkeypatch)
    write_queries(run, ["ok_a"], ".fna")
    screen_dir = install_addon(run)
    flow(run)
    saved = {path.name: path.read_bytes() for path in (screen_dir / "ok_a").iterdir()}
    refuse_to_process(monkeypatch)

    off = flow(run, resume=True, contig_screen=False)

    assert [row["screen_status"] for row in off] == ["disabled"]
    assert off[0]["screen_cellular_pct"] == "NaN"
    assert {
        path.name: path.read_bytes() for path in (screen_dir / "ok_a").iterdir()
    } == saved

    (screen_dir / "ok_a" / "args.json").unlink()
    on = flow(run, resume=True)

    assert [row["screen_status"] for row in on] == ["ok"]
    assert not (screen_dir / "ok_a" / "args.json").exists()


def test_resume_screens_again_after_the_addon_manifest_changed(
    run: RunConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    install_fake_screen(tmp_path, monkeypatch)
    write_queries(run, ["ok_a"], ".fna")
    screen_dir = install_addon(run)
    flow(run)
    (screen_dir / "ok_a" / "args.json").unlink()
    refuse_to_process(monkeypatch)
    (run.database_path / "contig_screen" / "manifest.json").write_text('{"v": 2}')

    rows = flow(run, resume=True)

    assert [row["screen_status"] for row in rows] == ["ok"]
    assert (screen_dir / "ok_a" / "args.json").exists()


def test_an_unreadable_manifest_never_costs_the_summary(
    run: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_queries(run, ["q1"], ".fna")
    write_queries(run, ["q2"], ".faa")

    def refuse(database: Path) -> str:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(layout, "manifest_digest", refuse)

    on = flow(run)
    off = flow(run, contig_screen=False)

    assert [(r["query"], r["screen_status"]) for r in on] == [
        ("q1", "failed"),
        ("q2", "not_applicable"),
    ]
    assert [r["screen_status"] for r in off] == ["disabled", "disabled"]
    assert all(r["taxonomy_majority"].startswith("processed_") for r in on + off)
