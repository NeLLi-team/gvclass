"""The run flow screens every query after the queries finish, before the summary."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest

from src.contig_screen import layout
from src.pipeline import parallel_runner
from src.pipeline.contig_screen_step import ScreenOptions
from tests.fake_screen import install_fake_screen
from tests.test_resume_sentinel import (
    _fake_query_processor,
    _fake_resume_validate,
    _seed_resume_complete_query,
)


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
def run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    query_dir = tmp_path / "queries"
    output = tmp_path / "out"
    query_dir.mkdir()
    output.mkdir()
    monkeypatch.setattr(
        parallel_runner,
        "validate_and_setup_task",
        _fake_resume_validate(query_dir, output),
    )
    monkeypatch.setattr(parallel_runner, "process_query_task", _fake_query_processor)
    return query_dir, output


def read_rows(summary_file: Path) -> list[dict[str, str]]:
    with summary_file.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_every_query_is_screened_before_the_final_summary_is_written(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    for name in ("q1", "q2"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
    spy = StepSpy(output)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    summary_file = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8
    )

    assert sorted(spy.queries) == ["q1", "q2"]
    assert not spy.summary_existed
    assert spy.database == output / "db"
    rows = read_rows(summary_file)
    assert [row["screen_status"] for row in rows] == ["ok", "ok"]
    assert [row["screen_flagged_pct"] for row in rows] == ["12.35", "12.35"]


def test_queries_skipped_by_resume_are_screened_too(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    for name in ("q1", "q2", "q3"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
    _seed_resume_complete_query(output, "q1", "skipped_q1")
    spy = StepSpy(output)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    summary_file = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8, resume=True
    )

    assert sorted(spy.queries) == ["q1", "q2", "q3"]
    assert sorted(spy.query_files) == ["q1", "q2", "q3"]
    assert spy.options is not None
    assert spy.options.resume is True
    assert [row["screen_status"] for row in read_rows(summary_file)] == ["ok"] * 3


def test_a_run_in_which_resume_skips_every_query_is_screened(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    for name in ("q1", "q2"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
        _seed_resume_complete_query(output, name, f"skipped_{name}")
    spy = StepSpy(output)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    summary_file = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8, resume=True
    )

    assert sorted(spy.queries) == ["q1", "q2"]
    assert not spy.summary_existed
    assert [row["screen_status"] for row in read_rows(summary_file)] == ["ok", "ok"]


def test_the_step_gets_the_screen_options_and_the_worker_plan_of_the_run(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    for name in ("q1", "q2"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
    spy = StepSpy(output)
    monkeypatch.setattr(parallel_runner, "add_screen_columns", spy)

    parallel_runner.gvclass_flow(
        str(query_dir),
        str(output),
        total_threads=8,
        contig_screen=False,
        contigs_mode=True,
    )

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
    run: tuple[Path, Path],
) -> None:
    query_dir, output = run
    (query_dir / "bin_a.fna").write_text(">c\nACGT\n")
    (query_dir / "bin_b.faa").write_text(">p\nM\n")

    summary_file = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8
    )

    rows = {row["query"]: row for row in read_rows(summary_file)}
    assert rows["bin_a"]["screen_status"] == "not_installed"
    assert rows["bin_b"]["screen_status"] == "not_applicable"
    assert rows["bin_a"]["screen_cellular_pct"] == "NaN"
    assert rows["bin_a"]["taxonomy_majority"] == "processed_bin_a"
    assert rows["bin_a"]["estimated_completeness"] == "91.25"


def test_a_screen_that_breaks_never_fails_the_run(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    (query_dir / "q1.fna").write_text(">c\nACGT\n")
    (output / "db" / "contig_screen").mkdir(parents=True)
    (output / "db" / "contig_screen" / "manifest.json").write_text("{}")

    def broken_run(*args: object, **kwargs: object) -> None:
        msg = "cannot start"
        raise OSError(msg)

    monkeypatch.setattr("src.pipeline.contig_screen_step.subprocess.Popen", broken_run)

    summary_file = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8
    )

    row = read_rows(summary_file)[0]
    assert row["screen_status"] == "failed"
    assert row["taxonomy_majority"] == "processed_q1"


def install_addon(output: Path) -> Path:
    """Put an add-on manifest where the flow looks for it; return the screen folder."""
    manifest = output / "db" / "contig_screen" / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("{}", encoding="utf-8")
    return output / "contig_screen"


def flow(query_dir: Path, output: Path, **options: Any) -> list[dict[str, str]]:
    summary = parallel_runner.gvclass_flow(
        str(query_dir), str(output), total_threads=8, **options
    )
    return read_rows(summary)


def complete_queries(output: Path, names: list[str]) -> None:
    """Leave the files that make resume skip the queries.

    The fake query processor does not record completion in ``run_status.json``, so
    the manifest of the earlier run goes; resume then reads the ``.SUCCESS`` files.
    """
    (output / "run_status.json").unlink(missing_ok=True)
    for name in names:
        _seed_resume_complete_query(output, name, f"skipped_{name}")


def refuse_to_process(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: object, **kwargs: object) -> None:
        msg = "resume must not process a completed query"
        raise AssertionError(msg)

    monkeypatch.setattr(parallel_runner, "process_query_task", refuse)


def test_resume_screens_skipped_queries_once_the_addon_is_installed(
    run: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    install_fake_screen(tmp_path, monkeypatch)
    for name in ("ok_a", "ok_b"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
    first = flow(query_dir, output)
    screen_dir = install_addon(output)
    complete_queries(output, ["ok_a", "ok_b"])
    refuse_to_process(monkeypatch)

    second = flow(query_dir, output, resume=True)

    assert [row["screen_status"] for row in first] == ["not_installed"] * 2
    assert [row["screen_status"] for row in second] == ["ok", "ok"]
    assert [row["screen_flagged_pct"] for row in second] == ["25.00", "25.00"]
    assert [row["taxonomy_majority"] for row in second] == [
        "skipped_ok_a",
        "skipped_ok_b",
    ]
    assert sorted(path.name for path in screen_dir.iterdir()) == ["ok_a", "ok_b"]


def test_resume_screens_a_skipped_query_and_a_processed_one_after_the_install(
    run: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    install_fake_screen(tmp_path, monkeypatch)
    for name in ("ok_a", "ok_b"):
        (query_dir / f"{name}.fna").write_text(">c\nACGT\n")
    install_addon(output)
    complete_queries(output, ["ok_a"])

    rows = flow(query_dir, output, resume=True)

    assert [row["taxonomy_majority"] for row in rows] == [
        "skipped_ok_a",
        "processed_ok_b",
    ]
    assert [row["screen_status"] for row in rows] == ["ok", "ok"]


def test_resume_without_the_screen_keeps_saved_results_and_reports_disabled(
    run: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    install_fake_screen(tmp_path, monkeypatch)
    (query_dir / "ok_a.fna").write_text(">c\nACGT\n")
    screen_dir = install_addon(output)
    flow(query_dir, output)
    saved = {path.name: path.read_bytes() for path in (screen_dir / "ok_a").iterdir()}
    complete_queries(output, ["ok_a"])
    refuse_to_process(monkeypatch)

    off = flow(query_dir, output, resume=True, contig_screen=False)

    assert [row["screen_status"] for row in off] == ["disabled"]
    assert off[0]["screen_cellular_pct"] == "NaN"
    assert {
        path.name: path.read_bytes() for path in (screen_dir / "ok_a").iterdir()
    } == saved

    (screen_dir / "ok_a" / "args.json").unlink()
    complete_queries(output, ["ok_a"])
    on = flow(query_dir, output, resume=True)

    assert [row["screen_status"] for row in on] == ["ok"]
    assert not (screen_dir / "ok_a" / "args.json").exists()


def test_resume_screens_again_after_the_addon_manifest_changed(
    run: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    install_fake_screen(tmp_path, monkeypatch)
    (query_dir / "ok_a.fna").write_text(">c\nACGT\n")
    screen_dir = install_addon(output)
    flow(query_dir, output)
    (screen_dir / "ok_a" / "args.json").unlink()
    complete_queries(output, ["ok_a"])
    refuse_to_process(monkeypatch)
    (output / "db" / "contig_screen" / "manifest.json").write_text('{"v": 2}')

    rows = flow(query_dir, output, resume=True)

    assert [row["screen_status"] for row in rows] == ["ok"]
    assert (screen_dir / "ok_a" / "args.json").exists()


def test_an_unreadable_manifest_never_costs_the_summary(
    run: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir, output = run
    (query_dir / "q1.fna").write_text(">c\nACGT\n")
    (query_dir / "q2.faa").write_text(">p\nM\n")

    def refuse(database: Path) -> str:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(layout, "manifest_digest", refuse)

    on = flow(query_dir, output)
    off = flow(query_dir, output, contig_screen=False)

    assert [(r["query"], r["screen_status"]) for r in on] == [
        ("q1", "failed"),
        ("q2", "not_applicable"),
    ]
    assert [r["screen_status"] for r in off] == ["disabled", "disabled"]
    assert all(r["taxonomy_majority"].startswith("processed_") for r in on + off)
