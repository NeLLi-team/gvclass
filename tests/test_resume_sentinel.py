"""The completion record, the resume rule and the private work directory.

Coverage:
* ``results.jsonl`` gets one fsynced line per finished query and ``run_status.json``
  is written twice per run.
* ``--resume`` skips a query only when the input, the result settings, the database
  and the software are the ones of its last ``complete`` record.
* The previous archive survives a failed archive step; a leftover work directory is
  deleted at query start; work directories cannot collide with output folders.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import tarfile
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock

import pytest

from src.pipeline import parallel_runner, query_processing_engine, run_status
from src.pipeline.parallel_runner import RunInputs
from src.pipeline.run_status import RESULTS_FILENAME, read_results
from src.utils import ValidationError
from src.utils.error_handling import ProcessingError
from tests.pipeline_helpers import (
    fake_flow,
    fake_setup,
    make_config,
    read_jsonl_lines,
    tar_member_text,
    write_queries,
)

LOGGER = logging.getLogger("test_resume_sentinel")
RECORD_FIELDS = {
    "query",
    "status",
    "error",
    "input_sha256",
    "settings",
    "database_version",
    "software_version",
    "finished_at",
    "archive",
    "summary",
}
EXTENDED_VALUES = {
    "mixed_evidence_contig_count": 3,
    "mixed_evidence_status": "review",
    "contig_length_basis": "exact_dna",
    "protein_estimated_contig_count": 2,
    "contamination_model_applicability": "applicable",
    "contamination_prediction_range_status": "above_training_range",
    "contamination_model_prediction_upper_bound": 0.8765432,
}


def read_rows(summary_file: Path) -> list[dict[str, str]]:
    with summary_file.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# results.jsonl and run_status.json.
# ---------------------------------------------------------------------------


def test_results_log_has_one_line_per_finished_query(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch)
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2", "q3"])

    parallel_runner.gvclass_flow(config)

    lines = read_jsonl_lines(config.output_dir / RESULTS_FILENAME)
    assert sorted(json.loads(line)["query"] for line in lines) == ["q1", "q2", "q3"]


def test_results_record_holds_what_resume_compares(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch, summaries={"q1": EXTENDED_VALUES})
    (tmp_path / "db").mkdir()
    (tmp_path / "db" / "DB_VERSION").write_text("v9.9.9\n")
    config = make_config(tmp_path)
    (query_file,) = write_queries(config, ["q1"])

    parallel_runner.gvclass_flow(config)

    (line,) = read_jsonl_lines(config.output_dir / RESULTS_FILENAME)
    record = json.loads(line)
    archive = config.output_dir / "q1.tar.gz"
    assert set(record) == RECORD_FIELDS
    assert record["status"] == "complete"
    assert record["error"] is None
    assert record["input_sha256"] == sha256_of(query_file)
    assert record["settings"] == config.result_settings()
    assert record["database_version"] == "v9.9.9"
    assert record["software_version"] == run_status.SOFTWARE_VERSION
    assert record["finished_at"].endswith("+00:00")
    assert record["archive"] == {
        "filename": "q1.tar.gz",
        "bytes": archive.stat().st_size,
        "sha256": sha256_of(archive),
    }
    assert record["summary"]["taxonomy_majority"] == "processed_q1"
    assert record["summary"]["estimated_completeness"] == "91.25"
    assert record["summary"]["mixed_evidence_status"] == "review"
    assert record["summary"]["contamination_model_prediction_upper_bound"] == "0.876543"


def test_failed_query_gets_a_failed_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch, failing={"q2": RuntimeError("simulated query failure")})
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2"])

    with pytest.raises(RuntimeError, match="failed for 1 of 2 queries"):
        parallel_runner.gvclass_flow(config)

    records = read_results(config.output_dir)
    assert records["q1"].status == "complete"
    assert records["q2"].status == "failed"
    assert records["q2"].error == "simulated query failure"
    assert records["q2"].archive is None
    assert records["q2"].summary is None


def test_rejected_input_gets_an_invalid_input_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch)
    config = make_config(tmp_path)
    good, bad = write_queries(config, ["good", "bad"])

    def setup_with_rejection(config_: object) -> RunInputs:
        inputs = fake_setup(config)
        return RunInputs(
            (good,), ((bad, "total length is below the minimum"),), inputs.database_path
        )

    monkeypatch.setattr(
        parallel_runner, "validate_and_setup_task", setup_with_rejection
    )

    with pytest.raises(RuntimeError, match="failed for 1 of 2 queries"):
        parallel_runner.gvclass_flow(config)

    records = read_results(config.output_dir)
    assert records["good"].status == "complete"
    assert records["bad"].status == "invalid_input"
    assert records["bad"].error == "total length is below the minimum"
    status = json.loads((config.output_dir / "run_status.json").read_text())
    assert status["queries"]["bad"]["status"] == "invalid_input"
    assert status["status"] == "failed"


def test_run_status_is_written_at_start_and_at_the_end_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch)
    config = make_config(tmp_path, threads=4, max_workers=2, threads_per_worker=2)
    write_queries(config, ["q1", "q2", "q3", "q4"])
    written: list[dict[str, object]] = []
    real_replace = os.replace

    def spy_replace(source: str | Path, target: str | Path) -> None:
        if Path(target).name == "run_status.json":
            written.append(json.loads(Path(source).read_text()))
        real_replace(source, target)

    monkeypatch.setattr(run_status.os, "replace", spy_replace)

    parallel_runner.gvclass_flow(config)

    assert [state["status"] for state in written] == ["running", "complete"]
    final = written[-1]
    assert final["software_version"] == run_status.SOFTWARE_VERSION
    assert final["settings"]["tree_method"] == config.tree_method
    assert {name: q["status"] for name, q in final["queries"].items()} == {
        f"q{i}": "complete" for i in range(1, 5)
    }
    run_log = (config.output_dir / "run.log").read_text()
    assert run_log.count("QUERY_START") == 4
    assert run_log.count("QUERY_COMPLETE") == 4


def test_a_crash_between_queries_leaves_a_readable_results_log(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch, failing={"q3": KeyboardInterrupt()})
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2", "q3"])

    with pytest.raises(KeyboardInterrupt):
        parallel_runner.gvclass_flow(config)

    lines = read_jsonl_lines(config.output_dir / RESULTS_FILENAME)
    assert [json.loads(line)["query"] for line in lines] == ["q1", "q2"]
    assert set(read_results(config.output_dir)) == {"q1", "q2"}
    status = json.loads((config.output_dir / "run_status.json").read_text())
    assert status["status"] == "failed"


def test_a_run_after_a_killed_run_keeps_the_new_records_readable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_flow(monkeypatch)
    config = make_config(tmp_path)
    write_queries(config, ["q1"])
    config.output_dir.mkdir()
    (config.output_dir / RESULTS_FILENAME).write_text('{"query": "old", "sta')

    parallel_runner.gvclass_flow(config)

    assert set(read_results(config.output_dir)) == {"q1"}


def test_read_results_keeps_the_last_record_of_a_query(tmp_path: Path) -> None:
    first = {**dict.fromkeys(RECORD_FIELDS), "query": "q1", "status": "failed"}
    second = {**first, "status": "complete", "error": None}
    (tmp_path / RESULTS_FILENAME).write_text(
        json.dumps(first) + "\n" + json.dumps(second) + "\n", encoding="utf-8"
    )

    assert read_results(tmp_path)["q1"].status == "complete"


def test_read_results_skips_a_torn_last_line(tmp_path: Path) -> None:
    first = {**dict.fromkeys(RECORD_FIELDS), "query": "q1", "status": "complete"}
    (tmp_path / RESULTS_FILENAME).write_text(
        json.dumps(first) + "\n" + '{"query": "q2", "sta', encoding="utf-8"
    )

    assert set(read_results(tmp_path)) == {"q1"}


def test_read_results_without_a_file_is_empty(tmp_path: Path) -> None:
    assert read_results(tmp_path) == {}


# ---------------------------------------------------------------------------
# Resume.
# ---------------------------------------------------------------------------


def run_twice(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    names: list[str],
    *,
    change: object = None,
    **second_overrides: object,
) -> list[str]:
    """Run, apply ``change``, resume; return the queries the second run processed."""
    calls: list[str] = []
    fake_flow(monkeypatch, calls)
    config = make_config(tmp_path)
    write_queries(config, names)
    parallel_runner.gvclass_flow(config)
    calls.clear()
    if callable(change):
        change(config)
    parallel_runner.gvclass_flow(replace(config, resume=True, **second_overrides))
    return calls


def test_resume_skips_an_unchanged_query(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"]) == []


def test_resume_reprocesses_after_the_input_bytes_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def edit_input(config: object) -> None:
        (tmp_path / "queries" / "q1.faa").write_text(">p\nMK\n")

    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"], change=edit_input) == ["q1"]


@pytest.mark.parametrize(
    "override",
    [
        {"tree_method": "iqtree"},
        {"iqtree_mode": "ufboot"},
        {"mode_fast": False},
        {"completeness_mode": "legacy"},
        {"sensitive_mode": False},
        {"species_tree": True},
        {"species_tree_trim": "none"},
    ],
)
def test_resume_reprocesses_after_a_result_setting_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, override: dict[str, object]
) -> None:
    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"], **override) == ["q1", "q2"]


@pytest.mark.parametrize(
    "override",
    [
        {"threads": 8, "max_workers": 2, "threads_per_worker": 4},
        {"contig_screen": False},
        {"contigs_mode": True},
    ],
)
def test_resume_skips_after_a_setting_that_leaves_results_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, override: dict[str, object]
) -> None:
    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"], **override) == []


def test_resume_reprocesses_after_the_database_version_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "db").mkdir()
    (tmp_path / "db" / "DB_VERSION").write_text("v1.0.0\n")

    def update_database(config: object) -> None:
        (tmp_path / "db" / "DB_VERSION").write_text("v2.0.0\n")

    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"], change=update_database) == [
        "q1",
        "q2",
    ]


def test_resume_reprocesses_after_the_software_version_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def upgrade(config: object) -> None:
        monkeypatch.setattr(run_status, "SOFTWARE_VERSION", "v99.0.0")

    assert run_twice(tmp_path, monkeypatch, ["q1", "q2"], change=upgrade) == [
        "q1",
        "q2",
    ]


def test_resume_reprocesses_a_query_whose_archive_is_gone_or_changed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def damage(config: object) -> None:
        (tmp_path / "out" / "q1.tar.gz").unlink()
        (tmp_path / "out" / "q2.tar.gz").write_bytes(b"short")

    assert run_twice(tmp_path, monkeypatch, ["q1", "q2", "q3"], change=damage) == [
        "q1",
        "q2",
    ]


def test_resume_reprocesses_a_query_whose_last_record_is_failed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []
    fake_flow(monkeypatch, calls, failing={"q2": RuntimeError("boom")})
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2"])
    with pytest.raises(RuntimeError, match="failed for 1 of 2 queries"):
        parallel_runner.gvclass_flow(config)
    calls.clear()
    fake_flow(monkeypatch, calls)

    parallel_runner.gvclass_flow(replace(config, resume=True))

    assert calls == ["q2"]
    assert read_results(config.output_dir)["q2"].status == "complete"


def test_resume_logs_the_field_that_differs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    def edit_input(config: object) -> None:
        (tmp_path / "queries" / "q1.faa").write_text(">p\nMK\n")

    with caplog.at_level(logging.INFO, logger="gvclass_runner"):
        run_twice(
            tmp_path,
            monkeypatch,
            ["q1"],
            change=edit_input,
            tree_method="iqtree",
        )

    assert any(
        "q1" in message and "input_sha256" in message and "tree_method" in message
        for message in caplog.messages
    )


def test_resume_summary_keeps_the_saved_rows_with_their_extended_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []
    fake_flow(
        monkeypatch,
        calls,
        summaries={"q1": {**EXTENDED_VALUES, "taxonomy_confidence": "strict_q1"}},
    )
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2"])
    parallel_runner.gvclass_flow(config)
    (config.query_path / "q3.faa").write_text(">p\nM\n")
    calls.clear()

    summary_file = parallel_runner.gvclass_flow(replace(config, resume=True))

    assert calls == ["q3"]
    rows = read_rows(summary_file)
    assert [row["query"] for row in rows] == ["q1", "q2", "q3"]
    assert rows[0]["taxonomy_majority"] == "processed_q1"
    assert rows[0]["taxonomy_confidence"] == "strict_q1"
    assert rows[0]["completeness_query_support"] == "supported"
    assert rows[0]["estimated_completeness"] == "91.25"
    assert rows[2]["taxonomy_majority"] == "processed_q3"
    extended = tar_member_text(
        config.output_dir / "gvclass_summary.extended.tar.gz",
        "gvclass_summary.extended.tsv",
    )
    q1_extended = next(
        row
        for row in csv.DictReader(extended.splitlines(), delimiter="\t")
        if row["query"] == "q1"
    )
    assert q1_extended["mixed_evidence_contig_count"] == "3"
    assert q1_extended["contig_length_basis"] == "exact_dna"
    assert q1_extended["contamination_model_prediction_upper_bound"] == "0.876543"


def test_resume_with_every_query_skipped_rebuilds_the_summary_and_records_skips(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []
    fake_flow(monkeypatch, calls)
    config = make_config(tmp_path)
    write_queries(config, ["q1", "q2"], suffix=".fna")
    parallel_runner.gvclass_flow(config)
    first_lines = read_jsonl_lines(config.output_dir / RESULTS_FILENAME)
    calls.clear()

    summary_file = parallel_runner.gvclass_flow(replace(config, resume=True))

    assert calls == []
    assert [row["taxonomy_majority"] for row in read_rows(summary_file)] == [
        "processed_q1",
        "processed_q2",
    ]
    assert read_jsonl_lines(config.output_dir / RESULTS_FILENAME) == first_lines
    status = json.loads((config.output_dir / "run_status.json").read_text())
    assert status["status"] == "complete"
    assert status["totals"] == {
        "queries": 2,
        "complete": 2,
        "failed": 0,
        "resume_skipped": 2,
    }
    assert all(query["resume_skipped"] for query in status["queries"].values())
    run_log = (config.output_dir / "run.log").read_text()
    assert "QUERY_RESUME_SKIPPED\tq1" in run_log
    assert "QUERY_RESUME_SKIPPED\tq2" in run_log


# ---------------------------------------------------------------------------
# The archive and the private work directory.
# ---------------------------------------------------------------------------


def make_work_directory(root: Path) -> Path:
    work = root / "work" / "q1"
    work.mkdir(parents=True)
    (work / "hello.txt").write_text("hello")
    return work


def test_archive_query_output_is_atomic_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An exception mid-tar publishes no ``.tar.gz`` and leaves no ``.part``."""

    output_base = tmp_path / "out"
    output_base.mkdir()
    work = make_work_directory(tmp_path)

    def exploding_tarfile_open(*args: object, **kwargs: object) -> None:
        raise RuntimeError("simulated tar failure")

    monkeypatch.setattr(tarfile, "open", exploding_tarfile_open)
    with pytest.raises(RuntimeError, match="simulated tar failure"):
        query_processing_engine._archive_query_output("q1", work, output_base, LOGGER)

    assert not (output_base / "q1.tar.gz").exists()
    assert not (output_base / "q1.tar.gz.part").exists()


def test_archive_query_output_publishes_a_valid_archive(tmp_path: Path) -> None:

    output_base = tmp_path / "out"
    output_base.mkdir()
    work = make_work_directory(tmp_path)

    query_processing_engine._archive_query_output("q1", work, output_base, LOGGER)

    assert tarfile.is_tarfile(output_base / "q1.tar.gz")
    assert tar_member_text(output_base / "q1.tar.gz", "q1/hello.txt") == "hello"
    assert not (output_base / "q1.tar.gz.part").exists()


def test_the_previous_archive_survives_a_failed_re_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:

    output_base = tmp_path / "out"
    output_base.mkdir()
    previous = output_base / "q1.tar.gz"
    previous.write_bytes(b"previous archive bytes")
    work = make_work_directory(tmp_path)

    def exploding_tarfile_open(*args: object, **kwargs: object) -> None:
        raise RuntimeError("simulated tar failure")

    monkeypatch.setattr(tarfile, "open", exploding_tarfile_open)
    with pytest.raises(RuntimeError, match="Post-processing failed for q1"):
        query_processing_engine._post_process_query("q1", work, output_base, LOGGER)

    assert previous.read_bytes() == b"previous archive bytes"
    assert work.is_dir()


def test_post_processing_replaces_the_previous_archive_and_removes_the_work_directory(
    tmp_path: Path,
) -> None:

    output_base = tmp_path / "out"
    output_base.mkdir()
    (output_base / "q1.tar.gz").write_bytes(b"previous archive bytes")
    work = make_work_directory(tmp_path)

    query_processing_engine._post_process_query("q1", work, output_base, LOGGER)

    assert tar_member_text(output_base / "q1.tar.gz", "q1/hello.txt") == "hello"
    assert not work.exists()


def test_initialize_query_gives_a_private_work_directory(tmp_path: Path) -> None:

    query = tmp_path / "sample.faa"
    query.write_text(">p\nM\n")

    name, work = query_processing_engine._initialize_query(query, tmp_path / "out")

    assert name == "sample"
    assert work == tmp_path / "out" / ".gvclass_work" / "sample"
    assert work.is_dir()


def test_initialize_query_deletes_a_leftover_work_directory(tmp_path: Path) -> None:

    query = tmp_path / "sample.faa"
    query.write_text(">p\nM\n")
    leftover = tmp_path / "out" / ".gvclass_work" / "sample"
    (leftover / "hmmout").mkdir(parents=True)
    (leftover / "hmmout" / "models.out.filtered").write_text("stale")

    _, work = query_processing_engine._initialize_query(query, tmp_path / "out")

    assert list(work.iterdir()) == []
    assert not list((tmp_path / "out").glob("sample.previous.*"))


def test_work_directories_do_not_collide_with_output_folders_or_other_queries(
    tmp_path: Path,
) -> None:

    output = tmp_path / "out"
    screen = output / "contig_screen" / "bin"
    screen.mkdir(parents=True)
    (screen / "keep.txt").write_text("screen result")
    names = ["contig_screen", "sample", "sample.previous.abc123"]
    works = []
    for name in names:
        query = tmp_path / f"{name}.faa"
        query.write_text(">p\nM\n")
        works.append(query_processing_engine._initialize_query(query, output)[1])

    assert len(set(works)) == 3
    assert all(work.parent == output / ".gvclass_work" for work in works)
    assert (screen / "keep.txt").read_text() == "screen result"


def test_initialize_query_refuses_a_work_directory_that_holds_the_input(
    tmp_path: Path,
) -> None:

    work = tmp_path / "out" / ".gvclass_work" / "sample"
    work.mkdir(parents=True)
    query = work / "sample.faa"
    query.write_text(">p\nM\n")

    with pytest.raises(ProcessingError, match="holds the input file"):
        query_processing_engine._initialize_query(query, tmp_path / "out")

    assert query.exists()


# ---------------------------------------------------------------------------
# Atomic DB install.
# ---------------------------------------------------------------------------


def test_atomic_swap_directory_rolls_back_on_failure(
    tmp_path: Path, monkeypatch
) -> None:
    from src.utils.database_manager import DatabaseManager

    db_path = tmp_path / "resources"
    staging_path = tmp_path / "resources.new"
    db_path.mkdir()
    (db_path / "existing.txt").write_text("old")
    staging_path.mkdir()
    (staging_path / "new.txt").write_text("new")

    # Force the second os.replace (staging -> db_path) to fail.
    import os as os_module

    real_replace = os_module.replace
    call_state = {"count": 0}

    def flaky_replace(src, dst):
        call_state["count"] += 1
        if call_state["count"] == 2:
            raise OSError("simulated rename failure")
        real_replace(src, dst)

    monkeypatch.setattr(os_module, "replace", flaky_replace)
    with pytest.raises(OSError, match="simulated rename failure"):
        DatabaseManager._atomic_swap_directory(staging_path, db_path)

    # Original db_path must have been restored.
    assert db_path.exists()
    assert (db_path / "existing.txt").exists()
    assert (db_path / "existing.txt").read_text() == "old"


def test_atomic_swap_directory_preserves_backup_when_restore_fails(
    tmp_path: Path, monkeypatch
) -> None:
    from src.utils.database_manager import DatabaseManager

    db_path = tmp_path / "resources"
    staging_path = tmp_path / "resources.new"
    db_path.mkdir()
    (db_path / "existing.txt").write_text("old")
    staging_path.mkdir()
    (staging_path / "new.txt").write_text("new")

    import os as os_module

    real_replace = os_module.replace
    call_state = {"count": 0}

    def double_fault_replace(src, dst):
        call_state["count"] += 1
        if call_state["count"] in {2, 3}:
            raise OSError(f"simulated rename failure {call_state['count']}")
        real_replace(src, dst)

    monkeypatch.setattr(os_module, "replace", double_fault_replace)
    with pytest.raises(OSError, match="simulated rename failure 2"):
        DatabaseManager._atomic_swap_directory(staging_path, db_path)

    assert not db_path.exists()
    (backup_path,) = tmp_path.glob(".resources.old.*")
    assert (backup_path / "existing.txt").read_text() == "old"
    assert (staging_path / "new.txt").read_text() == "new"


def test_atomic_swap_directory_replaces_previous_and_cleans_backup(
    tmp_path: Path,
) -> None:
    from src.utils.database_manager import DatabaseManager

    db_path = tmp_path / "resources"
    staging_path = tmp_path / "resources.new"
    db_path.mkdir()
    (db_path / "existing.txt").write_text("old")
    staging_path.mkdir()
    (staging_path / "new.txt").write_text("new")

    DatabaseManager._atomic_swap_directory(staging_path, db_path)

    assert db_path.exists()
    assert (db_path / "new.txt").read_text() == "new"
    assert not (db_path / "existing.txt").exists()
    assert not list(tmp_path.glob("*.old*"))
    assert not staging_path.exists()


# ---------------------------------------------------------------------------
# Input check.
# ---------------------------------------------------------------------------


def test_validate_and_setup_rejects_duplicate_stems_before_side_effects(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = make_config(tmp_path)
    config.query_path.mkdir()
    (config.query_path / "same.fna").write_text(">c\nACGT\n")
    (config.query_path / "same.faa").write_text(">p\nM\n")
    setup_database = Mock()
    monkeypatch.setattr(
        parallel_runner.DatabaseManager, "setup_database", setup_database
    )

    with pytest.raises(ValidationError, match=r"same: same\.faa, same\.fna"):
        parallel_runner.validate_and_setup_task(config)

    setup_database.assert_not_called()
    assert not config.output_dir.exists()
