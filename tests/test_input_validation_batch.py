"""One bad input file must not stop the batch; it is reported as a failed query."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from src.pipeline import parallel_runner
from src.pipeline.parallel_runner import validate_and_setup_task
from src.utils.error_handling import ValidationError
from src.utils.input_validation import RESERVED_QUERY_STEMS
from tests.test_resume_sentinel import _fake_query_processor, _fake_resume_validate

LONG_DNA = "ACGT" * 6000


@pytest.fixture
def query_dir(tmp_path: Path) -> Path:
    path = tmp_path / "queries"
    path.mkdir()
    return path


@pytest.fixture(autouse=True)
def database_setup(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Mock:
    setup_database = Mock(return_value=tmp_path / "db")
    monkeypatch.setattr(
        parallel_runner.DatabaseManager, "setup_database", setup_database
    )
    return setup_database


def _setup(query_dir: Path, tmp_path: Path, **kwargs: object) -> dict:
    return validate_and_setup_task(str(query_dir), str(tmp_path / "out"), **kwargs)


def test_setup_keeps_valid_files_and_rejects_invalid_ones(
    query_dir: Path, tmp_path: Path
) -> None:
    (query_dir / "good.fna").write_text(f">c\n{LONG_DNA}\n")
    (query_dir / "short.fna").write_text(">c\nACGT\n")
    (query_dir / "dup.faa").write_text(">p\nMKV\n>p\nMKV\n")

    config = _setup(query_dir, tmp_path)

    assert [path.name for path in config["query_files"]] == ["good.fna"]
    rejected = dict(config["rejected_inputs"])
    assert sorted(rejected) == ["dup", "short"]
    assert "minimum" in rejected["short"]
    assert "'p'" in rejected["dup"]


def test_setup_finds_an_uppercase_extension(query_dir: Path, tmp_path: Path) -> None:
    (query_dir / "X.FNA").write_text(f">c\n{LONG_DNA}\n")

    config = _setup(query_dir, tmp_path)

    assert [path.name for path in config["query_files"]] == ["X.FNA"]
    assert config["rejected_inputs"] == []


def test_setup_allow_short_still_accepts_short_files(
    query_dir: Path, tmp_path: Path
) -> None:
    (query_dir / "short.fna").write_text(">c\nACGT\n")

    config = _setup(query_dir, tmp_path, allow_short=True)

    assert [path.name for path in config["query_files"]] == ["short.fna"]
    assert config["rejected_inputs"] == []


def test_setup_without_any_valid_file_is_an_error(
    query_dir: Path, tmp_path: Path
) -> None:
    (query_dir / "short.fna").write_text(">c\nACGT\n")

    with pytest.raises(ValidationError, match="short"):
        _setup(query_dir, tmp_path)


def test_setup_without_query_files_is_an_error(query_dir: Path, tmp_path: Path) -> None:
    (query_dir / "notes.txt").write_text("x")

    with pytest.raises(ValidationError, match="No valid sequence files"):
        _setup(query_dir, tmp_path)


@pytest.mark.parametrize("stem", sorted(RESERVED_QUERY_STEMS))
def test_setup_rejects_a_reserved_query_name_for_the_whole_run(
    query_dir: Path, tmp_path: Path, database_setup: Mock, stem: str
) -> None:
    (query_dir / "good.fna").write_text(f">c\n{LONG_DNA}\n")
    (query_dir / f"{stem}.fna").write_text(f">c\n{LONG_DNA}\n")

    with pytest.raises(ValidationError, match=f"reserved.*{stem}"):
        _setup(query_dir, tmp_path)

    database_setup.assert_not_called()


def test_reserved_names_cover_the_pipeline_outputs() -> None:
    assert {"gvclass_summary", "gvclass_summary.extended", "_combined"} <= (
        RESERVED_QUERY_STEMS
    )


def test_flow_reports_a_rejected_input_as_a_failed_query(
    query_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "out"
    output.mkdir()
    (query_dir / "good.faa").write_text(">p\nM\n")
    fake_setup = _fake_resume_validate(query_dir, output)

    def setup_with_rejection(*args: object, **kwargs: object) -> dict:
        config = fake_setup(*args, **kwargs)
        config["rejected_inputs"] = [("bad", "Total sequence length is too short")]
        return config

    monkeypatch.setattr(
        parallel_runner, "validate_and_setup_task", setup_with_rejection
    )
    monkeypatch.setattr(parallel_runner, "process_query_task", _fake_query_processor)

    with pytest.raises(RuntimeError, match="failed for 1 of 2 queries"):
        parallel_runner.gvclass_flow(
            str(query_dir),
            str(output),
            total_threads=1,
            max_workers=1,
            threads_per_worker=1,
        )

    with (output / "gvclass_failed_queries.tsv").open(newline="") as handle:
        failed = list(csv.DictReader(handle, delimiter="\t"))
    assert failed == [
        {
            "query": "bad",
            "status": "invalid_input",
            "error": "Total sequence length is too short",
        }
    ]
    with (output / "gvclass_summary.tsv").open(newline="") as handle:
        assert [row["query"] for row in csv.DictReader(handle, delimiter="\t")] == [
            "good"
        ]
    status = json.loads((output / "run_status.json").read_text())
    assert status["status"] == "failed"
    assert status["queries"]["bad"]["status"] == "failed"
    assert status["queries"]["bad"]["error"] == "Total sequence length is too short"
