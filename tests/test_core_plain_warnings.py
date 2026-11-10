"""Keep low-level search warning text readable without decorative icons."""

import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from src.core import blast
from src.pipeline import query_processing_engine


@pytest.mark.parametrize(
    ("alignment_length", "warning"),
    [
        (0, "No hits found despite having query and reference sequences."),
        (6, "Suspicious patterns detected in BLAST results"),
        (12, "1/1 (100.0%) hits show perfect identity!"),
    ],
)
def test_blast_warnings_are_plain(
    alignment_length: int,
    warning: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Exercise warning branches using valid FASTA input and controlled search hits."""
    query = tmp_path / "query.faa"
    reference = tmp_path / "reference.faa"
    query.write_text(">query\nMKVLAAAMKVLA\n", encoding="utf-8")
    reference.write_text(">reference\nMKVLAAAMKVLA\n", encoding="utf-8")
    hit = SimpleNamespace(
        query_index=0,
        target_index=0,
        evalue=1e-10,
        score=100,
        result=SimpleNamespace(
            query_start=0,
            target_start=0,
            query_end=alignment_length - 1,
            target_end=alignment_length - 1,
            alignment="M" * alignment_length,
            identity=lambda: 1.0,
        ),
    )
    hits = [hit] if alignment_length else []
    monkeypatch.setattr(blast.pyswrd, "search", Mock(return_value=iter(hits)))

    with caplog.at_level(logging.WARNING, logger=blast.logger.name):
        blast.run_blastp(str(query), str(reference), str(tmp_path / "hits.m8"))

    assert any(warning in message for message in caplog.messages)
    assert "⚠" not in caplog.text
    assert "\x1b[" not in caplog.text


def test_post_processing_success_logs_are_plain(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Successful cleanup and archive messages retain text without checkmarks."""
    query = tmp_path / "query"
    query.mkdir()
    (query / "query_reformatted.fna").write_text(">q\nACGT\n", encoding="utf-8")
    (query / "result.txt").write_text("result", encoding="utf-8")
    logger = logging.getLogger("plain_post_processing")

    with caplog.at_level(logging.INFO, logger=logger.name):
        query_processing_engine._post_process_query("query", query, tmp_path, logger)

    assert (tmp_path / "query.tar.gz").is_file()
    assert not query.exists()
    assert "Removed reformatted file:" in caplog.text
    assert "Created archive:" in caplog.text
    assert "Removed original directory:" in caplog.text
    assert "Post-processing completed in" in caplog.text
    assert "✓" not in caplog.text


def test_post_processing_failure_log_is_plain(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed archive retains its ERROR diagnostic and raised exception."""
    query = tmp_path / "query"
    query.mkdir()
    logger = logging.getLogger("plain_post_processing")
    monkeypatch.setattr(
        query_processing_engine,
        "_archive_query_output",
        Mock(side_effect=OSError("archive unavailable")),
    )

    with (
        caplog.at_level(logging.ERROR, logger=logger.name),
        pytest.raises(RuntimeError, match="Post-processing failed for query"),
    ):
        query_processing_engine._post_process_query("query", query, tmp_path, logger)

    assert query.exists()
    assert caplog.records[0].levelno == logging.ERROR
    assert caplog.records[0].getMessage().startswith("Post-processing failed for query")
    assert "archive unavailable" in caplog.text
    assert "❌" not in caplog.text
