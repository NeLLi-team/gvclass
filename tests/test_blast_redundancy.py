"""Each marker is searched once and each input file is read once."""

import builtins
import io
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from src.core import blast, marker_processing
from src.core.blast import iter_blastp_records
from src.core.marker_processing import MarkerProcessor

QUERY_FAA = ">q1 first protein\nMKVLAAMKVL\n>q2\nMKVLAAMKVL\n"
REFERENCE_FAA = (
    ">r1 reference one\nAAAA\n"
    ">r2 reference two\nBBBB\n"
    ">r3 reference three\nCCCC\n"
    ">r4 reference four\nDDDD\n"
    ">r5 reference five\nEEEE\n"
)


def _search_hit(query_index: int, target_index: int, score: float) -> SimpleNamespace:
    return SimpleNamespace(
        query_index=query_index,
        target_index=target_index,
        evalue=1e-30,
        score=score,
        result=SimpleNamespace(
            identity=lambda: 0.9,
            query_start=0,
            query_end=9,
            target_start=0,
            target_end=9,
            alignment="M" * 10,
        ),
    )


@pytest.fixture
def opened_paths(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record the path of every file opened for reading."""
    paths: list[str] = []
    real_open = io.open

    def recording_open(file: object, *args: object, **kwargs: object) -> object:
        mode = str(args[0] if args else kwargs.get("mode", "r"))
        if isinstance(file, str | os.PathLike) and mode.startswith("r"):
            paths.append(os.fspath(file))
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", recording_open)
    monkeypatch.setattr(io, "open", recording_open)
    return paths


@pytest.fixture
def search_calls(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    """Replace pyswrd with a fake: q1 hits r1 to r4 and q2 hits only r5."""
    calls: list[dict[str, object]] = []

    def fake_search(
        queries: list[str], targets: list[str], **kwargs: object
    ) -> list[SimpleNamespace]:
        calls.append({"queries": queries, "targets": targets, **kwargs})
        return [
            _search_hit(0, 0, 90),
            _search_hit(0, 1, 80),
            _search_hit(0, 2, 70),
            _search_hit(0, 3, 60),
            _search_hit(1, 4, 10),
        ]

    monkeypatch.setattr(blast.pyswrd, "search", fake_search)
    return calls


@pytest.fixture
def marker_files(tmp_path: Path) -> tuple[MarkerProcessor, Path]:
    reference = tmp_path / "db" / "database" / "faa" / "OG1.faa"
    reference.parent.mkdir(parents=True)
    reference.write_text(REFERENCE_FAA, encoding="utf-8")
    query = tmp_path / "OG1.query.faa"
    query.write_text(QUERY_FAA, encoding="utf-8")
    processor = MarkerProcessor("OG1", tmp_path / "db", tmp_path / "out")
    return processor, query


def test_blast_and_merge_keeps_the_second_proteins_neighbour_in_the_merge(
    marker_files: tuple[MarkerProcessor, Path],
    search_calls: list[dict[str, object]],
) -> None:
    processor, query = marker_files

    merged = processor.blast_and_merge(query, max_hits=3, threads=1)

    assert merged.read_text(encoding="utf-8") == (
        ">q1 first protein\nMKVLAAMKVL\n"
        ">q2\nMKVLAAMKVL\n"
        ">r1 reference one\nAAAA\n"
        ">r5 reference five\nEEEE\n"
        ">r2 reference two\nBBBB\n"
    )


def test_blast_and_merge_searches_once_with_the_hit_cap(
    marker_files: tuple[MarkerProcessor, Path],
    search_calls: list[dict[str, object]],
) -> None:
    processor, query = marker_files

    processor.blast_and_merge(query, max_hits=3, threads=1)

    assert len(search_calls) == 1
    assert search_calls[0]["max_alignments"] == 3
    assert len(search_calls[0]["targets"]) == 5


def test_blast_and_merge_reads_the_reference_file_once(
    marker_files: tuple[MarkerProcessor, Path],
    search_calls: list[dict[str, object]],
    opened_paths: list[str],
    tmp_path: Path,
) -> None:
    processor, query = marker_files
    reference = tmp_path / "db" / "database" / "faa" / "OG1.faa"

    processor.blast_and_merge(query, max_hits=3, threads=1)

    assert opened_paths.count(str(reference)) == 1


def test_blast_and_merge_writes_the_m8_file_that_scoring_reads(
    marker_files: tuple[MarkerProcessor, Path],
    search_calls: list[dict[str, object]],
) -> None:
    processor, query = marker_files

    processor.blast_and_merge(query, max_hits=3, threads=1)

    records = iter_blastp_records(processor.blast_dir / "OG1.m8")
    assert [(record.query_id, record.subject_id) for record in records] == [
        ("q1", "r1"),
        ("q1", "r2"),
        ("q1", "r3"),
        ("q2", "r5"),
    ]


def test_blast_and_merge_rejects_a_missing_reference_set(
    tmp_path: Path,
    search_calls: list[dict[str, object]],
) -> None:
    query = tmp_path / "OG404.query.faa"
    query.write_text(QUERY_FAA, encoding="utf-8")
    processor = MarkerProcessor("OG404", tmp_path / "db", tmp_path / "out")

    with pytest.raises(FileNotFoundError, match="OG404"):
        processor.blast_and_merge(query, threads=1)

    assert search_calls == []
    assert not (processor.merged_dir / "OG404.faa").exists()


def test_build_tree_reads_the_alignment_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    opened_paths: list[str],
) -> None:
    processor = MarkerProcessor("OG1", tmp_path / "db", tmp_path / "out")
    alignment = tmp_path / "OG1.mafft01"
    alignment.write_text(">a\nMKT-\n>b\nMKTA\n>c\nMKSA\n", encoding="utf-8")
    monkeypatch.setattr(
        marker_processing, "run_veryfasttree", lambda _path, _threads: "(a,b,c);"
    )

    processor.build_tree(alignment, threads=1)

    assert opened_paths.count(str(alignment)) == 1
