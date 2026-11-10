"""The m8 file keeps pyswrd's identity and reports BLAST-style length and gaps."""

from pathlib import Path

import pytest

from src.core.blast import iter_blastp_records, run_blastp

QUERY = "MAGFLKVVQLLAKYGSKAVQWAWANKGKILDWLNAGQAIDWVVSKIKQILGIK"


def _search_one_target(tmp_path: Path, target: str) -> list[str]:
    query_faa = tmp_path / "query.faa"
    reference_faa = tmp_path / "reference.faa"
    m8 = tmp_path / "hits.m8"
    query_faa.write_text(f">query\n{QUERY}\n", encoding="utf-8")
    reference_faa.write_text(f">target\n{target}\n", encoding="utf-8")

    run_blastp(str(query_faa), str(reference_faa), str(m8), threads=1)

    lines = m8.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    return lines[0].split("\t")


def test_m8_length_includes_gap_columns_and_one_gap_opening_per_run(
    tmp_path: Path,
) -> None:
    target = QUERY[:20] + "GGG" + QUERY[20:]

    fields = _search_one_target(tmp_path, target)

    # query, target, identity, length, mismatches, gap openings, coordinates
    assert fields[:2] == ["query", "target"]
    assert fields[2:10] == ["100.0", "56", "0", "1", "1", "53", "1", "56"]


def test_m8_counts_each_separate_gap_as_one_opening(tmp_path: Path) -> None:
    target = QUERY[:12] + "PP" + QUERY[12:35] + "W" + QUERY[35:]

    fields = _search_one_target(tmp_path, target)

    assert fields[3] == "56"
    assert fields[5] == "2"


def test_m8_counts_mismatches_without_gaps(tmp_path: Path) -> None:
    target = QUERY[:10] + "A" + QUERY[11:]

    fields = _search_one_target(tmp_path, target)

    assert fields[2:6] == ["98.1", "53", "1", "0"]


def test_m8_identity_is_matches_over_aligned_pairs_not_over_alignment_length(
    tmp_path: Path,
) -> None:
    """Three gap columns leave identity at 100.0, which BLAST would put at 94.6."""
    target = QUERY[:20] + "GGG" + QUERY[20:]
    _search_one_target(tmp_path, target)

    (record,) = iter_blastp_records(tmp_path / "hits.m8")

    assert record.identity == pytest.approx(100.0)
