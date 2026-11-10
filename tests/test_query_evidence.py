"""Query evidence is parsed once per query; candidate types follow the reasons."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from src.core import marker_extraction
from src.core.summarize import paths_stats
from src.core.summarize.paths_stats import QueryEvidence, load_query_evidence
from src.core.summarize_full import FullSummarizer
from tests.conftest import RESOURCES_ROOT, skip_if_no_runtime_resources

TABLE_HEADER = "# target name\taccession\ttlen\tquery name\n"


def _write_hmmout(query_dir: Path, rows: list[str], counts: str = "") -> None:
    hmmout = query_dir / "hmmout"
    hmmout.mkdir(parents=True)
    (hmmout / "models.out.filtered").write_text(
        TABLE_HEADER + "".join(f"{row}\n" for row in rows), encoding="utf-8"
    )
    (hmmout / "models.counts").write_text(counts, encoding="utf-8")


def _row(protein: str, model: str, score: float) -> str:
    return f"{protein}\t-\t100\t{model}\t-\t50\t1e-20\t{score}\t0.0"


def test_load_query_evidence_reads_counts_hits_and_scores(tmp_path: Path) -> None:
    _write_hmmout(
        tmp_path,
        [_row("p1", "M1", 120.5), _row("p1", "M2", 80.0), _row("p2", "M1", 95.0)],
        counts="M1\t2\nM2\t1\n",
    )

    evidence = load_query_evidence(tmp_path)

    assert evidence.counts == {"M1": 2, "M2": 1}
    assert evidence.hits == {"M1": {"p1", "p2"}, "M2": {"p1"}}
    assert evidence.scores == {"p1": {"M1": 120.5, "M2": 80.0}, "p2": {"M1": 95.0}}


def test_load_query_evidence_keeps_the_best_score_of_a_protein_and_model(
    tmp_path: Path,
) -> None:
    _write_hmmout(tmp_path, [_row("p1", "M1", 50.0), _row("p1", "M1", 70.0)])

    assert load_query_evidence(tmp_path).scores == {"p1": {"M1": 70.0}}


def test_load_query_evidence_parses_the_filtered_table_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_hmmout(tmp_path, [_row("p1", "M1", 120.5)], counts="M1\t1\n")
    reads: list[Path] = []
    parse_scores = paths_stats.parse_hmm_scores

    def counting_parse(table: Path) -> dict[str, dict[str, float]]:
        reads.append(table)
        return parse_scores(table)

    monkeypatch.setattr(paths_stats, "parse_hmm_scores", counting_parse)

    load_query_evidence(tmp_path)

    assert reads == [tmp_path / "hmmout" / "models.out.filtered"]


def test_query_without_marker_tables_has_empty_evidence(tmp_path: Path) -> None:
    assert load_query_evidence(tmp_path) == QueryEvidence({}, {}, {})


def test_query_evidence_is_frozen() -> None:
    evidence = QueryEvidence({}, {}, {})

    with pytest.raises(FrozenInstanceError):
        evidence.counts = {"M1": 1}  # type: ignore[misc]  # frozen dataclass


@pytest.mark.parametrize(
    ("reason", "candidate_type"),
    [
        ("cellular_coherent_lineage", "cellular"),
        ("cellular_markers", "cellular"),
        ("cellular_hits", "cellular"),
        ("phage_markers", "phage"),
        ("phage_hits", "phage"),
        ("viral_mixture", "mixed_viral"),
        ("something_new", "uncertain"),
        ("", "uncertain"),
    ],
)
def test_candidate_type_follows_the_reason(reason: str, candidate_type: str) -> None:
    assert paths_stats._candidate_type_for_reason(reason) == candidate_type


@skip_if_no_runtime_resources()
def test_summarize_query_full_parses_the_marker_tables_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query_dir = tmp_path / "query"
    (query_dir / "query_faa").mkdir(parents=True)
    (query_dir / "query_faa" / "q.faa").write_text(">contig_1_1\nMKLVA\n")
    _write_hmmout(query_dir, [_row("contig_1_1", "COG0060", 150.0)], "COG0060\t1\n")
    scores_reads: list[Path] = []
    parse_scores = paths_stats.parse_hmm_scores

    def counting_parse(table: Path) -> dict[str, dict[str, float]]:
        scores_reads.append(table)
        return parse_scores(table)

    def forbidden_parse(_table: Path) -> dict[str, set[str]]:
        raise AssertionError("hits derive from the scores; no second parse")

    monkeypatch.setattr(paths_stats, "parse_hmm_scores", counting_parse)
    monkeypatch.setattr(marker_extraction, "parse_hmm_output", forbidden_parse)

    summarizer = FullSummarizer(RESOURCES_ROOT)

    result = summarizer.summarize_query_full("q", query_dir, {})

    assert scores_reads == [query_dir / "hmmout" / "models.out.filtered"]
    assert result["query"] == "q"


@skip_if_no_runtime_resources()
def test_summarizer_and_scorer_share_one_label_index() -> None:
    summarizer = FullSummarizer(RESOURCES_ROOT)

    assert summarizer.contamination_scorer.labels is summarizer.labels
