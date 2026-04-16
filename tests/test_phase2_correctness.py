"""Regression tests for Phase 2 of the v1.4.3 remediation plan.

Covers:
* 2.1 Taxonomy majority per marker (including the taxonomy_confidence
  column emitted by the summary writer).
* 2.2 Marker extraction keyed by record.id with record.description fallback.
* 2.3 Genetic-code deterministic tiebreak with an absolute margin over meta.
* 2.5 CLI resolver default for completeness_mode matches default_config.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace
from typing import Dict


# ---------------------------------------------------------------------------
# 2.1 Taxonomy majority per marker.
# ---------------------------------------------------------------------------


def _make_summarizer(tmp_path: Path):
    """Construct a FullSummarizer with a minimal on-disk database."""
    from src.core.summarize_full import FullSummarizer

    db = tmp_path / "db"
    db.mkdir()
    (db / "gvclassFeb26_labels.tsv").write_text("")
    (db / "order_completeness.tab").write_text(
        "Order\tOrthogroups\tAverage_Percent\tStd_Percent\n"
    )
    return FullSummarizer(db)


def test_per_marker_majority_requires_distinct_marker_support(tmp_path: Path) -> None:
    """A single paralog-heavy marker must not dominate the vote."""
    summarizer = _make_summarizer(tmp_path)

    paralog_counter: Counter = Counter()
    paralog_counter["NCLDV__Paraloplace"] = 10

    single_counter: Counter = Counter()
    single_counter["NCLDV__Imitervirales"] = 1

    per_marker = {
        "MCP_paralog": paralog_counter,
        "MARKER_A": Counter(single_counter),
        "MARKER_B": Counter(single_counter),
        "MARKER_C": Counter(single_counter),
    }

    winner, supporting, total = summarizer._per_marker_majority(per_marker)

    assert winner == "NCLDV__Imitervirales"
    assert supporting == 3
    assert total == 4


def test_build_consensus_taxonomies_emits_low_support_when_below_threshold(
    tmp_path: Path,
) -> None:
    summarizer = _make_summarizer(tmp_path)

    flat_counters: Dict[str, Counter] = {
        level: Counter() for level in summarizer.TAX_LEVELS
    }
    flat_counters["order"]["NCLDV__Imitervirales"] = 2

    per_marker: Dict[str, Dict[str, Counter]] = {
        level: defaultdict(Counter) for level in summarizer.TAX_LEVELS
    }
    per_marker["order"]["MARKER_A"]["NCLDV__Imitervirales"] = 1
    per_marker["order"]["MARKER_B"]["NCLDV__Imitervirales"] = 1

    majority, _, confidence = summarizer._build_consensus_taxonomies(
        flat_counters, per_marker, mode_fast=False
    )

    order_segment = [piece for piece in majority.split(";") if piece.startswith("o_")]
    assert order_segment == ["o_"]
    assert "low_support" in confidence


def test_build_consensus_taxonomies_fast_mode_relaxes_order_threshold(
    tmp_path: Path,
) -> None:
    summarizer = _make_summarizer(tmp_path)

    flat_counters: Dict[str, Counter] = {
        level: Counter() for level in summarizer.TAX_LEVELS
    }
    flat_counters["order"]["NCLDV__Imitervirales"] = 2

    per_marker: Dict[str, Dict[str, Counter]] = {
        level: defaultdict(Counter) for level in summarizer.TAX_LEVELS
    }
    per_marker["order"]["MARKER_A"]["NCLDV__Imitervirales"] = 1
    per_marker["order"]["MARKER_B"]["NCLDV__Imitervirales"] = 1

    majority, _, confidence = summarizer._build_consensus_taxonomies(
        flat_counters, per_marker, mode_fast=True
    )

    assert "o_Imitervirales" in majority
    assert "reduced_fastmode" in confidence


def test_taxonomy_confidence_column_present_in_summary_schemas() -> None:
    from src.pipeline.summary_writer import FINAL_SUMMARY_COLUMNS, LEGACY_SUMMARY_HEADERS

    assert "taxonomy_confidence" in FINAL_SUMMARY_COLUMNS
    assert "taxonomy_confidence" in LEGACY_SUMMARY_HEADERS


# ---------------------------------------------------------------------------
# 2.2 Marker extraction by record.id with description fallback.
# ---------------------------------------------------------------------------


def test_extract_marker_sequences_matches_by_record_id(tmp_path: Path) -> None:
    from src.core.marker_extraction import extract_marker_sequences

    query_faa = tmp_path / "q.faa"
    query_faa.write_text(
        ">prot_1 # start=1 # end=300 # annotation\n"
        "MKLILV\n"
        ">prot_2 some descriptor\n"
        "MTKAAG\n"
    )
    output_dir = tmp_path / "markers"

    marker_hits = {"GVOGm0022": {"prot_1"}, "GVOGm0054": {"prot_2"}}

    files = extract_marker_sequences(query_faa, marker_hits, output_dir)

    assert set(files.keys()) == {"GVOGm0022", "GVOGm0054"}
    for marker in files:
        contents = files[marker].read_text()
        assert contents.startswith(">prot_"), contents


# ---------------------------------------------------------------------------
# 2.3 Genetic-code tiebreak.
# ---------------------------------------------------------------------------


def _code_results(codes):
    """Build a minimal ``code_results`` dict keyed by integer genetic code."""
    return {
        code: {
            "complete_hits": metrics.get("complete_hits", 0),
            "avg_best_hit_score": metrics.get("avg_best_hit_score", 0.0),
            "coding_density": metrics.get("coding_density", 0.0),
        }
        for code, metrics in codes.items()
    }


def test_select_best_code_standard_beats_exotic_on_tie(tmp_path: Path) -> None:
    """With identical metrics, standard code 11 must beat exotic code 106 and
    the ``_CODE_PREFERENCE_RANK`` ordering must never rank 15 ahead of 4."""
    from src.core.genetic_code_optimizer import GeneticCodeOptimizer

    optimizer = GeneticCodeOptimizer(database_path=tmp_path, threads=1)
    results = _code_results({
        0: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
        4: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
        11: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
        15: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
        106: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
    })

    key_11 = optimizer._score_key(11, results[11])
    key_4 = optimizer._score_key(4, results[4])
    key_15 = optimizer._score_key(15, results[15])
    key_106 = optimizer._score_key(106, results[106])
    assert key_11 > key_4 > key_15
    assert key_4 > key_106

    best_code, _ = optimizer.select_best_code(results, list(results.keys()))
    assert best_code == 0


def test_select_best_code_switches_to_standard_with_strong_hits_improvement(
    tmp_path: Path,
) -> None:
    from src.core.genetic_code_optimizer import GeneticCodeOptimizer

    optimizer = GeneticCodeOptimizer(database_path=tmp_path, threads=1)
    results = _code_results({
        0: {"complete_hits": 3, "avg_best_hit_score": 80.0, "coding_density": 0.85},
        11: {"complete_hits": 5, "avg_best_hit_score": 90.0, "coding_density": 0.90},
    })
    best_code, best = optimizer.select_best_code(results, [0, 11])
    assert best_code == 11
    assert best["complete_hits"] == 5


def test_select_best_code_rejects_noise_driven_switch_to_exotic(
    tmp_path: Path,
) -> None:
    """A 3% score bump over meta must NOT flip to exotic code 106."""
    from src.core.genetic_code_optimizer import GeneticCodeOptimizer

    optimizer = GeneticCodeOptimizer(database_path=tmp_path, threads=1)
    results = _code_results({
        0: {"complete_hits": 5, "avg_best_hit_score": 100.0, "coding_density": 0.90},
        106: {"complete_hits": 5, "avg_best_hit_score": 103.0, "coding_density": 0.93},
    })
    best_code, _ = optimizer.select_best_code(results, [0, 106])
    assert best_code == 0


def test_select_best_code_code_106_no_longer_displaces_code_11(tmp_path: Path) -> None:
    """Codex-audit regression: code 11 must win when it has the most complete hits."""
    from src.core.genetic_code_optimizer import GeneticCodeOptimizer

    optimizer = GeneticCodeOptimizer(database_path=tmp_path, threads=1)
    results = _code_results({
        0: {"complete_hits": 2, "avg_best_hit_score": 60.0, "coding_density": 0.80},
        11: {"complete_hits": 8, "avg_best_hit_score": 110.0, "coding_density": 0.92},
        106: {"complete_hits": 7, "avg_best_hit_score": 110.0, "coding_density": 0.92},
    })
    best_code, _ = optimizer.select_best_code(results, [0, 11, 106])
    assert best_code == 11


# ---------------------------------------------------------------------------
# 2.5 Completeness-mode fallback.
# ---------------------------------------------------------------------------


def test_resolve_completeness_mode_falls_back_to_novelty_aware() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode=None)
    config = {"pipeline": {}}
    assert resolve_completeness_mode(args, config) == "novelty-aware"


def test_resolve_completeness_mode_honors_explicit_legacy() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode=None)
    config = {"pipeline": {"completeness_mode": "legacy"}}
    assert resolve_completeness_mode(args, config) == "legacy"


def test_resolve_completeness_mode_cli_flag_wins_over_config() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode="legacy")
    config = {"pipeline": {"completeness_mode": "novelty-aware"}}
    assert resolve_completeness_mode(args, config) == "legacy"
