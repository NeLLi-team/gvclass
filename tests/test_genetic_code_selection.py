"""M5: genetic-code selection adds a 5% coding-density override gate.

Rule (maintainer-confirmed): keep a non-meta code over `-p meta` iff
  complete_hits >= meta+2  OR  avg_best_hit_score >= meta*1.10
  OR (both fail) coding_density >= meta*1.05 ; else keep meta.
Only the top-ranked non-meta candidate is compared with meta.
"""

from __future__ import annotations

from pathlib import Path

from src.core.genetic_code_optimizer import GeneticCodeOptimizer


def test_genecalling_table_reports_best_hit_score(tmp_path):
    from src.pipeline.query_processing_engine import _write_genecalling_comparison_table

    (tmp_path / "stats").mkdir()
    metrics = {0: {"complete_hits": 2, "avg_score": 10.0, "avg_best_hit_score": 95.0}}
    _write_genecalling_comparison_table(tmp_path, "query", metrics)
    row = (tmp_path / "stats" / "query.genecalling.tab").read_text().splitlines()[1]
    assert row.split("\t")[-1] == "95.0"


def _metrics(complete_hits, score, density):
    return {
        "complete_hits": complete_hits,
        "avg_best_hit_score": score,
        "coding_density": density,
    }


def _select(code_results):
    opt = GeneticCodeOptimizer(Path("/tmp"))
    code, _ = opt.select_best_code(code_results, list(code_results))
    return code


def test_density_override_selects_nonmeta():
    # equal hits, score below 10% margin, density >= 5% over meta -> keep code 11.
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 105.0, 53.0),  # 53 >= 50*1.05 (52.5); 105 < 110
    }
    assert _select(code_results) == 11


def test_zero_score_and_tied_density_do_not_override_meta():
    for density in (0.0, 50.0):
        assert (
            _select({0: _metrics(0, 0.0, density), 11: _metrics(0, 0.0, density)}) == 0
        )

    assert _select({0: _metrics(0, 0.0, 0.0), 11: _metrics(0, 0.0, 5.0)}) == 11


def test_density_just_below_5pct_keeps_meta():
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 105.0, 52.4),  # 52.4 < 52.5; hits/score gates also fail
    }
    assert _select(code_results) == 0


def test_hits_plus_two_selects_nonmeta():
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(7, 100.0, 50.0),  # +2 hits
    }
    assert _select(code_results) == 11


def test_hits_plus_one_keeps_meta():
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(6, 100.0, 50.0),  # only +1 hits; score/density tie
    }
    assert _select(code_results) == 0


def test_score_margin_above_ten_percent_selects_nonmeta():
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 111.0, 50.0),  # clearly above 1.10x
    }
    assert _select(code_results) == 11


def test_score_just_below_ten_percent_keeps_meta():
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 109.9, 50.0),
    }
    assert _select(code_results) == 0


def test_density_does_not_override_when_meta_outranks():
    # meta has the best score (ranks #1), so it wins outright. A lower-ranked
    # code with higher density must NOT override a clearly-better meta
    # ("otherwise keep -p meta").
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 90.0, 55.0),  # higher density but lower score than meta
    }
    assert _select(code_results) == 0


def test_density_override_fires_for_top_ranked_candidate():
    # code 11 ranks #1 (higher score than meta) but within the 10% gate; its
    # >=5% higher density rescues it as the selected code.
    code_results = {
        0: _metrics(5, 100.0, 50.0),
        11: _metrics(5, 104.0, 53.0),  # score<110, density 53>=52.5
    }
    assert _select(code_results) == 11


def test_no_meta_keeps_best_ranked():
    code_results = {
        11: _metrics(5, 100.0, 50.0),
        1: _metrics(3, 90.0, 40.0),
    }
    assert _select(code_results) == 11
