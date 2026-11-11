"""Genetic-code selection needs marker evidence; coding density never decides.

A nonzero code replaces meta (code 0) only with at least one complete hit and
either two more complete hits than meta, or at least as many and a 10% higher
average best-hit score. Input under 20,000 bp is evaluated with meta alone.
"""

from __future__ import annotations

import logging
import random
from pathlib import Path

import pyhmmer
import pytest
from Bio.Data.CodonTable import standard_dna_table

from src.core import genetic_code_optimizer
from src.core.genetic_code_optimizer import (
    AssemblyStats,
    CodeEvaluation,
    GeneCallingResult,
    optimize_genetic_code,
    select_best_code,
)
from src.pipeline.query_processing_engine import (
    _find_query_gff_file,
    _write_genecalling_comparison_table,
    _write_selected_genecode_table,
)

ALL_CODES = [0, 1, 4, 6, 11, 15, 29, 106, 129]
MODEL_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLQDKPEAIRRLVDELAQLGVE"


def _evaluation(
    code: int, complete_hits: int, score: float, density: float = 50.0
) -> CodeEvaluation:
    return CodeEvaluation(
        code=code,
        gene_count=100,
        coding_density=density,
        complete_hits=complete_hits,
        avg_best_hit_score=score,
    )


def _select(*evaluations: CodeEvaluation) -> int:
    return select_best_code({item.code: item for item in evaluations})


def test_genecalling_table_reports_best_hit_score(tmp_path: Path) -> None:
    (tmp_path / "stats").mkdir()
    result = GeneCallingResult(
        best_code=0,
        evaluations={0: _evaluation(0, 2, 95.0)},
        assembly=AssemblyStats(contigs=1, length_bp=1000, gc_percent=40.0),
        faa=tmp_path / "q.faa",
        fna=tmp_path / "q.fna",
        gff=tmp_path / "q.gff",
    )

    _write_genecalling_comparison_table(tmp_path, "query", result)

    row = (tmp_path / "stats" / "query.genecalling.tab").read_text().splitlines()[1]
    assert row.split("\t") == ["0", "100", "1000", "50.000", "2", "95.0"]


def _result(tmp_path: Path, best_code: int) -> GeneCallingResult:
    return GeneCallingResult(
        best_code=best_code,
        evaluations={0: _evaluation(0, 2, 95.0), 11: _evaluation(11, 4, 99.0)},
        assembly=AssemblyStats(contigs=1, length_bp=1000, gc_percent=40.0),
        faa=tmp_path / "q.faa",
        fna=tmp_path / "q.fna",
        gff=tmp_path / "q.gff",
    )


@pytest.mark.parametrize(("best_code", "ttable"), [(0, "codemeta"), (11, "11")])
def test_selected_genecode_table_names_the_selected_code(
    tmp_path: Path, best_code: int, ttable: str
) -> None:
    (tmp_path / "stats").mkdir()

    _write_selected_genecode_table(tmp_path, "query", _result(tmp_path, best_code))

    lines = (tmp_path / "stats" / "query.stats.tab").read_text().splitlines()
    assert lines == [
        f"ttable\t{ttable}",
        "genes\t100",
        "bases\t1000",
        "coding_density\t50.000",
    ]


def test_query_gff_file_is_the_gene_calling_output_only(tmp_path: Path) -> None:
    gene_dir = tmp_path / "gene_calling" / "query_reformatted"
    gene_dir.mkdir(parents=True)
    (tmp_path / "unrelated.gff").write_text("unrelated\n")

    assert _find_query_gff_file(tmp_path, "query") is None

    gff = gene_dir / "query_reformatted.gff"
    gff.write_text("##gff-version  3\n")

    assert _find_query_gff_file(tmp_path, "query") == gff


def test_zero_hits_everywhere_selects_meta_whatever_the_coding_density() -> None:
    assert _select(_evaluation(0, 0, 0.0, 6.5), _evaluation(6, 0, 0.0, 72.0)) == 0


def test_one_more_complete_hit_keeps_meta() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 6, 100.0)) == 0


def test_two_more_complete_hits_select_the_nonmeta_code() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 7, 100.0)) == 11


def test_two_more_complete_hits_win_despite_a_lower_score() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 7, 60.0)) == 11


def test_ten_percent_score_lead_with_equal_hits_selects_the_nonmeta_code() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 5, 110.0)) == 11


def test_score_lead_below_ten_percent_keeps_meta() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 5, 109.9)) == 0


def test_score_lead_with_fewer_complete_hits_keeps_meta() -> None:
    assert _select(_evaluation(0, 5, 100.0), _evaluation(11, 4, 150.0)) == 0


def test_score_lead_without_a_complete_hit_keeps_meta() -> None:
    assert _select(_evaluation(0, 0, 10.0), _evaluation(11, 0, 50.0)) == 0


def test_first_complete_hit_selects_the_nonmeta_code_over_a_hitless_meta() -> None:
    assert _select(_evaluation(0, 0, 0.0), _evaluation(11, 1, 80.0)) == 11


def test_coding_density_never_overrides_meta() -> None:
    meta = _evaluation(0, 5, 100.0, density=50.0)
    other = _evaluation(11, 5, 105.0, density=90.0)

    assert _select(meta, other) == 0


def test_qualifying_codes_rank_by_hits_then_score_then_preference() -> None:
    meta = _evaluation(0, 2, 60.0)

    assert _select(meta, _evaluation(11, 8, 110.0), _evaluation(106, 7, 140.0)) == 11
    assert _select(meta, _evaluation(11, 8, 100.0), _evaluation(106, 8, 120.0)) == 106
    assert _select(meta, _evaluation(106, 8, 100.0), _evaluation(11, 8, 100.0)) == 11
    assert _select(meta, _evaluation(15, 8, 100.0), _evaluation(4, 8, 100.0)) == 4


def test_a_nonqualifying_code_never_beats_a_qualifying_one() -> None:
    meta = _evaluation(0, 5, 100.0)
    near_miss = _evaluation(11, 6, 100.0)
    qualifies = _evaluation(106, 7, 100.0)

    assert _select(meta, near_miss, qualifies) == 106


def _random_dna(length: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


def _write_fna(path: Path, contigs: dict[str, str]) -> Path:
    path.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in contigs.items()),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def tiny_hmm(tmp_path: Path) -> Path:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(
        name=b"tiny_model", sequence=MODEL_PROTEIN
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence, pyhmmer.plan7.Background(alphabet)
    )
    hmm_path = tmp_path / "tiny.hmm"
    with hmm_path.open("wb") as handle:
        hmm.write(handle)
    return hmm_path


def test_input_under_20kb_is_evaluated_with_meta_only(
    tmp_path: Path, tiny_hmm: Path, caplog: pytest.LogCaptureFixture
) -> None:
    query = _write_fna(tmp_path / "q.fna", {"q|c1": _random_dna(5000, 1)})

    with caplog.at_level(logging.INFO):
        result = optimize_genetic_code(
            query, tmp_path / "out", tiny_hmm, ALL_CODES, threads=4
        )

    assert result.best_code == 0
    assert set(result.evaluations) == {0}
    assert [
        r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING
    ] == []
    assert sum("meta mode only" in r.getMessage() for r in caplog.records) == 1


def test_random_dna_without_marker_hits_selects_meta(
    tmp_path: Path, tiny_hmm: Path
) -> None:
    query = _write_fna(tmp_path / "q.fna", {"q|c1": _random_dna(25000, 2)})

    result = optimize_genetic_code(
        query, tmp_path / "out", tiny_hmm, ALL_CODES, threads=2
    )

    assert set(result.evaluations) == set(ALL_CODES)
    assert result.best_code == 0


def test_an_error_in_one_code_evaluation_propagates(
    tmp_path: Path, tiny_hmm: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = _write_fna(tmp_path / "q.fna", {"q|c1": _random_dna(25000, 3)})
    call_genes = genetic_code_optimizer._call_genes

    def failing_call_genes(records: list, code: int):
        if code == 11:
            raise RuntimeError("injected failure for code 11")
        return call_genes(records, code)

    monkeypatch.setattr(genetic_code_optimizer, "_call_genes", failing_call_genes)

    with pytest.raises(RuntimeError, match="injected failure for code 11"):
        optimize_genetic_code(query, tmp_path / "out", tiny_hmm, [0, 11], threads=2)


def test_stats_file_states_the_selection_rule(tmp_path: Path, tiny_hmm: Path) -> None:
    query = _write_fna(tmp_path / "q.fna", {"q|c1": _random_dna(5000, 4)})

    result = optimize_genetic_code(query, tmp_path / "out", tiny_hmm, [0], threads=1)

    statsout = result.faa.with_suffix(".statsout").read_text(encoding="utf-8")
    assert "at least 2 more complete_hits" in statsout
    assert "10%" in statsout
    assert "coding_density" in statsout
    assert "meta only" in statsout


_CODON_BY_AMINO_ACID = {
    amino_acid: codon for codon, amino_acid in standard_dna_table.forward_table.items()
}
_MODEL_GENE = (
    "ATG" + "".join(_CODON_BY_AMINO_ACID[aa] for aa in MODEL_PROTEIN[1:]) + "TAA"
)


def _contig_with_model_genes(seed: int) -> str:
    flank = _random_dna(3000, seed)
    return flank + _MODEL_GENE + _random_dna(3000, seed + 100) + _MODEL_GENE + flank


def _output_files(result: GeneCallingResult) -> dict[str, bytes]:
    paths = (result.faa, result.fna, result.gff, result.faa.with_suffix(".statsout"))
    return {path.suffix: path.read_bytes() for path in paths}


def test_the_thread_count_changes_neither_the_selection_nor_any_output_file(
    tmp_path: Path, tiny_hmm: Path
) -> None:
    contigs = {f"q|c{seed}": _contig_with_model_genes(seed) for seed in range(1, 5)}
    query = _write_fna(tmp_path / "q.fna", contigs)

    single = optimize_genetic_code(
        query, tmp_path / "one", tiny_hmm, ALL_CODES, threads=1
    )
    several = optimize_genetic_code(
        query, tmp_path / "several", tiny_hmm, ALL_CODES, threads=4
    )

    assert set(single.evaluations) == set(ALL_CODES)
    assert single.evaluations[0].complete_hits > 0
    assert several.best_code == single.best_code
    assert several.evaluations == single.evaluations
    assert _output_files(several) == _output_files(single)
