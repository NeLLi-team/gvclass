"""Marker tables: raw rows, best row per protein and model, counts per model."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from src.core.hmm_search import (
    DomainHit,
    best_domain_per_pair,
    count_proteins_per_model,
    read_proteins,
    write_counts,
    write_domain_table,
)


def _domain(
    protein: str,
    model: str,
    score: float,
    domain_score: float,
    **changes: object,
) -> DomainHit:
    row = DomainHit(
        protein=protein,
        protein_length=100,
        model=model,
        model_accession="-",
        model_length=300,
        evalue=1e-20,
        score=score,
        bias=0.0,
        domain_index=1,
        domain_count=1,
        c_evalue=1e-20,
        i_evalue=1e-20,
        domain_score=domain_score,
        domain_bias=0.0,
        hmm_from=1,
        hmm_to=100,
        ali_from=1,
        ali_to=100,
        env_from=1,
        env_to=100,
    )
    return replace(row, **changes)


def _rows(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]


def test_read_proteins_forces_amino_alphabet(tmp_path: Path) -> None:
    query_file = tmp_path / "low_complexity.faa"
    query_file.write_text(">low_complexity_protein\n" + "A" * 64 + "\n")

    sequences = read_proteins(query_file)

    assert len(sequences) == 1
    assert str(sequences[0].alphabet) == "AA()"


def test_best_domain_per_pair_keeps_highest_scoring_row() -> None:
    domains = [
        _domain("prot1", "GVOGm0022", 50.0, 50.0),
        _domain("prot1", "GVOGm0022", 60.0, 55.0),
        _domain("prot2", "GVOGm0022", 55.0, 55.0),
        _domain("prot1", "GVOGm0054", 42.0, 42.0),
    ]

    best = best_domain_per_pair(domains)

    assert [(row.protein, row.model, row.score) for row in best] == [
        ("prot1", "GVOGm0022", 60.0),
        ("prot2", "GVOGm0022", 55.0),
        ("prot1", "GVOGm0054", 42.0),
    ]


def test_best_domain_per_pair_breaks_score_ties_by_domain_score() -> None:
    domains = [
        _domain("prot1", "M1", 50.0, 40.0),
        _domain("prot1", "M1", 50.0, 45.0),
    ]

    (best,) = best_domain_per_pair(domains)

    assert best.domain_score == 45.0


def test_best_domain_per_pair_keeps_first_row_on_full_tie() -> None:
    first = _domain("prot1", "M1", 50.0, 40.0, domain_index=1)
    second = _domain("prot1", "M1", 50.0, 40.0, domain_index=2)

    (best,) = best_domain_per_pair([first, second])

    assert best.domain_index == 1


def test_count_proteins_per_model_counts_each_protein_once() -> None:
    best = best_domain_per_pair(
        [
            _domain("prot1", "GVOGm0022", 50.0, 50.0),
            _domain("prot1", "GVOGm0022", 60.0, 60.0),
            _domain("prot2", "GVOGm0022", 55.0, 55.0),
            _domain("prot1", "GVOGm0054", 42.0, 42.0),
        ]
    )

    assert count_proteins_per_model(best) == {"GVOGm0022": 2, "GVOGm0054": 1}


def test_write_counts_sorts_models(tmp_path: Path) -> None:
    counts_file = tmp_path / "models.counts"

    write_counts({"GVOGm0054": 1, "GVOGm0022": 2}, counts_file)

    assert counts_file.read_text(encoding="utf-8") == "GVOGm0022\t2\nGVOGm0054\t1\n"


def test_write_domain_table_keeps_the_downstream_column_layout(
    tmp_path: Path,
) -> None:
    table = tmp_path / "models.out"
    row = _domain(
        "q|prot_1",
        "OG0001",
        58.4,
        57.9,
        protein_length=208,
        model_length=179,
        evalue=3.95e-17,
        bias=0.8,
        c_evalue=7.34e-20,
        i_evalue=3.95e-17,
        domain_bias=0.7,
        hmm_from=17,
        hmm_to=172,
        ali_from=9,
        ali_to=165,
        env_from=1,
        env_to=171,
    )

    write_domain_table([row], table)

    assert table.read_text(encoding="utf-8").startswith("# target name\t")
    (fields,) = _rows(table)
    assert len(fields) == 22
    assert fields[0] == "q|prot_1"
    assert fields[2] == "208"
    assert fields[3] == "OG0001"
    assert fields[5] == "179"
    assert fields[6] == "3.95e-17"
    assert fields[7] == "58.4"
    assert fields[13] == "57.9"
    assert fields[15:21] == ["17", "172", "9", "165", "1", "171"]
    assert fields[21] == "-"
