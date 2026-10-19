"""Tests of the copy-number, order and hit-agreement features of the screen."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.contig_screen.gv_order import (
    N_CLASSES,
    SKIP,
    Contig,
    IdentityCache,
    Priors,
    bin_features,
    composition,
    contig_from_genes,
    contig_taxonomy,
    pair_keys,
    taxonomy_features,
    top_code,
)

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"


def make_priors() -> Priors:
    return Priors(
        p_multi={"OG1": 0.01, "OG2": 0.5},
        default=0.1,
        mu=np.zeros(N_CLASSES + 1),
        sd=np.ones(N_CLASSES + 1),
        extra={
            "codon": (np.zeros(65), np.ones(65)),
            "struct": (np.zeros(4), np.ones(4)),
        },
    )


def make_contig(name: str, models: list[str], length: int = 10_000) -> Contig:
    rng = np.random.default_rng(len(name) + length)
    contig = Contig(
        name=name,
        genome="bin",
        length=length,
        models=models,
        proteins=[PROTEIN] * len(models),
        hashes=[f"h-{model}" for model in models],
    )
    contig.comp = composition("".join(rng.choice(list("ACGT"), size=500)))
    contig.codon = rng.random(65)
    contig.struct = rng.random(4)
    return contig


def gene_frame(rows: list[tuple[str, str, float | None]]) -> pd.DataFrame:
    """Rows are (model, class of the best hit, identity of the NCLDV hit)."""
    return pd.DataFrame(
        {
            "begin": range(1, len(rows) + 1),
            "end": range(100, 100 + len(rows)),
            "model": [row[0] for row in rows],
            "protein": PROTEIN,
            "prot_hash": [f"h{i}" for i in range(len(rows))],
            "tx_class": [row[1] for row in rows],
            "tx_family": ["Mimiviridae" if row[2] else "" for row in rows],
            "tx_genus": ["" for _ in rows],
            "tx_pident": [row[2] or np.nan for row in rows],
            "tx_bits": [100.0 if row[2] else np.nan for row in rows],
        }
    )


def test_a_marker_on_two_contigs_counts_as_a_duplicated_single_copy_marker() -> None:
    contigs = [make_contig("a", ["OG1", "OG2"]), make_contig("b", ["OG1"], 8_000)]

    table = bin_features(contigs, make_priors(), IdentityCache())

    first = table.iloc[0]
    assert first["n_dup_sc"] == 1
    assert first["dup_frac"] == pytest.approx(0.5)
    assert first["dup_weighted"] == pytest.approx(1.0)  # (-log10 0.01) / 2 markers
    assert table.iloc[1]["dup_frac"] == pytest.approx(1.0)


def test_a_marker_that_occurs_once_is_not_duplicated() -> None:
    contigs = [make_contig("a", ["OG1"]), make_contig("b", ["OG2"], 8_000)]

    table = bin_features(contigs, make_priors(), IdentityCache())

    assert table["n_dup_sc"].tolist() == [0, 0]
    assert table["dup_frac"].tolist() == [0.0, 0.0]


def test_pairs_of_markers_are_unordered_and_stop_after_the_skip_distance() -> None:
    models = [f"OG{i}" for i in range(SKIP + 3)]

    keys = pair_keys(models)

    assert ("OG0", "OG1") in keys
    assert ("OG0", f"OG{SKIP}") in keys
    assert ("OG0", f"OG{SKIP + 1}") not in keys
    assert pair_keys(["OG2", "OG1"]) == {("OG1", "OG2")}


def test_shared_neighbour_pairs_raise_pair_shared() -> None:
    contigs = [
        make_contig("a", ["OG1", "OG2"]),
        make_contig("b", ["OG2", "OG1"], 8_000),
    ]

    table = bin_features(contigs, make_priors(), IdentityCache())

    assert table["pair_shared"].tolist() == [1.0, 1.0]
    assert table["pair_shared_larger"].tolist() == [0.0, 1.0]


def test_contig_from_genes_keeps_marker_genes_inside_the_interval_in_order() -> None:
    genes = gene_frame([("OG2", "", None), ("none", "", None), ("OG1", "", None)])
    genes = genes.assign(strand=1, aa_len=50)

    contig = contig_from_genes(genes, "c", "bin", 0, 10_000)

    assert contig.models == ["OG2", "OG1"]
    assert contig.length == 10_000


def test_top_code_prefers_the_smaller_string_on_ties() -> None:
    assert top_code({"B": 2, "A": 2, "C": 1}) == "A"
    assert top_code({"A": 0}) == ""


def test_hit_agreement_flags_the_contig_whose_family_differs_from_the_bin() -> None:
    same = gene_frame([("OG1", "NCLDV", 80.0)] * 3)
    other = gene_frame([("OG1", "NCLDV", 80.0)] * 2)
    other["tx_family"] = "Phycodnaviridae"
    contigs = [make_contig(n, ["OG1"]) for n in ("a", "b", "c")]
    for contig, genes in zip(contigs, (same, same, other), strict=True):
        contig.tax = contig_taxonomy(genes)

    table = taxonomy_features(contigs)

    assert table["tx_family_differs"].tolist() == [0.0, 0.0, 1.0]
    assert table["tx_ncldv_frac"].tolist() == [1.0, 1.0, 1.0]


def test_contigs_without_hits_have_missing_agreement() -> None:
    genes = gene_frame([("OG1", "", None)])
    contigs = [make_contig("a", ["OG1"]), make_contig("b", ["OG2"])]
    for contig in contigs:
        contig.tax = contig_taxonomy(genes)

    table = taxonomy_features(contigs)

    assert table["tx_family_differs"].isna().all()
    assert table["tx_pident"].isna().all()


def test_the_taxonomy_columns_are_added_only_when_every_contig_has_hits() -> None:
    contigs = [make_contig("a", ["OG1"]), make_contig("b", ["OG1"], 8_000)]

    without = bin_features(contigs, make_priors(), IdentityCache())
    for contig in contigs:
        contig.tax = contig_taxonomy(gene_frame([("OG1", "NCLDV", 70.0)]))
    with_hits = bin_features(contigs, make_priors(), IdentityCache())

    assert "tx_hit_frac" not in without
    assert with_hits["tx_hit_frac"].tolist() == [1.0, 1.0]
