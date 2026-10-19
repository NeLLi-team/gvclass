"""Tests of the per-gene hit tokens."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.contig_screen import hit_tokens as ht


def index() -> pd.DataFrame:
    return pd.DataFrame({"contig_id": ["g|3|0", "g|3|1"], "n_genes": [3, 2]})


def best() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "qseqid": [
                "g|3|0|0",
                "g|3|0|0",
                "g|3|0|1",
                "g|3|1|0",
                "g|3|0|7",
                "x|9|0|0",
            ],
            "cls": ["BAC", "NCLDV", "CAUDO", "EUK", "BAC", "BAC"],
            "ref": [4, 9, 2, 6, 1, 1],
            "pident": [60.0, 40.0, 35.0, 90.0, 99.0, 99.0],
            "bitscore": [300.0, 100.0, 80.0, 500.0, 900.0, 900.0],
            "spa": [1.0, 0.5, 0.4, 2.5, 1.0, 1.0],
        }
    )


def test_build_genes_gives_one_token_per_gene_with_zero_tokens_for_no_hit() -> None:
    genes, offsets = ht.build_genes(best(), index())

    assert genes.shape == (5, len(ht.GENE_FEATURES))
    assert offsets.tolist() == [0, 3, 5]
    assert not genes[2].any()


def test_a_token_states_presence_identity_and_bit_score_per_class() -> None:
    genes, _ = ht.build_genes(best(), index())
    bac, ncldv = ht.CLASSES.index("BAC"), ht.CLASSES.index("NCLDV")

    assert genes[0, ht.HAS + bac] == 1.0
    assert genes[0, ht.HAS + ncldv] == 1.0
    assert np.isclose(genes[0, ht.IDENT + bac], 0.60)
    assert np.isclose(genes[0, ht.SPA + bac], 0.5)
    assert np.isclose(genes[0, ht.SPA + ncldv], 0.25)


def test_the_top_class_is_the_class_of_the_best_bit_score() -> None:
    genes, _ = ht.build_genes(best(), index())

    assert genes[0, ht.TOP + ht.CLASSES.index("BAC")] == 1.0
    assert genes[0, ht.TOP : ht.TOP + ht.N_CLASSES].sum() == 1.0
    assert genes[1, ht.TOP + ht.CLASSES.index("CAUDO")] == 1.0
    assert genes[0, ht.ANY_HIT] == 1.0


def test_the_margin_is_the_best_cellular_minus_the_best_viral_bit_score() -> None:
    genes, _ = ht.build_genes(best(), index())

    assert np.isclose(genes[0, ht.MARGIN], 0.5 - 0.25)
    assert np.isclose(genes[1, ht.MARGIN], -0.2)
    assert genes[2, ht.MARGIN] == 0.0


def test_hits_of_unknown_contigs_and_of_genes_beyond_the_count_are_dropped() -> None:
    genes, offsets = ht.build_genes(best(), index())

    assert genes[:, ht.ANY_HIT].sum() == 3.0
    assert len(genes) == offsets[-1]


def test_equal_bit_scores_are_resolved_by_identity_then_by_class_name() -> None:
    table = pd.DataFrame(
        {
            "qseqid": ["g|3|0|0"] * 3,
            "cls": ["EUK", "BAC", "ARC"],
            "ref": [1, 2, 3],
            "pident": [40.0, 40.0, 30.0],
            "bitscore": [100.0] * 3,
            "spa": [0.5] * 3,
        }
    )

    genes, _ = ht.build_genes(
        table, pd.DataFrame({"contig_id": ["g|3|0"], "n_genes": [1]})
    )

    assert genes[0, ht.TOP + ht.CLASSES.index("BAC")] == 1.0


def test_contig_tokens_follow_the_order_of_the_pieces_not_of_the_tables() -> None:
    pieces = pd.DataFrame(
        {"contig_id": ["b|5|0", "a|3|0", "c|10|0"], "length_kb": [5.0, 3.0, 10.0]}
    )
    table = pd.DataFrame({"contig_id": ["a|3|0", "b|5|0"], "n_genes": [2, 1]})
    hits = pd.DataFrame(
        {
            "qseqid": ["a|3|0|1", "b|5|0|0"],
            "cls": ["EUK", "BAC"],
            "ref": [1, 2],
            "pident": [50.0, 70.0],
            "bitscore": [100.0, 200.0],
            "spa": [1.0, 1.0],
        }
    )
    features = pd.DataFrame(
        {
            name: [10.0 * (k + 1), 20.0 * (k + 1)]
            for k, name in enumerate(ht.SCALARS[2:])
        }
    ).set_index(pd.Index(["a|3|0", "b|5|0"], name="contig_id"))

    genes, offsets, scalars = ht.contig_tokens(pieces, hits, table, features)

    assert offsets.tolist() == [0, 1, 3, 3]
    assert genes[0, ht.HAS + ht.CLASSES.index("BAC")] == 1.0
    assert genes[2, ht.HAS + ht.CLASSES.index("EUK")] == 1.0
    assert scalars[:, 0].tolist() == [5.0, 3.0, 10.0]
    assert scalars[:, 1].tolist() == [1.0, 2.0, 0.0]
    assert np.isnan(scalars[2, 2])


def test_window_tokens_sort_the_windows_by_contig_id() -> None:
    features = pd.DataFrame(
        {
            "contig_id": ["g|3|1", "g|3|0"],
            "genome_id": "g",
            "class": "BIN",
            "length_kb": [3, 3],
            **dict.fromkeys(ht.SCALARS[2:], 1.0),
        }
    )

    genes, offsets, scalars, pieces = ht.window_tokens(features, best(), index())

    assert pieces["contig_id"].tolist() == ["g|3|0", "g|3|1"]
    assert offsets.tolist() == [0, 3, 5]
    assert scalars.shape == (2, len(ht.SCALARS))
    assert genes[0, ht.ANY_HIT] == 1.0


def test_window_tokens_refuse_hits_that_match_no_gene_of_the_index() -> None:
    features = pd.DataFrame(
        {
            "contig_id": ["g|3|0"],
            "genome_id": "g",
            "class": "BIN",
            "length_kb": [3],
            **dict.fromkeys(ht.SCALARS[2:], 1.0),
        }
    )
    stray = best().assign(qseqid="other|1|0|0")

    with pytest.raises(ValueError, match="matches a gene"):
        ht.window_tokens(features, stray, index())


def mixed_token() -> np.ndarray:
    genes, _ = ht.build_genes(
        pd.DataFrame(
            {
                "qseqid": ["g|3|0|0"] * 3,
                "cls": ["BAC", "NCLDV", "CAUDO"],
                "ref": [1, 2, 3],
                "pident": [60.0, 40.0, 30.0],
                "bitscore": [300.0, 200.0, 100.0],
                "spa": [1.0, 0.8, 0.4],
            }
        ),
        pd.DataFrame({"contig_id": ["g|3|0"], "n_genes": [1]}),
    )
    return genes


def test_the_full_ablation_leaves_the_tokens_unchanged() -> None:
    genes = mixed_token()

    assert np.array_equal(ht.ablate_genes(genes, "full"), genes)


def test_merging_the_viral_classes_puts_them_into_one_slot() -> None:
    merged = ht.ablate_genes(mixed_token(), "viral_merged")
    slot, phage = ht.CLASSES.index("NCLDV"), ht.CLASSES.index("CAUDO")

    assert merged[0, ht.HAS + slot] == 1.0
    assert merged[0, ht.HAS + phage] == 0.0
    assert np.isclose(merged[0, ht.IDENT + slot], 0.40)


def test_unknown_ablation_is_an_error() -> None:
    with pytest.raises(ValueError, match="ablation"):
        ht.ablate_genes(mixed_token(), "nonsense")


def test_dropping_the_genes_without_a_hit_shortens_the_contigs() -> None:
    genes, offsets = ht.build_genes(best(), index())

    kept, cuts = ht.ablate(genes, offsets, "no_hit_dropped")

    assert len(kept) == 3
    assert cuts.tolist() == [0, 2, 3]
