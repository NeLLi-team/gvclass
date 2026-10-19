"""Tests of the protein embeddings appended to the hit tokens."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.contig_screen import gene_embeddings as ge

LONG = "ACDEFGHIKL" * 6
OTHER = "MNPQRSTVWY" * 6


class FakeEmbedder:
    """Embeds a protein as a vector that starts with its length."""

    def __init__(self) -> None:
        self.asked: list[list[str]] = []

    def embed(self, proteins: dict[str, str]) -> dict[str, np.ndarray]:
        self.asked.append(list(proteins))
        return {
            h: np.array([len(p), 0.5, 1.0] + [0.0] * (ge.EMBEDDING_DIM - 3))
            for h, p in proteins.items()
        }


def test_eligible_proteins_are_long_enough_and_have_no_unknown_residue() -> None:
    assert ge.eligible(LONG)
    assert not ge.eligible(LONG[:49])
    assert not ge.eligible("X" + LONG)


def test_protein_hash_is_the_sha1_of_the_sequence() -> None:
    assert ge.protein_hash("A") == "6dcd4ce23d88e2ee9568ba546c007c63d9131c1b"


def test_embed_proteins_embeds_each_distinct_eligible_protein_once_in_hash_order() -> (
    None
):
    embedder = FakeEmbedder()
    queries = [("w|0", LONG), ("w|1", OTHER), ("x|0", LONG), ("x|1", "SHORT")]

    names, matrix = ge.embed_proteins(queries, embedder)

    assert list(names) == sorted([ge.protein_hash(LONG), ge.protein_hash(OTHER)])
    assert embedder.asked == [list(names)]
    assert matrix.dtype == np.float16
    assert matrix.shape == (2, ge.EMBEDDING_DIM)


def test_embed_proteins_without_an_eligible_protein_gives_an_empty_matrix() -> None:
    names, matrix = ge.embed_proteins([("w|0", "SHORT")], FakeEmbedder())

    assert len(names) == 0
    assert matrix.shape == (0, ge.EMBEDDING_DIM)


def tokens() -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    genes = np.ones((3, 4), dtype=np.float32)
    offsets = np.array([0, 2, 3])
    pieces = pd.DataFrame({"contig_id": ["w1", "w2"]})
    return genes, offsets, pieces


def test_extend_tokens_appends_the_embedding_and_a_flag_to_every_token() -> None:
    genes, offsets, pieces = tokens()
    queries = [("w1|0", LONG), ("w1|1", "SHORT"), ("w2|0", OTHER)]

    out = ge.extend_tokens(queries, genes, offsets, pieces, FakeEmbedder())

    assert out.shape == (3, 4 + ge.EMBEDDING_DIM + 1)
    assert out.dtype == np.float32
    assert out[:, -1].tolist() == [1.0, 0.0, 1.0]
    assert out[0, 4] == len(LONG)
    assert not out[1, 4:-1].any()
    assert out[2, 4] == len(OTHER)


def test_extend_tokens_of_a_bin_without_an_eligible_protein_are_zero_padded() -> None:
    genes, offsets, pieces = tokens()
    queries = [("w1|0", "A"), ("w1|1", "A"), ("w2|0", "A")]

    out = ge.extend_tokens(queries, genes, offsets, pieces, FakeEmbedder())

    assert out.shape == (3, 4 + ge.EMBEDDING_DIM + 1)
    assert not out[:, 4:].any()


def test_extend_tokens_refuses_proteins_that_lack_a_token_gene() -> None:
    genes, offsets, pieces = tokens()

    with pytest.raises(ValueError, match="lack 1 of 3"):
        ge.extend_tokens(
            [("w1|0", LONG), ("w2|0", OTHER)], genes, offsets, pieces, FakeEmbedder()
        )


def test_extend_tokens_refuse_a_gene_beyond_the_genes_of_its_window() -> None:
    genes, offsets, pieces = tokens()

    with pytest.raises(ValueError, match="beyond its 2 tokens"):
        ge.extend_tokens([("w1|2", LONG)], genes, offsets, pieces, FakeEmbedder())


def test_extend_tokens_ignore_proteins_of_windows_that_are_not_scored() -> None:
    genes, offsets, pieces = tokens()
    queries = [("w1|0", LONG), ("w1|1", LONG), ("w2|0", OTHER), ("zz|0", LONG)]

    out = ge.extend_tokens(queries, genes, offsets, pieces, FakeEmbedder())

    assert out.shape[0] == 3
