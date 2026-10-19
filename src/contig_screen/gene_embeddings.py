"""Per-gene protein embeddings that extend the hit tokens of the set model.

The hit tokens say what the database knows about each gene; a gene without a hit is
the zero token. A language-model embedding of the protein itself tells something
about every gene, with or without a hit. ``embed_proteins`` embeds the distinct
proteins of the windows with the mean-pooled ESM-2 embedding of
``extra_features.Embedder`` (stored as float16, as in the stored shards);
``extend_tokens`` appends to every token the raw embedding of its protein and a flag
that says whether the gene has one (a protein under 50 residues or with an unknown
residue has none and gets zeros). The set model standardizes the embedding itself
(``set_model.EmbeddingScale``). The proteins are ``(<window id>|<gene>, sequence)``
pairs; a gene of the tokens that has no protein is an error.
"""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Iterator, Sequence
from typing import Protocol

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

MIN_PROTEIN_AA = 50
EMBEDDING_DIM = 320  # width of the mean-pooled embedding of ESM-2 t6 8M


class Embedder(Protocol):
    """A protein embedder such as ``extra_features.Embedder``."""

    def embed(self, proteins: dict[str, str]) -> dict[str, np.ndarray]:
        """Return the embedding of every protein, by hash."""
        ...


def protein_hash(sequence: str) -> str:
    """Return the key of a protein, the same as ``extra_features`` uses."""
    return hashlib.sha1(sequence.encode(), usedforsecurity=False).hexdigest()


def eligible(sequence: str) -> bool:
    """Return True if the protein gets an embedding."""
    return len(sequence) >= MIN_PROTEIN_AA and "X" not in sequence


def embed_proteins(
    queries: Sequence[tuple[str, str]], embedder: Embedder
) -> tuple[pd.Index, np.ndarray]:
    """Return the hashes and the embeddings (float16) of the eligible proteins.

    The proteins are the distinct eligible ones of ``queries``, in the order of
    their hashes. With none, the matrix has no row and ``EMBEDDING_DIM`` columns.
    """
    proteins = {}
    for _, sequence in queries:
        if eligible(sequence):
            proteins[protein_hash(sequence)] = sequence
    keys = sorted(proteins)
    if not keys:
        return pd.Index(keys), np.zeros((0, EMBEDDING_DIM), dtype=np.float16)
    vectors = embedder.embed({k: proteins[k] for k in keys})
    return pd.Index(keys), np.stack([vectors[k] for k in keys]).astype(np.float16)


def gene_rows(
    queries: Sequence[tuple[str, str]], offsets: np.ndarray, contig_ids: pd.Series
) -> Iterator[tuple[int, str]]:
    """Yield the token row and the protein of every gene of the chosen contigs."""
    piece_of = {c: i for i, c in enumerate(contig_ids)}
    for name, sequence in queries:
        contig, _, gene = name.rpartition("|")
        piece = piece_of.get(contig)
        if piece is None:
            continue
        n_genes = int(offsets[piece + 1] - offsets[piece])
        if int(gene) >= n_genes:
            msg = f"gene {gene} of {contig} is beyond its {n_genes} tokens"
            raise ValueError(msg)
        yield int(offsets[piece]) + int(gene), sequence


def extend_tokens(
    queries: Sequence[tuple[str, str]],
    genes: np.ndarray,
    offsets: np.ndarray,
    pieces: pd.DataFrame,
    embedder: Embedder,
) -> np.ndarray:
    """Return the tokens of all windows with the protein embedding appended."""
    names, matrix = embed_proteins(queries, embedder)
    embedding = np.zeros((len(genes), matrix.shape[1]), dtype=np.float32)
    present = np.zeros(len(genes), dtype=bool)
    seen = np.zeros(len(genes), dtype=bool)
    for row, sequence in gene_rows(queries, offsets, pieces["contig_id"]):
        seen[row] = True
        if eligible(sequence):
            embedding[row] = matrix[names.get_loc(protein_hash(sequence))]
            present[row] = True
    if not seen.all():
        msg = f"the proteins lack {int((~seen).sum())} of {len(genes)} token genes"
        raise ValueError(msg)
    logger.info(
        "%d windows, %d genes, %d with an embedding",
        len(pieces),
        len(genes),
        present.sum(),
    )
    extended = np.hstack([genes, embedding, present[:, None].astype(np.float32)])
    return extended.astype(np.float32)
