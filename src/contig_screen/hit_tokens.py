"""Per-gene hit tokens of windows, the input of the set model.

Every gene of a window becomes one token that states, for each reference class
(bacteria, archaea, eukaryotes, nucleocytoviruses, mirusviruses, tailed phages, other
viruses), whether the gene has a hit in that class, the identity of the best hit and
its bit score per amino acid; which class holds the best hit overall; and whether the
gene has any hit at all. A gene without a hit is the zero token with ``any_hit`` 0,
a gene with only viral hits has no cellular entry, and the three virus groups stay
apart. The tokens come from the best hits of ``ref_hits.reference_features``. A model
that pools the tokens of a window (the set model, ``set_model``) embeds the taxonomic
distribution of the hits together with identity and bit score, window by window.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

CLASSES = ("BAC", "ARC", "EUK", "NCLDV", "MIRUS", "CAUDO", "OTHERV")
CELLULAR = ("BAC", "ARC", "EUK")
VIRAL = ("NCLDV", "MIRUS", "CAUDO", "OTHERV")
N_CLASSES = len(CLASSES)
GENE_FEATURES = (
    [f"has_{c}" for c in CLASSES]
    + [f"id_{c}" for c in CLASSES]
    + [f"spa_{c}" for c in CLASSES]
    + [f"top_{c}" for c in CLASSES]
    + ["any_hit", "margin"]
)
HAS, IDENT, SPA, TOP = 0, N_CLASSES, 2 * N_CLASSES, 3 * N_CLASSES
ANY_HIT = 4 * N_CLASSES
MARGIN = ANY_HIT + 1
SPA_CAP = 2.0
SCALARS = (
    "length_kb",
    "n_genes",
    "s_genes_per_kb",
    "s_coding_density",
    "c_gc",
    "s_gene_len_mean",
    "s_strand_switch_rate",
    "s_gap_mean",
    "s_overlap_frac",
)
TOP_KEYS = ["bitscore", "pident", "cls", "ref"]
TOP_ASCENDING = [False, False, True, True]
ABLATIONS = (
    "full",
    "id_only",
    "bits_only",
    "class_only",
    "viral_merged",
    "no_hit_dropped",
    "tokens_only",
    "esm_only",
)


def gene_rows(best: pd.DataFrame, index: pd.DataFrame) -> pd.DataFrame:
    """Return for every row of ``best`` the contig row and the gene row of the table.

    Genes are numbered from 0 within a contig; rows of contigs that are not in the
    index, or of genes beyond ``n_genes``, are dropped.
    """
    n_genes = index["n_genes"].to_numpy(dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(n_genes)])[:-1]
    piece_of = pd.Series(np.arange(len(index)), index=index["contig_id"])
    parts = best["qseqid"].str.rpartition("|")
    piece = piece_of.reindex(parts[0]).to_numpy()
    gene = pd.to_numeric(parts[2], errors="coerce").to_numpy()
    keep = ~np.isnan(piece) & ~np.isnan(gene)
    out = best[keep].copy()
    piece = piece[keep].astype(np.int64)
    gene = gene[keep].astype(np.int64)
    inside = gene < n_genes[piece]
    out = out[inside]
    out["piece"] = piece[inside]
    out["row"] = starts[piece[inside]] + gene[inside]
    return out


def build_genes(
    best: pd.DataFrame, index: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    """Return the token matrix (genes by features) and the first gene of each contig."""
    n_genes = index["n_genes"].to_numpy(dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(n_genes)])
    genes = np.zeros((int(offsets[-1]), len(GENE_FEATURES)), dtype=np.float32)
    rows = gene_rows(best, index)
    class_of = {name: k for k, name in enumerate(CLASSES)}
    for name, k in class_of.items():
        part = rows[rows["cls"] == name]
        genes[part["row"], HAS + k] = 1.0
        genes[part["row"], IDENT + k] = part["pident"].to_numpy() / 100.0
        genes[part["row"], SPA + k] = (
            np.minimum(part["spa"].to_numpy(), SPA_CAP) / SPA_CAP
        )
    top = rows.sort_values(TOP_KEYS, ascending=TOP_ASCENDING, kind="stable")
    top = top.drop_duplicates("qseqid")
    top_class = top["cls"].map(class_of).to_numpy()
    genes[top["row"], TOP + top_class] = 1.0
    genes[top["row"], ANY_HIT] = 1.0
    cellular = [SPA + class_of[c] for c in CELLULAR]
    viral = [SPA + class_of[c] for c in VIRAL]
    genes[:, MARGIN] = genes[:, cellular].max(axis=1) - genes[:, viral].max(axis=1)
    return genes, offsets


def contig_tokens(
    pieces: pd.DataFrame,
    best: pd.DataFrame,
    index: pd.DataFrame,
    features: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return tokens, offsets and raw scalars of the contigs of ``pieces``.

    The order of the contigs is the order of ``pieces``.

    ``index`` has ``contig_id`` and ``n_genes`` (a contig that is missing has no
    genes); ``features`` is indexed by contig id and holds ``SCALARS[2:]``.
    """
    index = index.set_index("contig_id").reindex(pieces["contig_id"]).reset_index()
    index["n_genes"] = index["n_genes"].fillna(0).astype(int)
    genes, offsets = build_genes(best, index)
    scalars = pd.DataFrame(
        {
            "length_kb": pieces["length_kb"].to_numpy(),
            "n_genes": index["n_genes"].to_numpy(),
        }
    )
    for name in SCALARS[2:]:
        scalars[name] = features[name].reindex(pieces["contig_id"]).to_numpy()
    return genes, offsets, scalars[list(SCALARS)].to_numpy(dtype=np.float32)


def window_tokens(
    features: pd.DataFrame, best: pd.DataFrame, index: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """Return tokens, offsets, raw scalars and the window table of a bin.

    The windows are those of ``features`` sorted by contig id.

    Raises:
        ValueError: No best hit matches a gene of the index (wrong query names).
    """
    pieces = features[["contig_id", "genome_id", "class", "length_kb"]]
    pieces = pieces.sort_values("contig_id").assign(y=1).reset_index(drop=True)
    genes, offsets, scalars = contig_tokens(
        pieces, best, index, features.set_index("contig_id")
    )
    if len(best) and not genes[:, ANY_HIT].any():
        msg = (
            f"none of the {len(best)} best hits matches a gene of the index: "
            "check the query names"
        )
        raise ValueError(msg)
    logger.info(
        "%d windows, %d genes, %d with a hit",
        len(pieces),
        len(genes),
        int(genes[:, ANY_HIT].sum()),
    )
    return genes, offsets, scalars, pieces


def ablate_genes(genes: np.ndarray, kind: str) -> np.ndarray:
    """Return the tokens with part of their content removed, to ask what matters.

    ``id_only`` keeps the identity of the best hit per class and drops the bit score,
    ``bits_only`` the opposite, ``class_only`` keeps only which classes have a hit and
    which holds the best one, and ``viral_merged`` puts the four virus classes (the
    nucleocytoviruses, mirusviruses, tailed phages and other viruses) into one slot, the
    first of them, with the largest value of each kind. ``tokens_only`` zeroes the
    columns after the hit tokens (the protein embedding and its flag, if the directory
    has them) and ``esm_only`` zeroes the hit tokens.
    """
    out = genes.copy()
    if kind == "full":
        return out
    if kind == "tokens_only":
        out[:, len(GENE_FEATURES) :] = 0.0
        return out
    if kind == "esm_only":
        out[:, : len(GENE_FEATURES)] = 0.0
        return out
    if kind in {"id_only", "class_only"}:
        out[:, SPA : SPA + N_CLASSES] = 0.0
        out[:, MARGIN] = 0.0
    if kind in {"bits_only", "class_only"}:
        out[:, IDENT : IDENT + N_CLASSES] = 0.0
    if kind == "viral_merged":
        viral = [CLASSES.index(c) for c in VIRAL]
        for start in (HAS, IDENT, SPA, TOP):
            merged = out[:, [start + k for k in viral]].max(axis=1)
            out[:, [start + k for k in viral]] = 0.0
            out[:, start + viral[0]] = merged
    elif kind not in {"id_only", "bits_only", "class_only"}:
        msg = f"unknown ablation {kind!r}; choose one of {ABLATIONS}"
        raise ValueError(msg)
    return out


def drop_no_hit(
    genes: np.ndarray, offsets: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the tokens of genes with a hit and the offsets of the shorter contigs."""
    keep = genes[:, ANY_HIT] > 0
    owner = np.repeat(np.arange(len(offsets) - 1), np.diff(offsets))
    counts = np.bincount(owner[keep], minlength=len(offsets) - 1)
    return genes[keep], np.concatenate([[0], np.cumsum(counts)])


def ablate(
    genes: np.ndarray, offsets: np.ndarray, kind: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return tokens and offsets after the ablation ``kind`` (see ``ABLATIONS``).

    ``no_hit_dropped`` leaves out the genes without a hit, so that the pooled embedding
    cannot see the zero token; the other ablations change the content of the tokens.
    """
    if kind == "no_hit_dropped":
        return drop_no_hit(genes, offsets)
    return ablate_genes(genes, kind), offsets
