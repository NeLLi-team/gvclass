"""Per-window reference-hit features from DIAMOND searches of the reference.

``hits.tsv`` is DIAMOND outfmt 6 with ``qseqid sseqid pident length qlen slen evalue
bitscore``. A query protein is named ``<window id>|<gene>``. The best hit of every
gene to each class (BAC, ARC, EUK, NCLDV, MIRUS, CAUDO, OTHERV) is kept; a new bin
is scored against every reference subject, so no hold-out applies. Features per
window: the share of genes whose best hit lies in each class, the share by identity
bin, the mean best identity and bit score per amino acid to each class (0 without a
hit), the share with any hit, and the entropy of the best-hit classes. Equal bit
scores are resolved by identity, then by reference index or class name
(``CLASS_KEYS``, ``TOP_KEYS``).
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd

CLASSES = ("BAC", "ARC", "EUK", "NCLDV", "MIRUS", "CAUDO", "OTHERV")
CELLULAR = ("BAC", "ARC", "EUK")
VIRAL = ("NCLDV", "MIRUS", "CAUDO", "OTHERV")
ID_EDGES = (0.0, 35.0, 50.0, 70.0, 90.0, 100.1)
ID_NAMES = ("lt35", "35_50", "50_70", "70_90", "ge90")
FAMILY_TOKEN = re.compile(r"^f__|(viridae|aceae|idae)$")
COLUMNS = ("qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore")
BEST_COLUMNS = ["qseqid", "cls", "ref", "pident", "bitscore", "spa"]
CHUNK_ROWS = 4_000_000
# Tie rule (review 15). Hits are ordered by bit score; of equal scores the higher
# identity wins, then the lower reference index (best hit of one class) or the class
# name (best hit over the classes). A choice therefore depends only on the rows of its
# own query, never on the other rows of the table or on their order.
CLASS_KEYS = ["qseqid", "cls", "bitscore", "pident", "ref"]
CLASS_ASCENDING = [True, True, False, False, True]
TOP_KEYS = ["qseqid", "bitscore", "pident", "cls", "ref"]
TOP_ASCENDING = [True, False, False, True, True]
RANK_PREFIX = re.compile(r"^[a-z]_{1,2}(?=[A-Za-z0-9])")


def normalize_token(token: str) -> str:
    """Return a lineage token without its rank prefix (``f__X``, ``f_X`` to ``X``)."""
    return RANK_PREFIX.sub("", token.strip())


def family_token(lineage: str) -> str:
    """Return the last token of a lineage that names a family, or an empty string."""
    for token in reversed([t.strip() for t in str(lineage).split(";")]):
        if token and FAMILY_TOKEN.search(token):
            return token
    return ""


def read_hits(path: Path) -> Iterator[pd.DataFrame]:
    """Yield DIAMOND hit chunks (none for an empty file)."""
    if path.stat().st_size == 0:
        return
    frames = pd.read_csv(
        path,
        sep="\t",
        names=COLUMNS,
        header=None,
        chunksize=CHUNK_ROWS,
        dtype={"qseqid": str, "sseqid": str},
    )
    yield from frames


def load_reference(
    ref_genomes: Path, extra_genomes: Path, skip_classes: list[str]
) -> pd.DataFrame:
    """Return the reference table with the skipped classes emptied and extras added."""
    ref = pd.read_csv(ref_genomes, sep="\t", dtype={"lineage": str})
    ref.loc[ref["class"].isin(skip_classes), "class"] = ""
    extra = pd.read_csv(extra_genomes, sep="\t", dtype={"lineage": str})
    return pd.concat([ref, extra], ignore_index=True)


def best_per_class(chunk: pd.DataFrame, ref_class: np.ndarray) -> pd.DataFrame:
    """Return the best hit of every (query, class) pair of one chunk."""
    sub = chunk["sseqid"].str.extract(r"^r(\d+)_")[0].astype(int).to_numpy()
    keep = ref_class[sub] != ""
    frame = chunk[keep].copy()
    frame["ref"] = sub[keep]
    frame["cls"] = ref_class[sub[keep]]
    frame["spa"] = frame["bitscore"] / frame["qlen"].clip(lower=1)
    frame = frame.sort_values(CLASS_KEYS, ascending=CLASS_ASCENDING)
    return frame.drop_duplicates(["qseqid", "cls"])[BEST_COLUMNS]


def collect_best(paths: list[Path], ref_class: np.ndarray) -> pd.DataFrame:
    """Return the best hit per query and class over all hit files."""
    pieces = [
        best_per_class(chunk, ref_class) for path in paths for chunk in read_hits(path)
    ]
    if not pieces:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in BEST_COLUMNS})
    return pd.concat(pieces, ignore_index=True)


def coherence(top: pd.DataFrame, ref_name: np.ndarray, label: str) -> pd.DataFrame:
    """Return per contig the commonest name among the best hits that carry a name.

    The share divides by the number of best hits with a name; the count is the number
    of distinct names.
    """
    names = pd.Series(ref_name[top["ref"].to_numpy()], index=top.index)
    frame = pd.DataFrame({"contig_id": top["contig_id"], "name": names})
    frame = frame[frame["name"] != ""]
    counts = frame.groupby(["contig_id", "name"]).size()
    total = counts.groupby(level=0).sum()
    top_count = counts.groupby(level=0).max()
    return pd.DataFrame(
        {
            f"r_top_{label}_share": top_count / total,
            f"r_n_{label}": counts.groupby(level=0).size(),
        }
    )


def piece_features(
    best: pd.DataFrame,
    index: pd.DataFrame,
    ref_family: np.ndarray,
    ref_prefix: np.ndarray,
) -> pd.DataFrame:
    """Return the r_ features per piece."""
    if best.empty:
        best = pd.DataFrame(
            {
                "qseqid": pd.Series(dtype=str),
                "cls": pd.Series(dtype=str),
                "ref": pd.Series(dtype=int),
                "pident": pd.Series(dtype=float),
                "bitscore": pd.Series(dtype=float),
                "spa": pd.Series(dtype=float),
            }
        )
    best = best.sort_values(CLASS_KEYS, ascending=CLASS_ASCENDING)
    best = best.drop_duplicates(["qseqid", "cls"])
    best["contig_id"] = [q.rpartition("|")[0] for q in best["qseqid"]]
    top = best.sort_values(TOP_KEYS, ascending=TOP_ASCENDING).drop_duplicates("qseqid")
    out = pd.DataFrame(index=pd.Index(index["contig_id"], name="contig_id"))
    n_genes = index.set_index("contig_id")["n_genes"].clip(lower=1)
    out["r_frac_hit"] = (
        top.groupby("contig_id").size().reindex(out.index).fillna(0) / n_genes
    )
    share = pd.crosstab(top["contig_id"], top["cls"]).reindex(out.index, fill_value=0)
    share = share.reindex(columns=list(CLASSES), fill_value=0).div(n_genes, axis=0)
    for cls in CLASSES:
        out[f"r_share_{cls}"] = share[cls]
    out["r_frac_cell"] = share[list(CELLULAR)].sum(axis=1)
    out["r_frac_vir"] = share[list(VIRAL)].sum(axis=1)
    top = top.assign(
        ibin=pd.cut(top["pident"], bins=ID_EDGES, labels=ID_NAMES, right=False).astype(
            str
        )
    )
    for cls in CLASSES:
        part = top[top["cls"] == cls]
        counts = pd.crosstab(part["contig_id"], part["ibin"])
        counts = counts.reindex(out.index, fill_value=0).reindex(
            columns=list(ID_NAMES), fill_value=0
        )
        for name in ID_NAMES:
            out[f"r_icell_{cls}_{name}"] = counts[name] / n_genes
        mine = best[best["cls"] == cls].groupby("contig_id")
        out[f"r_idbest_{cls}"] = (
            mine["pident"].sum().reindex(out.index).fillna(0) / n_genes
        )
        out[f"r_bsmax_{cls}"] = mine["spa"].sum().reindex(out.index).fillna(0) / n_genes
    for label, names in (("family", ref_family), ("p4", ref_prefix)):
        part = coherence(top, names, label)
        for column in part.columns:
            out[column] = part[column].reindex(out.index)
    hit_share = share.div(share.sum(axis=1).replace(0, np.nan), axis=0)
    probs = hit_share.fillna(0.0).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(probs > 0, np.log2(probs), 0.0)
    out["r_entropy"] = -(probs * logp).sum(axis=1)
    return out.reset_index()


def best_table(best: pd.DataFrame) -> pd.DataFrame:
    """Return the best hit of every gene and class, one row each, in a fixed order.

    The table is what ``piece_features`` aggregates; it also feeds the set model,
    which reads the hits of single genes (``hit_tokens``). Genes without a hit have
    no row.
    """
    best = best.sort_values(CLASS_KEYS, ascending=CLASS_ASCENDING)
    best = best.drop_duplicates(["qseqid", "cls"])
    best = best.sort_values(["qseqid", "cls"], kind="stable")
    return best[BEST_COLUMNS]


def reference_features(
    ref: pd.DataFrame, hit_paths: list[Path], index: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the r_ features per window and the table of best hits per gene."""
    ref_class = np.full(int(ref["index"].max()) + 1, "", dtype=object)
    ref_class[ref["index"].to_numpy()] = ref["class"].to_numpy()
    best = collect_best(hit_paths, ref_class)
    ref_family = np.full(len(ref_class), "", dtype=object)
    ref_prefix = np.full(len(ref_class), "", dtype=object)
    for idx, lineage in zip(ref["index"], ref["lineage"].fillna(""), strict=True):
        ref_family[int(idx)] = family_token(lineage)
        tokens = [t.strip() for t in str(lineage).split(";") if t.strip()]
        ref_prefix[int(idx)] = ";".join(tokens[:4])
    return piece_features(best, index, ref_family, ref_prefix), best_table(best)
