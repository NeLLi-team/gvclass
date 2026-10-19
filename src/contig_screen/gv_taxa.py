"""Best catalog hit of every gene of a giant-virus bin.

The proteins of the bin are searched with DIAMOND against the catalog of reference
genomes (``ref7x``) without any hold-out, as a new bin is scored by ``cellular``.
Query ``bin|i`` is row ``i`` of the gene table. For each gene the table holds the class
of the best hit over the catalog classes, and the family and genus code, identity and
bit score of its best NCLDV hit.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.contig_screen import layout, ref_hits as rh
from src.contig_screen.cellular import diamond_binary, diamond_search

CATALOG_TARGETS = 100
QUERY_UNIT = "bin"
BEST_COLUMNS = ["qseqid", "cls", "ref", "pident", "bitscore", "spa"]
TABLE_COLUMNS = ["row", "tx_class", "tx_family", "tx_genus", "tx_pident", "tx_bits"]


def codes(lineage: str) -> tuple[str, str]:
    """Return the family and genus code of a catalog NCLDV lineage."""
    parts = [rh.normalize_token(p) for p in str(lineage).split(";")]
    return (parts[4] if len(parts) > 4 else ""), (parts[5] if len(parts) > 5 else "")


def fasta_text(unit: str, proteins: Sequence[str]) -> str:
    """Return the FASTA records ``>unit|row`` of the non-empty proteins."""
    return "".join(
        f">{unit}|{row}\n{protein}\n" for row, protein in enumerate(proteins) if protein
    )


def best_catalog_hits(hits: Path, ref: pd.DataFrame) -> pd.DataFrame:
    """Return the best hit of every (gene, class) pair."""
    if hits.stat().st_size == 0:
        return pd.DataFrame(columns=BEST_COLUMNS)
    ref_class = np.full(int(ref["index"].max()) + 1, "", dtype=object)
    ref_class[ref["index"].to_numpy()] = ref["class"].to_numpy()
    pieces = [rh.best_per_class(chunk, ref_class) for chunk in rh.read_hits(hits)]
    best = pd.concat(pieces, ignore_index=True)
    best = best.sort_values(rh.CLASS_KEYS, ascending=rh.CLASS_ASCENDING)
    return best.drop_duplicates(["qseqid", "cls"])


def gene_table(best: pd.DataFrame, lineage: dict[int, str]) -> pd.DataFrame:
    """Return per gene the best class and the codes of the best NCLDV hit."""
    if best.empty:
        return pd.DataFrame(columns=TABLE_COLUMNS)
    top = best.sort_values(rh.TOP_KEYS, ascending=rh.TOP_ASCENDING)
    top = top.drop_duplicates("qseqid")[["qseqid", "cls"]]
    ncldv = best[best["cls"] == "NCLDV"]
    by_ref = {i: codes(lineage[i]) for i in ncldv["ref"].unique()}
    ncldv = ncldv.assign(
        tx_family=ncldv["ref"].map(lambda i: by_ref[i][0]),
        tx_genus=ncldv["ref"].map(lambda i: by_ref[i][1]),
        tx_pident=ncldv["pident"],
        tx_bits=ncldv["bitscore"],
    )
    table = top.rename(columns={"cls": "tx_class"}).merge(
        ncldv[["qseqid", "tx_family", "tx_genus", "tx_pident", "tx_bits"]],
        on="qseqid",
        how="left",
    )
    table[["tx_family", "tx_genus"]] = table[["tx_family", "tx_genus"]].fillna("")
    table["row"] = table["qseqid"].str.rsplit("|", n=1).str[1].astype(int)
    return table.drop(columns="qseqid")[TABLE_COLUMNS]


def taxa_for_bin(
    proteins: Sequence[str],
    resources: layout.ScreenResources,
    work: Path,
    *,
    threads: int,
    diamond: Path | None = None,
) -> pd.DataFrame:
    """Search the proteins of one bin against the catalog; return its gene table.

    ``proteins`` has one entry per row of the gene table. A bin without proteins
    skips the search.
    """
    queries, hits = work / "proteins.faa", work / "hitsx.tsv"
    text = fasta_text(QUERY_UNIT, proteins)
    queries.write_text(text, encoding="utf-8")
    if text:
        database = resources.path(layout.REFX_DMND).with_suffix("")
        diamond_search(
            diamond_binary(diamond), database, queries, hits, CATALOG_TARGETS, threads
        )
    else:
        hits.write_text("", encoding="utf-8")
    ref = pd.read_csv(
        resources.path(layout.REFX_GENOMES_FILE), sep="\t", dtype={"lineage": str}
    )
    lineage = dict(zip(ref["index"], ref["lineage"].fillna(""), strict=True))
    return gene_table(best_catalog_hits(hits, ref), lineage)
