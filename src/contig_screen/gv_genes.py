"""Whole-bin gene calls with catalog model assignments for the giant-virus screen.

Calls genes on all contigs of 1 kb or more (pyrodigal, metagenomic mode, as
``extract_genes``) and searches the HMM catalog with pyhmmer (same E-value). The
result is one table row per gene: contig, order, position, strand, protein hash, best
model, score, and the protein sequence. The prototype wrote this table to a parquet
file; here it stays in memory.
"""

from __future__ import annotations

import gzip
import logging
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyrodigal

from src.contig_screen.extract_genes import gene_row, hmm_search, unique_proteins

logger = logging.getLogger(__name__)

MIN_CONTIG_BP = 1000
GENE_COLUMNS = [
    "contig",
    "gene_idx",
    "begin",
    "end",
    "strand",
    "aa_len",
    "prot_hash",
    "model",
    "score",
    "evalue",
    "protein",
]


@dataclass(frozen=True, slots=True)
class GeneCalls:
    """The gene table of a bin and the protein of every row (also without a model)."""

    table: pd.DataFrame
    proteins: list[str]


def read_fasta(path: Path) -> dict[str, str]:
    """Return contig name to upper-case sequence (plain or gzip)."""
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, list[str]] = {}
    name = ""
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name:
                records[name].append(line.strip().upper())
    return {key: "".join(parts) for key, parts in records.items()}


def call_genes(sequences: Mapping[str, str], hmm_path: Path, threads: int) -> GeneCalls:
    """Return the gene table of a bin, with the best catalog model of each gene.

    Genes without a model have an empty ``model``, score 0 and E-value 1.
    """
    started = time.monotonic()
    finder = pyrodigal.GeneFinder(meta=True)
    gene_rows: list[dict[str, object]] = []
    for name, sequence in sequences.items():
        if len(sequence) < MIN_CONTIG_BP:
            continue
        genes = finder.find_genes(sequence.encode("ascii"))
        gene_rows.extend(
            gene_row(name, index, gene) for index, gene in enumerate(genes)
        )
    hits = hmm_search(unique_proteins(gene_rows), hmm_path, threads)
    best = {
        key: max(found, key=lambda hit: hit[1]) for key, found in hits.items() if found
    }
    rows = [table_row(row, best.get(str(row["prot_hash"]))) for row in gene_rows]
    table = pd.DataFrame(rows, columns=GENE_COLUMNS)
    logger.info(
        "%d genes, %d with a model, %.0f s",
        len(table),
        int((table["model"] != "").sum()),
        time.monotonic() - started,
    )
    return GeneCalls(table, [str(row["protein"]) for row in gene_rows])


def table_row(
    gene: Mapping[str, object], model: tuple[str, float, float] | None
) -> dict[str, object]:
    """Return the table row of one gene and its best model (or None)."""
    return {
        "contig": gene["contig_id"],
        "gene_idx": gene["gene_idx"],
        "begin": gene["begin"],
        "end": gene["end"],
        "strand": gene["strand"],
        "aa_len": gene["aa_len"],
        "prot_hash": gene["prot_hash"],
        "model": model[0] if model else "",
        "score": model[1] if model else 0.0,
        "evalue": model[2] if model else 1.0,
        "protein": gene["protein"] if model else "",
    }
