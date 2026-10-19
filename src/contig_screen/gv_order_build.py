"""Contigs of a giant-virus bin with their marker genes, composition and hits.

Joins the gene table of ``gv_genes`` (and, when the model reads them, the best catalog
hits of ``gv_taxa``) to the sequences of the bin and returns one ``gv_order.Contig``
per contig of 5 kb or more. The synthetic-bin builder of the prototype is not ported.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from src.contig_screen.gv_order import (
    CODON_COLUMNS,
    MIN_FRAGMENT_BP,
    Contig,
    composition,
    contig_from_genes,
    gene_codon_counts,
)

TAXON_COLUMNS = ("tx_class", "tx_family", "tx_genus", "tx_pident", "tx_bits")
NUMERIC_TAXON_COLUMNS = ("tx_pident", "tx_bits")


def add_taxa(genes: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    """Return the gene table with the best catalog hit of every gene.

    ``taxa`` has one row per gene with a hit (``row``, the position in ``genes``, and
    the ``TAXON_COLUMNS``); genes without a hit get empty codes and missing numbers.
    """
    found = taxa.set_index("row").reindex(range(len(genes)))
    genes = genes.copy()
    for column in TAXON_COLUMNS:
        values = found[column].to_numpy()
        if column in NUMERIC_TAXON_COLUMNS:
            genes[column] = values.astype(float)
        else:
            genes[column] = pd.Series(values).fillna("").to_numpy()
    return genes


class Genome:
    """Sequence and marker genes of one bin."""

    def __init__(
        self,
        genome_id: str,
        sequences: Mapping[str, str],
        genes: pd.DataFrame,
        taxa: pd.DataFrame | None = None,
    ) -> None:
        """Count the codons of every gene (with hit taxonomy when ``taxa``)."""
        self.genome_id = genome_id
        self.sequences = dict(sequences)
        genes = genes.copy() if taxa is None else add_taxa(genes, taxa)
        counts = np.zeros((len(genes), 64), dtype=np.int32)
        for row, (name, begin, end, strand) in enumerate(
            zip(
                genes["contig"],
                genes["begin"],
                genes["end"],
                genes["strand"],
                strict=True,
            )
        ):
            if name in self.sequences:
                counts[row] = gene_codon_counts(
                    self.sequences[name], int(begin), int(end), int(strand)
                )
        genes[CODON_COLUMNS] = counts
        self.genes = dict(tuple(genes.groupby("contig")))

    def contigs(self) -> list[Contig]:
        """Return the whole contigs of 5 kb or more with composition."""
        out = []
        for name, sequence in self.sequences.items():
            if name not in self.genes or len(sequence) < MIN_FRAGMENT_BP:
                continue
            contig = contig_from_genes(
                self.genes[name], name, self.genome_id, 0, len(sequence)
            )
            contig.comp = composition(sequence)
            out.append(contig)
        return out
