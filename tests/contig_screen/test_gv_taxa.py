"""Tests of the best catalog hit per gene of the giant-virus screen."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.contig_screen import gv_taxa
from src.contig_screen.gv_order_build import add_taxa

PHYLUM = "d__Viruses;p__Nucleocytoviricota;c__Megaviricetes"
LINEAGE = {
    0: f"{PHYLUM};o__Imitervirales;f__Mimiviridae;g__Cafeteriavirus",
    1: "d__Bacteria",
    2: f"{PHYLUM};o__Algavirales;f__Phycodnaviridae;g__Chlorovirus",
}
REF = pd.DataFrame(
    {
        "index": [0, 1, 2, 3],
        "class": ["NCLDV", "BAC", "NCLDV", ""],
        "lineage": [*LINEAGE.values(), ""],
    }
)


def write_hits(path: Path, rows: list[tuple[str, str, float, float]]) -> None:
    lines = [
        f"{query}\t{subject}\t{pident}\t100\t200\t300\t1e-9\t{bits}"
        for query, subject, pident, bits in rows
    ]
    path.write_text("\n".join(lines) + "\n")


def test_codes_are_the_family_and_genus_without_rank_prefix() -> None:
    assert gv_taxa.codes(LINEAGE[0]) == ("Mimiviridae", "Cafeteriavirus")
    assert gv_taxa.codes("d__Viruses") == ("", "")


def test_fasta_text_numbers_the_rows_and_skips_empty_proteins() -> None:
    assert gv_taxa.fasta_text("bin", ["MKT", "", "MAA"]) == ">bin|0\nMKT\n>bin|2\nMAA\n"


def test_the_best_hit_per_gene_and_class_wins_and_unclassed_references_are_dropped(
    tmp_path: Path,
) -> None:
    hits = tmp_path / "hits.tsv"
    write_hits(
        hits,
        [
            ("bin|0", "r0_a", 40.0, 80.0),
            ("bin|0", "r2_a", 60.0, 120.0),
            ("bin|0", "r3_a", 99.0, 400.0),
        ],
    )

    best = gv_taxa.best_catalog_hits(hits, REF)

    assert best["ref"].tolist() == [2]
    assert best["bitscore"].tolist() == [120.0]


def test_an_empty_hit_file_gives_an_empty_table(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text("")

    assert gv_taxa.best_catalog_hits(hits, REF).empty


def test_gene_table_takes_the_best_class_and_the_ncldv_codes(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    write_hits(
        hits,
        [
            ("bin|0", "r0_a", 40.0, 80.0),
            ("bin|0", "r1_a", 90.0, 300.0),
            ("bin|3", "r2_a", 60.0, 120.0),
        ],
    )

    table = gv_taxa.gene_table(
        gv_taxa.best_catalog_hits(hits, REF),
        dict(zip(REF["index"], REF["lineage"], strict=True)),
    ).set_index("row")

    assert table.loc[0, "tx_class"] == "BAC"
    assert table.loc[0, "tx_family"] == "Mimiviridae"
    assert table.loc[0, "tx_pident"] == pytest.approx(40.0)
    assert table.loc[3, "tx_class"] == "NCLDV"
    assert table.loc[3, "tx_genus"] == "Chlorovirus"


def test_add_taxa_fills_genes_without_a_hit_with_empty_codes_and_missing_numbers() -> (
    None
):
    genes = pd.DataFrame({"contig": ["a", "a", "b"]})
    taxa = pd.DataFrame(
        {
            "row": [1],
            "tx_class": ["NCLDV"],
            "tx_family": ["Mimiviridae"],
            "tx_genus": [""],
            "tx_pident": [55.0],
            "tx_bits": [90.0],
        }
    )

    joined = add_taxa(genes, taxa)

    assert joined["tx_class"].tolist() == ["", "NCLDV", ""]
    assert joined["tx_pident"].isna().tolist() == [True, False, True]
    assert "tx_class" not in genes
