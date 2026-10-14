"""Verify custom gene-calling GFF records against Pyrodigal's native writer."""

from __future__ import annotations

import io
from pathlib import Path

import pyrodigal
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.core.genetic_code_optimizer import GeneticCodeOptimizer

EXAMPLE_GENOME = Path(__file__).parents[1] / "example" / "GVMAG-S-1096109-37.fna"


def _full_contig_id(record: SeqRecord) -> str:
    fields = record.id.split("|", maxsplit=1)
    return "|".join(field.replace("_reformatted", "") for field in fields)


def _native_gff_rows(records: list[SeqRecord]) -> list[list[str]]:
    finder = pyrodigal.GeneFinder(meta=True)
    rows: list[list[str]] = []
    for record in records:
        genes = finder.find_genes(str(record.seq))
        output = io.StringIO()
        genes.write_gff(output, sequence_id=_full_contig_id(record))
        rows.extend(
            line.split("\t")
            for line in output.getvalue().splitlines()
            if "\tCDS\t" in line
        )
    return rows


def _gff_rows(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if "\tCDS\t" in line
    ]


def _attribute(row: list[str], name: str) -> str:
    attributes = dict(
        field.split("=", maxsplit=1) for field in row[8].rstrip(";").split(";")
    )
    return attributes[name]


def _gff_signature(row: list[str]) -> tuple[str, str, str, str, str, str, str]:
    return (
        row[0],
        row[3],
        row[4],
        row[6],
        row[7],
        _attribute(row, "ID"),
        _attribute(row, "partial"),
    )


def test_custom_gff_matches_native_phase_for_complete_and_partial_genes(
    tmp_path: Path,
) -> None:
    records = list(SeqIO.parse(EXAMPLE_GENOME, "fasta"))
    optimizer = GeneticCodeOptimizer(Path("/unused"), threads=1)

    _, _, custom_gff, _, _ = optimizer.run_gene_calling_for_code(
        EXAMPLE_GENOME,
        0,
        tmp_path,
        records,
    )
    custom_rows = _gff_rows(custom_gff)
    native_rows = _native_gff_rows(records)

    assert [_gff_signature(row) for row in custom_rows] == [
        _gff_signature(row) for row in native_rows
    ]
    assert {(row[6], _attribute(row, "partial")) for row in custom_rows} >= {
        ("+", "00"),
        ("+", "01"),
        ("+", "10"),
        ("-", "00"),
        ("-", "01"),
        ("-", "10"),
    }
