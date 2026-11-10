"""Gene files come from pyrodigal's own writers and agree with each other.

A raw forward slice of the contig returns the sense strand for minus-strand
genes, so the nucleotide records must be strand-aware to translate to the
protein records. Protein, gene and GFF identifiers are ``<contig id>_<n>``.
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from src.core.genetic_code_optimizer import GeneCallingResult, optimize_genetic_code

MODEL_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLQDKPEAIRRLVDELAQLGVE"
FORWARD_ORF = "ATG" + "GCT" * 299 + "TAA"


@pytest.fixture
def tiny_hmm(tmp_path: Path) -> Path:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(
        name=b"tiny_model", sequence=MODEL_PROTEIN
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence, pyhmmer.plan7.Background(alphabet)
    )
    hmm_path = tmp_path / "tiny.hmm"
    with hmm_path.open("wb") as handle:
        hmm.write(handle)
    return hmm_path


def _call(tmp_path: Path, hmm_path: Path, contigs: dict[str, str]) -> GeneCallingResult:
    query = tmp_path / "q.fna"
    query.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in contigs.items()),
        encoding="utf-8",
    )
    return optimize_genetic_code(query, tmp_path / "out", hmm_path, [0], threads=1)


def _read_fasta(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(path), "fasta")}


def _gff_rows(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if "\tCDS\t" in line
    ]


def test_minus_strand_fna_translates_to_faa(tmp_path: Path, tiny_hmm: Path) -> None:
    contig = str(Seq(FORWARD_ORF).reverse_complement())

    result = _call(tmp_path, tiny_hmm, {"q|g": contig})

    faa = _read_fasta(result.faa)
    fna = _read_fasta(result.fna)
    assert faa, "pyrodigal found no genes in the synthetic ORF"
    assert [row[6] for row in _gff_rows(result.gff)] == ["-"] * len(faa)
    assert set(faa) == set(fna)
    for protein_id, protein in faa.items():
        translated = str(Seq(fna[protein_id]).translate(table=11)).rstrip("*")
        # Alternative start codons are translated to M by pyrodigal.
        assert translated[1:] == protein.rstrip("*")[1:]


def test_gene_identifiers_are_contig_id_and_gene_number(
    tmp_path: Path, tiny_hmm: Path
) -> None:
    contigs = {"q_reformatted|c1": FORWARD_ORF, "q_reformatted|c2": FORWARD_ORF}

    result = _call(tmp_path, tiny_hmm, contigs)

    expected = {"q_reformatted|c1_1", "q_reformatted|c2_1"}
    assert set(_read_fasta(result.faa)) == expected
    assert set(_read_fasta(result.fna)) == expected
    assert {row[0] for row in _gff_rows(result.gff)} == {
        "q_reformatted|c1",
        "q_reformatted|c2",
    }


def test_gff_uses_pyrodigal_confidence_attribute(
    tmp_path: Path, tiny_hmm: Path
) -> None:
    result = _call(tmp_path, tiny_hmm, {"q|g": FORWARD_ORF})

    (row,) = _gff_rows(result.gff)

    assert result.gff.read_text(encoding="utf-8").startswith("##gff-version")
    assert "conf=" in row[8]
    assert "rbs_spacer=" in row[8]
