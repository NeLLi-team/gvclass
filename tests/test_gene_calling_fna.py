"""M7: minus-strand gene .fna records must be reverse-complemented.

A raw forward slice of the contig returns the sense strand for minus-strand
genes, so the .fna would not translate to the corresponding .faa protein. The
optimizer must use pyrodigal's strand-aware ``gene.sequence()``.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.core.genetic_code_optimizer import GeneticCodeOptimizer


def _read_fasta(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(path), "fasta")}


def test_minus_strand_fna_translates_to_faa(tmp_path: Path) -> None:
    # A clean forward ORF (M + 299x Ala + stop), then reverse-complemented so
    # pyrodigal calls it on the MINUS strand of the contig.
    forward_orf = "ATG" + "GCT" * 299 + "TAA"
    contig = str(Seq(forward_orf).reverse_complement())
    records = [SeqRecord(Seq(contig), id="g", description="")]

    opt = GeneticCodeOptimizer(Path("/tmp"))
    faa_file, fna_file, gff_file, proteins, genes = opt.run_gene_calling_for_code(
        Path("unused.fna"), 0, tmp_path, records
    )
    assert proteins, "pyrodigal found no genes in the synthetic ORF"

    # At least one gene must be on the minus strand (else the test is vacuous).
    strands = [
        line.split("\t")[6]
        for line in gff_file.read_text().splitlines()
        if "\tCDS\t" in line
    ]
    assert "-" in strands, f"no minus-strand gene was called (strands={strands})"

    faa = _read_fasta(faa_file)
    fna = _read_fasta(fna_file)
    assert set(faa) == set(fna)

    for protein_id, prot in faa.items():
        nuc = fna[protein_id]
        translated = str(Seq(nuc).translate(table=11)).rstrip("*")
        # Compare ignoring the start residue (alternative start codons are
        # translated to M by pyrodigal but literally by Bio.Seq).
        assert translated[1:] == prot[1:], (
            f"{protein_id}: fna does not translate to faa "
            f"(translated={translated[:10]}..., faa={prot[:10]}...)"
        )
