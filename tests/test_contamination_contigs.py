"""Contig identity in the contamination scorer."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.core.contamination_scoring import ContaminationScorer
from src.core.labels import LabelIndex
from src.core.reformat import reformat_sequences
from src.core.summarize.paths_stats import QueryEvidence


def test_single_contig_keeps_one_identity_for_nucleotide_and_protein_records(
    tmp_path: Path,
) -> None:
    fna = tmp_path / "parent__ctg_reformatted.fna"
    records = reformat_sequences(
        [SeqRecord(Seq("ACGT" * 6000), id="ctg")], "parent__ctg.fna"
    )
    SeqIO.write(records, fna, "fasta")
    query_ids = [f"parent__ctg|ctg_{index}" for index in range(1, 6)]
    faa = tmp_path / "parent__ctg_reformatted.faa"
    faa.write_text(
        "".join(f">{name}\n{'M' * 100}\n" for name in query_ids), encoding="utf-8"
    )
    blast = tmp_path / "blastp_out"
    blast.mkdir()
    rows = [
        f"{name}\tBAC__ref|p1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100"
        for name in query_ids
    ]
    (blast / "OG1.m8").write_text("\n".join(rows) + "\n", encoding="utf-8")
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.labels = LabelIndex.from_taxonomy(
        {"BAC__ref": ["BAC", "", "", "OrderA", "FamA"]}
    )
    tree_neighbors = {"OG1": {name: {"BAC__ref|p1": 0.1} for name in query_ids}}

    result = scorer.collect_contig_features(
        fna,
        faa,
        scorer.scan_blast_dir(blast),
        QueryEvidence({}, {}, {}),
        "",
        tree_neighbors,
    )

    assert result["contig_count"] == 1
    assert result["suspicious_bp_fraction"] == 100.0
    assert result["cellular_coherent_bp_fraction"] == 100.0


def test_protein_to_contig_id_only_strips_numeric_cds_suffix() -> None:
    """Only a trailing numeric CDS suffix is stripped; other suffixes stay."""
    assert ContaminationScorer.protein_to_contig_id("query|contig_1_12") == "contig_1"
    assert ContaminationScorer.protein_to_contig_id("query|contig_A") == "contig_A"
    assert ContaminationScorer.protein_to_contig_id("contig_alpha") == "contig_alpha"
