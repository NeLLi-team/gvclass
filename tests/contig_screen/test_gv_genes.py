"""Tests of the gene calls and the catalog assignment of the giant-virus screen."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest

from src.contig_screen import gv_genes

STOPS = {"TAA", "TAG", "TGA"}


def orf_sequence(seed: int = 1, codons: int = 400) -> str:
    """Return a random open reading frame with flanks, about 1.4 kb long."""
    rng = np.random.default_rng(seed)
    triplets = ["".join(t) for t in rng.choice(list("ACGT"), size=(4 * codons, 3))]
    body = [t for t in triplets if t not in STOPS][:codons]
    flank = "".join(rng.choice(list("ACGT"), size=120))
    return flank + "ATG" + "".join(body) + "TAA" + flank


def test_read_fasta_upper_cases_and_joins_lines_in_plain_and_gzip_files(
    tmp_path: Path,
) -> None:
    text = ">c1 desc\nacgt\nAC\n>c2\nGG\n"
    plain, packed = tmp_path / "bin.fna", tmp_path / "bin.fna.gz"
    plain.write_text(text)
    with gzip.open(packed, "wt") as handle:
        handle.write(text)

    assert gv_genes.read_fasta(plain) == {"c1": "ACGTAC", "c2": "GG"}
    assert gv_genes.read_fasta(packed) == gv_genes.read_fasta(plain)


def test_call_genes_gives_one_row_per_gene_and_a_model_for_the_searched_protein(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    searched: list[dict[str, str]] = []

    def fake_search(
        proteins: dict[str, str], hmm_path: Path, threads: int
    ) -> dict[str, list[tuple[str, float, float]]]:
        searched.append(proteins)
        return {key: [("OG9", 30.0, 1e-9), ("OG3", 55.0, 1e-20)] for key in proteins}

    monkeypatch.setattr(gv_genes, "hmm_search", fake_search)
    sequences = {"long": orf_sequence(), "short": "ACGT" * 100}

    calls = gv_genes.call_genes(sequences, tmp_path / "x.hmm", 1)

    table = calls.table
    assert list(table.columns) == gv_genes.GENE_COLUMNS
    assert len(table) > 0
    assert set(table["contig"]) == {"long"}
    assert (table["model"] == "OG3").all()
    assert (table["score"] == 55.0).all()
    assert len(calls.proteins) == len(table)
    assert all(table["protein"] == calls.proteins)


def test_genes_without_a_model_keep_the_protein_only_in_the_protein_list(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(gv_genes, "hmm_search", lambda *args: {})

    calls = gv_genes.call_genes({"long": orf_sequence()}, tmp_path / "x.hmm", 1)

    assert (calls.table["model"] == "").all()
    assert (calls.table["protein"] == "").all()
    assert (calls.table["evalue"] == 1.0).all()
    assert all(len(p) > 0 for p in calls.proteins)


def test_call_genes_skips_contigs_shorter_than_one_kb(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(gv_genes, "hmm_search", lambda *args: {})

    calls = gv_genes.call_genes({"tiny": orf_sequence()[:900]}, tmp_path / "x.hmm", 1)

    assert calls.table.empty
    assert calls.proteins == []
