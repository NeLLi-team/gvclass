"""Tests of the windows, gene calls and searches of a bin."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

import pandas as pd
import pyhmmer
import pytest

from src.contig_screen import extract_genes as eg, layout
from tests.contig_screen.conftest import BundleFactory

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"


def random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def coding_dna(n_genes: int, seed: int) -> str:
    """Return a sequence of open reading frames separated by short spacers."""
    rng = random.Random(seed)
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
    sense = [c for c in codons if c not in {"TAA", "TAG", "TGA", "ATG"}]
    parts = []
    for _ in range(n_genes):
        body = "".join(rng.choice(sense) for _ in range(300))
        parts.append(
            "AGGAGG" + "ATATA" + "ATG" + body + "TAA" + random_dna(40, rng.random())
        )
    return "".join(parts)


def write_fasta(path: Path, contigs: dict[str, str]) -> dict[str, str]:
    path.write_text("".join(f">{n}\n{s}\n" for n, s in contigs.items()))
    return {"genome_id": "bin", "class": "BIN", "fasta_path": str(path)}


def test_tile_contigs_cuts_overlapping_windows_and_keeps_a_long_enough_tail(
    tmp_path: Path,
) -> None:
    genome = write_fasta(tmp_path / "b.fna", {"c1": random_dna(25_000, 1)})

    windows = eg.tile_contigs(genome, 10, 5)

    assert [w["source_start"] for w in windows] == [0, 5000, 10_000, 15_000, 20_000]
    assert [w["contig_id"] for w in windows][:2] == ["bin|c1|0", "bin|c1|5000"]
    assert [len(str(w["sequence"])) for w in windows] == [10_000] * 4 + [5000]
    assert windows[-1]["length_kb"] == 5


def test_tile_contigs_without_a_step_does_not_overlap(tmp_path: Path) -> None:
    genome = write_fasta(tmp_path / "b.fna", {"c1": random_dna(25_000, 1)})

    windows = eg.tile_contigs(genome, 10, None)

    assert [w["source_start"] for w in windows] == [0, 10_000, 20_000]


def test_tile_contigs_with_tile_zero_keeps_the_contigs(tmp_path: Path) -> None:
    genome = write_fasta(
        tmp_path / "b.fna", {"c1": random_dna(4000, 1), "short": random_dna(2999, 2)}
    )

    windows = eg.tile_contigs(genome, 0)

    assert [w["contig_id"] for w in windows] == ["bin|c1|0"]


def test_tile_contigs_skips_windows_that_are_mostly_unknown_bases(
    tmp_path: Path,
) -> None:
    sequence = random_dna(10_000, 1) + "N" * 10_000
    genome = write_fasta(tmp_path / "b.fna", {"c1": sequence})

    windows = eg.tile_contigs(genome, 10, None)

    assert [w["source_start"] for w in windows] == [0]


def test_read_fasta_reads_gzipped_files_in_upper_case(tmp_path: Path) -> None:
    path = tmp_path / "b.fna.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">c1 description\n" + "acgt" * 1000 + "\n")

    assert eg.read_fasta(str(path)) == {"c1": "ACGT" * 1000}


def test_call_genes_summarizes_every_window_and_numbers_its_genes(
    tmp_path: Path,
) -> None:
    genome = write_fasta(tmp_path / "b.fna", {"c1": coding_dna(8, 3)})
    windows = eg.tile_contigs(genome, 0)

    genes, rows = eg.call_genes(genome, windows)

    assert len(rows) == 1
    row = rows[0]
    assert row["n_genes"] == len(genes) > 0
    assert [g["gene_idx"] for g in genes] == list(range(len(genes)))
    assert 0 < int(str(row["coding_bp"])) <= int(str(row["length_bp"]))
    assert all(f"cod_{codon}" in row for codon in eg.CODONS)


def test_unique_proteins_keep_long_proteins_once() -> None:
    genes = [
        {"aa_len": 40, "prot_hash": "h1", "protein": "A" * 40},
        {"aa_len": 40, "prot_hash": "h1", "protein": "A" * 40},
        {"aa_len": 29, "prot_hash": "h2", "protein": "C" * 29},
    ]

    assert eg.unique_proteins(genes) == {"h1": "A" * 40}


def test_write_tables_adds_the_best_model_of_each_gene(tmp_path: Path) -> None:
    genes = [
        {"contig_id": "w", "prot_hash": "h1", "protein": "A", "gene_idx": 0},
        {"contig_id": "w", "prot_hash": "h2", "protein": "C", "gene_idx": 1},
    ]
    models = {"h1": [("M1", 10.0, 1e-9), ("M2", 30.0, 1e-20), ("M2", 5.0, 1e-3)]}

    eg.write_tables(tmp_path, "b", [{"contig_id": "w"}], genes, models, [])

    table = pd.read_parquet(tmp_path / "b.genes.parquet")
    assert table["best_model"].tolist() == ["M2", None]
    assert table["best_hmm_score"].iloc[0] == 30.0
    assert table["n_models"].tolist() == [2, 0]
    assert "protein" not in table
    assert pd.read_parquet(tmp_path / "b.hits.parquet").columns.tolist() == (
        eg.HIT_COLUMNS
    )


def write_catalog(path: Path) -> None:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(name=b"q", sequence=PROTEIN).digitize(
        alphabet
    )
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence, pyhmmer.plan7.Background(alphabet)
    )
    hmm.name = b"OG0001"
    with path.open("wb") as handle:
        hmm.write(handle)


def test_hmm_search_returns_the_models_that_hit_each_protein(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.hmm"
    write_catalog(catalog)

    found = eg.hmm_search({"h1": PROTEIN, "h2": "W" * 60}, catalog, 1)

    assert list(found) == ["h1"]
    assert found["h1"][0][0] == "OG0001"


def test_hmm_search_of_no_protein_is_empty(tmp_path: Path) -> None:
    assert eg.hmm_search({}, tmp_path / "missing.hmm", 1) == {}


def test_reference_search_finds_the_subjects_of_the_marker_group(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    refs = layout.MARKER_REFS_DIR
    resources = make_resources(
        tmp_path,
        {
            layout.MODEL_TO_GROUP_FILE: b'{"OG0001": "GRP"}',
            f"{refs}/GRP.faa": f">sub1\n{PROTEIN}\n>sub2\n{PROTEIN[::-1]}\n".encode(),
        },
    )
    models = {"h1": [("OG0001", 50.0, 1e-9)], "h2": [("OTHER", 50.0, 1e-9)]}

    rows = eg.reference_search({"h1": PROTEIN, "h2": PROTEIN}, models, resources, 1)

    assert {r["prot_hash"] for r in rows} == {"h1"}
    best = rows[0]
    assert (best["marker"], best["rank"], best["subject"]) == ("GRP", 1, "sub1")
    assert best["identity"] == pytest.approx(100.0)


def test_extract_bin_refuses_a_bin_without_a_contig_of_3_kb(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    genome = write_fasta(tmp_path / "b.fna", {"c1": random_dna(2000, 1)})
    resources = make_resources(tmp_path / "bundle", {"x.txt": b"x"})

    with pytest.raises(ValueError, match="no contig of at least 3 kb"):
        eg.extract_bin(genome, resources, tmp_path, tile_kb=10, step_kb=5, threads=1)
