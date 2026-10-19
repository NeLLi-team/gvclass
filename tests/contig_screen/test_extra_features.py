"""Tests of the protein sample and the window embeddings (no model needed)."""

from __future__ import annotations

import random

import numpy as np
import pandas as pd

from src.contig_screen import extra_features as ef
from src.contig_screen.extract_genes import SEED


def protein(seed: int, length: int = 80) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(length))


def gene(contig: str, index: int, sequence: str) -> dict[str, object]:
    return {"contig_id": contig, "gene_idx": index, "protein": sequence}


class FakeEmbedder:
    """Embeds a protein as its length repeated; records what it was asked."""

    def __init__(self) -> None:
        self.asked: list[list[str]] = []

    def embed(self, proteins: dict[str, str]) -> dict[str, np.ndarray]:
        self.asked.append(list(proteins))
        return {
            h: np.full(3, float(len(p)), dtype=np.float32) for h, p in proteins.items()
        }


def test_contig_proteins_keep_long_proteins_without_unknown_residues_once() -> None:
    genes = [
        gene("w", 0, protein(1)),
        gene("w", 1, protein(1)),
        gene("w", 2, protein(2, 49)),
        gene("w", 3, "A" * 30 + "X" + "A" * 30),
    ]

    found = ef.contig_proteins(genes, 8, random.Random(0))

    assert list(found.values()) == [protein(1)]


def test_contig_proteins_sample_at_most_the_limit_from_the_sorted_hashes() -> None:
    genes = [gene("w", i, protein(i)) for i in range(20)]

    first = ef.contig_proteins(genes, 8, random.Random(5))
    again = ef.contig_proteins(genes, 8, random.Random(5))
    other = ef.contig_proteins(genes, 8, random.Random(6))

    assert len(first) == 8
    assert list(first) == list(again)
    assert list(first) != list(other)


def test_window_embeddings_average_the_embeddings_of_the_sampled_proteins() -> None:
    contigs = [{"contig_id": "w1"}, {"contig_id": "w2"}]
    genes = [gene("w1", 0, protein(1, 60)), gene("w1", 1, protein(2, 100))]

    out = ef.window_embeddings("bin", contigs, genes, FakeEmbedder(), 8)

    assert out["contig_id"].tolist() == ["w1", "w2"]
    assert out.loc[0, ["e_0", "e_1", "e_2"]].tolist() == [80.0, 80.0, 80.0]
    assert out.loc[1, ["e_0", "e_1", "e_2"]].isna().all()


def test_window_embeddings_sample_depends_on_the_bin_name_and_window_order() -> None:
    contigs = [{"contig_id": "w1"}, {"contig_id": "w2"}]
    genes = [gene("w1", i, protein(i)) for i in range(20)]
    genes += [gene("w2", i, protein(100 + i)) for i in range(20)]

    def sampled(name: str, order: list[dict[str, object]]) -> set[str]:
        embedder = FakeEmbedder()
        ef.window_embeddings(name, order, genes, embedder, 8)
        return set(embedder.asked[0])

    assert sampled("bin", contigs) == sampled("bin", contigs)
    assert sampled("bin", contigs) != sampled("renamed", contigs)
    assert sampled("bin", contigs) != sampled("bin", contigs[::-1])


def test_the_sample_is_seeded_as_the_prototype_seeds_it() -> None:
    genes = [gene("w1", i, protein(i)) for i in range(20)]
    expected = ef.contig_proteins(genes, 8, random.Random(f"{SEED}|esm|bin"))
    embedder = FakeEmbedder()

    ef.window_embeddings("bin", [{"contig_id": "w1"}], genes, embedder, 8)

    assert embedder.asked[0] == list(expected)


def test_window_embeddings_of_a_bin_without_proteins_have_only_the_window_ids() -> None:
    out = ef.window_embeddings("bin", [{"contig_id": "w1"}], [], FakeEmbedder(), 8)

    assert isinstance(out, pd.DataFrame)
    assert out.columns.tolist() == ["contig_id"]


def toy_bin(names: list[str], seqs: list[str], genome: str) -> tuple[dict, list]:
    """Return contig name to sequence and the 10 kb windows (step 5 kb) of a toy bin."""
    sequences = dict(zip(names, seqs, strict=True))
    windows = [
        {
            "contig_id": f"{genome}|{name}|{start}",
            "source_contig": name,
            "source_start": start,
        }
        for name, seq in sequences.items()
        for start in range(0, len(seq), 5000)
    ]
    return sequences, windows


class HashEmbedder:
    """Embeds a protein as a number that depends on its content; records requests."""

    def __init__(self) -> None:
        self.asked: list[list[str]] = []

    def embed(self, proteins: dict[str, str]) -> dict[str, np.ndarray]:
        self.asked.append(list(proteins))
        return {h: np.full(2, int(h[:6], 16) / 1e6, dtype=np.float32) for h in proteins}


def window_genes(
    windows: list[dict[str, object]], sequences: dict[str, str]
) -> list[dict[str, object]]:
    """Return 20 proteins per window, determined by the content of the window."""
    genes = []
    for window in windows:
        base = ef.sequence_hash(sequences[str(window["source_contig"])])
        key = f"{base}|{window['source_start']}"
        genes += [
            gene(
                str(window["contig_id"]),
                i,
                protein(int(ef.sequence_hash(f"{key}|{i}")[:8], 16)),
            )
            for i in range(20)
        ]
    return genes


SEQS = ["ACGT" * 4000, "TTGCA" * 3200, "GGATC" * 3300]


def sampled_by_content(
    names: list[str], seqs: list[str], genome: str, sample_seed: str | None
) -> dict[tuple[str, int], set[str]]:
    """Return the sampled protein hashes of each window, keyed by content."""
    sequences, windows = toy_bin(names, seqs, genome)
    by_name = {n: ef.sequence_hash(s) for n, s in sequences.items()}
    seed, ordered = ef.draw_plan(windows, sequences, sample_seed=sample_seed)
    genes = window_genes(windows, sequences)
    out = ef.window_embeddings(seed, ordered, genes, HashEmbedder(), 8)
    kept = {}
    for window in windows:
        row = out[out["contig_id"] == window["contig_id"]].iloc[0]
        kept[(by_name[str(window["source_contig"])], int(window["source_start"]))] = {
            round(float(row["e_0"]), 6)
        }
    return kept


def test_renaming_the_file_and_its_contigs_does_not_change_the_sample() -> None:
    first = sampled_by_content(["c1", "c2", "c3"], SEQS, "binA", None)
    second = sampled_by_content(["x", "y", "z"], SEQS, "binB", None)

    assert first == second


def test_reversing_the_contig_order_does_not_change_the_sample() -> None:
    first = sampled_by_content(["c1", "c2", "c3"], SEQS, "bin", None)
    second = sampled_by_content(["c3", "c2", "c1"], SEQS[::-1], "bin", None)

    assert first == second


def test_two_different_bins_have_different_seeds() -> None:
    a, _ = toy_bin(["c1", "c2"], SEQS[:2], "bin")
    b, _ = toy_bin(["c1", "c2"], [SEQS[0], SEQS[2]], "bin")

    assert ef.content_seed(a) != ef.content_seed(b)


def test_the_content_seed_ignores_contig_names_and_order() -> None:
    assert ef.content_seed({"a": "AC", "b": "GT"}) == ef.content_seed(
        {"q": "GT", "p": "AC"}
    )


def test_a_given_sample_seed_reproduces_the_prototype_draw() -> None:
    sequences, windows = toy_bin(["c1", "c2"], SEQS[:2], "bin")
    genes = window_genes(windows, sequences)

    seed, ordered = ef.draw_plan(windows, sequences, sample_seed="NAME")
    embedder = HashEmbedder()
    ef.window_embeddings(seed, ordered, genes, embedder, 8)

    rng = random.Random(f"{SEED}|esm|NAME")
    expected: dict[str, str] = {}
    for window in windows:
        mine = [g for g in genes if g["contig_id"] == window["contig_id"]]
        expected.update(ef.contig_proteins(mine, 8, rng))
    assert seed == "NAME"
    assert ordered == windows
    assert set(embedder.asked[0]) == set(expected)
