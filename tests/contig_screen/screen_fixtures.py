"""Small deterministic bins whose features the screen models were trained on.

The expected values in ``data/feature_baseline.npz`` were captured from the code
before the order-independence and efficiency changes. Input without ties and without
duplicate sequences must keep producing them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.contig_screen import context_stage as cs, extra_features as ef
from src.contig_screen.gv_order import (
    N_CLASSES,
    Contig,
    IdentityCache,
    Priors,
    bin_features,
)

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"
AMINO = "ACDEFGHIKLMNPQRSTVWY"
COMP = ["c_gc", "c_aa_A", "e_0", "e_1"]
SEQUENCES = {"c1": "ACGT" * 4000, "c2": "TTGCA" * 3200, "c3": "GGATC" * 3300}


def context_windows() -> pd.DataFrame:
    """Return 14 windows of two bins with distinct scores and positions."""
    rng = np.random.default_rng(11)
    n = 14
    frame = pd.DataFrame(
        {
            "bin_id": ["a"] * 8 + ["b"] * 6,
            "p": rng.uniform(0.02, 0.98, n),
            "length_bp": rng.integers(8_000, 10_001, n),
            "c_gc": rng.uniform(0.3, 0.7, n),
            "c_aa_A": rng.uniform(0.05, 0.15, n),
            "e_0": rng.normal(size=n),
            "e_1": rng.normal(size=n),
        }
    )
    for name in cs.OWN:
        frame[name] = rng.uniform(0.0, 1.0, n)
    return frame


def context_features(windows: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return the stage 2 context features of ``context_windows``."""
    frame = context_windows() if windows is None else windows
    mu = np.array([0.5, 0.1, 0.0, 0.0])
    sd = np.array([0.1, 0.03, 1.0, 1.0])
    return cs.context_table(frame, COMP, mu, sd)


def order_contigs() -> list[Contig]:
    """Return five contigs whose marker partners and distances have no ties."""
    rng = np.random.default_rng(5)
    layout = {
        "a": (["OG1", "OG2", "OG3", "OG4", "OG5"], 120_000),
        "b": (["OG1", "OG2", "OG3"], 80_000),
        "c": (["OG4", "OG6"], 30_000),
        "d": (["OG5", "OG7", "OG8"], 21_000),
        "e": (["OG9"], 9_000),
    }
    contigs = []
    for name, (models, length) in layout.items():
        contig = Contig(
            name=name,
            genome="bin",
            length=length,
            models=models,
            proteins=[PROTEIN] * len(models),
            hashes=[f"h-{model}" for model in models],
        )
        contig.comp = np.append(
            rng.dirichlet(np.ones(N_CLASSES)), rng.uniform(0.3, 0.7)
        )
        contig.codon = rng.random(65)
        contig.struct = rng.random(4)
        contigs.append(contig)
    return contigs


def order_priors() -> Priors:
    """Return priors with distinct weights per marker."""
    return Priors(
        p_multi={"OG1": 0.01, "OG2": 0.5, "OG3": 0.05, "OG4": 0.02},
        default=0.1,
        mu=np.full(N_CLASSES + 1, 0.02),
        sd=np.full(N_CLASSES + 1, 0.05),
        extra={
            "codon": (np.full(65, 0.1), np.full(65, 0.2)),
            "struct": (np.full(4, 0.5), np.full(4, 0.4)),
        },
    )


def order_features(contigs: list[Contig] | None = None) -> pd.DataFrame:
    """Return the giant-virus order features of ``order_contigs``."""
    chosen = order_contigs() if contigs is None else contigs
    return bin_features(chosen, order_priors(), IdentityCache())


def protein(seed: int, length: int = 80) -> str:
    """Return a reproducible protein."""
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list(AMINO), size=length))


class HashEmbedder:
    """Embeds a protein as numbers that depend on its content; records requests."""

    def __init__(self) -> None:
        self.asked: list[list[str]] = []

    def embed(self, proteins: dict[str, str]) -> dict[str, np.ndarray]:
        self.asked.append(list(proteins))
        return {
            h: np.array([int(h[:6], 16) / 1e6, int(h[6:12], 16) / 1e6])
            for h in proteins
        }


def toy_windows(
    sequences: dict[str, str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Return the 10 kb windows (step 5 kb) of the bin and 20 proteins per window.

    The proteins depend only on the content of a window (sequence hash and start).
    """
    windows: list[dict[str, object]] = []
    genes: list[dict[str, object]] = []
    for name, sequence in sequences.items():
        base = ef.sequence_hash(sequence)
        for start in range(0, len(sequence), 5000):
            window = {
                "contig_id": f"bin|{name}|{start}",
                "source_contig": name,
                "source_start": start,
            }
            windows.append(window)
            genes += [
                {
                    "contig_id": window["contig_id"],
                    "gene_idx": i,
                    "protein": protein(
                        int(ef.sequence_hash(f"{base}|{start}|{i}")[:8], 16)
                    ),
                }
                for i in range(20)
            ]
    return windows, genes


def embedding_features(sequences: dict[str, str] | None = None) -> pd.DataFrame:
    """Return the window embeddings of a toy bin (indexed by window id)."""
    bin_sequences = SEQUENCES if sequences is None else sequences
    windows, genes = toy_windows(bin_sequences)
    seed, ordered = ef.draw_plan(windows, bin_sequences)
    out = ef.window_embeddings(seed, ordered, genes, HashEmbedder(), 8)
    return out.set_index("contig_id").sort_index()


def call_windows() -> tuple[pd.DataFrame, dict[str, int]]:
    """Return overlapping window scores for three contigs and their lengths."""
    rng = np.random.default_rng(3)
    lengths = {"x": 52_300, "y": 31_000, "z": 17_450}
    rows = []
    for contig, size in lengths.items():
        for start in range(0, size, 5000):
            piece = min(10_000, size - start)
            if piece >= 3000:
                rows.append((contig, start, piece, rng.uniform()))
    frame = pd.DataFrame(
        rows, columns=["source_contig", "source_start", "length_bp", "p"]
    )
    return frame, lengths
