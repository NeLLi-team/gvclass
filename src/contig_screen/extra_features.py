"""Mean ESM-2 protein embedding of the windows of a bin (the ``e_`` features).

Per window the output holds the mean ESM-2 embedding of up to ``esm_genes`` called
proteins. The proteins are those of the gene calls of ``extract_genes``; the sample
is seeded with a digest of the sequences of the bin and drawn window by window in
the order of the window content, so it does not depend on the file name or on the
order of the contigs (``draw_plan``). ``sample_seed`` reproduces the prototype, which
seeded with the file stem and drew in file order.
"""

from __future__ import annotations

import hashlib
import logging
import random
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.contig_screen.extract_genes import SEED

logger = logging.getLogger(__name__)

MIN_PROTEIN_AA = 50
MAX_PROTEIN_AA = 600


class Embedder:
    """Mean-pooled ESM-2 embeddings with a per-process cache keyed by protein hash."""

    def __init__(self, model_dir: Path, threads: int, pooling: str = "legacy") -> None:
        """Load the model from a local directory once; ``pooling``: legacy or mean.

        The directory is read without network access.
        """
        import torch  # noqa: PLC0415  # optional dependency, see torch_stack
        from transformers import AutoModel, AutoTokenizer  # noqa: PLC0415

        torch.set_num_threads(threads)
        self.torch = torch
        self.pooling = pooling
        self.tokenizer = AutoTokenizer.from_pretrained(
            str(model_dir), local_files_only=True
        )
        self.model = AutoModel.from_pretrained(
            str(model_dir), local_files_only=True
        ).eval()
        self.cache: dict[str, np.ndarray] = {}

    def embed(self, proteins: dict[str, str], batch: int = 16) -> dict[str, np.ndarray]:
        """Return hash to embedding for the proteins, using the cache."""
        todo = sorted(
            (h for h in proteins if h not in self.cache),
            key=lambda h: len(proteins[h]),
        )
        for start in range(0, len(todo), batch):
            keys = todo[start : start + batch]
            tokens = self.tokenizer(
                [proteins[k][:MAX_PROTEIN_AA] for k in keys],
                return_tensors="pt",
                padding=True,
            )
            with self.torch.no_grad():
                hidden = self.model(**tokens).last_hidden_state
            # The mask is a view of the attention mask, so clearing the start token
            # shortens `lengths` by one and the second clear lands on the last
            # residue, not on the end token. The stored embeddings were computed
            # this way (mean of residues 1 to L-1 and the end token); the code is
            # kept as it is so that they stay reproducible.
            if self.pooling == "mean":
                lengths = tokens["attention_mask"].sum(dim=1)
                mask = tokens["attention_mask"].clone().unsqueeze(-1)
                mask[:, 0] = 0
                for row, length in enumerate(lengths):
                    mask[row, length - 1] = 0
            else:
                mask = tokens["attention_mask"].unsqueeze(-1)
                mask[:, 0] = 0
                lengths = tokens["attention_mask"].sum(dim=1)
                for row, length in enumerate(lengths):
                    mask[row, length - 1] = 0
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
            for key, vector in zip(keys, pooled.numpy(), strict=True):
                self.cache[key] = vector
        return {h: self.cache[h] for h in proteins}


def sequence_hash(sequence: str) -> str:
    """Return the SHA-256 of a sequence, the name-free identity of a contig."""
    return hashlib.sha256(sequence.encode()).hexdigest()


def content_seed(sequences: Mapping[str, str]) -> str:
    """Return a digest of the scored contigs that ignores their names and order."""
    hashes = sorted(sequence_hash(sequence) for sequence in sequences.values())
    return hashlib.sha256("\n".join(hashes).encode()).hexdigest()


def draw_plan(
    windows: Sequence[dict[str, object]],
    sequences: Mapping[str, str],
    *,
    sample_seed: str | None = None,
) -> tuple[str, list[dict[str, object]]]:
    """Return the seed text and the windows in the order in which they are drawn.

    By default the seed is ``content_seed`` of the contigs of 3 kb or more and the
    windows are sorted by (sequence hash of the contig, start). With ``sample_seed``
    (the file stem of the prototype) the seed is that name and the windows keep the
    file order, which reproduces the numbers of the research prototype.
    """
    if sample_seed is not None:
        return sample_seed, list(windows)
    hashes = {name: sequence_hash(sequence) for name, sequence in sequences.items()}
    ordered = sorted(
        windows,
        key=lambda w: (hashes[str(w["source_contig"])], int(str(w["source_start"]))),
    )
    return content_seed(sequences), ordered


def contig_proteins(
    genes: list[dict[str, object]], limit: int, rng: random.Random
) -> dict[str, str]:
    """Return hash to protein for up to ``limit`` called proteins of a window.

    ``genes`` are the gene rows of the window, in gene order.
    """
    found = {}
    for gene in genes:
        protein = str(gene["protein"])
        if len(protein) >= MIN_PROTEIN_AA and "X" not in protein:
            digest = hashlib.sha1(protein.encode(), usedforsecurity=False).hexdigest()
            found[digest] = protein
    keys = sorted(found)
    if len(keys) > limit:
        keys = rng.sample(keys, limit)
    return {k: found[k] for k in keys}


def window_embeddings(
    seed: str,
    contigs: list[dict[str, object]],
    genes: list[dict[str, object]],
    embedder: Embedder,
    esm_genes: int,
) -> pd.DataFrame:
    """Return ``contig_id`` and the ``e_`` columns of every window.

    A window without an eligible protein has no embedding (missing values). The
    windows are visited in the order given and share one random generator seeded
    with ``seed`` (see ``draw_plan``).
    """
    by_contig: dict[str, list[dict[str, object]]] = {}
    for gene in genes:
        by_contig.setdefault(str(gene["contig_id"]), []).append(gene)
    rng = random.Random(f"{SEED}|esm|{seed}")
    rows: list[dict[str, object]] = [
        {"contig_id": contig["contig_id"]} for contig in contigs
    ]
    wanted = [
        contig_proteins(by_contig.get(str(contig["contig_id"]), []), esm_genes, rng)
        for contig in contigs
    ]
    every = {h: p for found in wanted for h, p in found.items()}
    vectors = embedder.embed(every)
    for record, found in zip(rows, wanted, strict=True):
        keys = list(found)
        if keys:
            mean = np.mean([vectors[k] for k in keys], axis=0)
            record.update({f"e_{i}": float(v) for i, v in enumerate(mean)})
    return pd.DataFrame(rows)
