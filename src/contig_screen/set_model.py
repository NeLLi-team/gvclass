"""Set model over per-gene hit tokens: scoring with the fitted models.

Every gene of a window is a token (class of the hits, identity, bit score per amino
acid, the zero token for a gene without a hit, and the protein embedding); a shared
network embeds each token, the embeddings of a window are pooled by mean and by
maximum, and a second network scores the pool together with a few window scalars
(length, gene density; with the variant ``hits_struct`` also gene length, strand
switches, coding density and GC). The score does not depend on the order of the genes.
The fitted models are ``model_<variant>_k<k>.pt`` (state dicts, loaded with
``weights_only``) and ``scaler_<variant>.json``, as written by the training code of
the prototype. This module imports torch; import it only after
``torch_stack.require_torch``.
"""

from __future__ import annotations

import io
import json
from dataclasses import dataclass

import numpy as np
import torch
from torch import nn

from src.contig_screen import hit_tokens as ht, layout

STRUCT_COLUMNS = list(range(2, len(ht.SCALARS)))
EMBEDDING_START = len(ht.GENE_FEATURES)


@dataclass(frozen=True)
class Scaler:
    """Fill values, means and standard deviations of the scalar inputs.

    A missing value (NaN) is replaced by the mean of the non-missing values of the
    rows that the scaler was fitted on (0 for a column that is missing everywhere),
    so that a window is scored the same alone and in any batch.
    """

    fill: np.ndarray
    mean: np.ndarray
    std: np.ndarray

    def transform(self, values: np.ndarray) -> np.ndarray:
        """Return the matrix with missing values filled and the columns standardized."""
        filled = np.where(np.isnan(values), self.fill, values)
        return ((filled - self.mean) / self.std).astype(np.float32)


@dataclass(frozen=True)
class EmbeddingScale:
    """Mean and deviation of the protein embedding of the genes that have one.

    A token matrix with embeddings (``gene_embeddings.extend_tokens``) holds the raw
    embedding after the hit tokens and, in the last column, a flag that says whether
    the gene has one. Genes without an embedding stay zero.
    """

    mean: np.ndarray
    std: np.ndarray

    def transform(self, genes: np.ndarray) -> np.ndarray:
        """Return the tokens with the embedding standardized (zero without a flag)."""
        out = genes.copy()
        present = genes[:, -1] > 0
        block = (genes[:, EMBEDDING_START:-1] - self.mean) / self.std
        out[:, EMBEDDING_START:-1] = np.where(present[:, None], block, 0.0)
        return out


def set_threads(threads: int) -> None:
    """Set the number of torch threads of the scoring."""
    torch.set_num_threads(threads)


def scalar_matrix(raw: np.ndarray, variant: str) -> np.ndarray:
    """Return the scalar inputs of a variant from the raw ``SCALARS`` columns.

    Missing values stay missing here; the scaler fills them.
    """
    length = np.log10(np.maximum(raw[:, 0], 1e-3))
    density = raw[:, 1] / np.maximum(raw[:, 0], 1e-3)
    columns = [length, density]
    if variant == "hits_struct":
        columns += list(raw[:, STRUCT_COLUMNS].T)
    return np.stack(columns, axis=1).astype(np.float32)


def batch_genes(
    genes: np.ndarray, offsets: np.ndarray, rows: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the tokens of the chosen contigs and the batch position of each token."""
    starts = offsets[rows]
    counts = offsets[rows + 1] - starts
    total = int(counts.sum())
    base = np.repeat(starts - (np.cumsum(counts) - counts), counts)
    index = base + np.arange(total)
    return genes[index], np.repeat(np.arange(len(rows)), counts)


class SetNet(nn.Module):
    """Token network, mean and max pooling over the tokens of a window, score head."""

    def __init__(self, n_features: int, n_scalars: int, hidden: int = 64) -> None:
        """Build the token and head networks for the given input widths."""
        super().__init__()
        self.hidden = hidden
        self.gene = nn.Sequential(
            nn.Linear(n_features, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
        )
        self.head = nn.Sequential(
            nn.Linear(2 * hidden + n_scalars, hidden),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden, 1),
        )

    def forward(
        self,
        genes: torch.Tensor,
        piece: torch.Tensor,
        n_pieces: int,
        scalars: torch.Tensor,
    ) -> torch.Tensor:
        """Return the logit of every window of the batch."""
        embedded = self.gene(genes)
        counts = torch.zeros(n_pieces).index_add_(0, piece, torch.ones(len(piece)))
        total = torch.zeros(n_pieces, self.hidden).index_add_(0, piece, embedded)
        mean = total / counts.clamp(min=1).unsqueeze(1)
        spread = piece.unsqueeze(1).expand_as(embedded)
        peak = torch.zeros(n_pieces, self.hidden).scatter_reduce(
            0, spread, embedded, reduce="amax", include_self=False
        )
        return self.head(torch.cat([mean, peak, scalars], dim=1)).squeeze(1)


def predict(
    model: SetNet,
    genes: np.ndarray,
    offsets: np.ndarray,
    scalars: np.ndarray,
    rows: np.ndarray,
    batch: int = 2048,
) -> np.ndarray:
    """Return the logits of the chosen contigs."""
    model.eval()
    out = []
    with torch.no_grad():
        for start in range(0, len(rows), batch):
            chosen = rows[start : start + batch]
            tokens, piece = batch_genes(genes, offsets, chosen)
            logits = model(
                torch.from_numpy(tokens),
                torch.from_numpy(piece),
                len(chosen),
                torch.from_numpy(scalars[chosen]),
            )
            out.append(logits.numpy())
    return np.concatenate(out) if out else np.zeros(0, dtype=np.float32)


def mean_logits(
    models: list[SetNet],
    genes: np.ndarray,
    offsets: np.ndarray,
    scalars: np.ndarray,
    rows: np.ndarray,
) -> np.ndarray:
    """Return the mean over models of the logits of the chosen contigs."""
    return np.mean([predict(m, genes, offsets, scalars, rows) for m in models], axis=0)


def load_fitted(
    resources: layout.ScreenResources, variant: str
) -> tuple[list[SetNet], Scaler, dict]:
    """Return the saved models, scaler and settings of a variant.

    Every file is read once, checked against the manifest, and deserialized from
    those bytes (``torch.load`` with ``weights_only``).
    """
    base = layout.SET_MODEL_DIR
    meta = json.loads(resources.read_bytes(f"{base}/scaler_{variant}.json"))
    scaler = Scaler(
        np.array(meta["fill"], dtype=np.float32),
        np.array(meta["mean"], dtype=np.float32),
        np.array(meta["std"], dtype=np.float32),
    )
    models = []
    for k in range(meta["n_models"]):
        raw = resources.read_bytes(f"{base}/model_{variant}_k{k}.pt")
        model = SetNet(meta["n_features"], meta["n_scalars"], meta["hidden"])
        model.load_state_dict(torch.load(io.BytesIO(raw), weights_only=True))
        models.append(model.eval())
    return models, scaler, meta


def score_contigs(
    resources: layout.ScreenResources,
    variant: str,
    genes: np.ndarray,
    offsets: np.ndarray,
    raw: np.ndarray,
) -> np.ndarray:
    """Return the logits of all windows of a token table for a saved variant."""
    models, scaler, meta = load_fitted(resources, variant)
    tokens, cuts = ht.ablate(genes, offsets, meta["ablation"])
    if "emb_mean" in meta:
        embedding_scale = EmbeddingScale(
            np.array(meta["emb_mean"], dtype=np.float32),
            np.array(meta["emb_std"], dtype=np.float32),
        )
        tokens = embedding_scale.transform(tokens)
    scaled = scaler.transform(scalar_matrix(raw, meta["scalar_set"]))
    rows = np.arange(len(offsets) - 1)
    return mean_logits(models, tokens, cuts, scaled, rows)
