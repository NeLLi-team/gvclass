"""Second stage: score each window in the context of the bin that holds it.

Stage 1 scores every window alone. A giant-virus bin is mostly one genome, so a
window that differs in composition from the bulk of its bin, or that sits next to
windows that stage 1 calls cellular, is more likely foreign. The bin-context
features are the stage 1 score, the length, the distance of the window to the bulk
composition of the other windows of its bin (leave one out), to the score-weighted
viral and cellular centroids, the mean score of its nearest windows in composition
space, and bin totals. The windows of one tiling of a bin form a context
(``tilings``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

OWN = (
    "c_gc",
    "s_coding_density",
    "s_genes_per_kb",
    "s_gene_len_mean",
    "h_frac_nohit",
    "h_share_cellular",
    "h_share_giant",
    "h_share_PHAGE",
)
K_NEAREST = 5


def logit(p: np.ndarray) -> np.ndarray:
    """Return clipped log odds."""
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.clip(np.log(p / (1 - p)), -12.0, 12.0)


def loo_centroid_distance(u: np.ndarray, weight: np.ndarray) -> np.ndarray:
    """Return the distance of each row to the weighted centroid of the other rows."""
    total = weight.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        rest = total - weight
        centroid = ((u * weight[:, None]).sum(axis=0) - u * weight[:, None]) / rest[
            :, None
        ]
        distance = np.linalg.norm(u - centroid, axis=1)
    distance[~np.isfinite(distance) | (rest <= 1e-9)] = np.nan
    return distance


def nearest_means(distance: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the mean score and mean distance of the nearest windows of each window.

    Windows at equal distance are taken in the order of their score, so the choice
    depends on the content of the windows and not on their row order. ``distance``
    is the full square matrix; its diagonal is ignored.
    """
    n = len(distance)
    distance = distance.copy()
    np.fill_diagonal(distance, np.inf)
    scores = np.broadcast_to(z, distance.shape)
    nearest = np.lexsort((scores, distance), axis=1)[:, : min(K_NEAREST, n - 1)]
    rows = np.arange(n)[:, None]
    return z[nearest].mean(axis=1), distance[rows, nearest].mean(axis=1)


def bin_context(
    part: pd.DataFrame, comp: list[str], mu: np.ndarray, sd: np.ndarray
) -> pd.DataFrame:
    """Return the context features of the windows of one bin (rows keep their order)."""
    n = len(part)
    both = np.nan_to_num((part[comp].to_numpy(dtype=float) - mu) / sd)
    is_embedding = np.array([c.startswith("e_") for c in comp], dtype=bool)
    u = both[:, ~is_embedding]
    length = part["length_bp"].to_numpy(dtype=float)
    p = part["p"].to_numpy(dtype=float)
    z = logit(p)
    total = length.sum()
    out = pd.DataFrame(index=part.index)
    out["z"] = z
    out["log_len"] = np.log10(length)
    for name in OWN:
        out[name] = part[name].to_numpy(dtype=float)
    out["n_win"] = n
    out["bin_bp"] = total
    out["frac_cell_bp"] = float((length * (p > 0.5)).sum() / total)
    out["mean_z_bp"] = float((length * z).sum() / total)
    out["d_bulk"] = loo_centroid_distance(u, length)
    out["d_viral"] = loo_centroid_distance(u, length * (1 - p) ** 2)
    out["d_cell"] = loo_centroid_distance(u, length * p**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["d_ratio"] = out["d_viral"] / (out["d_viral"] + out["d_cell"])
    if n > 1:
        out["knn_z"], out["knn_d"] = nearest_means(pairwise_distances(u), z)
        out["z_minus_knn"] = z - out["knn_z"]
        # tied scores share their mean rank: the row order of a bin must not decide it
        ranks = pd.Series(z).rank(method="average").to_numpy() - 1
        out["z_rank"] = ranks / max(n - 1, 1)
    else:
        for name in ("knn_z", "knn_d", "z_minus_knn", "z_rank"):
            out[name] = np.nan
    if is_embedding.any():
        extra = embedding_context(both[:, is_embedding], length, p, z)
        out = out.join(extra.set_axis(out.index))
    return out


def embedding_context(
    v: np.ndarray, length: np.ndarray, p: np.ndarray, z: np.ndarray
) -> pd.DataFrame:
    """Return distance and neighbour features in protein embedding space."""
    n = len(v)
    out = pd.DataFrame(index=range(n))
    out["d_bulk_e"] = loo_centroid_distance(v, length)
    out["d_viral_e"] = loo_centroid_distance(v, length * (1 - p) ** 2)
    out["d_cell_e"] = loo_centroid_distance(v, length * p**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["d_ratio_e"] = out["d_viral_e"] / (out["d_viral_e"] + out["d_cell_e"])
    if n > 1:
        out["knn_z_e"], out["knn_d_e"] = nearest_means(pairwise_distances(v), z)
        out["z_minus_knn_e"] = z - out["knn_z_e"]
    else:
        for name in ("knn_z_e", "knn_d_e", "z_minus_knn_e"):
            out[name] = np.nan
    return out


def context_table(
    windows: pd.DataFrame, comp: list[str], mu: np.ndarray, sd: np.ndarray
) -> pd.DataFrame:
    """Return context features for every window; ``windows`` needs a ``bin_id``."""
    parts = [
        bin_context(part, comp, mu, sd)
        for _, part in windows.groupby("bin_id", sort=False)
    ]
    return pd.concat(parts).loc[windows.index]


def tilings(frame: pd.DataFrame, tile_kb: int, step_kb: int) -> pd.Series:
    """Return the tiling of each window (windows of one tiling do not overlap)."""
    if not step_kb or step_kb >= tile_kb:
        return pd.Series(0, index=frame.index)
    phases = -(-tile_kb // step_kb)
    return (frame["source_start"] // (step_kb * 1000)) % phases
