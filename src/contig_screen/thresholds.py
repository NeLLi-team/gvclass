"""Call base pairs from overlapping window scores."""

from __future__ import annotations

import numpy as np
import pandas as pd


def window_segments(
    rows: pd.DataFrame, size: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the edges of the stretches between window ends and their window stats.

    ``rows`` are the windows of the contig (``source_start``, ``length_bp``, ``p``);
    a window past ``size`` is cut. Stretch ``i`` spans ``edges[i]`` to
    ``edges[i + 1]``; the window count and the mean score are constant within it.
    """
    start = rows["source_start"].to_numpy(dtype=int)
    end = np.minimum(start + rows["length_bp"].to_numpy(dtype=int), size)
    p = rows["p"].to_numpy(dtype=float)
    edges = np.unique(np.concatenate([start, end]))
    first, last = np.searchsorted(edges, start), np.searchsorted(edges, end)
    total = np.zeros(len(edges))
    count = np.zeros(len(edges))
    np.add.at(total, first, p)
    np.add.at(total, last, -p)
    np.add.at(count, first, 1.0)
    np.add.at(count, last, -1.0)
    covered = np.cumsum(count)[:-1]
    mean = np.cumsum(total)[:-1] / np.maximum(covered, 1.0)
    return edges, covered, mean


def called_bp(
    windows: pd.DataFrame, lengths: dict[str, int], tau: float
) -> dict[str, int]:
    """Return per source contig the number of base pairs called from the windows.

    A base pair is called when the mean P(cellular) of the windows that cover it
    reaches ``tau``; a base pair that no window covers is not called.
    """
    calls: dict[str, int] = {}
    for contig, rows in windows.groupby("source_contig"):
        edges, covered, mean = window_segments(rows, int(lengths[str(contig)]))
        called = (covered > 0.5) & (mean >= tau)
        calls[str(contig)] = int(np.diff(edges)[called].sum())
    return calls
