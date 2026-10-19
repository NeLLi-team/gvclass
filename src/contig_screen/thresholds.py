"""Call base pairs from overlapping window scores."""

from __future__ import annotations

import numpy as np
import pandas as pd


def window_cover(rows: pd.DataFrame, size: int) -> tuple[np.ndarray, np.ndarray]:
    """Return, per base pair of one contig, the window count and the mean score.

    ``rows`` are the windows of the contig (``source_start``, ``length_bp``, ``p``).
    """
    start = rows["source_start"].to_numpy(dtype=int)
    end = np.minimum(start + rows["length_bp"].to_numpy(dtype=int), size)
    p = rows["p"].to_numpy(dtype=float)
    total = np.zeros(size + 1)
    count = np.zeros(size + 1)
    np.add.at(total, start, p)
    np.add.at(total, end, -p)
    np.add.at(count, start, 1.0)
    np.add.at(count, end, -1.0)
    covered = np.cumsum(count)[:size]
    mean = np.cumsum(total)[:size] / np.maximum(covered, 1.0)
    return covered, mean


def sliding_calls(
    windows: pd.DataFrame, lengths: dict[str, int], tau: float
) -> dict[str, np.ndarray]:
    """Return per source contig the call of every base pair from overlapping windows.

    A base pair is called when the mean P(cellular) of the windows that cover it
    reaches ``tau``; a base pair that no window covers is not called.
    """
    calls: dict[str, np.ndarray] = {}
    for contig, rows in windows.groupby("source_contig"):
        covered, mean = window_cover(rows, int(lengths[str(contig)]))
        calls[str(contig)] = (covered > 0.5) & (mean >= tau)
    return calls
