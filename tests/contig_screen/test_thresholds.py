"""Tests of the base-pair calls from overlapping window scores."""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.contig_screen.thresholds import sliding_calls, window_cover


def windows(starts: list[int], scores: list[float], length: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_contig": "c1",
            "source_start": starts,
            "length_bp": length,
            "p": scores,
        }
    )


def test_window_cover_counts_windows_and_averages_their_scores() -> None:
    frame = windows([0, 5], [0.2, 0.8], 10)

    covered, mean = window_cover(frame, 15)

    assert covered.tolist() == [1] * 5 + [2] * 5 + [1] * 5
    assert mean[:5].tolist() == [0.2] * 5
    assert mean[5:10].tolist() == [0.5] * 5
    assert mean[10:].tolist() == [0.8] * 5


def test_a_window_past_the_end_of_the_contig_is_cut() -> None:
    covered, _ = window_cover(windows([5], [0.5], 10), 8)

    assert covered.tolist() == [0] * 5 + [1] * 3


def test_a_base_pair_is_called_when_the_mean_of_its_windows_reaches_the_threshold() -> (
    None
):
    frame = windows([0, 5], [0.2, 0.8], 10)

    calls = sliding_calls(frame, {"c1": 15}, 0.5)

    assert calls["c1"].tolist() == [False] * 5 + [True] * 10


def test_a_base_pair_that_no_window_covers_is_not_called() -> None:
    calls = sliding_calls(windows([0], [1.0], 5), {"c1": 12}, 0.1)

    assert calls["c1"].sum() == 5
    assert not calls["c1"][5:].any()
    assert isinstance(calls["c1"], np.ndarray)
