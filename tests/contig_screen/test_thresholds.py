"""Tests of the base-pair calls from overlapping window scores."""

from __future__ import annotations

import pandas as pd

from src.contig_screen.thresholds import called_bp, window_segments


def windows(starts: list[int], scores: list[float], length: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_contig": "c1",
            "source_start": starts,
            "length_bp": length,
            "p": scores,
        }
    )


def test_window_segments_count_windows_and_average_their_scores() -> None:
    frame = windows([0, 5], [0.2, 0.8], 10)

    edges, covered, mean = window_segments(frame, 15)

    assert edges.tolist() == [0, 5, 10, 15]
    assert covered.tolist() == [1, 2, 1]
    assert mean.tolist() == [0.2, 0.5, 0.8]


def test_a_window_past_the_end_of_the_contig_is_cut() -> None:
    edges, covered, _ = window_segments(windows([5], [0.5], 10), 8)

    assert edges.tolist() == [5, 8]
    assert covered.tolist() == [1]


def test_a_base_pair_is_called_when_the_mean_of_its_windows_reaches_the_threshold() -> (
    None
):
    frame = windows([0, 5], [0.2, 0.8], 10)

    calls = called_bp(frame, {"c1": 15}, 0.5)

    assert calls == {"c1": 10}


def test_a_base_pair_that_no_window_covers_is_not_called() -> None:
    calls = called_bp(windows([0], [1.0], 5), {"c1": 12}, 0.1)

    assert calls == {"c1": 5}


def test_a_gap_between_windows_is_not_called() -> None:
    calls = called_bp(windows([0, 20], [1.0, 1.0], 5), {"c1": 30}, 0.5)

    assert calls == {"c1": 10}
