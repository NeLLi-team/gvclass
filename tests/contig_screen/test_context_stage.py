"""Tests of the bin-context features of stage 2."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.contig_screen import context_stage as cs

COMP = ["c_gc", "c_aa_A"]


def bin_windows(scores: list[float], gc: list[float]) -> pd.DataFrame:
    n = len(scores)
    frame = pd.DataFrame(
        {
            "p": scores,
            "length_bp": [10_000] * n,
            "c_gc": gc,
            "c_aa_A": np.linspace(0.05, 0.15, n),
        }
    )
    for name in cs.OWN:
        if name not in frame:
            frame[name] = 0.5
    return frame


def test_tilings_of_overlapping_windows_group_those_that_do_not_overlap() -> None:
    frame = pd.DataFrame({"source_start": [0, 5000, 10_000, 15_000, 20_000]})

    phase = cs.tilings(frame, 10, 5)

    assert phase.tolist() == [0, 1, 0, 1, 0]


@pytest.mark.parametrize(("tile_kb", "step_kb"), [(10, 0), (10, 10), (10, 20)])
def test_tilings_without_overlap_are_one_tiling(tile_kb: int, step_kb: int) -> None:
    frame = pd.DataFrame({"source_start": [0, 10_000, 20_000]})

    assert cs.tilings(frame, tile_kb, step_kb).tolist() == [0, 0, 0]


def test_logit_is_clipped_at_both_ends() -> None:
    values = cs.logit(np.array([0.0, 0.5, 1.0]))

    assert values.tolist() == pytest.approx([-12.0, 0.0, 12.0])


def test_the_distance_to_the_other_windows_leaves_the_window_itself_out() -> None:
    u = np.array([[0.0], [2.0], [4.0]])

    distance = cs.loo_centroid_distance(u, np.ones(3))

    assert distance.tolist() == pytest.approx([3.0, 0.0, 3.0])


def test_a_lone_window_has_no_distance_and_no_neighbours() -> None:
    frame = bin_windows([0.1], [0.4])

    out = cs.bin_context(frame, COMP, np.zeros(2), np.ones(2))

    assert np.isnan(out.loc[0, "d_bulk"])
    assert np.isnan(out.loc[0, "knn_z"])
    assert out.loc[0, "n_win"] == 1


def test_bin_context_ranks_tied_scores_alike_whatever_the_row_order() -> None:
    frame = bin_windows([0.2, 0.2, 0.9, 0.2], [0.4, 0.5, 0.6, 0.7])
    shuffled = frame.iloc[[3, 0, 2, 1]]

    first = cs.bin_context(frame, COMP, np.zeros(2), np.ones(2))
    second = cs.bin_context(shuffled, COMP, np.zeros(2), np.ones(2))

    assert len(set(first.loc[[0, 1, 3], "z_rank"])) == 1
    assert second.loc[frame.index, "z_rank"].tolist() == first["z_rank"].tolist()


def test_context_table_keeps_the_window_order_across_bins() -> None:
    frame = pd.concat(
        [
            bin_windows([0.1, 0.9], [0.4, 0.6]).assign(bin_id="b1"),
            bin_windows([0.5, 0.5, 0.5], [0.3, 0.5, 0.7]).assign(bin_id="b2"),
        ],
        ignore_index=True,
    )

    table = cs.context_table(frame, COMP, np.zeros(2), np.ones(2))

    assert table.index.tolist() == [0, 1, 2, 3, 4]
    assert table["n_win"].tolist() == [2, 2, 3, 3, 3]


def test_a_bin_with_embedding_columns_gets_features_in_embedding_space() -> None:
    frame = bin_windows([0.1, 0.5, 0.9], [0.4, 0.5, 0.6])
    frame["e_0"] = [0.0, 1.0, 2.0]

    out = cs.bin_context(frame, [*COMP, "e_0"], np.zeros(3), np.ones(3))

    assert "d_bulk_e" in out
    assert "knn_z_e" in out


def tied_windows() -> pd.DataFrame:
    """Return seven windows with the same composition and different scores."""
    frame = bin_windows([0.05, 0.9, 0.3, 0.6, 0.15, 0.75, 0.45], [0.5] * 7)
    frame["c_aa_A"] = 0.1
    frame["e_0"] = 1.0
    return frame


@pytest.mark.parametrize("comp", [COMP, [*COMP, "e_0"]])
def test_bin_context_with_tied_neighbours_does_not_depend_on_the_row_order(
    comp: list[str],
) -> None:
    frame = tied_windows()
    shuffled = frame.iloc[[4, 6, 1, 0, 5, 3, 2]]
    mu, sd = np.zeros(len(comp)), np.ones(len(comp))

    first = cs.bin_context(frame, comp, mu, sd)
    second = cs.bin_context(shuffled, comp, mu, sd)

    assert np.allclose(
        second.loc[frame.index].to_numpy(dtype=float),
        first.to_numpy(dtype=float),
        equal_nan=True,
    )


def test_tied_neighbours_are_taken_in_the_order_of_their_score() -> None:
    distance = np.zeros((7, 7))
    z = np.array([5.0, 1.0, 4.0, 2.0, 7.0, 3.0, 6.0])

    knn_z, knn_d = cs.nearest_means(distance, z)

    assert knn_z[0] == pytest.approx(np.mean([1.0, 2.0, 3.0, 4.0, 6.0]))
    assert knn_d.tolist() == [0.0] * 7
