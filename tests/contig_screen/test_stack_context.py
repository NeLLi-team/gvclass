"""Tests of the stage 1 stack of the boosting probability and the set logit."""

from __future__ import annotations

import numpy as np
import pytest

from src.contig_screen.stack_context import stack_scores


def test_stack_scores_are_the_sigmoid_of_the_mean_of_the_two_logits() -> None:
    gbm = np.array([0.2, 0.9])
    logit = np.array([-1.0, 2.0])

    scores = stack_scores(gbm, logit)

    mixed = 0.5 * (np.log(gbm / (1 - gbm)) + logit)
    assert scores.tolist() == pytest.approx((1 / (1 + np.exp(-mixed))).tolist())


def test_stack_scores_clip_probabilities_of_zero_and_one_at_one_in_a_million() -> None:
    scores = stack_scores(np.array([0.0, 1.0]), np.array([0.0, 0.0]))

    edge = 0.5 * np.log(1e-6 / (1 - 1e-6))
    assert np.isfinite(scores).all()
    assert scores[0] == pytest.approx(1 / (1 + np.exp(-edge)))
    assert scores[1] == pytest.approx(1 / (1 + np.exp(edge)))


def test_a_confident_set_model_can_move_a_neutral_boosting_score() -> None:
    scores = stack_scores(np.array([0.5]), np.array([6.0]))

    assert scores[0] == pytest.approx(1 / (1 + np.exp(-3.0)))
