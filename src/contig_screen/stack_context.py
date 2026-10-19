"""Stage 1 of the screen: the stack of the boosting model and the set model.

The stage 1 score of a window is the sigmoid of the mean of two logits: that of the
boosting probability (clipped at 1e-6 from either end) and that of the set model.
Stage 2 (``context_stage``) takes this score as its input.
"""

from __future__ import annotations

import numpy as np

PROBABILITY_CLIP = 1e-6


def stack_scores(gbm: np.ndarray, logit: np.ndarray) -> np.ndarray:
    """Return the sigmoid of the mean of the logit of a probability and a logit."""
    clipped = np.clip(gbm, PROBABILITY_CLIP, 1 - PROBABILITY_CLIP)
    mixed = 0.5 * (np.log(clipped / (1 - clipped)) + logit)
    return 1 / (1 + np.exp(-mixed))
