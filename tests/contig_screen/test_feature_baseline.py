"""Features of bins without ties or duplicate sequences stay what the models expect.

``data/feature_baseline.npz`` holds the output of the code before contigs were
made order independent (neighbour ties, partner ties, identical sequences).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.contig_screen.thresholds import called_bp
from tests.contig_screen import screen_fixtures as fx

BASELINE = Path(__file__).parent / "data" / "feature_baseline.npz"


@pytest.fixture(scope="module")
def baseline() -> dict[str, np.ndarray]:
    with np.load(BASELINE) as stored:
        return {name: stored[name] for name in stored.files}


def test_context_features_match_the_baseline(baseline: dict[str, np.ndarray]) -> None:
    table = fx.context_features()

    assert table.columns.tolist() == baseline["context_columns"].tolist()
    assert np.allclose(
        table.to_numpy(dtype=float), baseline["context"], rtol=1e-9, equal_nan=True
    )


def test_order_features_match_the_baseline(baseline: dict[str, np.ndarray]) -> None:
    table = fx.order_features().drop(columns=["contig", "genome"])

    assert table.columns.tolist() == baseline["order_columns"].tolist()
    assert np.allclose(
        table.to_numpy(dtype=float), baseline["order"], rtol=1e-9, equal_nan=True
    )


def test_window_embeddings_match_the_baseline(
    baseline: dict[str, np.ndarray],
) -> None:
    table = fx.embedding_features()

    assert np.allclose(
        table.to_numpy(dtype=float), baseline["embedding"], rtol=1e-9, equal_nan=True
    )


@pytest.mark.parametrize("tau", [0.3, 0.5, 0.7])
def test_called_base_pairs_match_the_baseline(
    baseline: dict[str, np.ndarray], tau: float
) -> None:
    windows, lengths = fx.call_windows()

    called = called_bp(windows, lengths, tau)

    assert [called[contig] for contig in lengths] == baseline[f"called_{tau}"].tolist()
