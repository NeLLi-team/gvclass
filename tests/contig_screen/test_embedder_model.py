"""The local ESM-2 model reproduces the embeddings of release 2.1.1 (dbaa765).

``data/esm.npz`` was captured with the code of dbaa765: the legacy and the mean
pooled embeddings of 8 proteins, and the mean pooled embeddings of 40 proteins
(those 8 and 32 more) from a fresh embedder, with one thread. The tests need the
runtime resources in ``<repo>/resources`` and skip without them.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest

from src.contig_screen import extra_features as ef
from src.contig_screen.layout import ESM_DIR
from tests.contig_screen.screen_fixtures import protein

MODEL = Path(__file__).parents[2] / "resources" / ESM_DIR
REFERENCE = Path(__file__).parent / "data" / "esm.npz"

pytestmark = pytest.mark.skipif(
    not (MODEL / "config.json").is_file(), reason="the runtime resources are missing"
)


def proteins() -> dict[str, str]:
    """Return 40 proteins of different lengths, by hash, in hash order."""
    chosen = {}
    for i in range(40):
        sequence = protein(100 + i, 60 + (i * 37) % 150)
        digest = hashlib.sha1(sequence.encode(), usedforsecurity=False).hexdigest()
        chosen[digest] = sequence
    return {key: chosen[key] for key in sorted(chosen)}


def sampled(everything: dict[str, str]) -> dict[str, str]:
    """Return every fifth protein, as the sample of the windows would."""
    return {key: everything[key] for key in list(everything)[::5]}


def stacked(found: dict[str, np.ndarray], keys: list[str]) -> np.ndarray:
    """Return the embeddings in the order of ``keys``."""
    return np.stack([found[key] for key in keys])


@pytest.fixture(scope="module")
def loaded() -> ef.Embedder:
    return ef.Embedder(MODEL, 1, "mean")


@pytest.fixture(scope="module")
def reference() -> dict[str, np.ndarray]:
    with np.load(REFERENCE) as stored:
        return {name: stored[name] for name in stored.files}


def test_a_second_embedder_with_the_weights_of_the_first_has_its_own_cache(
    loaded: ef.Embedder, reference: dict[str, np.ndarray]
) -> None:
    everything = proteins()
    sample = sampled(everything)
    first = loaded.with_pooling("mean")
    first.embed(sample)

    second = first.with_pooling("mean")
    found = second.embed(everything)

    assert second.model is first.model
    assert len(first.cache) == len(sample)
    assert np.allclose(
        stacked(found, list(everything)), reference["mean_full"], rtol=0, atol=1e-6
    )


def test_a_legacy_embedder_with_the_weights_of_a_mean_embedder_pools_as_before(
    loaded: ef.Embedder, reference: dict[str, np.ndarray]
) -> None:
    sample = sampled(proteins())

    legacy = loaded.with_pooling("legacy").embed(sample)

    assert np.allclose(
        stacked(legacy, list(sample)), reference["legacy_sampled"], rtol=0, atol=1e-6
    )
    assert not np.allclose(reference["legacy_sampled"], reference["mean_sampled"])


def test_a_mean_embedder_matches_the_released_embeddings_of_the_sample(
    loaded: ef.Embedder, reference: dict[str, np.ndarray]
) -> None:
    sample = sampled(proteins())

    mean = loaded.with_pooling("mean").embed(sample)

    assert np.allclose(
        stacked(mean, list(sample)), reference["mean_sampled"], rtol=0, atol=1e-6
    )
