"""Tests of the error for a missing torch or transformers."""

from __future__ import annotations

import importlib.util
from collections.abc import Callable

import pytest
import threadpoolctl
from sklearn.ensemble import (
    HistGradientBoostingClassifier,  # noqa: F401  # loads OpenMP
)

from src.contig_screen import torch_stack


def only_missing(*names: str) -> Callable[[str], object | None]:
    """Return a stand-in for ``find_spec`` that finds every package but ``names``."""
    return lambda name: None if name in names else object()


def test_require_torch_names_both_missing_packages_and_the_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        importlib.util, "find_spec", only_missing("torch", "transformers")
    )

    with pytest.raises(torch_stack.MissingTorchError) as error:
        torch_stack.require_torch()

    message = str(error.value)
    assert "torch and transformers" in message
    assert "pixi run -e" not in message
    assert "pixi environment of this GVClass release" in message


def test_require_torch_names_only_the_package_that_is_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", only_missing("transformers"))

    with pytest.raises(torch_stack.MissingTorchError, match="needs transformers,"):
        torch_stack.require_torch()


def test_require_torch_passes_when_both_packages_are_found(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", only_missing())

    torch_stack.require_torch()


def test_limit_threads_caps_the_openmp_and_blas_pools_inside_the_block() -> None:
    with torch_stack.limit_threads(1):
        inside = [pool["num_threads"] for pool in threadpoolctl.threadpool_info()]

    assert inside
    assert set(inside) == {1}
