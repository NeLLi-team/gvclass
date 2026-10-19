"""Tests of the model bundle and threshold loading of the giant-virus screen."""

from __future__ import annotations

import io
import sys
import types
from pathlib import Path

import joblib
import numpy as np
import pytest
from sklearn.ensemble import HistGradientBoostingClassifier

from src.contig_screen import gv_model, layout
from src.contig_screen.gv_order import Priors
from tests.contig_screen.conftest import BundleFactory

COLUMNS = ["dup_frac", "pair_shared"]


def make_priors() -> Priors:
    return Priors({"OG1": 0.02}, 0.1, np.zeros(3), np.ones(3), {})


def make_model() -> HistGradientBoostingClassifier:
    features = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]])
    model = HistGradientBoostingClassifier(min_samples_leaf=1, max_iter=5)
    return model.fit(features, [0, 1, 0, 1])


def dump(document: object) -> bytes:
    buffer = io.BytesIO()
    joblib.dump(document, buffer)
    return buffer.getvalue()


def prototype_bytes(monkeypatch: pytest.MonkeyPatch) -> bytes:
    """Return a model file with priors named ``gv_order.Priors`` (prototype)."""
    module = types.ModuleType("gv_order")
    monkeypatch.setattr(Priors, "__module__", "gv_order")
    module.Priors = Priors  # type: ignore[attr-defined]  # a stand-in module
    with monkeypatch.context() as patch:
        patch.setitem(sys.modules, "gv_order", module)
        return dump(
            {"model": make_model(), "columns": COLUMNS, "priors": make_priors()}
        )


def test_the_prototype_priors_class_is_mapped_to_the_ported_class(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, make_resources: BundleFactory
) -> None:
    raw = prototype_bytes(monkeypatch)
    monkeypatch.undo()
    assert "gv_order" not in sys.modules
    resources = make_resources(tmp_path, {layout.GV_MODEL_FILE: raw})

    bundle = gv_model.load_bundle(resources)

    assert type(bundle.priors) is Priors
    assert bundle.priors.p_multi == {"OG1": 0.02}
    assert bundle.columns == COLUMNS
    assert bundle.model.predict_proba(np.array([[1.0, 1.0]])).shape == (1, 2)


def test_a_model_without_priors_is_refused(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    raw = dump({"model": make_model(), "columns": COLUMNS, "priors": None})
    resources = make_resources(tmp_path, {layout.GV_MODEL_FILE: raw})

    with pytest.raises(ValueError, match="no priors"):
        gv_model.load_bundle(resources)


def test_a_model_file_that_differs_from_the_manifest_is_not_loaded(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    raw = dump({"model": make_model(), "columns": COLUMNS, "priors": make_priors()})
    resources = make_resources(tmp_path, {layout.GV_MODEL_FILE: raw})
    (tmp_path / layout.GV_MODEL_FILE).write_bytes(raw + b" ")

    with pytest.raises(layout.IntegrityError, match="mismatch"):
        gv_model.load_bundle(resources)


def test_thresholds_are_read_by_their_clean_contig_share(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    text = b'{"0.01": 0.74299, "0.05": 0.34855}'
    resources = make_resources(tmp_path, {layout.GV_THRESHOLDS_FILE: text})

    assert gv_model.load_thresholds(resources) == {"0.01": 0.74299, "0.05": 0.34855}
