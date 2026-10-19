"""Load the model bundle and the thresholds of the giant-virus screen.

``model.joblib`` is a dict with the fitted ``model``, its feature ``columns`` and the
single-copy ``priors`` (a ``Priors`` object of the prototype module ``gv_order``).
The file is read once, checked against the manifest, and unpickled from those bytes;
the prototype's class name is mapped to ``gv_order.Priors`` here, so no module named
``gv_order`` has to exist.
"""

from __future__ import annotations

import io
import json
from dataclasses import dataclass
from typing import Any

from joblib.numpy_pickle import NumpyUnpickler

from src.contig_screen import layout
from src.contig_screen.gv_order import Priors

OPERATING_POINT = "0.01"
PROTOTYPE_CLASSES = {("gv_order", "Priors"): Priors}


@dataclass(frozen=True, slots=True)
class Bundle:
    """The fitted classifier with its feature columns and priors."""

    model: Any  # a scikit-learn estimator; scikit-learn ships no types
    columns: list[str]
    priors: Priors


class BundleUnpickler(NumpyUnpickler):
    """The joblib unpickler, with the prototype's ``gv_order.Priors`` mapped."""

    def find_class(self, module: str, name: str) -> Any:  # noqa: ANN401  # pickle API
        """Return the class for ``module.name``, mapping the prototype's classes."""
        return PROTOTYPE_CLASSES.get((module, name)) or super().find_class(module, name)


def load_bundle(resources: layout.ScreenResources) -> Bundle:
    """Return the model, columns and priors of the verified model file.

    Raises:
        ValueError: The file has no priors (a model that needs a retrain).
    """
    raw = resources.read_bytes(layout.GV_MODEL_FILE)
    document = BundleUnpickler(layout.GV_MODEL_FILE, io.BytesIO(raw), True).load()
    if document.get("priors") is None:
        msg = f"{layout.GV_MODEL_FILE} has no priors; retrain the model"
        raise ValueError(msg)
    return Bundle(document["model"], list(document["columns"]), document["priors"])


def load_thresholds(resources: layout.ScreenResources) -> dict[str, float]:
    """Return the thresholds by the share of clean contigs called (``"0.01"``...)."""
    document = json.loads(resources.read_bytes(layout.GV_THRESHOLDS_FILE))
    return {key: float(value) for key, value in document.items()}
