"""Regression tests for the Phase 1.4 SHA-256 gate on the contamination model.

v1.5.0 moved the bundled model out of ``src/bundled_models/`` into the runtime
resources bundle at ``resources/contamination/model.joblib``. The scorer now
resolves the path via ``database_path / CONTAMINATION_MODEL_FILE`` instead of
``Path(__file__).resolve().parents[1] / "bundled_models"`` so that model
rotation is a resources-tarball concern, not a Python-package concern.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

RESOURCES_DIR = Path(__file__).resolve().parents[1] / "resources"


def _bundled_model_path() -> Path:
    from src.core.contamination_scoring import CONTAMINATION_MODEL_FILE

    return RESOURCES_DIR / CONTAMINATION_MODEL_FILE


def test_bundled_model_sha256_matches_repo_tracked_constant() -> None:
    """The code-level constant must match the on-disk digest in the resources
    tree; drift means the pipeline refuses to load the model at runtime."""
    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256
    from src.utils.common import sha256_file

    model_path = _bundled_model_path()
    assert model_path.exists(), f"Contamination model missing at {model_path}"
    assert sha256_file(model_path) == CONTAMINATION_MODEL_SHA256


def test_bundled_model_card_advertises_matching_sha256() -> None:
    """The model card YAML (beside the .joblib) must advertise the same digest."""
    import yaml

    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256

    card_path = _bundled_model_path().with_suffix(".yaml")
    assert card_path.exists(), f"Model card missing at {card_path}"
    card = yaml.safe_load(card_path.read_text())
    assert card["sha256"] == CONTAMINATION_MODEL_SHA256


def test_load_ml_model_raises_on_sha256_mismatch_and_never_reaches_joblib(
    tmp_path: Path,
) -> None:
    """Any tampering with the bundle must fail-closed before ``joblib.load``."""
    import shutil

    from src.core.contamination_scoring import (
        CONTAMINATION_MODEL_FILE,
        ContaminationScorer,
    )

    # Point the scorer at a tmp database with a tampered model copy.
    fake_db = tmp_path / "resources"
    (fake_db / "contamination").mkdir(parents=True)
    shutil.copyfile(
        _bundled_model_path(), fake_db / CONTAMINATION_MODEL_FILE
    )
    # Overwrite the last byte to change the digest without touching the path.
    target = fake_db / CONTAMINATION_MODEL_FILE
    target.write_bytes(target.read_bytes() + b"X")

    with patch("joblib.load") as joblib_load:
        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            ContaminationScorer(fake_db)
    joblib_load.assert_not_called()


def test_load_ml_model_raises_when_bundled_file_missing(tmp_path: Path) -> None:
    """If the model is physically missing under resources/contamination/, the
    scorer must fail loudly rather than silently degrading."""
    from src.core.contamination_scoring import ContaminationScorer

    # tmp_path has no contamination/model.joblib → must raise.
    with pytest.raises(RuntimeError, match="Contamination model not found"):
        ContaminationScorer(tmp_path)


def test_sensitive_mode_scorer_still_loads_model(tmp_path: Path) -> None:
    """Sensitive mode must NOT skip the load. Construction against the real
    resources tree should succeed with ``ml_available=True``."""
    from src.core.contamination_scoring import ContaminationScorer

    scorer = ContaminationScorer(RESOURCES_DIR, sensitive_mode=True)
    assert scorer.sensitive_mode is True
    assert scorer.ml_available is True
    assert scorer.ml_model is not None


def test_rotate_contamination_model_script_prints_current_digest() -> None:
    """The rotation helper must print the on-disk digest so maintainers can
    update the code constant after retraining."""
    import runpy
    import sys
    from io import StringIO

    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256

    script_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "rotate_contamination_model.py"
    )
    original_stdout, original_argv = sys.stdout, sys.argv
    sys.stdout = captured = StringIO()
    sys.argv = ["rotate_contamination_model.py", str(_bundled_model_path())]
    try:
        with pytest.raises(SystemExit) as exc_info:
            runpy.run_path(str(script_path), run_name="__main__")
    finally:
        sys.stdout = original_stdout
        sys.argv = original_argv

    assert exc_info.value.code == 0
    assert CONTAMINATION_MODEL_SHA256 in captured.getvalue()
