"""Regression tests for the Phase 1.4 SHA-256 gate on the bundled
contamination model.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from unittest.mock import patch

import pytest


def _bundled_model_path() -> Path:
    import src.core.contamination_scoring as mod

    return Path(mod.__file__).resolve().parents[1] / "bundled_models" / mod.CONTAMINATION_MODEL_FILE


def test_bundled_model_sha256_matches_repo_tracked_constant() -> None:
    """The constant must always match the bundled model that ships with the
    source tree; drift between the two means the pipeline refuses to load the
    model at runtime."""
    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256
    from src.utils.common import sha256_file

    model_path = _bundled_model_path()
    assert model_path.exists(), "Bundled contamination model is missing"
    assert sha256_file(model_path) == CONTAMINATION_MODEL_SHA256


def test_bundled_model_card_advertises_matching_sha256() -> None:
    """The model card YAML must advertise the same digest as the code constant."""
    import yaml

    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256

    card_path = _bundled_model_path().with_suffix(".yaml")
    assert card_path.exists(), "Model card contamination_model.yaml is missing"
    card = yaml.safe_load(card_path.read_text())
    assert card["sha256"] == CONTAMINATION_MODEL_SHA256


def test_load_ml_model_raises_on_sha256_mismatch(tmp_path: Path) -> None:
    """Any tampering with the bundle must fail-closed before ``joblib.load``."""
    from src.core.contamination_scoring import ContaminationScorer

    # Point the scorer at the bundled model but patch the expected digest to
    # a value that does not match the on-disk file.
    with patch(
        "src.core.contamination_scoring.CONTAMINATION_MODEL_SHA256",
        "deadbeef" * 8,
    ):
        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            ContaminationScorer(tmp_path)


def test_load_ml_model_raises_when_bundled_file_missing(
    tmp_path: Path, monkeypatch
) -> None:
    """If the bundled model is physically missing, construction must fail
    loudly rather than silently degrading to ``ml_available=False``."""
    import src.core.contamination_scoring as scoring

    # Redirect the bundled-model lookup to a path that does not exist.
    missing_root = tmp_path / "no_bundled_models"

    original_init = scoring.ContaminationScorer.__init__

    def patched_init(self, database_path, sensitive_mode=False):
        # Ensure the bundled path resolver walks up to a nonexistent tree.
        monkeypatch.setattr(
            scoring.Path,
            "__file__",
            str(missing_root / "contamination_scoring.py"),
            raising=False,
        )
        original_init(self, database_path, sensitive_mode)

    # Instead of monkeypatching __file__ (which does not propagate via
    # parents[1]), rename the bundled file temporarily.
    real_path = _bundled_model_path()
    backup = real_path.with_suffix(".joblib.test_backup")
    real_path.rename(backup)
    try:
        with pytest.raises(RuntimeError, match="Contamination model not found"):
            scoring.ContaminationScorer(tmp_path)
    finally:
        backup.rename(real_path)


def test_sensitive_mode_scorer_skips_model_load(tmp_path: Path, monkeypatch) -> None:
    """When sensitive_mode=True, the SHA-256 gate and ``joblib.load`` must not
    execute at all. Removing the bundled file would normally raise; under
    sensitive_mode construction succeeds with ``ml_available=False``."""
    from src.core.contamination_scoring import ContaminationScorer

    real_path = _bundled_model_path()
    backup = real_path.with_suffix(".joblib.test_backup")
    real_path.rename(backup)
    try:
        scorer = ContaminationScorer(tmp_path, sensitive_mode=True)
    finally:
        backup.rename(real_path)

    assert scorer.sensitive_mode is True
    assert scorer.ml_available is False
    assert scorer.ml_model is None


def test_rotate_contamination_model_script_prints_current_digest() -> None:
    """The rotation helper must print the on-disk digest of the bundled file
    so maintainers can update the code constant after retraining."""
    import runpy
    import sys
    from io import StringIO

    from src.core.contamination_scoring import CONTAMINATION_MODEL_SHA256

    script_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "rotate_contamination_model.py"
    )
    # Run the module-as-main in-process to inspect stdout without spawning
    # a subprocess (keeps the test fast and deterministic).
    original_stdout, original_argv = sys.stdout, sys.argv
    sys.stdout = captured = StringIO()
    sys.argv = ["rotate_contamination_model.py"]
    try:
        with pytest.raises(SystemExit) as exc_info:
            runpy.run_path(str(script_path), run_name="__main__")
    finally:
        sys.stdout = original_stdout
        sys.argv = original_argv

    assert exc_info.value.code == 0
    out = captured.getvalue()
    assert CONTAMINATION_MODEL_SHA256 in out
