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


def test_load_ml_model_raises_on_sha256_mismatch_and_never_reaches_joblib(
    tmp_path: Path,
) -> None:
    """Any tampering with the bundle must fail-closed before ``joblib.load``.

    Also asserts that ``joblib.load`` is never called when the digest check
    fails — closes a TOCTOU regression where a pre-hash/post-load swap could
    still reach the pickle deserializer.
    """
    import src.core.contamination_scoring as scoring

    with patch(
        "src.core.contamination_scoring.CONTAMINATION_MODEL_SHA256",
        "deadbeef" * 8,
    ), patch("joblib.load") as joblib_load:
        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            scoring.ContaminationScorer(tmp_path)
    # Critical invariant: the untrusted bytes must never reach joblib.load.
    joblib_load.assert_not_called()


def _redirect_bundled_model_to(tmp_dir: Path) -> Path:
    """Copy the bundled contamination model into a sibling ``bundled_models``
    directory under ``tmp_dir`` and return the faux
    ``contamination_scoring.py`` file whose ``parents[1]`` points at
    ``tmp_dir``.

    The scorer resolves its bundled-model path via
    ``Path(__file__).resolve().parents[1] / "bundled_models" / CONTAMINATION_MODEL_FILE``.
    Pointing ``__file__`` at a file two directories deep lets each test
    work on its own isolated tree without renaming the shared repo artefact
    (which would race with concurrent pytest workers).
    """
    core_dir = tmp_dir / "src" / "core"
    bundled_dir = tmp_dir / "src" / "bundled_models"
    core_dir.mkdir(parents=True)
    bundled_dir.mkdir(parents=True)
    fake_scoring_module = core_dir / "contamination_scoring.py"
    fake_scoring_module.write_text("# test redirect target\n")
    return fake_scoring_module


def test_load_ml_model_raises_when_bundled_file_missing(tmp_path: Path) -> None:
    """If the bundled model is physically missing, construction must fail
    loudly rather than silently degrading to ``ml_available=False``.

    Uses an isolated ``tmp_path`` tree with no bundled model present so the
    test never touches the shared repo artefact.
    """
    import src.core.contamination_scoring as scoring

    fake_module = _redirect_bundled_model_to(tmp_path)
    # Intentionally do NOT drop a joblib into bundled_dir.
    with patch.object(scoring, "__file__", str(fake_module)):
        with pytest.raises(RuntimeError, match="Contamination model not found"):
            scoring.ContaminationScorer(tmp_path)


def test_sensitive_mode_scorer_skips_model_load(tmp_path: Path) -> None:
    """When sensitive_mode=True, the SHA-256 gate and ``joblib.load`` must not
    execute at all. Even with no bundled file present, construction must
    succeed with ``ml_available=False``.

    Uses an isolated ``tmp_path`` tree so the test does not rename or
    otherwise perturb the shared repo artefact.
    """
    import src.core.contamination_scoring as scoring

    fake_module = _redirect_bundled_model_to(tmp_path)
    with patch.object(scoring, "__file__", str(fake_module)), patch(
        "joblib.load"
    ) as joblib_load:
        scorer = scoring.ContaminationScorer(tmp_path, sensitive_mode=True)

    assert scorer.sensitive_mode is True
    assert scorer.ml_available is False
    assert scorer.ml_model is None
    joblib_load.assert_not_called()


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
