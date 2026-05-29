"""M29: SHA-256 gate on the completeness model pickle.

Mirrors tests/test_contamination_model_sha256.py. A digest mismatch must RAISE
(refuse to deserialize an untrusted pickle); a genuinely missing model must
soft-degrade to available=False (the completeness model is optional, unlike the
mandatory contamination model).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from src.core import novelty_completeness as nc
from src.core.novelty_completeness import NoveltyAwareCompletenessScorer

REPO_ROOT = Path(__file__).resolve().parents[1]
SHIPPED_MODEL = REPO_ROOT / "resources" / "completeness" / "model.joblib"


def _seed_completeness_dir(db: Path, model_bytes: bytes | None) -> None:
    comp = db / "completeness"
    comp.mkdir(parents=True, exist_ok=True)
    # Content does not need to be valid: the digest gate runs before parsing.
    (comp / "config.json").write_text("{}")
    (comp / "tiers.tsv").write_text("group_level\tgroup_name\torder_name\ttier\tmarker_name\n")
    (comp / "baselines.tsv").write_text(
        "group_level\tgroup_name\torder_name\tbaseline_mean\tbaseline_median\tn_refs\n"
    )
    (comp / "model_metadata.tsv").write_text("order_name\ttest_r2\n")
    if model_bytes is not None:
        (comp / "model.joblib").write_bytes(model_bytes)


def test_tampered_completeness_model_raises(tmp_path: Path) -> None:
    _seed_completeness_dir(tmp_path, model_bytes=b"tampered-pickle-bytes")
    with pytest.raises(RuntimeError, match="SHA-256"):
        NoveltyAwareCompletenessScorer(tmp_path)


def test_missing_completeness_model_soft_degrades(tmp_path: Path) -> None:
    # All other files present, but no model.joblib -> optional, soft-degrade.
    _seed_completeness_dir(tmp_path, model_bytes=None)
    scorer = NoveltyAwareCompletenessScorer(tmp_path)
    assert scorer.available is False


def test_completeness_model_sha256_constant_is_set() -> None:
    assert isinstance(nc.COMPLETENESS_MODEL_SHA256, str)
    assert len(nc.COMPLETENESS_MODEL_SHA256) == 64


@pytest.mark.requires_db
def test_shipped_completeness_model_matches_constant() -> None:
    if not SHIPPED_MODEL.exists():
        pytest.skip("completeness model not present (runtime bundle absent)")
    digest = hashlib.sha256(SHIPPED_MODEL.read_bytes()).hexdigest()
    assert digest == nc.COMPLETENESS_MODEL_SHA256
