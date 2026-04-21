import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

RESOURCES_ROOT = REPO_ROOT / "resources"

_CONTAMINATION_MODEL = RESOURCES_ROOT / "contamination" / "model.joblib"


def skip_if_no_runtime_resources():
    """Return a pytest.mark.skipif that fires when the runtime resources bundle
    is not installed locally.

    Tests that construct :class:`ContaminationScorer` or :class:`FullSummarizer`
    need the bundled contamination model on disk (v1.5.0 ships it in
    ``resources/contamination/model.joblib``, downloaded via ``pixi run setup-db``
    rather than committed to the repo). On a fresh CI checkout or a dev machine
    that has not yet fetched the resources tarball, that file is missing and
    these tests should be skipped rather than failed.
    """
    return pytest.mark.skipif(
        not _CONTAMINATION_MODEL.exists(),
        reason=(
            "Runtime resources not installed at "
            f"{_CONTAMINATION_MODEL}; run `pixi run setup-db` to fetch the "
            "v1.5.0 resources bundle so these tests can execute."
        ),
    )


def stage_db_resources(db_path: Path, *, markers: bool = True) -> None:
    """Populate a temp database dir with the runtime resources v1.5.0 needs.

    The production scorer resolves the contamination model via
    ``database_path / "contamination" / "model.joblib"`` and the novelty
    scorer resolves the completeness bundle under
    ``database_path / "completeness"``. For tests that construct a scorer
    against a tmp directory, we symlink those two subtrees (and the
    ``labels.tsv`` file) to the real resources tree so the SHA-256 gate
    on the bundled model passes without copying the 7 MB joblib into
    every tmp dir.

    Callers stage their own ``labels.tsv`` or marker fixtures before or
    after calling this helper; the symlinks only cover the dirs whose
    contents tests do not typically customise.
    """
    db_path.mkdir(parents=True, exist_ok=True)
    for name in ("contamination", "completeness"):
        target = RESOURCES_ROOT / name
        if target.exists():
            link = db_path / name
            if link.exists() or link.is_symlink():
                continue
            link.symlink_to(target)
    if markers:
        markers_target = RESOURCES_ROOT / "markers"
        link = db_path / "markers"
        if markers_target.exists() and not (link.exists() or link.is_symlink()):
            link.symlink_to(markers_target)
