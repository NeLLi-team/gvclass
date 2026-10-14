"""Shared fixtures for tests that need the bundled reference models."""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

RESOURCES_ROOT = REPO_ROOT / "resources"

_CONTAMINATION_MODEL = RESOURCES_ROOT / "contamination" / "model.joblib"

TEST_REFERENCE_LABEL = (
    "NCLDV__test_reference\tViruses|Nucleocytoviricota|Megaviricetes|"
    "Imitervirales|Mimiviridae|Mimivirus|test_reference\n"
)


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
    """Stage bundled models and a nonempty synthetic label for scorer tests.

    Model directories are symlinked so byte-hash checks use the bundled files.
    Existing labels are preserved. Callers can replace the synthetic label with
    their own fixture. Marker resources are linked only when requested.
    """
    db_path.mkdir(parents=True, exist_ok=True)
    labels_file = db_path / "labels.tsv"
    if not labels_file.exists():
        labels_file.write_text(TEST_REFERENCE_LABEL, encoding="utf-8")
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
