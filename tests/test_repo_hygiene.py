"""Repository hygiene checks for public-facing release metadata."""

from __future__ import annotations

import re
import subprocess
import tomllib
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def _git_ls_files(path: str) -> str:
    return subprocess.check_output(
        ["git", "ls-files", "--stage", "--", path],
        cwd=REPO_ROOT,
        text=True,
    ).strip()


def _software_version() -> str:
    version_file = REPO_ROOT / "src" / "__version__.py"
    match = re.search(
        r'^__version__\s*=\s*"([^"]+)"',
        version_file.read_text(),
        flags=re.MULTILINE,
    )
    assert match is not None
    return match.group(1)


def test_readme_local_links_resolve_to_tracked_files() -> None:
    readme = (REPO_ROOT / "README.md").read_text()
    targets = re.findall(r"\[[^\]]+\]\(([^)]+)\)", readme)

    missing_or_untracked = []
    for target in targets:
        if "://" in target or target.startswith("#"):
            continue
        path = target.split("#", 1)[0]
        if not path:
            continue
        if not _git_ls_files(path):
            missing_or_untracked.append(target)

    assert missing_or_untracked == []


def test_readme_uses_runtime_resource_model_path() -> None:
    readme = (REPO_ROOT / "README.md").read_text()

    assert "src/bundled_models" not in readme
    assert "resources/contamination/model.joblib" in readme


def test_release_version_matches_user_facing_wrappers() -> None:
    version = _software_version()
    pixi = tomllib.loads((REPO_ROOT / "pixi.toml").read_text())

    assert pixi["workspace"]["version"] == version
    for path in [
        "gvclass-a",
        "src/bin/gvclass_apptainer.py",
        "src/bin/gvclass_apptainer.sh",
        "src/bin/gvclass_docker.sh",
        "src/bin/gvclass_shifter.sh",
    ]:
        text = (REPO_ROOT / path).read_text()
        assert f"gvclass:{version}" in text
        assert "gvclass:1.4.3" not in text


def test_gitmodules_does_not_reference_missing_submodules() -> None:
    gitmodules = REPO_ROOT / ".gitmodules"
    if not gitmodules.exists() or not gitmodules.read_text().strip():
        return

    for path in re.findall(r"^\s*path\s*=\s*(\S+)", gitmodules.read_text(), re.M):
        stage = _git_ls_files(path)
        assert stage.startswith("160000 "), path


def test_ci_enforces_static_analysis_and_hygiene() -> None:
    ci = yaml.safe_load((REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text())
    jobs = ci["jobs"]
    maintenance = (
        REPO_ROOT / ".github" / "workflows" / "repository-maintenance.yml"
    ).read_text()

    assert "repo-hygiene" in jobs
    for job in jobs.values():
        for step in job.get("steps", []):
            assert step.get("continue-on-error") is not True
    assert "schedule:" in maintenance
    assert "workflow_dispatch:" in maintenance
    assert "Stale branch report" in maintenance


def test_private_and_runtime_workspaces_stay_ignored() -> None:
    gitignore = (REPO_ROOT / ".gitignore").read_text().splitlines()

    for pattern in [
        "resources/",
        "docs/",
        "tasks/",
        ".codexloop/",
        ".codex",
        "benchmarking/contamination/v1_5_0/",
    ]:
        assert pattern in gitignore
