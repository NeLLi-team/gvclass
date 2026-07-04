"""Repository hygiene checks for public-facing release metadata."""

from __future__ import annotations

import re
import shlex
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
    # pyproject.toml is the single version provenance (see src/__version__.py).
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    return str(pyproject["project"]["version"])


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
    expected_image_ref = re.compile(
        rf"(?:gvclass/)?gvclass:{re.escape(version)}"
    )
    expected_apptainer_dev_ref = re.compile(
        r"gvclass_v2\.0-dev\.sif|gvclass:v2\.0-dev"
    )
    stale_image_ref = re.compile(
        rf"(?:gvclass/)?gvclass:(?!{re.escape(version)}\b)\d+\.\d+\.\d+"
    )

    assert pixi["workspace"]["version"] == version
    for path in [
        "gvclass-a",
        "src/bin/gvclass_apptainer.py",
        "src/bin/gvclass_apptainer.sh",
    ]:
        text = (REPO_ROOT / path).read_text()
        assert expected_apptainer_dev_ref.search(text) or "gvclass-a" in text, path
        assert not stale_image_ref.search(text), path

    for path in ["src/bin/gvclass_docker.sh", "src/bin/gvclass_shifter.sh"]:
        text = (REPO_ROOT / path).read_text()
        assert expected_image_ref.search(text), path
        assert not stale_image_ref.search(text), path


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


def test_no_pixi_run_gvclass_in_container_assets() -> None:
    """H5: the image ENTRYPOINT already runs `./gvclass`; a `pixi run gvclass`
    prefix in compose/scripts would be forwarded as bogus argv."""
    for rel in [
        "containers/docker/docker-compose.yml",
        "containers/docker/build_containers.sh",
        "containers/docker/Dockerfile",
    ]:
        text = (REPO_ROOT / rel).read_text()
        assert "pixi run gvclass" not in text, rel


def test_compose_command_is_valid_gvclass_argv() -> None:
    """H5: the effective container command must start with the input path
    (`/data`), i.e. valid `./gvclass` argv, not a launcher prefix."""
    compose = yaml.safe_load(
        (REPO_ROOT / "containers" / "docker" / "docker-compose.yml").read_text()
    )
    command = compose["services"]["gvclass"]["command"]
    tokens = command if isinstance(command, list) else shlex.split(command)
    assert tokens, "compose command is empty"
    assert tokens[0] == "/data", tokens
    assert "pixi" not in tokens, tokens


def test_private_and_runtime_workspaces_stay_ignored() -> None:
    gitignore = (REPO_ROOT / ".gitignore").read_text().splitlines()

    for pattern in [
        "resources/",
        "docs/handoffs/",
        "docs/plans/",
        "docs/quality_metrics.md",
        "tasks/",
        ".codexloop/",
        ".codex",
        "benchmarking/contamination/v1_5_0/",
    ]:
        assert pattern in gitignore
