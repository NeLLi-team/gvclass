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


def test_release_version_matches_wrappers_and_container_assets() -> None:
    version = _software_version()
    pixi = tomllib.loads((REPO_ROOT / "pixi.toml").read_text())
    expected_apptainer_ref = re.compile(rf"gvclass:{re.escape(version)}")
    stale_image_ref = re.compile(
        rf"(?:gvclass/)?gvclass:(?!{re.escape(version)}\b)\d+\.\d+\.\d+"
    )

    assert pixi["workspace"]["version"] == version
    assert f"gvclass:{version}" in pixi["tasks"]["build-docker"]
    assert f"gvclass_{version}.sif" in pixi["tasks"]["build-apptainer"]
    for path in [
        "gvclass",
        "src/bin/gvclass_apptainer.py",
        "containers/apptainer/gvclass.def",
        "containers/docker/Dockerfile",
        "containers/docker/docker-compose.yml",
        "containers/build.sh",
    ]:
        text = (REPO_ROOT / path).read_text()
        assert expected_apptainer_ref.search(text) or f"Version {version}" in text or f"GVCLASS_VERSION={version}" in text, path
        assert not stale_image_ref.search(text), path

    for path in [
        "containers/build.sh",
        "containers/docker/Dockerfile",
        "containers/docker/docker-compose.yml",
        "containers/docker/build_containers.sh",
        "containers/apptainer/gvclass.def",
    ]:
        text = (REPO_ROOT / path).read_text()
        assert version in text, path
        assert not stale_image_ref.search(text), path


def test_apptainer_wrapper_copies_stay_in_sync() -> None:
    """``gvclass`` embeds the standalone Apptainer runner from ``src/bin/gvclass_apptainer.py``.

    ``gvclass`` can be downloaded standalone via ``wget`` without a repo checkout,
    while ``src/bin/gvclass_apptainer.py`` is the importable module inside ``src/``.
    Verify both files declare the same default image URI, resource cache version,
    and container helper definitions.
    """
    wrapper = (REPO_ROOT / "gvclass").read_text()
    module = (REPO_ROOT / "src" / "bin" / "gvclass_apptainer.py").read_text()

    for line in [
        f'DEFAULT_IMAGE = "library://nelligroup-jgi/gvclass/gvclass:{_software_version()}"',
        f'SYLABS_IMAGE = "library://nelligroup-jgi/gvclass/gvclass:{_software_version()}"',
        'Path.home() / ".cache" / "gvclass" / "resource-cache" / "v2.1.0"',
        "def ensure_local_image(image: str) -> str:",
        "def run_container(",
    ]:
        assert line in wrapper, f"Missing {line!r} in gvclass"
        assert line in module, f"Missing {line!r} in src/bin/gvclass_apptainer.py"


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

    container_job = jobs["container-build"]
    assert set(container_job["needs"]) >= {
        "lint-and-test",
        "repo-hygiene",
        "golden-file",
    }
    steps_by_name = {step["name"]: step for step in container_job["steps"]}
    tag_step = steps_by_name["Verify release tag matches software version"]
    assert "GITHUB_REF_NAME" in tag_step["run"]
    assert "pyproject.toml" in tag_step["run"]
    identity_step = steps_by_name["Verify SIF software, resource, and model identities"]
    assert "DB_VERSION" in identity_step["run"]
    assert 'database_version}" != "v2.1.0"' in identity_step["run"]
    assert "CONTAMINATION_MODEL_SHA256" in identity_step["run"]
    assert "COMPLETENESS_MODEL_SHA256" in identity_step["run"]
    assert "create_experimental_quality_scorer" in identity_step["run"]
    assert steps_by_name["Generate SBOM"]["with"]["image"] == "./gvclass.sif"
    artifact_step = steps_by_name["Upload SIF and SBOM artifacts"]
    assert "gvclass.sif" in artifact_step["with"]["path"]
    step_names = list(steps_by_name)
    assert step_names.index("Upload SIF and SBOM artifacts") < step_names.index(
        "Publish SIF to Sylabs library"
    )


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


def test_container_builds_use_pinned_base_and_repo_compose_context() -> None:
    dockerfile = (REPO_ROOT / "containers" / "docker" / "Dockerfile").read_text()
    definition = (REPO_ROOT / "containers" / "apptainer" / "gvclass.def").read_text()
    compose = yaml.safe_load(
        (REPO_ROOT / "containers" / "docker" / "docker-compose.yml").read_text()
    )
    pinned_base = re.compile(r"ghcr\.io/prefix-dev/pixi@sha256:[0-9a-f]{64}")

    docker_base = pinned_base.search(dockerfile)
    definition_base = pinned_base.search(definition)
    assert docker_base is not None
    assert definition_base is not None
    assert docker_base.group() == definition_base.group()
    assert "RUN pixi install --frozen" in dockerfile
    assert "    pixi install --frozen" in definition
    assert compose["services"]["gvclass"]["build"] == {
        "context": "../..",
        "dockerfile": "containers/docker/Dockerfile",
    }


def test_container_runtime_copy_excludes_development_and_host_files() -> None:
    dockerfile = (REPO_ROOT / "containers" / "docker" / "Dockerfile").read_text()
    definition = (REPO_ROOT / "containers" / "apptainer" / "gvclass.def").read_text()

    assert "COPY example /opt/gvclass/example" in dockerfile
    assert "    example /opt/gvclass/example" in definition
    for excluded in [
        "COPY src /opt/gvclass/src",
        "    src /opt/gvclass/src",
        "src/test",
        "src/bin/build_novelty_completeness_resources.py",
        "src/bin/install_cli.py",
        "src/bin/gvclass_apptainer.py",
        "src/bin/gvclass_apptainer.sh",
        "src/bin/gvclass_docker.sh",
        "src/bin/gvclass_shifter.sh",
        "images /opt/gvclass/images",
        "containers /opt/gvclass/containers",
        "gvclass-a /opt/gvclass/gvclass-a",
    ]:
        assert excluded not in dockerfile, excluded
        assert excluded not in definition, excluded


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
