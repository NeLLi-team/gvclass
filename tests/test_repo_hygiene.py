"""Repository hygiene checks for public-facing release metadata."""

from __future__ import annotations

import re
import shlex
import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest
import tomllib
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
    stale_image_ref = re.compile(
        rf"(?:gvclass/)?gvclass:(?!{re.escape(version)}\b)\d+\.\d+\.\d+"
    )

    for path in [
        "pyproject.toml",
        "pixi.toml",
        "gvclass",
        "containers/build.sh",
        "containers/apptainer/gvclass.def",
        "containers/docker/Dockerfile",
        "containers/docker/docker-compose.yml",
    ]:
        text = (REPO_ROOT / path).read_text()
        assert version in text, path
        assert not stale_image_ref.search(text), path

    pixi = tomllib.loads((REPO_ROOT / "pixi.toml").read_text())
    assert pixi["workspace"]["version"] == version
    assert f"gvclass:{version}" in pixi["tasks"]["build-docker"]
    assert f"gvclass_{version}.sif" in pixi["tasks"]["build-apptainer"]


def test_tracked_files_hold_no_lbl_address_but_the_licensing_office() -> None:
    local_part = re.compile(
        r"(?<![\w.+-])([\w.+-]+)@(?:[\w-]+\.)*lbl\.gov", re.IGNORECASE
    )
    tracked = subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=REPO_ROOT, text=True
    ).split("\0")

    found = {}
    for name in filter(None, tracked):
        path = REPO_ROOT / name
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for match in local_part.findall(text):
            found.setdefault(name, set()).add(match.lower())

    assert found == {"LICENCE": {"ipo"}}


def test_package_metadata_names_authors_without_addresses() -> None:
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    pixi = tomllib.loads((REPO_ROOT / "pixi.toml").read_text())

    assert pyproject["project"]["authors"] == [{"name": "Frederik Schulz"}]
    assert pixi["workspace"]["authors"] == ["Frederik Schulz"]
    assert "click" not in pixi["dependencies"]
    assert "click" not in pixi["feature"]["dev"]["dependencies"]


def test_containers_ignore_host_python_packages() -> None:
    definition = (REPO_ROOT / "containers" / "apptainer" / "gvclass.def").read_text()
    dockerfile = (REPO_ROOT / "containers" / "docker" / "Dockerfile").read_text()

    assert "export PYTHONPATH=/opt/gvclass\n" in definition
    assert "export PYTHONNOUSERSITE=1\n" in definition
    assert 'ENV PYTHONPATH="/opt/gvclass"\n' in dockerfile
    assert 'ENV PYTHONNOUSERSITE="1"\n' in dockerfile
    assert "maintainer" not in dockerfile
    assert (
        'LABEL org.opencontainers.image.source="https://github.com/NeLLi-team/gvclass"'
    ) in dockerfile


def _apptainer_output_snippet() -> str:
    """Return the shell block of the runscript that adds ``-o /output``."""
    text = (REPO_ROOT / "containers" / "apptainer" / "gvclass.def").read_text()
    match = re.search(r"^ *has_output=0\n.*?\n *fi\n", text, re.DOTALL | re.MULTILINE)
    assert match, "the runscript has no has_output block"
    return match.group(0)


def _docker_output_snippet() -> str:
    """Return the shell block of the Docker entry point that adds ``-o /output``."""
    lines = (
        (REPO_ROOT / "containers" / "docker" / "Dockerfile").read_text().splitlines()
    )
    start = next(i for i, line in enumerate(lines) if "'has_output=0'" in line)
    end = next(
        i for i, line in enumerate(lines[start:], start) if '$has_output" = 0' in line
    )
    return "\n".join(
        line.strip().removesuffix("\\").strip().strip("'")
        for line in lines[start : end + 1]
    )


def _arguments_after(shell: str, snippet: str, args: list[str]) -> list[str]:
    """Run ``snippet`` with ``args`` as its positional parameters; return them after."""
    script = snippet + '\nprintf "%s\\n" "$@"'
    done = subprocess.run(
        [shell, "-c", script, "entry-point", *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return done.stdout.splitlines()


OUTPUT_OPTION_FORMS = [
    ["-o", "/results"],
    ["-o/results"],
    ["--output-dir", "/results"],
    ["--output-dir=/results"],
]


ENTRY_POINTS = pytest.mark.parametrize(
    ("shell", "snippet"),
    [
        pytest.param("sh", _apptainer_output_snippet, id="apptainer"),
        pytest.param("bash", _docker_output_snippet, id="docker"),
    ],
)


@ENTRY_POINTS
def test_container_entry_point_sends_results_to_the_output_mount_by_default(
    shell: str, snippet: Callable[[], str]
) -> None:
    args = ["/data/input", "-t", "8"]

    assert _arguments_after(shell, snippet(), args) == [*args, "-o", "/output"]


@ENTRY_POINTS
@pytest.mark.parametrize("option", OUTPUT_OPTION_FORMS)
def test_container_entry_point_keeps_every_form_of_the_output_option(
    shell: str, snippet: Callable[[], str], option: list[str]
) -> None:
    args = ["/data/input", *option]

    assert _arguments_after(shell, snippet(), args) == args


def test_legacy_container_build_script_is_removed() -> None:
    assert not (REPO_ROOT / "containers" / "docker" / "build_containers.sh").exists()


def test_container_mode_has_one_implementation_in_the_launcher() -> None:
    """The container mode lives only in ``gvclass``, which users download alone."""
    launcher = (REPO_ROOT / "gvclass").read_text()

    assert not (REPO_ROOT / "src" / "bin" / "gvclass_apptainer.py").exists()
    for line in [
        f'DEFAULT_IMAGE = "oras://ghcr.io/nelli-team/gvclass:{_software_version()}"',
        "def ensure_local_image(image: str) -> str:",
        "def container_main(",
    ]:
        assert line in launcher, f"Missing {line!r} in gvclass"


def test_gitmodules_does_not_reference_missing_submodules() -> None:
    gitmodules = REPO_ROOT / ".gitmodules"
    if not gitmodules.exists() or not gitmodules.read_text().strip():
        return

    for path in re.findall(
        r"^\s*path\s*=\s*(\S+)", gitmodules.read_text(), re.MULTILINE
    ):
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
    assert steps_by_name["Generate SBOM"]["with"]["image"] == "./gvclass.sif"
    artifact_step = steps_by_name["Upload SIF and SBOM artifacts"]
    assert "gvclass.sif" in artifact_step["with"]["path"]
    step_names = list(steps_by_name)
    assert step_names.index("Upload SIF and SBOM artifacts") < step_names.index(
        "Publish SIF to GHCR"
    )


def test_ci_workflow_runs_checked_in_scripts_and_the_launcher() -> None:
    text = (REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text()
    workflow = yaml.safe_load(text)

    assert workflow["permissions"] == {"contents": "read"}
    for script in set(re.findall(r"scripts/ci/[\w.]+\.py", text)):
        assert (REPO_ROOT / script).is_file(), script
    for removed in [
        "gvclass_apptainer",
        "src.bin.gvclass_cli",
        "stable_columns = [",
        '["started_at"]',
        '["completed_at"]',
        "--cluster",
        ".SUCCESS",
    ]:
        assert removed not in text, removed
    heredocs = re.findall(r"<<'(\w+)'\n(.*?)\n\s*\1\n", text, re.DOTALL)
    for _, body in heredocs:
        assert body.count("\n") < 30, "long Python programs belong in scripts/ci"


def test_no_pixi_run_gvclass_in_container_assets() -> None:
    """H5: the image ENTRYPOINT already runs `./gvclass`; a `pixi run gvclass`
    prefix in compose/scripts would be forwarded as bogus argv."""
    for rel in [
        "containers/docker/docker-compose.yml",
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
