"""The release workflow publishes the SIF to GHCR and checks an anonymous pull."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

WORKFLOW = Path(__file__).resolve().parents[1] / ".github" / "workflows" / "ci.yml"
ADDRESS = "oras://ghcr.io/nelli-team/gvclass"
TAG_PUSH = "github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')"


@pytest.fixture(scope="module")
def workflow() -> dict[str, Any]:
    return yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def steps(workflow: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {step["name"]: step for step in workflow["jobs"]["container-build"]["steps"]}


def test_only_the_container_job_may_write_packages(workflow: dict[str, Any]) -> None:
    # The workflow default is read-only; only the container job adds a write grant.
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["jobs"]["container-build"]["permissions"] == {
        "contents": "read",
        "packages": "write",
    }
    for name, job in workflow["jobs"].items():
        if name != "container-build":
            assert "packages" not in job.get("permissions", {}), name


@pytest.mark.parametrize(
    "step_name", ["Publish SIF to GHCR", "Verify anonymous pull of the published SIF"]
)
@pytest.mark.parametrize(
    ("event", "ref", "expected"),
    [
        ("push", "refs/tags/v2.1.1", True),
        ("push", "refs/heads/main", False),
        ("workflow_dispatch", "refs/heads/main", False),
        ("workflow_dispatch", "refs/tags/v2.1.1", False),
    ],
)
def test_publication_requires_a_tag_push(
    steps: dict[str, dict[str, Any]],
    step_name: str,
    event: str,
    ref: str,
    expected: bool,
) -> None:
    """Check both release-step predicates against real event/ref combinations."""
    predicates = {
        "github.event_name == 'push'": event == "push",
        "startsWith(github.ref, 'refs/tags/v')": ref.startswith("refs/tags/v"),
    }
    terms = steps[step_name]["if"].split("&&")
    assert all(predicates[term.strip()] for term in terms) is expected


def test_package_writing_job_does_not_persist_checkout_credentials(
    workflow: dict[str, Any],
) -> None:
    """Keep the package-writing token out of the checkout used by image tests."""
    checkout = next(
        step
        for step in workflow["jobs"]["container-build"]["steps"]
        if step.get("uses", "").startswith("actions/checkout@")
    )
    assert checkout.get("with", {}).get("persist-credentials") is False


def test_nothing_of_sylabs_remains_in_the_workflow() -> None:
    text = WORKFLOW.read_text(encoding="utf-8").lower()

    assert "sylabs" not in text
    assert "library://" not in text


def test_the_publish_step_logs_in_and_pushes_the_version_and_latest(
    steps: dict[str, dict[str, Any]],
) -> None:
    step = steps["Publish SIF to GHCR"]
    script = step["run"]

    assert step["if"] == TAG_PUSH
    assert step["env"] == {"GITHUB_TOKEN": "${{ secrets.GITHUB_TOKEN }}"}
    assert "--password-stdin" in script
    assert '--username "${GITHUB_ACTOR}"' in script
    assert "apptainer registry login" in script
    assert "oras://ghcr.io\n" in script or "oras://ghcr.io " in script
    assert f'apptainer push gvclass.sif "{ADDRESS}:${{version}}"' in script
    assert f'apptainer push gvclass.sif "{ADDRESS}:latest"' in script


def test_the_publish_follows_the_artifact_upload_and_precedes_its_check(
    steps: dict[str, dict[str, Any]],
) -> None:
    names = list(steps)

    assert names.index("Upload SIF and SBOM artifacts") < names.index(
        "Publish SIF to GHCR"
    )
    assert names.index("Publish SIF to GHCR") + 1 == names.index(
        "Verify anonymous pull of the published SIF"
    )


def test_the_check_pulls_anonymously_and_compares_the_bytes_of_the_built_file(
    steps: dict[str, dict[str, Any]],
) -> None:
    step = steps["Verify anonymous pull of the published SIF"]
    script = step["run"]

    assert step["if"] == TAG_PUSH
    assert "registry login" not in script
    assert "apptainer registry logout oras://ghcr.io" in script
    assert "APPTAINER_CONFIGDIR" in script
    assert "APPTAINER_CACHEDIR" in script
    assert f'"{ADDRESS}:${{version}}"' in script
    assert script.count("sha256sum") == 2
    assert "gvclass.sif" in script
    assert "private" in script


def test_the_failure_of_the_anonymous_pull_stays_visible(
    steps: dict[str, dict[str, Any]],
) -> None:
    step = steps["Verify anonymous pull of the published SIF"]

    assert "continue-on-error" not in step
    assert "|| true" not in step["run"]
    assert 'exit "${status}"' in step["run"]


FAKE_APPTAINER = """#!/bin/sh
printf '%s|%s\\n' "$*" "${APPTAINER_CONFIGDIR:-}" >> "$FAKE_LOG"
if [ "$1 $2" = "registry login" ]; then cat > "$FAKE_STDIN"; fi
if [ "$1" = pull ]; then
  [ "${FAKE_PULL_STATUS:-0}" = 0 ] || exit "$FAKE_PULL_STATUS"
  cp "$FAKE_PULLED" "$2"
fi
"""
TOKEN = "ghs_not_a_real_token"  # noqa: S105  # a test value


class Runner:
    """Run a workflow step script with bash and a fake apptainer, as the runner does."""

    def __init__(self, tmp_path: Path) -> None:
        self.tmp = tmp_path
        self.work = tmp_path / "work"
        self.work.mkdir()
        (self.work / "gvclass.sif").write_bytes(b"built sif")
        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        (fake_bin / "apptainer").write_text(FAKE_APPTAINER, encoding="utf-8")
        (fake_bin / "apptainer").chmod(0o755)
        (tmp_path / "temp").mkdir()
        (tmp_path / "pulled.sif").write_bytes(b"built sif")
        self.log = tmp_path / "calls.log"
        self.stdin = tmp_path / "login.stdin"
        self.env = {
            "PATH": f"{fake_bin}:/usr/bin:/bin",
            "HOME": str(tmp_path / "home"),
            "RUNNER_TEMP": str(tmp_path / "temp"),
            "GITHUB_REF_NAME": "v2.1.1",
            "GITHUB_ACTOR": "octo-cat",
            "GITHUB_TOKEN": TOKEN,
            "FAKE_LOG": str(self.log),
            "FAKE_STDIN": str(self.stdin),
            "FAKE_PULLED": str(tmp_path / "pulled.sif"),
        }

    def run(self, script: str, **env: str) -> subprocess.CompletedProcess[str]:
        file = self.tmp / "step.sh"
        file.write_text(script, encoding="utf-8")
        return subprocess.run(
            ["bash", "-e", str(file)],
            capture_output=True,
            text=True,
            check=False,
            cwd=self.work,
            env={**self.env, **env},
        )

    def calls(self) -> list[tuple[str, str]]:
        lines = self.log.read_text(encoding="utf-8").splitlines()
        return [tuple(line.split("|", 1)) for line in lines]  # type: ignore[misc]


@pytest.fixture
def runner(tmp_path: Path) -> Runner:
    return Runner(tmp_path)


def test_the_publish_script_sends_the_token_on_stdin_and_pushes_both_tags(
    steps: dict[str, dict[str, Any]], runner: Runner
) -> None:
    done = runner.run(steps["Publish SIF to GHCR"]["run"])

    assert done.returncode == 0, done.stderr
    assert [call for call, _ in runner.calls()] == [
        "registry login --username octo-cat --password-stdin oras://ghcr.io",
        f"push gvclass.sif {ADDRESS}:2.1.1",
        f"push gvclass.sif {ADDRESS}:latest",
        "registry logout oras://ghcr.io",
    ]
    assert runner.stdin.read_text(encoding="utf-8") == TOKEN
    assert TOKEN not in runner.log.read_text(encoding="utf-8")
    assert TOKEN not in done.stdout + done.stderr


def verify(steps: dict[str, dict[str, Any]], runner: Runner, **env: str):
    return runner.run(steps["Verify anonymous pull of the published SIF"]["run"], **env)


def test_the_check_passes_when_the_anonymous_pull_has_the_bytes_of_the_built_file(
    steps: dict[str, dict[str, Any]], runner: Runner
) -> None:
    done = verify(steps, runner)

    assert done.returncode == 0, done.stderr
    [logout, pull] = runner.calls()
    assert logout[0] == "registry logout oras://ghcr.io"
    command, config_dir = pull
    assert command.startswith("pull ")
    assert command.endswith(f" {ADDRESS}:2.1.1")
    assert command.split()[1].startswith(str(runner.tmp / "temp"))
    assert config_dir.startswith(str(runner.tmp / "temp"))
    assert not (runner.tmp / "home").exists()


def test_the_check_fails_when_the_pulled_bytes_differ(
    steps: dict[str, dict[str, Any]], runner: Runner
) -> None:
    (runner.tmp / "pulled.sif").write_bytes(b"another sif")

    done = verify(steps, runner)

    assert done.returncode != 0
    assert "built " in done.stdout
    assert "pulled " in done.stdout


def test_the_check_fails_with_the_exit_code_of_a_refused_pull_and_names_the_cause(
    steps: dict[str, dict[str, Any]], runner: Runner
) -> None:
    done = verify(steps, runner, FAKE_PULL_STATUS="7")

    assert done.returncode == 7
    assert "private" in done.stderr
    assert "visibility" in done.stderr
