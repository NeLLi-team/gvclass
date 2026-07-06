"""GVClass software version, sourced from the ``pyproject.toml`` provenance.

The version is declared once in ``pyproject.toml`` (``[project] version``). It is
read from there for source/dev runs (the live source of truth) and falls back to
the installed package metadata when the source tree is unavailable, so the version
is never hardcoded in the Python code.
"""

from pathlib import Path
from typing import Optional


def _version_from_pyproject() -> Optional[str]:
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    if not pyproject.is_file():
        return None
    try:
        import tomllib

        with pyproject.open("rb") as handle:
            return str(tomllib.load(handle)["project"]["version"])
    except Exception:
        return None


def _version_from_metadata() -> Optional[str]:
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("gvclass")
        except PackageNotFoundError:
            return None
    except ImportError:
        return None


__version__ = _version_from_pyproject() or _version_from_metadata() or "0.0.0"
