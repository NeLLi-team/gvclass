"""Load the ``gvclass`` launcher script as a module, for the tests."""

from __future__ import annotations

import importlib.machinery
import importlib.util
from pathlib import Path
from types import ModuleType

LAUNCHER_PATH = Path(__file__).resolve().parents[1] / "gvclass"


def load_launcher(path: Path = LAUNCHER_PATH) -> ModuleType:
    """Return a fresh module of the launcher script, so a test can patch it freely."""
    loader = importlib.machinery.SourceFileLoader("gvclass_launcher", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module
