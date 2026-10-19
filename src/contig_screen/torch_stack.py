"""The optional deep-learning packages of the screen, and the error when absent."""

from __future__ import annotations

import importlib.util
from collections.abc import Iterator
from contextlib import contextmanager

REQUIRED_PACKAGES = ("torch", "transformers")


class MissingTorchError(RuntimeError):
    """torch or transformers is not installed in the running environment."""


def require_torch() -> None:
    """Raise a clear error when torch or transformers cannot be imported.

    The check does not import the packages, so it is cheap and runs before any
    long step of the screen.
    """
    missing = [
        name for name in REQUIRED_PACKAGES if importlib.util.find_spec(name) is None
    ]
    if missing:
        msg = (
            f"the cellular contig screen needs {' and '.join(missing)}, which "
            f"{'is' if len(missing) == 1 else 'are'} not installed here; run it in "
            "the pixi environment of this GVClass release (`pixi run gvclass ...`), "
            "which includes them"
        )
        raise MissingTorchError(msg)


@contextmanager
def limit_threads(threads: int) -> Iterator[None]:
    """Cap the OpenMP and BLAS pools (scikit-learn, numpy) at ``threads``."""
    from threadpoolctl import threadpool_limits  # noqa: PLC0415

    with threadpool_limits(limits=threads):
        yield
