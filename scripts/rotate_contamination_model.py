#!/usr/bin/env python3
"""Print the SHA-256 digest of the bundled contamination model.

Use this helper after retraining to update the ``CONTAMINATION_MODEL_SHA256``
constant in ``src/core/contamination_scoring.py`` and the ``sha256:`` field
in ``src/bundled_models/contamination_model.yaml``. The pipeline refuses to
load the pickle at runtime unless both values match the on-disk digest, so
rotation goes through code review by design.

Usage
-----

    python scripts/rotate_contamination_model.py
    # or
    python scripts/rotate_contamination_model.py /path/to/model.joblib
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path


DEFAULT_MODEL_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "bundled_models"
    / "contamination_model.joblib"
)


def compute_sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    """Stream the file in 1 MiB chunks to compute its SHA-256 digest."""
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> int:
    if len(sys.argv) == 1:
        model_path = DEFAULT_MODEL_PATH
    elif len(sys.argv) == 2:
        model_path = Path(sys.argv[1]).expanduser().resolve()
    else:
        sys.stderr.write(
            "usage: rotate_contamination_model.py [path/to/model.joblib]\n"
        )
        return 2

    if not model_path.exists():
        sys.stderr.write(f"error: {model_path} not found\n")
        return 1

    digest = compute_sha256(model_path)
    size = model_path.stat().st_size
    print(f"path:   {model_path}")
    print(f"bytes:  {size}")
    print(f"sha256: {digest}")
    print()
    print("Next steps:")
    print("  1. Update CONTAMINATION_MODEL_SHA256 in src/core/contamination_scoring.py")
    print("  2. Update the sha256: field in src/bundled_models/contamination_model.yaml")
    print("  3. Run pixi run test and commit both changes together.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
