#!/usr/bin/env python3
"""
Install GVClass CLI entrypoints into ~/bin.
"""

from pathlib import Path
import shutil
import sys


def main() -> int:
    repo_dir = Path(__file__).resolve().parents[2]
    bin_dir = Path.home() / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    config_dir = Path.home() / ".config" / "gvclass"
    config_dir.mkdir(parents=True, exist_ok=True)
    repo_path_file = config_dir / "repo_path"

    targets = ["gvclass", "gvclass-a"]
    installed = []
    for name in targets:
        src = repo_dir / name
        if not src.exists():
            print(f"Skipping {name}: not found at {src}")
            continue
        dst = bin_dir / name
        shutil.copy2(src, dst)
        dst.chmod(src.stat().st_mode)
        installed.append(dst)

    if installed:
        repo_path_file.write_text(str(repo_dir) + "\n")

    if not installed:
        print("No CLI scripts installed.")
        return 1

    print("Installed:")
    for path in installed:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
