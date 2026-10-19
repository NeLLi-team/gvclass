"""Split a release archive into parts for upload, with the SHA-256 sidecars.

Usage: split_archive.py <archive.tar.gz> <out dir> [--part-mib 128]

Writes ``<name>.partNNN`` (consecutive numbers from 000), ``<name>.sha256`` (the joined
archive) and ``<name>.parts.sha256`` (every part), in ``sha256sum`` format. The
installer (``DatabaseManager``) downloads the parts, checks each, joins them in order
and checks the joined archive. The upload to a record that limits the files per record
(Zenodo: 100) and the length of one request needs parts of 100 MiB or more.
"""

from __future__ import annotations

import argparse
import hashlib
from collections.abc import Sequence
from pathlib import Path

CHUNK = 8 * 1024 * 1024


def split_archive(archive: Path, out_dir: Path, *, part_bytes: int) -> list[Path]:
    """Write the parts and sidecars of ``archive`` to a new ``out_dir``."""
    if part_bytes <= 0:
        msg = f"part size must be positive, got {part_bytes}"
        raise ValueError(msg)
    out_dir.mkdir(parents=True, exist_ok=False)
    whole = hashlib.sha256()
    parts: list[Path] = []
    lines: list[str] = []
    with archive.open("rb") as source:
        while True:
            part = out_dir / f"{archive.name}.part{len(parts):03d}"
            digest = hashlib.sha256()
            written = 0
            with part.open("wb") as handle:
                while written < part_bytes:
                    chunk = source.read(min(CHUNK, part_bytes - written))
                    if not chunk:
                        break
                    handle.write(chunk)
                    digest.update(chunk)
                    whole.update(chunk)
                    written += len(chunk)
            if written == 0:
                part.unlink()
                break
            parts.append(part)
            lines.append(f"{digest.hexdigest()}  {part.name}\n")
    (out_dir / f"{archive.name}.sha256").write_text(
        f"{whole.hexdigest()}  {archive.name}\n", encoding="utf-8"
    )
    (out_dir / f"{archive.name}.parts.sha256").write_text(
        "".join(lines), encoding="utf-8"
    )
    return parts


def main(argv: Sequence[str] | None = None) -> int:
    """Split the archive and print the number of parts."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("archive", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--part-mib", type=int, default=128)
    args = parser.parse_args(argv)
    parts = split_archive(
        args.archive, args.out_dir, part_bytes=args.part_mib * 1024 * 1024
    )
    total = sum(p.stat().st_size for p in parts)
    print(f"{len(parts)} parts, {total} bytes, in {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
