#!/usr/bin/env python3
"""Wrapper to run the published GVClass Apptainer image with minimal typing."""

import argparse
import subprocess
from pathlib import Path

DEFAULT_IMAGE = "library://nelligroup-jgi/gvclass/gvclass:1.2.1"


def run_container(
    query: Path, output: Path, threads: int, image: str, sensitive_mode: bool = False
) -> int:
    query_abs = query.resolve()
    output_abs = output.resolve()

    if not query_abs.exists():
        raise SystemExit(f"Input directory not found: {query_abs}")

    output_abs.mkdir(parents=True, exist_ok=True)

    cmd = [
        "apptainer",
        "run",
        "-B",
        f"{query_abs}:/input",
        "-B",
        f"{output_abs}:/output",
        image,
        "/input",
        "-o",
        "/output",
        "--threads",
        str(threads),
    ]
    if sensitive_mode:
        cmd.append("--sensitive")

    return subprocess.call(cmd)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convenience wrapper for running GVClass via Apptainer"
    )
    parser.add_argument(
        "query_dir", help="Directory containing FASTA files to classify"
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help="Directory for GVClass results (default: <query_dir>_results)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="Total threads to allocate (default: 16)",
    )
    parser.add_argument(
        "--image",
        default=DEFAULT_IMAGE,
        help="Apptainer image URI or path (default: %(default)s)",
    )
    parser.add_argument(
        "--sensitive",
        action="store_true",
        help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
    )

    args = parser.parse_args()

    query_path = Path(args.query_dir)
    if args.output_dir:
        output_path = Path(args.output_dir)
    else:
        output_path = query_path.with_name(f"{query_path.name}_results")

    rc = run_container(
        query_path, output_path, args.threads, args.image, args.sensitive
    )
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
