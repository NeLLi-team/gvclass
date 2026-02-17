#!/usr/bin/env python3
"""Wrapper to run the published GVClass Apptainer image with minimal typing."""

import argparse
import re
import subprocess
from pathlib import Path

DEFAULT_IMAGE = "library://nelligroup-jgi/gvclass/gvclass:1.2.2"
PUBLIC_LIBRARY_URL = "https://library.sylabs.io"
IMAGE_CACHE_DIR = Path.home() / ".cache" / "gvclass" / "images"


def is_library_uri(image: str) -> bool:
    return image.startswith("library://")


def cache_filename_for_library_image(image: str) -> str:
    slug = image.split("://", 1)[1]
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", slug).strip("_")
    if not slug:
        slug = "gvclass_image"
    return f"{slug}.sif"


def ensure_local_image(image: str) -> str:
    if not is_library_uri(image):
        return image

    IMAGE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    local_image = IMAGE_CACHE_DIR / cache_filename_for_library_image(image)
    if local_image.exists():
        return str(local_image)

    pull_cmd = [
        "apptainer",
        "pull",
        "--library",
        PUBLIC_LIBRARY_URL,
        "-F",
        str(local_image),
        image,
    ]
    if subprocess.call(pull_cmd) != 0:
        raise SystemExit(f"Failed to pull Apptainer image: {image}")
    return str(local_image)


def run_container(
    query_path: Path,
    output: Path,
    threads: int,
    image: str,
    tree_method: str = None,
    mode_fast: bool = None,
    max_workers: int = None,
    sensitive_mode: bool = False,
    contigs_mode: bool = False,
) -> int:
    query_abs = query_path.resolve()
    output_abs = output.resolve()

    if not query_abs.exists():
        raise SystemExit(f"Input path not found: {query_abs}")
    if query_abs.is_file() and not contigs_mode:
        raise SystemExit(
            f"Input is a file: {query_abs}. Use --contigs for per-contig mode."
        )
    if contigs_mode and query_abs.is_file():
        if query_abs.suffix.lower() not in {".fna", ".fa", ".fasta"}:
            raise SystemExit(
                "--contigs requires an FNA/FA/FASTA file when input is a file."
            )

    output_abs.mkdir(parents=True, exist_ok=True)
    resolved_image = ensure_local_image(image)

    if query_abs.is_file():
        bind_source = query_abs.parent
        query_in_container = f"/input/{query_abs.name}"
    else:
        bind_source = query_abs
        query_in_container = "/input"

    cmd = [
        "apptainer",
        "run",
        "-B",
        f"{bind_source}:/input",
        "-B",
        f"{output_abs}:/output",
        resolved_image,
        query_in_container,
        "-o",
        "/output",
        "--threads",
        str(threads),
    ]
    if tree_method:
        cmd.extend(["--tree-method", tree_method])
    if mode_fast is not None:
        cmd.append("--mode-fast" if mode_fast else "--extended")
    if max_workers:
        cmd.extend(["-j", str(max_workers)])
    if sensitive_mode:
        cmd.append("--sensitive")
    if contigs_mode:
        cmd.append("--contigs")

    return subprocess.call(cmd)


def _add_base_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "query_path",
        help="Input directory containing FASTA files, or single FNA file with --contigs",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help="Output directory (legacy positional form; optional)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir_flag",
        default=None,
        help="Output directory for GVClass results (default: <query>_results)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="Total threads to allocate (default: 16)",
    )
    parser.add_argument(
        "-j",
        "--max-workers",
        type=int,
        default=None,
        help="Maximum parallel workers (default: auto-calculated)",
    )


def _add_mode_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tree-method",
        choices=["fasttree", "iqtree"],
        default=None,
        help="Tree building method (default: fasttree)",
    )
    parser.add_argument(
        "--mode-fast",
        action="store_true",
        default=None,
        help="Enable fast mode - skip order-level markers for faster analysis",
    )
    parser.add_argument(
        "-e",
        "--extended",
        action="store_false",
        dest="mode_fast",
        help="Extended mode: build trees for all markers (slower but more comprehensive)",
    )
    parser.add_argument(
        "--sensitive",
        action="store_true",
        help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
    )
    parser.add_argument(
        "-C",
        "--contigs",
        action="store_true",
        help="Per-contig mode: split input FNA(s) and classify each contig separately",
    )
    parser.add_argument(
        "--image",
        default=DEFAULT_IMAGE,
        help=(
            "Apptainer image URI or path. library:// URIs are auto-pulled via "
            f"{PUBLIC_LIBRARY_URL} (default: %(default)s)"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convenience wrapper for running GVClass via Apptainer"
    )
    _add_base_args(parser)
    _add_mode_args(parser)
    return parser


def resolve_output_path(
    query_path: Path, output_dir_positional: str | None, output_dir_flag: str | None
) -> Path:
    if output_dir_positional and output_dir_flag:
        raise ValueError("Use either positional output_dir or -o/--output-dir, not both.")

    resolved_output = output_dir_flag or output_dir_positional
    if resolved_output:
        return Path(resolved_output)

    query_name = query_path.stem if query_path.is_file() else query_path.name
    return query_path.with_name(f"{query_name}_results")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    query_path = Path(args.query_path)
    try:
        output_path = resolve_output_path(
            query_path, args.output_dir, args.output_dir_flag
        )
    except ValueError as exc:
        parser.error(str(exc))

    rc = run_container(
        query_path,
        output_path,
        args.threads,
        args.image,
        tree_method=args.tree_method,
        mode_fast=args.mode_fast,
        max_workers=args.max_workers,
        sensitive_mode=args.sensitive,
        contigs_mode=args.contigs,
    )
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
