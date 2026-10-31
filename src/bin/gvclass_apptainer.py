#!/usr/bin/env python3
"""Wrapper to run the published GVClass Apptainer image with minimal typing."""

import argparse
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.parse import urlparse

DEFAULT_IMAGE = "oras://ghcr.io/nelli-team/gvclass:2.1.1"
IMAGE_CACHE_DIR = Path.home() / ".cache" / "gvclass" / "images"
DEFAULT_RESOURCE_CACHE_DIR = (
    Path.home() / ".cache" / "gvclass" / "resource-cache" / "v2.1.0"
)
CONTAINER_RESOURCE_CACHE_DIR = "/resource-cache"
CONTAINER_DATABASE_DIR = "/opt/gvclass/resources"


def is_remote_image_uri(image: str) -> bool:
    return image.startswith(
        ("oras://", "docker://", "library://", "http://", "https://")
    )


def cache_filename_for_image(image: str) -> str:
    if image.startswith(("http://", "https://")):
        slug = Path(urlparse(image).path).name
    else:
        slug = image.split("://", 1)[1]
    # Two addresses can share a readable name: one basename on two hosts, one path under
    # oras:// and docker://, or names that differ only in punctuation. The digest of the
    # whole address keeps their caches apart.
    prefix = hashlib.sha256(image.encode("utf-8")).hexdigest()[:12] + "_"
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", slug).strip("_")
    if not slug:
        slug = "gvclass_image"
    if not slug.endswith(".sif"):
        slug = f"{slug}.sif"
    return prefix + slug


def validate_cached_image(path: Path) -> None:
    """Reject non-SIF cache entries and objects beyond the downloaded file."""
    if path.is_symlink() or not path.is_file():
        raise SystemExit(f"Invalid Apptainer image: {path}")
    try:
        inspected = subprocess.run(
            ["apptainer", "sif", "list", str(path)],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise SystemExit(
            f"Could not inspect Apptainer image: {path}: {error}"
        ) from error
    objects = re.findall(
        r"^\s*\d+\s*\|[^|]*\|[^|]*\|\s*(\d+)-(\d+)\s*\|([^\n]+)$",
        inspected.stdout,
        re.MULTILINE,
    )
    size = path.stat().st_size
    if (
        inspected.returncode != 0
        or not objects
        or any(not 0 <= int(start) <= int(end) <= size for start, end, _ in objects)
        or not any(
            "FS (" in kind and "/*System/" in kind and int(start) < int(end)
            for start, end, kind in objects
        )
    ):
        raise SystemExit(f"Invalid Apptainer image: {path}")


def ensure_local_image(image: str) -> str:
    if not is_remote_image_uri(image):
        return image

    IMAGE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    local_image = IMAGE_CACHE_DIR / cache_filename_for_image(image)
    if local_image.exists() or local_image.is_symlink():
        validate_cached_image(local_image)
        return str(local_image)

    with tempfile.TemporaryDirectory(prefix=".pull-", dir=IMAGE_CACHE_DIR) as staging:
        staged_image = Path(staging) / "image.sif"
        pull_cmd = ["apptainer", "pull", "-F", str(staged_image), image]
        if subprocess.call(pull_cmd) != 0:
            raise SystemExit(f"Failed to pull Apptainer image: {image}")
        validate_cached_image(staged_image)
        try:
            # A hard link installs atomically without replacing a concurrent entry.
            local_image.hardlink_to(staged_image)
        except FileExistsError:
            validate_cached_image(local_image)
        except OSError as error:
            raise SystemExit(
                f"Could not install cached Apptainer image: {local_image}: {error}"
            ) from error
    return str(local_image)


def _addon_mount(resource_cache: Path) -> list[str]:
    """Return the ``-B`` option that mounts the screen add-on, or nothing.

    The root is the one that ``gvclass screen`` uses: ``GVCLASS_DB`` if set, else the
    resource cache. The add-on is mounted read-only into the database of the image,
    where the screen step looks for it; without one, ``screen_status`` reads
    ``not_installed``.
    """
    root = Path(os.environ.get("GVCLASS_DB", "").strip() or resource_cache)
    addon = root.expanduser().resolve() / "contig_screen"
    if not (addon / "manifest.json").is_file():
        return []
    return ["-B", f"{addon}:{CONTAINER_DATABASE_DIR}/contig_screen:ro"]


def run_container(
    query_path: Path,
    output: Path,
    threads: int,
    image: str,
    resource_cache: Path,
    tree_method: str | None = None,
    mode_fast: bool | None = None,
    max_workers: int | None = None,
    sensitive_mode: bool = False,
    no_sensitive: bool = False,
    contigs_mode: bool = False,
    min_length: int | None = None,
    contigs_min_length: int | None = None,
    threads_per_worker: int | None = None,
    iqtree_mode: str | None = None,
    completeness_mode: str | None = None,
    species_tree: bool = False,
    species_tree_combined: bool = False,
    species_tree_trim: str | None = None,
    resume: bool = False,
    allow_short: bool = False,
    plain_output: bool = False,
    verbose: bool = False,
    *,
    no_screen: bool = False,
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
        if query_abs.suffix.lower() not in {".fna", ".fa", ".fasta", ".fas"}:
            raise SystemExit(
                "--contigs requires an FNA/FA/FASTA/FAS file when input is a file."
            )

    output_abs.mkdir(parents=True, exist_ok=True)
    resource_cache_abs = resource_cache.expanduser().resolve()
    resource_cache_abs.mkdir(parents=True, exist_ok=True)
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
        f"{bind_source}:/input:ro",
        "-B",
        f"{output_abs}:/output",
        "-B",
        f"{resource_cache_abs}:{CONTAINER_RESOURCE_CACHE_DIR}",
        *_addon_mount(resource_cache_abs),
        "--env",
        f"GVCLASS_RESOURCE_CACHE={CONTAINER_RESOURCE_CACHE_DIR}",
        # Apptainer passes the host environment; a host GVCLASS_DB would otherwise
        # replace the database of the image, which holds the mounted add-on.
        "--env",
        f"GVCLASS_DB={CONTAINER_DATABASE_DIR}",
        resolved_image,
        query_in_container,
        "-o",
        "/output",
        "--threads",
        str(threads),
    ]

    # Add optional flags
    if tree_method:
        cmd.extend(["--tree-method", tree_method])
    if iqtree_mode:
        cmd.extend(["--iqtree-mode", iqtree_mode])
    if mode_fast is not None:
        cmd.append("--mode-fast" if mode_fast else "--extended")
    if max_workers:
        cmd.extend(["-j", str(max_workers)])
    if threads_per_worker:
        cmd.extend(["--threads-per-worker", str(threads_per_worker)])
    # Same precedence as the native CLI: --no-sensitive wins over --sensitive.
    if no_sensitive:
        cmd.append("--no-sensitive")
    elif sensitive_mode:
        cmd.append("--sensitive")
    if completeness_mode:
        cmd.extend(["--completeness-mode", completeness_mode])
    if contigs_mode:
        cmd.append("--contigs")
    if min_length is not None:
        cmd.extend(["--min-length", str(min_length)])
    if contigs_min_length is not None:
        cmd.extend(["--contigs-min-length", str(contigs_min_length)])
    if species_tree:
        cmd.append("--species-tree")
    if species_tree_combined:
        cmd.append("--species-tree-combined")
    if species_tree_trim:
        cmd.extend(["--species-tree-trim", species_tree_trim])
    if resume:
        cmd.append("--resume")
    if allow_short:
        cmd.append("--allow-short")
    if plain_output:
        cmd.append("--plain-output")
    if verbose:
        cmd.append("--verbose")
    if no_screen:
        cmd.append("--no-screen")

    return subprocess.call(cmd)


def _add_base_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "query_path",
        nargs="?",
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
    parser.add_argument(
        "--threads-per-worker",
        type=int,
        default=None,
        help="Threads per worker (default: auto-calculated)",
    )


def _add_mode_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--apptainer",
        action="store_true",
        help="Run GVClass via Apptainer container",
    )
    parser.add_argument(
        "--tree-method",
        choices=["veryfasttree", "iqtree", "fasttree"],
        default=None,
        help="Tree building method (default: veryfasttree)",
    )
    parser.add_argument(
        "--iqtree-mode",
        choices=["fast", "ufboot"],
        default=None,
        help="IQ-TREE species-tree search mode (default: fast)",
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
        action="store_true",
        help="Extended mode: build trees for all markers (slower but more comprehensive)",
    )
    parser.add_argument(
        "--sensitive",
        action="store_true",
        help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
    )
    parser.add_argument(
        "--no-sensitive",
        action="store_true",
        help="Disable sensitive HMM mode (use model-defined GA cutoffs when available)",
    )
    parser.add_argument(
        "--completeness-mode",
        choices=["legacy", "novelty-aware"],
        default=None,
        help="Completeness estimator to surface as the primary estimate",
    )
    parser.add_argument(
        "-C",
        "--contigs",
        action="store_true",
        help="Per-contig mode: split input FNA(s) and classify each contig separately",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Minimum total nucleotide length (bp) for each input FNA",
    )
    parser.add_argument(
        "--contigs-min-length",
        type=int,
        default=None,
        help="Minimum contig length (bp) when splitting inputs with --contigs",
    )
    parser.add_argument(
        "--species-tree",
        action="store_true",
        help="Build one supermatrix species tree per query",
    )
    parser.add_argument(
        "--species-tree-combined",
        action="store_true",
        help="Implies --species-tree; also build a combined species tree across queries",
    )
    parser.add_argument(
        "--species-tree-trim",
        choices=["witchi", "pytrimal", "none"],
        default=None,
        help="Supermatrix column trimming for --species-tree",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous run, skipping completed queries",
    )
    parser.add_argument(
        "--allow-short",
        action="store_true",
        help="Bypass the nucleotide length floor",
    )
    parser.add_argument(
        "--plain-output",
        action="store_true",
        help="Disable emojis and ANSI colors in CLI output",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output",
    )
    parser.add_argument(
        "--no-screen",
        action="store_true",
        help="Skip the contig screen (screen_status is 'disabled')",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version information",
    )
    parser.add_argument(
        "--image",
        default=DEFAULT_IMAGE,
        help=(
            "Apptainer image address or path. oras://, docker://, library:// and "
            "https:// addresses are pulled and cached locally; a library:// address "
            "needs a library endpoint that Apptainer knows (default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--resource-cache-dir",
        default=str(DEFAULT_RESOURCE_CACHE_DIR),
        help=(
            "Host directory for Parquet materialization cache (default: %(default)s)"
        ),
    )


def build_parser(*, add_help: bool = True) -> argparse.ArgumentParser:
    """Build the classification parser, optionally without its help action."""
    parser = argparse.ArgumentParser(
        description="Run GVClass via Apptainer container", add_help=add_help
    )
    _add_base_args(parser)
    _add_mode_args(parser)
    return parser


def resolve_output_path(
    query_path: Path, output_dir_positional: str | None, output_dir_flag: str | None
) -> Path:
    if output_dir_positional and output_dir_flag:
        raise ValueError(
            "Use either positional output_dir or -o/--output-dir, not both."
        )

    resolved_output = output_dir_flag or output_dir_positional
    if resolved_output:
        return Path(resolved_output)

    resolved_query = query_path.resolve()
    query_name = (
        resolved_query.stem if resolved_query.is_file() else resolved_query.name
    )
    return resolved_query.with_name(f"{query_name}_results")


def _build_screen_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run contig screens via Apptainer")
    parser.add_argument("command", choices=["screen"], help="Run contig screens")
    parser.add_argument("input", nargs="?", type=Path, help="FASTA file or folder")
    parser.add_argument(
        "--install", action="store_true", help="Install screen resources"
    )
    parser.add_argument("-o", "--output-dir", type=Path, help="Default: <input>_screen")
    parser.add_argument("-d", "--database", help="Host database (default: GVCLASS_DB)")
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="Threads (default: 4)"
    )
    parser.add_argument("--apptainer", action="store_true", help="Use Apptainer")
    parser.add_argument(
        "--image", default=DEFAULT_IMAGE, help="Image URI or local path"
    )
    parser.add_argument(
        "--resource-cache-dir",
        type=Path,
        default=DEFAULT_RESOURCE_CACHE_DIR,
        help="Host cache for the optional screen resources",
    )
    return parser


def _screen_database_bind(args: argparse.Namespace) -> tuple[str, str]:
    host_database = args.database or os.environ.get("GVCLASS_DB", "").strip()
    root = Path(host_database or args.resource_cache_dir).expanduser().resolve()
    if args.install:
        root.mkdir(parents=True, exist_ok=True)
        return f"{root}:/screen-resources", "/screen-resources"
    if not (root / "contig_screen" / "manifest.json").is_file():
        raise SystemExit(
            f"No contig-screen resources in {root}; run: gvclass screen --install"
        )
    if host_database:
        return f"{root}:/screen-resources:ro", "/screen-resources"
    return (
        f"{root}/contig_screen:/opt/gvclass/resources/contig_screen:ro",
        "/opt/gvclass/resources",
    )


def _screen_main(argv: list[str]) -> int:
    parser = _build_screen_parser()
    args = parser.parse_intermixed_args(argv)
    if not args.install and args.input is None:
        parser.error("an input bin or folder is needed (or --install)")
    database_bind, database = _screen_database_bind(args)
    cmd = ["apptainer", "run", "-B", database_bind]
    if args.install:
        cmd.extend(
            [ensure_local_image(args.image), "screen", "--install", "-d", database]
        )
        return subprocess.call(cmd)
    source = args.input.expanduser().resolve()
    if not source.exists():
        parser.error(f"Input path not found: {source}")
    output = args.output_dir or source.with_name(source.name + "_screen")
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    bind_source = source.parent if source.is_file() else source
    container_input = f"/input/{source.name}" if source.is_file() else "/input"
    cmd.extend(
        [
            "-B",
            f"{bind_source}:/input:ro",
            "-B",
            f"{output}:/output",
            ensure_local_image(args.image),
            "screen",
            container_input,
            "-o",
            "/output",
            "-d",
            database,
            "-t",
            str(args.threads),
        ]
    )
    return subprocess.call(cmd)


def main(argv: list[str] | None = None) -> int:
    """Dispatch screening or classification and return the container exit code."""
    given = sys.argv[1:] if argv is None else list(argv)
    probe = build_parser(add_help=False)
    probe.add_argument("-d", "--database")
    probe.add_argument("--install", action="store_true")
    detected, _ = probe.parse_known_args(given)
    if detected.query_path == "screen":
        return _screen_main(given)
    parser = build_parser()
    args = parser.parse_args(given)
    if args.version:
        print(f"GVClass Apptainer launcher (default image: {DEFAULT_IMAGE})")
        return 0
    if not args.query_path:
        if given:
            parser.error("the query path is required")
        parser.print_help()
        return 0

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
        args.image or DEFAULT_IMAGE,
        Path(args.resource_cache_dir or DEFAULT_RESOURCE_CACHE_DIR),
        tree_method=args.tree_method,
        mode_fast=False if args.extended else args.mode_fast,
        max_workers=args.max_workers,
        sensitive_mode=args.sensitive,
        no_sensitive=args.no_sensitive,
        contigs_mode=args.contigs,
        min_length=args.min_length,
        contigs_min_length=args.contigs_min_length,
        threads_per_worker=args.threads_per_worker,
        iqtree_mode=args.iqtree_mode,
        completeness_mode=args.completeness_mode,
        species_tree=args.species_tree,
        species_tree_combined=args.species_tree_combined,
        species_tree_trim=args.species_tree_trim,
        resume=args.resume,
        allow_short=args.allow_short,
        plain_output=args.plain_output,
        verbose=args.verbose,
        no_screen=args.no_screen,
    )
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
