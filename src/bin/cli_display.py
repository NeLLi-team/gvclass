"""Display helpers for the GVClass CLI."""

from __future__ import annotations

import sys
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.bin.gvclass_cli import ContigInput, WorkerPlan

PLAIN_BANNER_LINES = (
    "    ================================================================",
    "    GVClass - Giant Virus Classification Tool",
    "    ================================================================",
)


def _build_color_banner_lines(software_version: str) -> list[str]:
    reset = "\033[0m"
    c = ["\033[95m", "\033[95m", "\033[91m", "\033[31m", "\033[31m", "\033[38;5;208m", "\033[93m"]
    return [
        f"  {c[0]}██████╗{reset} {c[1]}██╗   ██╗{reset}{c[2]} ██████╗{reset}{c[3]}██╗{reset}      {c[4]}█████╗{reset} {c[5]}███████╗{reset}{c[6]}███████╗{reset}",
        f" {c[0]}██╔════╝{reset} {c[1]}██║   ██║{reset}{c[2]}██╔════╝{reset}{c[3]}██║{reset}     {c[4]}██╔══██╗{reset}{c[5]}██╔════╝{reset}{c[6]}██╔════╝{reset}",
        f" {c[0]}██║  ███╗{reset}{c[1]}██║   ██║{reset}{c[2]}██║{reset}     {c[3]}██║{reset}     {c[4]}███████║{reset}{c[5]}███████╗{reset}{c[6]}███████╗{reset}",
        f" {c[0]}██║   ██║{reset}{c[1]}╚██╗ ██╔╝{reset}{c[2]}██║{reset}     {c[3]}██║{reset}     {c[4]}██╔══██║{reset}{c[5]}╚════██║{reset}{c[6]}╚════██║{reset}",
        f" {c[0]}╚██████╔╝{reset} {c[1]}╚████╔╝{reset} {c[2]}╚██████╗{reset}{c[3]}███████╗{reset}{c[4]}██║  ██║{reset}{c[5]}███████║{reset}{c[6]}███████║{reset}",
        f"  {c[0]}╚═════╝{reset}   {c[1]}╚═══╝{reset}   {c[2]}╚═════╝{reset}{c[3]}╚══════╝{reset}{c[4]}╚═╝  ╚═╝{reset}{c[5]}╚══════╝{reset}{c[6]}╚══════╝{reset}",
        f"         Giant Virus Classification Tool {software_version}",
    ]


def print_banner(plain_output: bool, software_version: str) -> None:
    print()
    if plain_output or not sys.stdout.isatty():
        for line in PLAIN_BANNER_LINES:
            print(line)
        return

    ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
    width = 60
    print(f"    ╔{'═' * width}╗")
    for line in _build_color_banner_lines(software_version):
        padding = max(0, width - len(ANSI_RE.sub("", line)))
        print(f"    ║{line}{' ' * padding}║")
    print(f"    ╚{'═' * width}╝")

def print_configuration(
    args,
    query_dir_abs: Path,
    output_dir: Path,
    database: Path,
    n_queries: int,
    n_input_files: int,
    n_skipped: int,
    mode_fast: bool,
    sensitive_mode: bool,
    threads: int,
    workers: "WorkerPlan",
    contigs: "ContigInput",
    input_min_length: int,
    contigs_min_length: int,
    validation_min_length: int,
) -> None:
    print()
    print("=" * 60)
    print("                    Pipeline Configuration")
    print("=" * 60)
    print(f"Config file: {args.config}")
    if contigs.enabled and contigs.input_file:
        print(f"Input file: {contigs.input_file}")
        print(f"Contigs mode: ENABLED ({n_queries} contigs from 1 file)")
        print(f"Contigs min length filter: {contigs_min_length} bp")
        print(f"Validation min length: {validation_min_length} bp")
    elif contigs.enabled and contigs.input_dir:
        print(f"Input directory: {contigs.input_dir}")
        print(f"Contigs mode: ENABLED ({n_queries} contigs from {n_input_files} files)")
        print(f"Contigs min length filter: {contigs_min_length} bp")
        print(f"Validation min length: {validation_min_length} bp")
    else:
        print(f"Query directory: {query_dir_abs}")
        print(f"Input min length: {input_min_length} bp")

    threads_used = workers.workers * workers.threads_per_worker
    print(f"Output directory: {output_dir}")
    print(f"Database: {database}")
    print(
        f"Threads: {threads_used} used / {threads} requested "
        f"(Workers: {workers.workers} × {workers.threads_per_worker} threads)"
    )
    print(f"Tree method: {args.tree_method if args.tree_method else 'veryfasttree'}")
    print(f"Fast mode: {mode_fast}")
    print(f"Sensitive mode: {sensitive_mode}")
    if args.resume:
        print("Resume mode: ENABLED")
        if n_skipped > 0:
            print(f"  - Skipping {n_skipped} completed queries")
            print(f"  - Processing {n_queries - n_skipped} remaining queries")
        else:
            print("  - No completed queries found")
    else:
        print("Resume mode: DISABLED")
    print(f"Total queries: {n_queries}")
    print("=" * 60)
