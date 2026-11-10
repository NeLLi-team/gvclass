"""Display helpers for the GVClass CLI."""

from __future__ import annotations

import sys
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.bin.gvclass_cli import PipelineContext

PLAIN_BANNER_LINES = (
    "    ================================================================",
    "    GVClass - Giant Virus Classification Tool",
    "    ================================================================",
)


def _build_color_banner_lines(software_version: str) -> list[str]:
    reset = "\033[0m"
    c = [
        "\033[95m",
        "\033[95m",
        "\033[91m",
        "\033[31m",
        "\033[31m",
        "\033[38;5;208m",
        "\033[93m",
    ]
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


def print_configuration(context: PipelineContext) -> None:
    """Print the resolved settings and the worker plan of the run."""
    config = context.config
    contigs = context.contigs
    print()
    print("=" * 60)
    print("                    Pipeline Configuration")
    print("=" * 60)
    print(f"Config file: {context.config_file or 'built-in defaults'}")
    if contigs.enabled and contigs.input_file:
        print(f"Input file: {contigs.input_file}")
        print(f"Contigs mode: ENABLED ({context.n_queries} contigs from 1 file)")
        print(f"Contigs min length filter: {config.contigs_min_length} bp")
    elif contigs.enabled and contigs.input_dir:
        print(f"Input directory: {contigs.input_dir}")
        print(
            f"Contigs mode: ENABLED ({context.n_queries} contigs "
            f"from {context.n_input_files} files)"
        )
        print(f"Contigs min length filter: {config.contigs_min_length} bp")
    else:
        print(f"Query directory: {config.query_path}")
        print(f"Input min length: {config.min_sequence_length} bp")

    threads_used = context.workers * context.threads_per_worker
    print(f"Output directory: {config.output_dir}")
    print(f"Database: {config.database_path}")
    print(
        f"Threads: {threads_used} used / {config.threads} requested "
        f"(Workers: {context.workers} × {context.threads_per_worker} threads)"  # noqa: RUF001  # the banner shows the multiplication sign
    )
    print(f"Tree method: {config.tree_method}")
    print(f"Fast mode: {config.mode_fast}")
    print(f"Sensitive mode: {config.sensitive_mode}")
    if config.resume:
        print("Resume mode: ENABLED")
        if context.n_skipped > 0:
            print(f"  - Skipping {context.n_skipped} completed queries")
            print(f"  - Processing {context.n_remaining} remaining queries")
        else:
            print("  - No completed queries found")
    else:
        print("Resume mode: DISABLED")
    print(f"Total queries: {context.n_queries}")
    print("=" * 60)
