"""Run the two contig screens on the bins of a folder: ``gvclass screen``.

Usage: gvclass screen <bin.fna | folder> [-o OUT] [-d DATABASE] [-t THREADS]
       gvclass screen --install [-d DATABASE]

For every bin the cellular screen (bacterial, archaeal and eukaryotic sequence) and the
giant-virus screen (contigs of another giant virus) write ``<bin>.cellular.tsv`` and
``<bin>.giant_virus.tsv``; their printed summaries go to ``<bin>.summary.txt``.
``contig_screen_summary.tsv`` has one row per bin with the status of each screen. The
resources are in ``<database>/contig_screen/``; ``--install`` downloads them (checked
against a SHA-256) next to the 2.1.0 database. A call is a review flag, not a
contamination estimate.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
from collections.abc import Callable, Sequence
from pathlib import Path

import pandas as pd

from src.contig_screen import cellular, giant_virus
from src.utils.database_manager import DatabaseManager

DATABASE_PATH_ENV = "GVCLASS_DB"
FASTA_SUFFIXES = (
    ".fna",
    ".fa",
    ".fasta",
    ".fas",
    ".fna.gz",
    ".fa.gz",
    ".fasta.gz",
    ".fas.gz",
)


def run_cellular(argv: list[str]) -> int:
    """Run the cellular screen with command-line arguments."""
    return cellular.main(argv)


def run_giant_virus(argv: list[str]) -> int:
    """Run the giant-virus screen with command-line arguments."""
    return giant_virus.main(argv)


def bundle_present(database: Path) -> bool:
    """Return True if the contig-screen resources are installed in ``database``."""
    return (database / "contig_screen" / "manifest.json").is_file()


def list_bins(source: Path) -> list[Path]:
    """Return the FASTA files of a folder, sorted, or the file itself."""
    if source.is_file():
        return [source]
    return sorted(
        path
        for path in source.iterdir()
        if path.is_file() and path.name.lower().endswith(FASTA_SUFFIXES)
    )


def bin_name(path: Path) -> str:
    """Return the file name of a bin without its FASTA (and .gz) suffix."""
    name = path.name
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def name_clashes(bins: list[Path]) -> list[list[str]]:
    """Return the groups of bins whose files get the same output name."""
    by_name: dict[str, list[str]] = {}
    for path in bins:
        by_name.setdefault(bin_name(path), []).append(path.name)
    return [names for names in by_name.values() if len(names) > 1]


def default_output_dir(source: Path) -> Path:
    """Return ``<input>_screen`` beside the input (also for ``.``)."""
    resolved = source.resolve()
    return resolved.with_name(resolved.name + "_screen")


def resolve_database(option: str | None, repo_dir: Path) -> Path:
    """Return the bundle root: the option, ``GVCLASS_DB``, else ``resources``."""
    setting = option or os.environ.get(DATABASE_PATH_ENV, "").strip() or "resources"
    path = Path(setting).expanduser()
    return (path if path.is_absolute() else repo_dir / path).resolve()


def run_screen(runner: Callable[[list[str]], int], argv: list[str]) -> tuple[int, str]:
    """Run one screen, return its exit code and what it printed."""
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stdout(buffer):
            code = runner(argv)
    except Exception as err:  # noqa: BLE001  # one bad bin must not stop the batch
        return 1, f"{buffer.getvalue()}error: {type(err).__name__}: {err}\n"
    return code, buffer.getvalue()


def screen_bin(
    path: Path, out_dir: Path, database: Path, threads: int
) -> dict[str, str]:
    """Screen one bin and return its row of the summary table."""
    name = bin_name(path)
    common = ["--resources", str(database), "--threads", str(threads)]
    cellular_table = out_dir / f"{name}.cellular.tsv"
    giant_table = out_dir / f"{name}.giant_virus.tsv"
    cellular_table.unlink(missing_ok=True)
    giant_table.unlink(missing_ok=True)
    code, text = run_screen(
        run_cellular, [str(path), "--out", str(cellular_table), *common]
    )
    cellular_ok = code == 0 and cellular_table.exists()
    summary = [text]
    virus_args = [str(path), "--out", str(giant_table), *common]
    if cellular_ok:
        virus_args += ["--cellular", str(cellular_table)]
    code_virus, text_virus = run_screen(run_giant_virus, virus_args)
    summary.append(text_virus)
    (out_dir / f"{name}.summary.txt").write_text("".join(summary), encoding="utf-8")
    return {
        "bin": name,
        "cellular_status": "ok" if cellular_ok else "failed",
        "giant_virus_status": "ok" if code_virus == 0 else "failed",
    }


def install_bundle(database: Path) -> int:
    """Download and install the contig-screen resources; return the exit code."""
    source = DatabaseManager.CONTIG_SCREEN_SOURCE
    if source is None:
        print("error: no download source is configured for this build", file=sys.stderr)
        return 1
    database.mkdir(parents=True, exist_ok=True)
    DatabaseManager.install_addon(database, "contig_screen", source)
    print(f"installed contig_screen in {database}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Return the command line parser."""
    parser = argparse.ArgumentParser(prog="gvclass screen", description=__doc__)
    parser.add_argument(
        "input", type=Path, nargs="?", help="a bin (FASTA) or a folder of bins"
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="download the contig-screen resources into the database",
    )
    parser.add_argument("-o", "--output-dir", type=Path, help="default: <input>_screen")
    parser.add_argument("-d", "--database", help="resource bundle (default GVCLASS_DB)")
    parser.add_argument("-t", "--threads", type=int, default=4)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Screen every bin; return 1 if no bin was found or a screen failed."""
    parser = build_parser()
    args = parser.parse_args(argv)
    repo_dir = Path(__file__).resolve().parents[2]
    database = resolve_database(args.database, repo_dir)
    if args.install:
        return install_bundle(database)
    if args.input is None:
        parser.error("an input bin or folder is needed (or --install)")
    if not bundle_present(database):
        print(
            f"error: no contig-screen resources in {database}; "
            "run: gvclass screen --install",
            file=sys.stderr,
        )
        return 1
    bins = list_bins(args.input)
    if not bins:
        print(f"error: no FASTA file in {args.input}", file=sys.stderr)
        return 1
    clashes = name_clashes(bins)
    if clashes:
        print(
            "error: these files would write to the same output name: "
            + "; ".join(", ".join(group) for group in clashes),
            file=sys.stderr,
        )
        return 1
    out_dir = args.output_dir or default_output_dir(args.input)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [screen_bin(path, out_dir, database, args.threads) for path in bins]
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "contig_screen_summary.tsv", sep="\t", index=False)
    failed = (table[["cellular_status", "giant_virus_status"]] == "failed").any().any()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
