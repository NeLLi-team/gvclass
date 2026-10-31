"""Run the two contig screens on the bins of a folder: ``gvclass screen``.

Usage: gvclass screen <bin.fna | folder> [-o OUT] [-d DATABASE] [-t THREADS]
       gvclass screen --install [-d DATABASE]

For every bin the cellular screen (bacterial, archaeal and eukaryotic sequence) and the
giant-virus screen (contigs of another giant virus) write ``<bin>.cellular.tsv`` and
``<bin>.giant_virus.tsv``; their printed summaries go to ``<bin>.summary.txt``.
``contig_screen_summary.tsv`` has one row per bin with the status of each screen, the
base pairs they call and the SHA-256 of the add-on manifest and of the bin file. The
add-on is in ``<database>/contig_screen/``; ``--install`` downloads it (checked
against a SHA-256) into the database, which must be writable. In the Apptainer
image the launcher mounts a host directory there. A call is a review flag, not a
contamination estimate.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import pandas as pd

from src.contig_screen import cellular, extract_genes, giant_virus, layout
from src.contig_screen.layout import bundle_present
from src.contig_screen.tables import read_table
from src.utils.common import sha256_file
from src.utils.database_manager import DatabaseManager

DATABASE_PATH_ENV = "GVCLASS_DB"
OK = "ok"
FAILED = "failed"
NOT_SCORED = "not_scored"
SUMMARY_FILE = "contig_screen_summary.tsv"
BP_COLUMNS = ["input_bp", "cellular_bp", "giant_virus_bp", "flagged_bp"]
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


@dataclass(frozen=True, slots=True)
class Addon:
    """What the summary needs of the add-on: its manifest digest and the cellular cut.

    The cut is read on demand, so that a broken ``thresholds.json`` fails the
    cellular screen of a bin and not the giant-virus screen.
    """

    manifest_sha256: str
    stage2: Callable[[], float]


def load_addon(database: Path) -> Addon:
    """Open the add-on after its manifest check; return its digest and cut reader."""
    resources = layout.ScreenResources.open(database)
    return Addon(
        resources.manifest_sha256,
        partial(_read_stage2, resources),
    )


def _read_stage2(resources: layout.ScreenResources) -> float:
    return layout.load_thresholds(resources)["stage2"]


def run_cellular_screen(
    path: Path, table: Path, common: list[str], lengths: dict[str, int]
) -> tuple[str, str]:
    """Run the cellular screen on one bin; return its status and printed text.

    A bin with no contig of 3 kb or more has nothing to score, so the screen does not
    run: that is not a failure.
    """
    longest = max(lengths.values(), default=0)
    if longest < extract_genes.MIN_SOURCE_CONTIG:
        text = (
            f"no contig of {extract_genes.MIN_SOURCE_CONTIG} bp or more: not scored\n"
        )
        return NOT_SCORED, text
    code, text = run_screen(run_cellular, [str(path), "--out", str(table), *common])
    return (OK if code == 0 and table.exists() else FAILED), text


def run_giant_virus_screen(
    path: Path, table: Path, common: list[str], cellular_table: Path | None
) -> tuple[str, str]:
    """Run the giant-virus screen on one bin; return its status and printed text.

    The screen exits cleanly without a table for a bin with fewer than two scored
    contigs: that is ``not_scored``.
    """
    argv = [str(path), "--out", str(table), *common]
    if cellular_table is not None:
        argv += ["--cellular", str(cellular_table)]
    code, text = run_screen(run_giant_virus, argv)
    if code != 0:
        return FAILED, text
    return (OK if table.exists() else NOT_SCORED), text


def percent(part: float, whole: int) -> float:
    """Return ``part`` as a percentage of ``whole``, NaN for an empty input."""
    return 100.0 * part / whole if whole else float("nan")


def bin_totals(
    lengths: dict[str, int],
    cellular_table: Path | None,
    giant_table: Path | None,
    stage2: float,
) -> dict[str, float]:
    """Return the base pairs and percentages the screens call; NaN where unscored.

    ``cellular_table`` and ``giant_table`` are the tables of the screens that
    scored the bin, else None. The union needs both.
    """
    nan = float("nan")
    input_bp = sum(lengths.values())
    called = None
    if giant_table is not None:
        called = giant_virus.called_contigs(read_table(giant_table, ["contig"]))
    union = None
    if cellular_table is not None:
        foreign = set(called["contig"]) if called is not None else set()
        union = giant_virus.union_totals(cellular_table, stage2, lengths, foreign)
    cellular_bp = nan if union is None else union.cellular_bp
    giant_bp = nan if called is None else float(called["length_bp"].sum())
    flagged_bp = nan if union is None or called is None else union.union_bp
    return {
        "input_bp": input_bp,
        "cellular_bp": cellular_bp,
        "cellular_pct": percent(cellular_bp, input_bp),
        "giant_virus_bp": giant_bp,
        "giant_virus_pct": percent(giant_bp, input_bp),
        "flagged_bp": flagged_bp,
        "flagged_pct": percent(flagged_bp, input_bp),
    }


def run_cellular_result(
    path: Path,
    table: Path,
    common: list[str],
    lengths: dict[str, int],
    addon: Addon,
) -> tuple[str, str, float]:
    """Run the cellular screen and read its cut; return status, text and the cut.

    Both belong to the cellular screen: a ``thresholds.json`` that cannot be read
    fails this screen only.
    """
    status, text = run_cellular_screen(path, table, common, lengths)
    if status != OK:
        return status, text, float("nan")
    try:
        return OK, text, addon.stage2()
    except (layout.IntegrityError, KeyError, ValueError) as err:
        return FAILED, f"{text}error: cellular thresholds: {err}\n", float("nan")


def empty_totals() -> dict[str, float]:
    """Return the numbers of a bin that no screen scored: all missing."""
    nan = float("nan")
    return dict.fromkeys(
        [
            "input_bp",
            "cellular_bp",
            "cellular_pct",
            "giant_virus_bp",
            "giant_virus_pct",
            "flagged_bp",
            "flagged_pct",
        ],
        nan,
    )


def screen_bin(
    path: Path, out_dir: Path, database: Path, threads: int, addon: Addon
) -> dict[str, object]:
    """Screen one bin and return its row of the summary table.

    A bin that cannot be read or totalled gets a row with both screens ``failed``
    and no numbers, so the other bins of the batch still run.
    """
    name = bin_name(path)
    try:
        return screen_bin_unguarded(path, name, out_dir, database, threads, addon)
    except Exception as err:  # noqa: BLE001  # one bad bin must not stop the batch
        message = f"error: {type(err).__name__}: {err}\n"
        return failed_row(path, out_dir, message, addon.manifest_sha256)


def failed_row(
    path: Path, out_dir: Path, message: str, manifest_sha256: str
) -> dict[str, object]:
    """Return the row of a bin that no screen scored: both ``failed``, no numbers.

    ``message`` goes to ``<bin>.summary.txt``.
    """
    name = bin_name(path)
    with contextlib.suppress(OSError):
        (out_dir / f"{name}.summary.txt").write_text(message, encoding="utf-8")
    return {
        "bin": name,
        "cellular_status": FAILED,
        "giant_virus_status": FAILED,
        **empty_totals(),
        "manifest_sha256": manifest_sha256,
        "input_sha256": input_digest(path),
    }


def input_digest(path: Path) -> str:
    """Return the SHA-256 of the bin file, empty if it cannot be read."""
    try:
        return sha256_file(path)
    except OSError:
        return ""


def screen_bin_unguarded(
    path: Path,
    name: str,
    out_dir: Path,
    database: Path,
    threads: int,
    addon: Addon,
) -> dict[str, object]:
    """Screen one bin; any error in reading or totalling it propagates."""
    common = ["--resources", str(database), "--threads", str(threads)]
    cellular_table = out_dir / f"{name}.cellular.tsv"
    giant_table = out_dir / f"{name}.giant_virus.tsv"
    cellular_table.unlink(missing_ok=True)
    giant_table.unlink(missing_ok=True)
    fingerprint = sha256_file(path)  # before any screen reads the bin
    lengths = cellular.contig_lengths(path)
    cellular_status, text, stage2 = run_cellular_result(
        path, cellular_table, common, lengths, addon
    )
    cellular_ok = cellular_status == OK
    giant_status, text_virus = run_giant_virus_screen(
        path, giant_table, common, cellular_table if cellular_ok else None
    )
    (out_dir / f"{name}.summary.txt").write_text(text + text_virus, encoding="utf-8")
    totals = bin_totals(
        lengths,
        cellular_table if cellular_ok else None,
        giant_table if giant_status == OK else None,
        stage2,
    )
    return {
        "bin": name,
        "cellular_status": cellular_status,
        "giant_virus_status": giant_status,
        **totals,
        "manifest_sha256": addon.manifest_sha256,
        "input_sha256": fingerprint,
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


def write_summary(rows: list[dict[str, object]], out_dir: Path) -> pd.DataFrame:
    """Write ``contig_screen_summary.tsv``; missing numbers read ``NaN``.

    The table goes to a file beside it and is renamed into place, so a reader sees
    the whole table or the one of the earlier run, never a torn one.
    """
    table = pd.DataFrame(rows).astype(dict.fromkeys(BP_COLUMNS, "Int64"))
    partial_file = out_dir / f"{SUMMARY_FILE}.part"
    try:
        table.to_csv(
            partial_file, sep="\t", index=False, na_rep="NaN", float_format="%.2f"
        )
        partial_file.replace(out_dir / SUMMARY_FILE)
    finally:
        partial_file.unlink(missing_ok=True)
    return table


def build_parser() -> argparse.ArgumentParser:
    """Return the command line parser."""
    parser = argparse.ArgumentParser(prog="gvclass screen", description=__doc__)
    parser.add_argument(
        "input", type=Path, nargs="?", help="a bin (FASTA) or a folder of bins"
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="download the add-on into the database",
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
    try:
        addon = load_addon(database)
    except (layout.IntegrityError, KeyError, ValueError) as err:
        # No screen runs with a rejected add-on; each bin says why it has no result.
        print(f"error: {err}", file=sys.stderr)
        rows = [failed_row(path, out_dir, f"error: {err}\n", "") for path in bins]
    else:
        rows = [screen_bin(p, out_dir, database, args.threads, addon) for p in bins]
    table = write_summary(rows, out_dir)
    failed = (table[["cellular_status", "giant_virus_status"]] == "failed").any().any()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
