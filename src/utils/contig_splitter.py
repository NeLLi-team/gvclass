"""Split multi-contig FNA files into individual contig files."""

import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO


class ContigIdCollisionError(ValueError):
    """Raised when two or more contig IDs sanitize to the same filename."""


def _plan_output_filename(record, prefix: str | None) -> str:
    safe_id = sanitize_contig_id(record.id)
    if prefix:
        return f"{prefix}__{safe_id}.fna"
    return f"{safe_id}.fna"


def _detect_collisions(records, min_length: int, prefix: str | None) -> List:
    """Two-pass collision detection.

    Iterate all records once to build a ``filename -> [contig_id]`` map. If
    any bucket has more than one contig, raise before a single file is
    written — the split stays transactional and no partial output is left
    on disk on failure (Codex-audit requirement for Phase 3.1).
    """
    filename_to_ids: Dict[str, List[str]] = defaultdict(list)
    retained = []
    for record in records:
        if len(record.seq) < min_length:
            continue
        retained.append(record)
        filename_to_ids[_plan_output_filename(record, prefix)].append(record.id)

    collisions = {
        filename: ids for filename, ids in filename_to_ids.items() if len(ids) > 1
    }
    if collisions:
        pretty = "; ".join(
            f"{filename} <- {ids}" for filename, ids in sorted(collisions.items())
        )
        raise ContigIdCollisionError(
            f"Sanitized contig filenames collide: {pretty}. "
            "Rename the conflicting contigs in the source FASTA or use "
            "--contigs with pre-disambiguated IDs."
        )
    return retained


def split_contigs(
    input_fna: Path,
    output_dir: Path | None = None,
    min_length: int = 0,
    prefix: str | None = None,
) -> tuple[Path, int]:
    """Split a multi-contig FNA file into individual contig files.

    The split is transactional: a first pass validates that no two contig
    IDs sanitize to the same filename, and only then do writes begin. A
    collision therefore never leaves a partially-populated output dir
    behind.

    Args:
        input_fna: Path to input FNA file with multiple contigs.
        output_dir: Directory for contig files. If None, creates temp dir.
        min_length: Minimum contig length to include (default: 0, no filter).
        prefix: Optional prefix for output filenames (e.g., source filename).

    Returns:
        Tuple of (output_directory_path, number_of_contigs_written).

    Raises:
        FileNotFoundError: If input file doesn't exist.
        ContigIdCollisionError: If two contig IDs sanitize to the same filename.
        ValueError: If no contigs found or all filtered out.
    """
    input_fna = Path(input_fna)
    if not input_fna.exists():
        raise FileNotFoundError(f"Input file not found: {input_fna}")

    # Pass 1: collect retained records and validate there are no collisions.
    # SeqIO.parse is a generator, so we materialize it once here.
    all_records = list(SeqIO.parse(input_fna, "fasta"))
    retained = _detect_collisions(all_records, min_length=min_length, prefix=prefix)

    if not retained:
        raise ValueError(
            f"No contigs found in {input_fna} "
            f"(min_length filter: {min_length} bp)"
        )

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="gvclass_contigs_"))
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Pass 2: actually write the files. At this point collisions are
    # impossible by construction.
    contig_count = 0
    for record in retained:
        filename = _plan_output_filename(record, prefix)
        contig_file = output_dir / filename
        with open(contig_file, "w") as f:
            SeqIO.write(record, f, "fasta")
        contig_count += 1

    return output_dir, contig_count


def split_contigs_from_directory(
    input_dir: Path,
    output_dir: Path | None = None,
    min_length: int = 0,
) -> tuple[Path, int, int]:
    """Split all FNA files in a directory into individual contig files.

    Each contig filename is prefixed with its source filename to ensure
    uniqueness when processing multiple input files.

    Args:
        input_dir: Directory containing FNA files.
        output_dir: Directory for contig files. If None, creates temp dir.
        min_length: Minimum contig length to include (default: 0, no filter).

    Returns:
        Tuple of (output_directory_path, total_contigs, files_processed).

    Raises:
        FileNotFoundError: If input directory doesn't exist.
        ValueError: If no FNA files found or no contigs extracted.
    """
    input_dir = Path(input_dir)
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise ValueError(f"Input path is not a directory: {input_dir}")

    # Find all FNA files first so we can plan the full split before touching
    # the output directory. The Codex-audit requirement is that a collision
    # anywhere in the directory leaves no partial files on disk.
    fna_extensions = {".fna", ".fa", ".fasta"}
    fna_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in fna_extensions
    )
    if not fna_files:
        raise ValueError(f"No FNA files found in {input_dir}")

    # Plan phase: collect every record we intend to write and check for
    # cross-file filename collisions before any output is materialized.
    planned: list[tuple[str, object]] = []
    filename_to_sources: Dict[str, List[str]] = defaultdict(list)
    files_with_data = 0
    for fna_file in fna_files:
        prefix = sanitize_contig_id(fna_file.stem)
        file_records = list(SeqIO.parse(fna_file, "fasta"))
        retained = _detect_collisions(
            file_records, min_length=min_length, prefix=prefix
        )
        if not retained:
            continue
        files_with_data += 1
        for record in retained:
            filename = _plan_output_filename(record, prefix)
            filename_to_sources[filename].append(f"{fna_file.name}:{record.id}")
            planned.append((filename, record))

    cross_file_collisions = {
        filename: sources
        for filename, sources in filename_to_sources.items()
        if len(sources) > 1
    }
    if cross_file_collisions:
        pretty = "; ".join(
            f"{filename} <- {sources}"
            for filename, sources in sorted(cross_file_collisions.items())
        )
        raise ContigIdCollisionError(
            f"Sanitized contig filenames collide across input files: {pretty}. "
            "Rename the conflicting contigs before running --contigs on this "
            "directory."
        )

    if not planned:
        raise ValueError(
            f"No contigs extracted from {len(fna_files)} files in {input_dir} "
            f"(min_length filter: {min_length} bp)"
        )

    # Write phase: only now do we create the output directory and write
    # records. An exception during writes still leaves the caller with a
    # well-defined (albeit partially-populated) directory, but by
    # construction we've validated there are no collisions to worry about.
    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="gvclass_contigs_"))
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    for filename, record in planned:
        contig_file = output_dir / filename
        with open(contig_file, "w") as f:
            SeqIO.write(record, f, "fasta")

    return output_dir, len(planned), files_with_data


def sanitize_contig_id(contig_id: str) -> str:
    """Sanitize contig ID for use as filename.

    Replaces characters that are problematic in filenames with underscores.

    Args:
        contig_id: Original contig ID from FASTA header.

    Returns:
        Sanitized string safe for use as filename.
    """
    # Characters to replace (filesystem-unfriendly)
    replacements = {
        "/": "_",
        "\\": "_",
        ":": "_",
        "*": "_",
        "?": "_",
        '"': "_",
        "<": "_",
        ">": "_",
        "|": "_",
        " ": "_",
    }

    result = contig_id
    for char, replacement in replacements.items():
        result = result.replace(char, replacement)

    # Collapse multiple underscores
    while "__" in result:
        result = result.replace("__", "_")

    # Remove leading/trailing underscores
    result = result.strip("_")

    # Ensure non-empty
    if not result:
        result = "contig"

    return result
