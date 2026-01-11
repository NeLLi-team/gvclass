"""Split multi-contig FNA files into individual contig files."""

import tempfile
from pathlib import Path

from Bio import SeqIO


def split_contigs(
    input_fna: Path,
    output_dir: Path | None = None,
    min_length: int = 0,
    prefix: str | None = None,
) -> tuple[Path, int]:
    """Split a multi-contig FNA file into individual contig files.

    Args:
        input_fna: Path to input FNA file with multiple contigs.
        output_dir: Directory for contig files. If None, creates temp dir.
        min_length: Minimum contig length to include (default: 0, no filter).
        prefix: Optional prefix for output filenames (e.g., source filename).

    Returns:
        Tuple of (output_directory_path, number_of_contigs_written).

    Raises:
        FileNotFoundError: If input file doesn't exist.
        ValueError: If no contigs found or all filtered out.
    """
    input_fna = Path(input_fna)
    if not input_fna.exists():
        raise FileNotFoundError(f"Input file not found: {input_fna}")

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="gvclass_contigs_"))
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    contig_count = 0
    for record in SeqIO.parse(input_fna, "fasta"):
        if len(record.seq) < min_length:
            continue

        # Sanitize contig ID for filename (replace problematic characters)
        safe_id = sanitize_contig_id(record.id)

        # Add prefix if provided (for multi-file uniqueness)
        if prefix:
            filename = f"{prefix}__{safe_id}.fna"
        else:
            filename = f"{safe_id}.fna"

        contig_file = output_dir / filename

        with open(contig_file, "w") as f:
            SeqIO.write(record, f, "fasta")

        contig_count += 1

    if contig_count == 0:
        raise ValueError(
            f"No contigs found in {input_fna} "
            f"(min_length filter: {min_length} bp)"
        )

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

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="gvclass_contigs_"))
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Find all FNA files
    fna_extensions = {".fna", ".fa", ".fasta"}
    fna_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in fna_extensions
    )

    if not fna_files:
        raise ValueError(f"No FNA files found in {input_dir}")

    total_contigs = 0
    files_processed = 0

    for fna_file in fna_files:
        # Use filename (without extension) as prefix for uniqueness
        prefix = sanitize_contig_id(fna_file.stem)

        try:
            _, n_contigs = split_contigs(
                fna_file,
                output_dir=output_dir,
                min_length=min_length,
                prefix=prefix,
            )
            total_contigs += n_contigs
            files_processed += 1
        except ValueError:
            # Skip files with no valid contigs
            continue

    if total_contigs == 0:
        raise ValueError(
            f"No contigs extracted from {files_processed} files in {input_dir} "
            f"(min_length filter: {min_length} bp)"
        )

    return output_dir, total_contigs, files_processed


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
