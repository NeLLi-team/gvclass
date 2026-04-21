"""
Extract marker-specific hits from HMM search results.

This module implements the functionality of the extract_qhits rule from the Snakemake workflow.
"""

from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, Iterable, Set, List, Tuple
from Bio import SeqIO

from src.utils import setup_logging

logger = setup_logging(__name__)


def parse_hmm_output(hmm_file: Path) -> Dict[str, Set[str]]:
    """
    Parse HMM output file and extract hits for each marker.

    Args:
        hmm_file: Path to filtered HMM output file

    Returns:
        Dictionary mapping marker names to sets of query sequence IDs
    """
    marker_hits = defaultdict(set)

    try:
        with open(hmm_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split()
                if len(fields) >= 4:
                    query_id = fields[0]
                    marker_name = fields[3]
                    marker_hits[marker_name].add(query_id)

        logger.info(f"Found hits for {len(marker_hits)} markers")
        return dict(marker_hits)

    except Exception as e:
        logger.error(f"Error parsing HMM output: {e}")
        raise


def extract_marker_sequences(
    query_faa: Path, marker_hits: Dict[str, Set[str]], output_dir: Path
) -> Dict[str, Path]:
    """
    Extract sequences for each marker and write to separate files.

    Args:
        query_faa: Path to query protein sequences
        marker_hits: Dictionary of marker -> sequence IDs
        output_dir: Directory to write marker-specific FASTA files

    Returns:
        Dictionary mapping marker names to output file paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all sequences into memory for efficient lookup.
    #
    # HMM target names (``fields[0]`` in parse_hmm_output) are Biopython's
    # ``record.id`` — the first whitespace-delimited token of the header.
    # The previous index was keyed by ``record.description`` (the full
    # header), which silently missed every lookup when headers carried
    # annotation after the ID (e.g. "prot_1 # start # end ..." from
    # pyrodigal). Index by both so either form resolves, preferring
    # ``record.id`` for the canonical match.
    logger.info(f"Loading sequences from {query_faa}")
    sequences_by_id: Dict[str, Any] = {}
    sequences_by_description: Dict[str, Any] = {}
    for record in SeqIO.parse(str(query_faa), "fasta"):
        sequences_by_id[record.id] = record
        sequences_by_description[record.description] = record

    # Extract sequences for each marker
    marker_files = {}

    for marker, hit_ids in marker_hits.items():
        output_file = output_dir / f"{marker}.faa"

        # Extract sequences that hit this marker
        marker_sequences = []
        for hit_id in hit_ids:
            record = sequences_by_id.get(hit_id) or sequences_by_description.get(hit_id)
            if record is not None:
                marker_sequences.append(record)
            else:
                logger.warning(f"Sequence {hit_id} not found for marker {marker}")

        # Write sequences to marker-specific file
        if marker_sequences:
            SeqIO.write(marker_sequences, str(output_file), "fasta")
            marker_files[marker] = output_file
            logger.info(f"Wrote {len(marker_sequences)} sequences for marker {marker}")
        else:
            logger.warning(f"No sequences found for marker {marker}")

    return marker_files


def extract_marker_hits(
    hmm_output: Path,
    query_faa: Path,
    output_dir: Path,
    min_hits: int = 1,
    mode_fast: bool = False,
) -> List[Tuple[str, Path]]:
    """
    Main function to extract marker-specific hits from HMM search results.

    Args:
        hmm_output: Path to filtered HMM output
        query_faa: Path to query protein sequences
        output_dir: Directory for marker-specific outputs
        min_hits: Minimum number of hits required to process a marker
        mode_fast: If True, skip OG markers (for tree building)

    Returns:
        List of tuples (marker_name, fasta_file_path) for markers with sufficient hits
    """
    logger.info(f"Starting marker extraction from {hmm_output}")
    logger.info(f"Query FAA file: {query_faa}, exists: {query_faa.exists()}")
    logger.info(
        f"Mode fast: {mode_fast} (OG markers will {'be skipped' if mode_fast else 'be included'})"
    )

    # Parse HMM output
    marker_hits = parse_hmm_output(hmm_output)
    logger.info(f"Parsed {len(marker_hits)} markers from HMM output")

    # Filter markers by minimum hits
    filtered_markers = {
        marker: hits for marker, hits in marker_hits.items() if len(hits) >= min_hits
    }

    # Skip OG markers in fast mode
    if mode_fast:
        og_markers_before = len(filtered_markers)
        filtered_markers = {
            marker: hits
            for marker, hits in filtered_markers.items()
            if not marker.startswith("OG")
        }
        og_markers_skipped = og_markers_before - len(filtered_markers)
        if og_markers_skipped > 0:
            logger.info(f"Skipped {og_markers_skipped} OG markers in fast mode")

    logger.info(f"Processing {len(filtered_markers)} markers with >= {min_hits} hits")

    # Extract sequences for each marker
    marker_files = extract_marker_sequences(query_faa, filtered_markers, output_dir)
    logger.info(f"Extracted sequences for {len(marker_files)} markers")

    # Return list of marker-file pairs
    return [(marker, path) for marker, path in marker_files.items()]


def get_marker_database_path(marker: str, database_path: Path) -> Path:
    """
    Get the database file path for a specific marker.

    Args:
        marker: Marker name
        database_path: Base database directory

    Returns:
        Path to marker-specific FASTA file
    """
    faa_path = database_path / "database" / "faa" / f"{marker}.faa"

    # Check if file exists
    if not faa_path.exists():
        # Try alternative location
        alt_path = database_path / "faa" / f"{marker}.faa"
        if alt_path.exists():
            return alt_path
        else:
            raise FileNotFoundError(
                f"Database file not found for marker {marker}: {faa_path}"
            )

    return faa_path


def count_marker_hits(hmm_output: Path) -> Dict[str, int]:
    """
    Count the number of hits for each marker.

    Args:
        hmm_output: Path to filtered HMM output

    Returns:
        Dictionary mapping marker names to hit counts
    """
    marker_hits = parse_hmm_output(hmm_output)
    return {marker: len(hits) for marker, hits in marker_hits.items()}


def count_unique_proteins_for_markers(
    marker_hits: Dict[str, Set[str]], markers: Iterable[str]
) -> int:
    """Count unique proteins matched by any marker in a category."""
    proteins: Set[str] = set()
    for marker in markers:
        proteins.update(marker_hits.get(marker, set()))
    return len(proteins)


def count_unique_proteins_for_prefixes(
    marker_hits: Dict[str, Set[str]], prefixes: Iterable[str]
) -> int:
    """Count unique proteins matched by markers sharing one of the prefixes."""
    prefix_tuple = tuple(prefixes)
    proteins: Set[str] = set()
    for marker, hits in marker_hits.items():
        if marker.startswith(prefix_tuple):
            proteins.update(hits)
    return len(proteins)


def count_unique_proteins_by_category_models(
    marker_hits: Dict[str, Set[str]], category_models: Dict[str, Iterable[str]]
) -> Dict[str, int]:
    """Count unique proteins for each category defined by explicit marker names."""
    return {
        category: count_unique_proteins_for_markers(marker_hits, models)
        for category, models in category_models.items()
    }


def count_unique_proteins_by_category_prefixes(
    marker_hits: Dict[str, Set[str]], category_prefixes: Dict[str, str]
) -> Dict[str, int]:
    """Count unique proteins for each category defined by a marker prefix."""
    return {
        category: count_unique_proteins_for_prefixes(marker_hits, [prefix])
        for category, prefix in category_prefixes.items()
    }
