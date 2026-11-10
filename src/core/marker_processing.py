"""
Process individual markers through BLAST, alignment, and tree building.

This module handles the complete processing pipeline for a single marker.
"""

from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple

import pytrimal
from Bio import SeqIO

from src.core.alignment import IQTREE_MODEL, align_sequences_pyfamsa, run_veryfasttree
from src.core.blast import read_fasta, search_proteins, select_references
from src.utils import setup_logging
from src.utils.resource_store import ResourceStore

logger = setup_logging(__name__)


class MarkerProcessor:
    """Handle processing of a single marker through the complete pipeline."""

    def __init__(self, marker: str, database_path: Path, output_dir: Path):
        """
        Initialize marker processor.

        Args:
            marker: Marker name
            database_path: Path to database directory
            output_dir: Base output directory
        """
        self.marker = marker
        self.database_path = database_path
        self.output_dir = output_dir

        # Create marker-specific output directories
        self.blast_dir = output_dir / "blastp_out"
        self.merged_dir = output_dir / "query_hits_merged_faa"
        self.align_dir = output_dir / "queryrefs_aligned"
        self.tree_dir = output_dir / "queryrefs_genetrees"

        for dir_path in [
            self.blast_dir,
            self.merged_dir,
            self.align_dir,
            self.tree_dir,
        ]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def get_database_path(self) -> Path:
        """Return the marker's reference FASTA.

        Raises:
            FileNotFoundError: If the marker has no reference set.
        """
        return ResourceStore(self.database_path).marker_faa_path(self.marker)

    def blast_and_merge(
        self, query_hits_faa: Path, max_hits: int = 100, threads: int = 4
    ) -> Path:
        """Search the marker's reference set and merge the chosen references.

        Every query protein gets its own best hits among the references
        (see :func:`src.core.blast.select_references`). The merged FASTA holds
        the query sequences first, then the references in selection order. The
        reference file is read once. The search is written to ``blastp_out``
        as ``<marker>.m8``, which contamination scoring reads.

        Args:
            query_hits_faa: Path to query sequences that hit this marker
            max_hits: Maximum number of hits per query and of references kept
            threads: Number of threads for the search

        Returns:
            Path to merged FASTA file

        Raises:
            FileNotFoundError: If the marker has no reference set.
        """
        logger.info("Processing BLAST for marker %s", self.marker)
        merged_out = self.merged_dir / f"{self.marker}.faa"

        queries = read_fasta(query_hits_faa)
        references = read_fasta(self.get_database_path())
        hits_by_query = search_proteins(
            queries,
            references,
            self.blast_dir / f"{self.marker}.m8",
            threads=threads,
            top_per_query=max_hits,
        )
        selected_ids = select_references(hits_by_query, max_hits)
        logger.info(
            "Found %d unique reference sequences for %s", len(selected_ids), self.marker
        )

        references_by_id = {record.id: record for record in references}
        with merged_out.open("w", encoding="utf-8") as out_f:
            for record in queries:
                out_f.write(f">{record.title}\n{record.sequence}\n")
            for subject_id in selected_ids:
                record = references_by_id[subject_id]
                out_f.write(f">{record.title}\n{record.sequence}\n")

        logger.info(
            "Merged %d reference sequences with query for %s",
            len(selected_ids),
            self.marker,
        )
        return merged_out

    def align_and_trim(self, merged_faa: Path, threads: int = 4) -> Tuple[Path, Path]:
        """
        Align sequences and trim alignment.

        Args:
            merged_faa: Path to merged FASTA file
            threads: Number of threads for alignment

        Returns:
            Tuple of (alignment_path, trimmed_alignment_path)
        """
        logger.info(f"Aligning sequences for marker {self.marker}")

        # Output files
        align_out = self.align_dir / f"{self.marker}.mafft"
        trim_out = self.align_dir / f"{self.marker}.mafft01"

        # Read sequences
        sequences = list(SeqIO.parse(str(merged_faa), "fasta"))

        if len(sequences) < 3:
            logger.warning(
                f"Too few sequences ({len(sequences)}) for alignment of {self.marker}"
            )
            # Just copy sequences as "alignment"
            SeqIO.write(sequences, str(align_out), "fasta")
            SeqIO.write(sequences, str(trim_out), "fasta")
            return align_out, trim_out

        # Align with FAMSA via the shared helper (identical params to the
        # historical inline call) and write the MSA in the same format as before.
        aligned = align_sequences_pyfamsa(sequences, threads)

        with open(align_out, "w") as f:
            for seq_id, seq_str in aligned:
                f.write(f">{seq_id}\n{seq_str}\n")

        # Trim alignment
        cleaned_align: Optional[Path] = None
        try:
            # First clean the alignment by replacing stop codons with gaps
            # pytrimal doesn't handle stop codons (*) well
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            ) as tmp:
                with open(align_out) as inf:
                    for line in inf:
                        if not line.startswith(">"):
                            # Replace stop codons with gaps
                            line = line.replace("*", "-")
                        tmp.write(line)
                cleaned_align = Path(tmp.name)

            trimmed_aln = pytrimal.Alignment.load(str(cleaned_align))

            # Use automated trimming
            trimmer = pytrimal.AutomaticTrimmer(method="automated1")
            trimmed = trimmer.trim(trimmed_aln)

            # Write trimmed alignment
            trimmed.dump(str(trim_out))

        except Exception as e:
            logger.warning(
                f"Trimming failed for {self.marker}: {e}. Using untrimmed alignment."
            )
            import shutil

            shutil.copy(align_out, trim_out)
        finally:
            if cleaned_align is not None:
                cleaned_align.unlink(missing_ok=True)

        return align_out, trim_out

    def build_tree(
        self,
        trimmed_alignment: Path,
        tree_method: str = "veryfasttree",
        threads: int = 4,
    ) -> Optional[Path]:
        """
        Build phylogenetic tree from alignment.

        Args:
            trimmed_alignment: Path to trimmed alignment
            tree_method: Tree building method (fasttree or iqtree)
            threads: Number of threads

        Returns:
            Path to tree file
        """
        logger.info(f"Building tree for marker {self.marker} using {tree_method}")

        tree_out = self.tree_dir / f"{self.marker}.treefile"
        for stale_output in (
            tree_out,
            tree_out.with_suffix(".SKIPPED"),
            tree_out.with_suffix(".ERROR"),
            tree_out.with_suffix(".TIMEOUT"),
        ):
            stale_output.unlink(missing_ok=True)

        # Short marker sets are intentionally copied without alignment, so the
        # equal-width check applies only to sets that go on to a tree.
        records = list(SeqIO.parse(trimmed_alignment, "fasta"))
        sequence_count = len(records)

        # Critical validation for tree building
        if sequence_count < 3:
            logger.warning(
                "Too few sequences (%d) for tree building of %s; skipping",
                sequence_count,
                self.marker,
            )

            # Create a SKIP file to indicate tree was not built
            skip_file = tree_out.with_suffix(".SKIPPED")
            with open(skip_file, "w") as f:
                f.write(f"TREE BUILDING SKIPPED AT {datetime.now()}\n")
                f.write(f"Marker: {self.marker}\n")
                f.write(f"Sequences: {sequence_count}\n")
                f.write("Minimum required: 3\n")
                f.write(
                    "Phylogenetic analysis not possible with insufficient sequences.\n"
                )

            # Return None to indicate no tree was built
            return None

        if len({len(record.seq) for record in records}) != 1:
            raise ValueError(
                f"Sequences in {trimmed_alignment} must all be the same length "
                f"to build the tree of {self.marker}"
            )

        if tree_method in ("veryfasttree", "fasttree"):
            # Use VeryFastTree
            try:
                # Forward `threads` so the per-marker worker budget is
                # honored; the previous call was single-threaded regardless
                # of how many threads the pipeline allocated per marker.
                tree_result = run_veryfasttree(trimmed_alignment, threads)

                with open(tree_out, "w") as f:
                    f.write(tree_result)

                # Validate the tree was created successfully
                if not tree_out.exists() or tree_out.stat().st_size == 0:
                    raise RuntimeError(
                        f"Tree file empty or not created for {self.marker}"
                    )

                logger.info(
                    "Successfully built tree for %s with %d sequences",
                    self.marker,
                    sequence_count,
                )

            except Exception as e:
                error_msg = f"ERROR: VeryFastTree failed for {self.marker}: {e}"
                logger.error(error_msg)

                # Create an error file to document the failure
                error_file = tree_out.with_suffix(".ERROR")
                with open(error_file, "w") as f:
                    f.write(f"TREE BUILDING FAILED AT {datetime.now()}\n")
                    f.write(f"Marker: {self.marker}\n")
                    f.write(f"Method: {tree_method}\n")
                    f.write(f"Sequences: {sequence_count}\n")
                    f.write(f"Error: {str(e)}\n")

                raise

        else:
            # Use IQ-TREE
            import subprocess

            cmd = [
                "iqtree",
                "-fast",  # Fast mode must come early to set optimization parameters
                "-s",
                str(trimmed_alignment),
                "-m",
                IQTREE_MODEL,
                "-nt",
                str(threads),
                "-pre",
                str(self.tree_dir / self.marker),
                "-quiet",
                "-redo",  # overwrite the outputs of an earlier run
            ]

            logger.info(f"Running IQ-TREE with command: {' '.join(cmd)}")

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=600,  # 10 minute timeout
                )
            except subprocess.TimeoutExpired:
                error_msg = (
                    f"ERROR: IQ-TREE timed out after 10 minutes for {self.marker}"
                )
                logger.error(error_msg)

                # Create an error file to document the timeout
                error_file = tree_out.with_suffix(".TIMEOUT")
                with open(error_file, "w") as f:
                    f.write(f"TREE BUILDING TIMEOUT AT {datetime.now()}\n")
                    f.write(f"Marker: {self.marker}\n")
                    f.write(f"Method: {tree_method}\n")
                    f.write(f"Sequences: {sequence_count}\n")
                    f.write("Timeout: 600 seconds\n")
                    f.write(f"Command: {' '.join(cmd)}\n")

                raise RuntimeError(f"IQ-TREE timed out for {self.marker}")
            except Exception as e:
                error_msg = f"ERROR: Failed to run IQ-TREE for {self.marker}: {e}"
                logger.error(error_msg)
                raise RuntimeError(f"Failed to run IQ-TREE: {e}")

            if result.returncode != 0:
                error_msg = f"ERROR: IQ-TREE failed for {self.marker}: {result.stderr}"
                logger.error(error_msg)

                # Create an error file to document the failure
                error_file = tree_out.with_suffix(".ERROR")
                with open(error_file, "w") as f:
                    f.write(f"TREE BUILDING FAILED AT {datetime.now()}\n")
                    f.write(f"Marker: {self.marker}\n")
                    f.write(f"Method: {tree_method}\n")
                    f.write(f"Sequences: {sequence_count}\n")
                    f.write(f"Error: {result.stderr}\n")
                    f.write(f"Command: {' '.join(cmd)}\n")

                raise RuntimeError(f"IQ-TREE failed: {result.stderr}")

            # Validate the tree was created successfully
            if not tree_out.exists() or tree_out.stat().st_size == 0:
                raise RuntimeError(f"Tree file empty or not created for {self.marker}")

            logger.info(
                "Successfully built tree for %s with %d sequences using IQ-TREE",
                self.marker,
                sequence_count,
            )

        return tree_out

    def process_marker(
        self,
        query_hits_faa: Path,
        max_blast_hits: int = 100,
        tree_method: str = "veryfasttree",
        threads: int = 4,
    ) -> Dict[str, Path]:
        """
        Process a single marker through the complete pipeline.

        Args:
            query_hits_faa: Path to query sequences that hit this marker
            max_blast_hits: Maximum number of BLAST hits
            tree_method: Tree building method
            threads: Number of threads

        Returns:
            Dictionary with paths to all output files
        """
        logger.info(f"Processing marker {self.marker}")

        try:
            # BLAST and merge
            merged_faa = self.blast_and_merge(query_hits_faa, max_blast_hits, threads)

            # Align and trim
            align_path, trim_path = self.align_and_trim(merged_faa, threads)

            # Build tree (may return None if insufficient sequences)
            tree_path = self.build_tree(trim_path, tree_method, threads)

            result = {
                "marker": self.marker,
                "query_hits": query_hits_faa,
                "merged_faa": merged_faa,
                "alignment": align_path,
                "trimmed_alignment": trim_path,
            }

            # Only include tree if it was actually built
            if tree_path is not None:
                result["tree"] = tree_path
            else:
                result["tree_skipped"] = True
                logger.info(
                    f"Tree building was skipped for {self.marker} due to insufficient sequences"
                )

            return result

        except Exception as e:
            logger.error(f"Error processing marker {self.marker}: {e}")
            raise
