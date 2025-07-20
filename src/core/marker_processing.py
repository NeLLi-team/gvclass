"""
Process individual markers through BLAST, alignment, and tree building.

This module handles the complete processing pipeline for a single marker.
"""

from pathlib import Path
from typing import Dict, Tuple
from Bio import SeqIO, AlignIO
import pyfamsa
import pytrimal
import veryfasttree

from src.core.blast import run_blastp, parse_blastp
from src.utils import setup_logging

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
        
        for dir_path in [self.blast_dir, self.merged_dir, self.align_dir, self.tree_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
            
    def get_database_path(self) -> Path:
        """Get the path to the marker-specific database."""
        faa_path = self.database_path / "database" / "faa" / f"{self.marker}.faa"
        
        if not faa_path.exists():
            # Try alternative location
            alt_path = self.database_path / "faa" / f"{self.marker}.faa"
            if alt_path.exists():
                return alt_path
                
        return faa_path
        
    def blast_and_merge(
        self, 
        query_hits_faa: Path,
        max_hits: int = 100,
        threads: int = 4
    ) -> Path:
        """
        Run BLAST against marker database and merge with query sequences.
        
        Args:
            query_hits_faa: Path to query sequences that hit this marker
            max_hits: Maximum number of BLAST hits to include
            threads: Number of threads for BLAST
            
        Returns:
            Path to merged FASTA file
        """
        logger.info(f"Processing BLAST for marker {self.marker}")
        
        # Output files
        blast_out = self.blast_dir / f"{self.marker}.m8"
        merged_out = self.merged_dir / f"{self.marker}.faa"
        
        # Get reference database
        ref_faa = self.get_database_path()
        
        if not ref_faa.exists():
            logger.warning(f"Reference database not found for {self.marker}: {ref_faa}")
            # If no reference database, just copy query sequences
            import shutil
            shutil.copy(query_hits_faa, merged_out)
            return merged_out
            
        # Run BLAST
        run_blastp(
            queryfaa=str(query_hits_faa),
            refdb=str(ref_faa),
            blastpout=str(blast_out),
            threads=threads
        )
        
        # Parse BLAST results and extract top hits
        blast_hits = parse_blastp(str(blast_out))
        
        # Get top hit IDs (up to max_hits)
        # parse_blastp returns a flat list of best hits
        top_hit_ids = set(blast_hits[:max_hits])
                
        logger.info(f"Found {len(top_hit_ids)} unique reference sequences for {self.marker}")
        
        # Merge query sequences with top reference hits
        with open(merged_out, 'w') as out_f:
            # Write query sequences
            for record in SeqIO.parse(str(query_hits_faa), "fasta"):
                SeqIO.write(record, out_f, "fasta")
                
            # Write reference sequences
            ref_count = 0
            for record in SeqIO.parse(str(ref_faa), "fasta"):
                if record.id in top_hit_ids or record.description in top_hit_ids:
                    SeqIO.write(record, out_f, "fasta")
                    ref_count += 1
                    
        logger.info(f"Merged {ref_count} reference sequences with query for {self.marker}")
        
        return merged_out
        
    def align_and_trim(
        self,
        merged_faa: Path,
        threads: int = 4
    ) -> Tuple[Path, Path]:
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
            logger.warning(f"Too few sequences ({len(sequences)}) for alignment of {self.marker}")
            # Just copy sequences as "alignment"
            SeqIO.write(sequences, str(align_out), "fasta")
            SeqIO.write(sequences, str(trim_out), "fasta")
            return align_out, trim_out
            
        # Convert to pyfamsa format
        famsa_seqs = []
        for seq in sequences:
            famsa_seq = pyfamsa.Sequence(seq.id.encode(), str(seq.seq).encode())
            famsa_seqs.append(famsa_seq)
            
        # Perform alignment
        aligner = pyfamsa.Aligner(threads=threads)
        msa = aligner.align(famsa_seqs)
        
        # Write alignment
        with open(align_out, 'w') as f:
            for seq in msa:
                f.write(f">{seq.id.decode()}\n{seq.sequence.decode()}\n")
                
        # Trim alignment
        try:
            trimmed_aln = pytrimal.Alignment.load(str(align_out))
            
            # Use automated trimming
            trimmer = pytrimal.AutomaticTrimmer(method='automated1')
            trimmed = trimmer.trim(trimmed_aln)
            
            # Write trimmed alignment
            trimmed.dump(str(trim_out))
            
        except Exception as e:
            logger.warning(f"Trimming failed for {self.marker}: {e}. Using untrimmed alignment.")
            import shutil
            shutil.copy(align_out, trim_out)
            
        return align_out, trim_out
        
    def build_tree(
        self,
        trimmed_alignment: Path,
        tree_method: str = "fasttree",
        threads: int = 4
    ) -> Path:
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
        
        # Check if alignment has enough sequences
        aln_records = list(AlignIO.read(str(trimmed_alignment), "fasta"))
        if len(aln_records) < 3:
            logger.warning(f"Too few sequences ({len(aln_records)}) for tree building of {self.marker}")
            # Create a simple tree
            with open(tree_out, 'w') as f:
                if len(aln_records) == 2:
                    f.write(f"({aln_records[0].id}:0.1,{aln_records[1].id}:0.1);\n")
                elif len(aln_records) == 1:
                    f.write(f"({aln_records[0].id}:0.0);\n")
                else:
                    f.write("();\n")
            return tree_out
            
        if tree_method == "fasttree":
            # Use VeryFastTree
            try:
                tree_result = veryfasttree.run(
                    alignment=str(trimmed_alignment),
                    quiet=True
                )
                
                with open(tree_out, 'w') as f:
                    f.write(tree_result)
                    
            except Exception as e:
                logger.error(f"VeryFastTree failed for {self.marker}: {e}")
                raise
                
        else:
            # Use IQ-TREE
            import subprocess
            cmd = [
                "iqtree",
                "-s", str(trimmed_alignment),
                "-m", "TEST",
                "-mset", "LG,WAG,JTT",
                "-nt", str(threads),
                "-pre", str(self.tree_dir / self.marker),
                "-quiet"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"IQ-TREE failed for {self.marker}: {result.stderr}")
                raise RuntimeError(f"IQ-TREE failed: {result.stderr}")
                
        return tree_out
        
    def process_marker(
        self,
        query_hits_faa: Path,
        max_blast_hits: int = 100,
        tree_method: str = "fasttree",
        threads: int = 4
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
            
            # Build tree
            tree_path = self.build_tree(trim_path, tree_method, threads)
            
            return {
                "marker": self.marker,
                "query_hits": query_hits_faa,
                "merged_faa": merged_faa,
                "alignment": align_path,
                "trimmed_alignment": trim_path,
                "tree": tree_path
            }
            
        except Exception as e:
            logger.error(f"Error processing marker {self.marker}: {e}")
            raise