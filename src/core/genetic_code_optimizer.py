"""
Genetic code optimization module for GVClass.

This module implements the same logic as opgecall.py for selecting the optimal
genetic code based on HMM search results, coding density, and other metrics.

Uses thread-based parallelism for better Dask compatibility.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple
import shutil
from Bio import SeqIO
import pyrodigal  # Using tomasbruna's fork with codes 106, 129
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import pandas as pd

from src.core.gene_calling import run_hmmsearch, calc_stats

logger = logging.getLogger(__name__)


class GeneticCodeOptimizer:
    """Optimize genetic code selection based on multiple metrics using thread-based parallelism."""
    
    def __init__(self, database_path: Path, threads: int = None):
        """
        Initialize with database path and thread configuration.
        
        Args:
            database_path: Path to database
            threads: Total threads available (defaults to CPU count)
        """
        self.database_path = database_path
        self.hmm_file = str(database_path / "models" / "combined.hmm")
        self.threads = threads or multiprocessing.cpu_count()
        
    def run_gene_calling_for_code(
        self, 
        query_file: Path,
        code: int,
        temp_dir: Path,
        records: List
    ) -> Tuple[Path, Path, List[str], List[str]]:
        """
        Run gene calling for a specific genetic code.
        
        Returns:
            Tuple of (faa_file, fna_file, proteins_list, genes_list)
        """
        faa_file = temp_dir / f"code_{code}.faa"
        fna_file = temp_dir / f"code_{code}.fna"
        proteins = []
        genes = []
        
        # tomasbruna's fork supports codes 106 and 129 natively!
        if code == 0:
            # Metagenomic mode
            orf_finder = pyrodigal.GeneFinder(meta=True)
        else:
            # Single genome mode with specific genetic code
            orf_finder = pyrodigal.GeneFinder(meta=False)
            # Prepare training sequence
            training_seq = self._prepare_training_sequence(records)
            orf_finder.train(training_seq, force_nonsd=True, translation_table=code)
        
        # Process all records
        for record in records:
            if code == 0:
                genes_found = orf_finder.find_genes(str(record.seq))
            else:
                genes_found = orf_finder.find_genes(bytes(str(record.seq), 'utf-8'))
                
            for i, gene in enumerate(genes_found):
                protein_id = f"{record.id}_CDS_{i+1}"
                if code == 0:
                    protein_seq = gene.translate()
                else:
                    protein_seq = gene.translate(translation_table=code)
                proteins.append(f">{protein_id}\n{protein_seq}")
                genes.append(f">{protein_id}\n{str(record.seq)[gene.begin-1:gene.end]}")
        
        # Write output files
        with open(faa_file, 'w') as f:
            f.write('\n'.join(proteins))
        with open(fna_file, 'w') as f:
            f.write('\n'.join(genes))
        
        logger.info(f"Gene calling for code {code}: found {len(proteins)} proteins")
            
        return faa_file, fna_file, proteins, genes
    
    def _prepare_training_sequence(self, records: List) -> bytes:
        """Prepare concatenated training sequence for pyrodigal."""
        sequences = []
        for i, record in enumerate(records):
            if i > 0:
                # Add stop codons in all frames between sequences
                sequences.append("TTAATTAATTAA")
            sequences.append(str(record.seq))
        if len(sequences) > 1:
            sequences.append("TTAATTAATTAA")
        return bytes("".join(sequences), 'utf-8')
    
    def evaluate_genetic_code(
        self,
        query_file: Path,
        faa_file: Path,
        code: int
    ) -> Dict[str, float]:
        """
        Evaluate a genetic code by running HMM search and calculating metrics.
        
        Returns:
            Dictionary with all metrics
        """
        # Run HMM search with detailed metrics
        hit_count, avg_score, profile_hits, avg_best_hit_score, avg_coverage, complete_hits, avg_complete_score, avg_proteins_per_profile = run_hmmsearch(
            str(faa_file), 
            self.hmm_file, 
            completeCutoff=0.66
        )
        
        # Calculate coding density and other stats
        stats = calc_stats(str(query_file), str(faa_file))
        coding_density = stats[4] if len(stats) > 4 else 0.0
        
        return {
            'code': code,
            'hit_count': hit_count,
            'avg_score': avg_score,
            'profile_hits': profile_hits,
            'avg_best_hit_score': avg_best_hit_score,
            'avg_coverage': avg_coverage,
            'complete_hits': complete_hits,
            'avg_complete_score': avg_complete_score,
            'avg_proteins_per_profile': avg_proteins_per_profile,
            'coding_density': coding_density,
            'gene_count': stats[3] if len(stats) > 3 else 0,
            'faa_file': str(faa_file),
            # Store all stats for later use
            'contigs': stats[0] if len(stats) > 0 else 0,
            'LENbp': stats[1] if len(stats) > 1 else 0,
            'GCperc': stats[2] if len(stats) > 2 else 0.0
        }
    
    def select_best_code(
        self,
        code_results: Dict[int, Dict],
        genetic_codes: List[int]
    ) -> Tuple[int, Dict]:
        """
        Select the best genetic code using the same logic as opgecall.py.
        
        Returns:
            Tuple of (best_code, best_metrics)
        """
        # Get meta coding density if available
        meta_coding_density = 0.0
        if 0 in code_results:
            meta_coding_density = code_results[0]['coding_density']
        
        # Apply selection logic from opgecall.py
        if meta_coding_density > 0 and 0 in code_results:
            # Start with meta as baseline
            best_code = 0
            best_metrics = code_results[0]
            
            # Compare each code against current best
            for code in genetic_codes:
                if code == 0 or code not in code_results:
                    continue
                
                metrics = code_results[code]
                
                # Selection criteria (in order of priority):
                # 1. More complete hits
                # 2. Same complete hits but at least 5% higher average best hit score than meta
                # 3. Same complete hits but at least 5% better coding density than meta
                if (metrics['complete_hits'] > best_metrics['complete_hits'] or
                    (metrics['complete_hits'] == best_metrics['complete_hits'] and 
                     metrics['avg_best_hit_score'] > code_results[0]['avg_best_hit_score'] * 1.05) or
                    (metrics['complete_hits'] == best_metrics['complete_hits'] and 
                     metrics['coding_density'] > meta_coding_density * 1.05)):
                    
                    best_code = code
                    best_metrics = metrics
                    
        else:
            # No meta baseline - select by complete hits, then score, then coding density
            best_code = max(
                code_results.keys(),
                key=lambda c: (
                    code_results[c]['complete_hits'],
                    code_results[c]['avg_best_hit_score'],
                    code_results[c]['coding_density']
                )
            )
            best_metrics = code_results[best_code]
        
        return best_code, best_metrics
    
    def optimize_genetic_code(
        self,
        query_file: Path,
        output_dir: Path,
        genetic_codes: List[int] = [0, 1, 4, 6, 11, 15, 29, 106, 129]
    ) -> Tuple[int, Dict[str, Path], Dict, Dict[int, Dict]]:
        """
        Run genetic code optimization with thread-based parallelism.
        
        Returns:
            Tuple of (best_code, output_paths, metrics, all_code_results)
        """
        # Create temp directory
        temp_dir = output_dir / "gene_calling" / "temp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Read input sequences
        records = list(SeqIO.parse(str(query_file), "fasta"))
        logger.info(f"Processing {len(records)} sequences with {len(genetic_codes)} genetic codes")
        logger.info(f"Using threaded processing with {self.threads} total threads")
        
        # Calculate threads for parallel execution
        # Use fewer workers to leave threads for HMM search within each code test
        n_workers = min(len(genetic_codes), max(1, self.threads // 2))
        
        logger.info(f"Thread strategy: {n_workers} parallel genetic code tests")
        
        # Process codes in parallel using threads
        code_results = {}
        
        def test_single_code(code):
            """Test a single genetic code."""
            logger.info(f"Testing genetic code {code}")
            
            try:
                # Run gene calling
                faa_file, fna_file, proteins, genes = self.run_gene_calling_for_code(
                    query_file, code, temp_dir, records
                )
                
                # Skip if no genes found
                if not proteins:
                    logger.warning(f"No genes found with code {code}")
                    return None
                
                # Evaluate the code
                metrics = self.evaluate_genetic_code(query_file, faa_file, code)
                metrics['fna_file'] = str(fna_file)
                
                # Add total_bases from the records
                total_bases = sum(len(rec.seq) for rec in records)
                metrics['total_bases'] = total_bases
                
                logger.info(
                    f"Code {code}: complete_hits={metrics['complete_hits']}, "
                    f"avg_score={metrics['avg_best_hit_score']:.2f}, "
                    f"coding_density={metrics['coding_density']:.2f}%"
                )
                
                return (code, metrics)
                
            except Exception as e:
                logger.error(f"Failed to process code {code}: {e}")
                import traceback
                traceback.print_exc()
                return None
        
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            # Submit all tasks
            futures = {
                executor.submit(test_single_code, code): code
                for code in genetic_codes
            }
            
            # Collect results as they complete
            for future in as_completed(futures):
                result = future.result()
                if result:
                    code, metrics = result
                    code_results[code] = metrics
        
        if not code_results:
            raise RuntimeError("No genetic codes produced valid results")
        
        logger.info(f"Tested {len(code_results)} genetic codes successfully")
        
        # Select best code
        best_code, best_metrics = self.select_best_code(code_results, genetic_codes)
        
        logger.info(
            f"Selected genetic code {best_code} with complete_hits={best_metrics['complete_hits']}, "
            f"avg_score={best_metrics['avg_best_hit_score']:.2f}, "
            f"coding_density={best_metrics['coding_density']:.2f}%"
        )
        
        # Prepare output files
        prefix = query_file.stem
        gene_dir = output_dir / "gene_calling" / prefix
        gene_dir.mkdir(parents=True, exist_ok=True)
        
        outputs = {
            "fna": gene_dir / f"{prefix}.fna",
            "faa": gene_dir / f"{prefix}.faa",
            "gvogout": gene_dir / f"{prefix}.gvogout",
            "scoreout": gene_dir / f"{prefix}.scoreout",
            "statsout": gene_dir / f"{prefix}.statsout"
        }
        
        # Copy best results
        shutil.copy(code_results[best_code]['faa_file'], outputs["faa"])
        shutil.copy(code_results[best_code]['fna_file'], outputs["fna"])
        
        # Create stats file
        stats_df = pd.DataFrame([{
            'query': query_file.stem,
            'genetic_code': best_code,
            'complete_hits': best_metrics['complete_hits'],
            'avg_best_hit_score': best_metrics['avg_best_hit_score'],
            'coding_density': best_metrics['coding_density'],
            'gene_count': best_metrics['gene_count'],
            'ttable': 'codemeta',  # Always show codemeta as requested
            'contigs': best_metrics.get('contigs', 1),
            'LENbp': best_metrics.get('LENbp', 0),
            'GCperc': best_metrics.get('GCperc', 0.0),
            'genecount': best_metrics['gene_count'],
            'CODINGperc': best_metrics['coding_density']
        }])
        stats_df.to_csv(outputs["statsout"], sep='\t', index=False)
        
        # Log output paths
        logger.info(f"Final outputs: faa={outputs['faa']}, fna={outputs['fna']}")
        logger.info(f"Gene count in final faa: {best_metrics['gene_count']}")
        
        # Clean up temp directory
        shutil.rmtree(temp_dir)
        
        return best_code, outputs, best_metrics, code_results