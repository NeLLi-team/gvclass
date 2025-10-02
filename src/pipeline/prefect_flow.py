"""
GVClass pipeline properly implemented with Prefect 2 and Dask.

This module provides a true Prefect-based workflow with:
- Proper @flow and @task decorators
- Dynamic DaskTaskRunner configuration
- Parallel query processing
- Automatic worker/thread distribution
- Support for local and cluster execution
"""

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import logging
from datetime import timedelta
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

from prefect import flow, task, get_run_logger
from prefect_dask import DaskTaskRunner

# Import core functionality
from src.core.reformat import reformat_sequences, calculate_stats
from src.core.tree_analysis import TreeAnalyzer
from src.utils import InputValidator
from src.utils.database_manager import DatabaseManager
from src.core.marker_extraction import extract_marker_hits, get_marker_database_path
from src.core.marker_processing import MarkerProcessor
from src.core.hmm_search import run_pyhmmer_search_with_filtering
from src.core.summarize_full import FullSummarizer
from src.core.genetic_code_optimizer import GeneticCodeOptimizer
from src.utils.error_handling import ProcessingError

from Bio import SeqIO


# Configure task defaults
task_config = {
    "retries": 2,
    "retry_delay_seconds": 30,
    "persist_result": True,
    "cache_expiration": timedelta(hours=24),
}


@task(name="validate_and_setup", **task_config)
def validate_and_setup_task(
    query_dir: str, output_dir: str, database_path: Optional[str] = None
) -> Dict[str, Any]:
    """Validate inputs and setup directories."""
    logger = get_run_logger()
    logger.info(f"Validating inputs - Query: {query_dir}, Output: {output_dir}")

    # Validate directories
    query_path = InputValidator.validate_query_directory(query_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Setup database
    db_path = DatabaseManager.setup_database(database_path)
    logger.info(f"Using database at: {db_path}")

    # Get query files
    query_files = []
    for ext in [".fna", ".faa"]:
        query_files.extend(query_path.glob(f"*{ext}"))

    logger.info(f"Found {len(query_files)} query files")

    return {
        "query_path": query_path,
        "output_path": output_path,
        "database_path": db_path,
        "query_files": query_files,
    }


@task(name="process_query", **task_config)
def process_query_task(
    query_file: Path,
    output_base: Path,
    database_path: Path,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    threads: int = 4,
) -> Dict[str, Any]:
    """
    Process a single query through the entire pipeline.
    This is the main unit of parallelization.
    """
    logger = get_run_logger()
    query_name = query_file.stem
    logger.info(f"Processing query: {query_name} with {threads} threads")

    # Create output directory for this query
    query_output_dir = output_base / query_name

    # Robust directory creation to handle race conditions
    # Try multiple times with small delays to avoid permission conflicts
    import time

    max_attempts = 5
    for attempt in range(max_attempts):
        try:
            query_output_dir.mkdir(parents=True, exist_ok=True)
            break
        except PermissionError as e:
            if attempt < max_attempts - 1:
                time.sleep(0.1 * (attempt + 1))  # Exponential backoff
                continue
            else:
                raise e

    # Step 1: Reformat input
    if query_file.suffix == ".fna":
        # Nucleotide sequences - need gene calling
        logger.info(f"Reformatting nucleotide sequences: {query_name}")

        # Reformat
        reformatted_file = query_output_dir / f"{query_name}_reformatted.fna"
        stats_file = query_output_dir / "stats" / f"{query_name}_stats.tsv"
        stats_file.parent.mkdir(exist_ok=True)

        records = list(SeqIO.parse(str(query_file), "fasta"))
        reformatted_records = reformat_sequences(records, str(query_file))
        SeqIO.write(reformatted_records, str(reformatted_file), "fasta")

        stats_df = calculate_stats(records, str(query_file))
        stats_df.to_csv(str(stats_file), sep="\t", index=False)

        # Step 2: Gene calling with genetic code optimization
        logger.info(f"Optimizing genetic code for: {query_name}")
        # Use parallel optimizer with allocated threads
        optimizer = GeneticCodeOptimizer(database_path=database_path, threads=threads)

        # Optimize genetic code in parallel
        best_code, outputs, best_metrics, all_code_results = (
            optimizer.optimize_genetic_code(
                query_file=reformatted_file,
                output_dir=query_output_dir,
                genetic_codes=genetic_codes,
            )
        )

        # Write genetic code comparison table
        stats_dir = query_output_dir / "stats"
        comparison_file = stats_dir / f"{query_name}.genecalling.tab"
        with open(comparison_file, "w") as f:
            f.write(
                "ttable\tGenes\tBases\tCodingDensity\tcomplete_bestHits\tavg_bestHit_score\n"
            )
            # Convert all_code_results dict to list format
            results_list = []
            for code, metrics in all_code_results.items():
                results_list.append(
                    {
                        "code": code,
                        "genes": metrics.get("gene_count", 0),
                        "bases": metrics.get("total_bases", 0),
                        "coding_density": metrics.get("coding_density", 0.0),
                        "complete_hits": metrics.get("complete_hits", 0),
                        "avg_best_hit_score": metrics.get("avg_score", 0.0),
                    }
                )

            sorted_results = sorted(
                results_list,
                key=lambda x: (-x["complete_hits"], -x["avg_best_hit_score"]),
            )
            for result in sorted_results:
                ttable = str(result["code"])
                f.write(
                    f"{ttable}\t{result['genes']}\t{result['bases']}\t"
                    f"{result['coding_density']:.3f}\t{result['complete_hits']}\t"
                    f"{result['avg_best_hit_score']:.1f}\n"
                )

        # Write stats file
        stats_file = stats_dir / f"{query_name}.stats.tab"
        ttable_str = "codemeta" if best_code == 0 else str(best_code)
        with open(stats_file, "w") as f:
            f.write(f"ttable\t{ttable_str}\n")
            if best_code in all_code_results:
                selected = all_code_results[best_code]
                f.write(f"genes\t{selected.get('gene_count', 0)}\n")
                f.write(f"bases\t{selected.get('total_bases', 0)}\n")
                f.write(f"coding_density\t{selected.get('coding_density', 0.0):.3f}\n")

        protein_file = Path(outputs["faa"])

    else:
        # Protein sequences - use directly
        logger.info(f"Using protein sequences directly: {query_name}")
        protein_file = query_file
        best_code = None

        # Reformat protein sequences
        reformatted_file = query_output_dir / "query_faa" / f"{query_name}.faa"
        reformatted_file.parent.mkdir(exist_ok=True, parents=True)

        records = list(SeqIO.parse(str(query_file), "fasta"))
        reformatted_records = reformat_sequences(records, str(query_file))
        SeqIO.write(reformatted_records, str(reformatted_file), "fasta")
        protein_file = reformatted_file

        # Create stats file for protein input
        stats_dir = query_output_dir / "stats"
        stats_dir.mkdir(exist_ok=True)
        stats_file = stats_dir / f"{query_name}.stats.tab"

        # Count proteins and estimate genome length
        protein_count = len(reformatted_records)
        total_aa_length = sum(len(rec.seq) for rec in reformatted_records)
        # Estimate genome length: AA length * 3 * 1.1 (assuming ~90% coding density)
        estimated_genome_length = int(total_aa_length * 3 * 1.1)

        with open(stats_file, "w") as f:
            f.write("ttable\tno_fna\n")
            f.write(f"genes\t{protein_count}\n")
            f.write(f"bases\t{estimated_genome_length}\n")
            f.write("coding_density\t0.0\n")

        # Also create a basic stats TSV for compatibility
        stats_tsv = stats_dir / f"{query_name}_stats.tsv"
        with open(stats_tsv, "w") as f:
            f.write("contigs\tLENbp\tGCperc\tgenecount\tCODINGperc\tttable\n")
            f.write(
                f"0\t{estimated_genome_length}\t0.0\t{protein_count}\t0.0\tno_fna\n"
            )

    # Step 3: HMM search
    logger.info(f"Running HMM search: {query_name}")
    hmmout_dir = query_output_dir / "hmmout"
    hmmout_dir.mkdir(exist_ok=True)

    # Run HMM search with multiple HMM files
    hmm_files = [
        str(database_path / "models" / "busco69.hmm"),
        str(database_path / "models" / "genomad_hmm_v1.7_selection20.hmm"),
        str(database_path / "models" / "GVOG9.hmm"),
        str(database_path / "models" / "mirus.hmm"),
        str(database_path / "models" / "MRYA.hmm"),
        str(database_path / "models" / "OGv2_order.hmm"),
        str(database_path / "models" / "UNI56.hmm"),
    ]

    # Filter to only existing files
    hmm_files = [f for f in hmm_files if Path(f).exists()]

    if not hmm_files:
        # Fallback to single combined.hmm if it exists
        combined_hmm = database_path / "models" / "combined.hmm"
        if combined_hmm.exists():
            hmm_files = [str(combined_hmm)]
            get_run_logger().info("Using fallback combined.hmm file")
        else:
            raise ProcessingError("No HMM model files found", step="hmm_search")

    cutoffs_file = None  # No longer needed - cutoffs are extracted from HMM headers

    # Output files matching expected names
    models_out = hmmout_dir / "models.out"
    models_out_filtered = hmmout_dir / "models.out.filtered"
    models_counts = hmmout_dir / "models.counts"
    models_score = hmmout_dir / "models.score"

    # Check if we should use multi-file search
    if len(hmm_files) > 1:
        # Use multi-file HMM search
        from src.core.hmm_search_multi import run_pyhmmer_search_multi

        run_pyhmmer_search_multi(
            hmm_files=hmm_files,
            query_file=str(protein_file),
            output_file=str(models_out),
            threads=threads,
        )

        # For multi-file search, filtering happens during search
        # Copy output to filtered output
        shutil.copy2(str(models_out), str(models_out_filtered))

        # Generate counts file from output
        from src.core.hmm_search import generate_counts_from_output

        generate_counts_from_output(str(models_out_filtered), str(models_counts))

    else:
        # Use single-file HMM search (original method)
        run_pyhmmer_search_with_filtering(
            query_file=str(protein_file),
            hmm_file=hmm_files[0],
            output_file=str(models_out),
            filtered_output=str(models_out_filtered),
            cutoffs_file=str(cutoffs_file),
            counts_file=str(models_counts),
            score_file=str(models_score),
            threads=threads,
        )

    # Load HMM results
    from src.core.marker_extraction import parse_hmm_output

    hmm_results = parse_hmm_output(str(models_out_filtered))

    # Step 4: BLAST search
    logger.info(f"Running BLAST search: {query_name}")
    blast_dir = query_output_dir / "blastp_out"
    blast_dir.mkdir(exist_ok=True)

    # Run BLAST for each marker hit
    blast_results = {}
    if hmm_results:
        from src.core.blast import run_blastp

        for marker in hmm_results:
            ref_faa = database_path / "database" / "alignment" / f"{marker}.ref.faa"
            if ref_faa.exists():
                blast_out = blast_dir / f"{query_name}.{marker}.blastpout"
                run_blastp(
                    queryfaa=str(protein_file),
                    refdb=str(ref_faa),
                    blastpout=str(blast_out),
                    threads=threads,
                )
                blast_results[marker] = str(blast_out)

    # Step 5: Extract marker hits and process markers for tree building
    logger.info(f"Extracting and processing markers: {query_name}")

    # Extract marker hits from HMM results
    query_hits_dir = query_output_dir / "query_hits_faa"
    query_hits_dir.mkdir(exist_ok=True)

    marker_files = []
    if models_out_filtered.exists():
        marker_files = extract_marker_hits(
            hmm_output=models_out_filtered,
            query_faa=protein_file,
            output_dir=query_hits_dir,
            mode_fast=mode_fast,
        )

    # Process each marker (BLAST, align, tree) in parallel
    marker_results = {}
    tree_nn_results = {}

    if marker_files:
        # Calculate optimal thread distribution for marker processing
        # Strategy: Balance between parallel markers and threads per marker
        num_markers = len(marker_files)

        # If we have more threads than markers, give each marker multiple threads
        if threads >= num_markers * 2:
            threads_per_marker = min(threads // num_markers, 4)  # Cap at 4
            max_parallel = num_markers
        else:
            # Otherwise, prioritize running more markers in parallel
            threads_per_marker = 1
            max_parallel = min(threads, num_markers)

        logger.info(
            f"Processing {len(marker_files)} markers: {max_parallel} parallel × {threads_per_marker} threads each"
        )

        def process_single_marker(marker_data):
            """Process a single marker - used for parallel execution."""
            marker, marker_faa = marker_data
            logger.debug(f"Processing marker {marker}, file: {marker_faa}")

            if not marker_faa.exists():
                logger.warning(f"Marker file does not exist: {marker_faa}")
                return None, None

            # Check if marker database exists
            try:
                marker_db = get_marker_database_path(marker, database_path)
                if not marker_db.exists():
                    logger.warning(f"Marker database does not exist: {marker_db}")
                    return None, None
            except FileNotFoundError:
                logger.warning(f"Marker database not found for {marker}")
                return None, None

            # Process this marker
            processor = MarkerProcessor(marker, database_path, query_output_dir)
            result = processor.process_marker(
                marker_faa,
                max_blast_hits=100,
                tree_method=tree_method,
                threads=threads_per_marker,  # Dynamic thread allocation based on available resources
            )

            logger.debug(f"Marker {marker} processing result: {result is not None}")
            if result and result.get("tree"):
                logger.info(
                    f"Successfully processed marker {marker} with tree file: {result['tree']}"
                )
                return marker, result
            else:
                logger.warning(
                    f"Marker {marker} processing did not produce a tree file"
                )
                return None, None

        # Use ThreadPoolExecutor to process markers in parallel
        # Use calculated max_parallel to avoid oversubscription
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            # Submit all marker processing tasks
            futures = [
                executor.submit(process_single_marker, marker_data)
                for marker_data in marker_files
            ]

            # Collect results as they complete
            for future in as_completed(futures):
                try:
                    marker, result = future.result()
                    if marker and result:
                        marker_results[marker] = result
                        logger.info(f"Completed marker: {marker}")
                except Exception as e:
                    logger.error(f"Error processing marker: {e}")

    # Get nearest neighbors from trees
    logger.info(f"Marker results: {len(marker_results)} markers processed")
    if marker_results:
        tree_dir = query_output_dir / "queryrefs_genetrees"
        logger.info(f"Looking for tree files in {tree_dir}")
        logger.info(f"Tree directory exists: {tree_dir.exists()}")

        if tree_dir.exists():
            tree_files = list(tree_dir.glob("*.treefile"))
            logger.info(f"Found {len(tree_files)} tree files")
            if tree_files:
                # Show first few tree files
                for i, tf in enumerate(tree_files[:5]):
                    logger.info(f"  Tree file {i+1}: {tf.name}")

                labels_file = database_path / "gvclassSeptember25_labels.tsv"
                logger.info(
                    f"Labels file: {labels_file}, exists: {labels_file.exists()}"
                )

                try:
                    analyzer = TreeAnalyzer(labels_file)
                    logger.info(
                        f"TreeAnalyzer initialized with labels file: {labels_file}"
                    )

                    # Ensure stats directory exists
                    stats_dir = query_output_dir / "stats"
                    stats_dir.mkdir(exist_ok=True, parents=True)
                    logger.info(f"Created stats directory: {stats_dir}")

                    tree_nn_file = stats_dir / f"{query_name}.tree_nn"
                    logger.info(f"Will write tree_nn results to {tree_nn_file}")

                    tree_nn_results = analyzer.process_marker_trees(
                        tree_dir, query_name, tree_nn_file
                    )
                    logger.info(
                        f"Tree analysis complete. Results: {len(tree_nn_results)} markers"
                    )
                except Exception as e:
                    logger.error(f"Tree analysis failed: {e}")
                    import traceback

                    logger.error(f"Traceback: {traceback.format_exc()}")
                    tree_nn_results = {}

                # Check if file was created
                if tree_nn_file.exists():
                    logger.info(f"Tree_nn file created successfully: {tree_nn_file}")
                    logger.info(
                        f"Tree_nn file size: {tree_nn_file.stat().st_size} bytes"
                    )
                else:
                    logger.error(f"Tree_nn file was NOT created: {tree_nn_file}")
            else:
                logger.warning("Tree directory exists but no .treefile files found")
        else:
            logger.warning(f"Tree directory does not exist: {tree_dir}")

    # Step 6: Use FullSummarizer to generate complete results
    logger.info(f"Generating full summary: {query_name}")

    try:
        # Initialize the full summarizer
        summarizer = FullSummarizer(database_path)

        # Generate the full summary with all metrics
        summary_data = summarizer.summarize_query_full(
            query_id=query_name,
            query_output_dir=query_output_dir,
            tree_nn_results=tree_nn_results,
            mode_fast=mode_fast,
        )
        logger.info(f"Summary generated successfully for {query_name}")
    except Exception as e:
        logger.error(f"Failed to generate summary for {query_name}: {e}")
        import traceback

        logger.error(traceback.format_exc())
        # Create minimal summary data to allow pipeline to continue
        summary_data = {
            "query": query_name,
            "length": "",
            "genecount": "",
            "coding_density": "",
            "ttable": "",
            "classification": "Error during classification",
        }

    # Step 7: Create output directories and files
    query_faa_dir = query_output_dir / "query_faa"
    query_fna_dir = query_output_dir / "query_fna"
    query_gff_dir = query_output_dir / "query_gff"

    for dir_path in [query_faa_dir, query_fna_dir, query_gff_dir]:
        dir_path.mkdir(exist_ok=True, parents=True)

    # Copy final files
    if query_file.suffix == ".fna" and best_code is not None:
        # Copy protein file
        if protein_file.exists():
            shutil.copy2(protein_file, query_faa_dir / f"{query_name}.faa")

        # Copy nucleotide file
        if reformatted_file.exists():
            shutil.copy2(reformatted_file, query_fna_dir / f"{query_name}.fna")

        # Copy GFF file if exists
        # First try the main gene_calling directory
        gff_file = (
            query_output_dir
            / "gene_calling"
            / f"{query_name}_reformatted"
            / f"{query_name}_reformatted.gff"
        )
        if gff_file.exists():
            shutil.copy2(gff_file, query_gff_dir / f"{query_name}.gff")
        else:
            # Try alternative location
            gff_file = query_output_dir / "gene_calling" / f"{query_name}.gff"
            if gff_file.exists():
                shutil.copy2(gff_file, query_gff_dir / f"{query_name}.gff")
            else:
                # Try looking in output directory for GFF files
                for gff in query_output_dir.glob("**/*.gff"):
                    shutil.copy2(gff, query_gff_dir / f"{query_name}.gff")
                    break
                else:
                    logger.warning(f"No GFF file found for {query_name}")

    # Create individual summary file with full data
    individual_summary = query_output_dir / f"{query_name}.summary.tab"

    try:
        logger.info(f"Writing summary file to {individual_summary}")
        # Write comprehensive summary file with LEGACY headers for backward compatibility
        with open(individual_summary, "w") as f:
            # Headers matching the OLD format exactly (minus taxonomy_strict)
            headers = [
                "query",
                "taxonomy_majority",
                "species",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "domain",
                "avgdist",
                "order_dup",
                "order_completeness",
                "gvog4_unique",
                "gvog8_unique",
                "gvog8_total",
                "gvog8_dup",
                "mcp_total",
                "mirus_unique",
                "mirus_total",
                "mirus_dup",
                "mrya_unique",
                "mrya_total",
                "phage_unique",
                "phage_total",
                "cellular_unique",
                "cellular_total",
                "cellular_dup",
                "contigs",
                "LENbp",
                "GCperc",
                "genecount",
                "CODINGperc",
                "ttable",
                "weighted_order_completeness",  # Keep the new metric at the end
            ]
            f.write("\t".join(headers) + "\n")

            # Map from summary_data keys (returned by FullSummarizer) to header names
            # Updated for LEGACY format - raw counts not percentages
            key_mapping = {
                "query": "query",
                "taxonomy_majority": "taxonomy_majority",
                "species": "species",
                "genus": "genus",
                "family": "family",
                "order": "order",
                "class": "class",
                "phylum": "phylum",
                "domain": "domain",
                "avgdist": "avgdist",
                "order_dup": "order_dup",
                "order_completeness": "order_completeness",  # Old % based metric
                "gvog4_unique": "gvog4_unique",  # Raw count
                "gvog8_unique": "gvog8_unique",  # Raw count
                "gvog8_total": "gvog8_total",
                "gvog8_dup": "gvog8_dup",
                "mcp_total": "mcp_total",
                "mirus_unique": "mirus_unique",  # Raw count
                "mirus_total": "mirus_total",
                "mirus_dup": "mirus_dup",
                "mrya_unique": "mrya_unique",  # Raw count
                "mrya_total": "mrya_total",
                "phage_unique": "phage_unique",  # Raw count
                "phage_total": "phage_total",
                "cellular_unique": "cellular_unique",
                "cellular_total": "cellular_total",
                "cellular_dup": "cellular_dup",
                "contigs": "contigs",
                "LENbp": "LENbp",
                "GCperc": "GCperc",
                "genecount": "genecount",
                "CODINGperc": "CODINGperc",
                "ttable": "ttable",
                "order_weighted_completeness": "weighted_order_completeness",  # Keep new metric
            }

            # Write data row from summary_data
            row_data = []
            for header in headers:
                # Find the corresponding key in summary_data
                for data_key, header_name in key_mapping.items():
                    if header_name == header:
                        value = summary_data.get(data_key, "")
                        row_data.append(str(value))
                        break
                else:
                    # If no mapping found, try direct match
                    row_data.append(str(summary_data.get(header, "")))

            f.write("\t".join(row_data) + "\n")

        logger.info(f"Summary file written successfully: {individual_summary}")
    except Exception as e:
        logger.error(f"Failed to write summary file for {query_name}: {e}")
        import traceback

        logger.error(traceback.format_exc())

    # Post-processing: Create tar.gz and copy summary.tab
    logger.info(f"Post-processing results for query {query_name}")
    logger.debug(f"Query output directory: {query_output_dir}")
    logger.debug(f"Output base directory: {output_base}")

    post_process_start = time.time()
    try:
        import tarfile

        # Copy summary.tab file to output root - it's in the query directory root
        summary_tab_file = query_output_dir / f"{query_name}.summary.tab"
        logger.debug(f"Looking for summary file: {summary_tab_file}")
        if summary_tab_file.exists():
            file_size = summary_tab_file.stat().st_size
            logger.debug(f"Found summary file, size: {file_size} bytes")
            dest_summary = output_base / f"{query_name}.summary.tab"
            shutil.copy2(summary_tab_file, dest_summary)
            logger.info(f"✓ Copied summary.tab to {dest_summary} ({file_size} bytes)")
        else:
            # Try in stats directory as backup
            summary_tab_file_stats = (
                query_output_dir / "stats" / f"{query_name}.summary.tab"
            )
            if summary_tab_file_stats.exists():
                dest_summary = output_base / f"{query_name}.summary.tab"
                shutil.copy2(summary_tab_file_stats, dest_summary)
                logger.info(f"Copied summary.tab from stats dir to {dest_summary}")
            else:
                logger.warning(
                    f"Summary tab file not found in {summary_tab_file} or {summary_tab_file_stats}"
                )

        # tree_nn file stays in the stats directory (inside the tar.gz archive)
        # No longer copying to output root to keep it cleaner
        tree_nn_file = query_output_dir / "stats" / f"{query_name}.tree_nn"
        if tree_nn_file.exists():
            file_size = tree_nn_file.stat().st_size
            logger.debug(f"Found tree_nn file in stats directory: {file_size} bytes")

        # Create tar.gz archive of the query directory
        tar_file = output_base / f"{query_name}.tar.gz"
        logger.debug(f"Creating archive: {tar_file}")

        # Count files to be archived
        file_count = sum(1 for _ in query_output_dir.rglob("*") if _.is_file())
        logger.debug(f"Archiving {file_count} files from {query_output_dir}")

        with tarfile.open(tar_file, "w:gz") as tar:
            tar.add(query_output_dir, arcname=query_name)

        archive_size = tar_file.stat().st_size
        logger.info(
            f"✓ Created archive: {tar_file} ({archive_size:,} bytes, {file_count} files)"
        )

        # Clean up the reformatted file from the output root
        reformatted_file = output_base / f"{query_name}_reformatted.fna"
        if reformatted_file.exists():
            reformatted_file.unlink()
            logger.info(f"✓ Removed reformatted file: {reformatted_file}")

        # Remove the original directory
        logger.debug(f"Removing original directory: {query_output_dir}")
        shutil.rmtree(query_output_dir)
        logger.info(f"✓ Removed original directory: {query_output_dir}")

        post_process_time = time.time() - post_process_start
        logger.info(f"✓ Post-processing completed in {post_process_time:.2f} seconds")

    except Exception as e:
        post_process_time = time.time() - post_process_start
        logger.error(
            f"❌ Post-processing failed for {query_name} after {post_process_time:.2f} seconds: {e}"
        )
        import traceback

        logger.error(f"Traceback: {traceback.format_exc()}")

        # Log what was partially completed
        logger.debug("Checking partial outputs:")
        if (output_base / f"{query_name}.summary.tab").exists():
            logger.debug("  - Summary file was copied")
        if (output_base / f"{query_name}.tar.gz").exists():
            logger.debug("  - Archive was created")
        if query_output_dir.exists():
            logger.debug("  - Original directory still exists")

        # Re-raise to ensure the error is visible
        raise RuntimeError(f"Post-processing failed for {query_name}: {e}") from e

    # Return complete results
    return {
        "query": query_name,
        "status": "complete",
        "summary_data": summary_data,  # Full data from FullSummarizer
        "genetic_code": best_code if best_code is not None else "N/A",
        "output_dir": str(
            output_base
        ),  # Return base output dir since query dir is removed
    }


@task(name="create_final_summary", **task_config)
def create_final_summary_task(results: List[Dict[str, Any]], output_dir: Path) -> Path:
    """Create the final summary file combining all results."""
    logger = get_run_logger()
    logger.info(f"Creating final summary for {len(results)} queries")

    summary_file = output_dir / "gvclass_summary.tsv"

    # Define the expected column order
    columns = [
        "query",
        "taxonomy_majority",
        "taxonomy_strict",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "domain",
        "avgdist",
        "order_dup",
        "order_completeness",
        "order_weighted_completeness",
        "order_confidence_score",
        "gvog4_unique",
        "gvog8_unique",
        "gvog8_total",
        "gvog8_dup",
        "mcp_total",
        "mirus_unique",
        "mirus_total",
        "mirus_dup",
        "mrya_unique",
        "mrya_total",
        "phage_unique",
        "phage_total",
        "cellular_unique",
        "cellular_total",
        "cellular_dup",
        "contigs",
        "LENbp",
        "GCperc",
        "genecount",
        "CODINGperc",
        "ttable",
    ]

    # Write the summary file
    with open(summary_file, "w") as f:
        # Write header
        f.write("\t".join(columns) + "\n")

        # Write results
        for result in sorted(results, key=lambda x: x["query"]):
            if result["status"] == "complete" and "summary_data" in result:
                # Get the full summary data
                summary_data = result["summary_data"]

                # Write row with all columns
                row_data = []
                for col in columns:
                    value = summary_data.get(col, "")
                    # Format float values
                    if isinstance(value, float):
                        if col in [
                            "avgdist",
                            "order_dup",
                            "gvog8_dup",
                            "mirus_dup",
                            "cellular_dup",
                        ]:
                            value = f"{value:.2f}"
                        elif col in [
                            "order_completeness",
                            "order_weighted_completeness",
                            "order_confidence_score",
                            "GCperc",
                            "CODINGperc",
                        ]:
                            value = f"{value:.2f}"
                        else:
                            value = f"{value:.0f}"
                    row_data.append(str(value))

                f.write("\t".join(row_data) + "\n")
            else:
                # Failed query - write minimal data
                row_data = [result["query"]] + [""] * (len(columns) - 1)
                f.write("\t".join(row_data) + "\n")

    logger.info(f"Summary written to: {summary_file}")

    # Post-processing is already done in process_query_task
    # tree_nn files are kept inside the tar.gz archives to keep output directory clean
    logger.info("All query processing complete")

    logger.info("Post-processing complete")
    return summary_file


def calculate_optimal_workers(
    n_queries: int,
    total_threads: int,
    min_threads_per_worker: int = 4,
    max_workers: Optional[int] = None,
) -> Tuple[int, int]:
    """
    Calculate optimal number of workers and threads per worker.

    Strategy:
    - Each worker needs at least min_threads_per_worker threads
    - Maximize parallelism while maintaining thread efficiency
    - For few queries with many threads, use fewer workers with more threads
    - For many queries with limited threads, use more workers with fewer threads
    """
    if max_workers is None:
        max_workers = n_queries

    # Calculate possible configurations
    max_possible_workers = min(
        n_queries,  # No more workers than queries
        total_threads // min_threads_per_worker,  # Respect minimum threads
        max_workers,  # Respect user limit
    )

    # For small numbers of queries, prefer fewer workers with more threads
    if n_queries <= 4 and total_threads >= n_queries * 8:
        n_workers = n_queries
        threads_per_worker = total_threads // n_workers
    else:
        # For larger datasets, be conservative with workers to limit memory usage
        # Target 4 threads per worker as the sweet spot
        target_threads_per_worker = 4

        # Calculate workers based on available threads
        n_workers = total_threads // target_threads_per_worker

        # Apply constraints
        n_workers = min(n_workers, max_possible_workers)

        # For very large datasets, further limit workers to prevent memory issues
        if n_queries > 100:
            # Limit to 10-20 workers max for large datasets
            n_workers = min(n_workers, 20)

        n_workers = max(1, n_workers)
        threads_per_worker = total_threads // n_workers

    return max(1, n_workers), max(min_threads_per_worker, threads_per_worker)


@flow(
    name="gvclass-pipeline",
    description="GVClass pipeline with dynamic Dask execution",
    persist_result=True,
)
def gvclass_flow(
    query_dir: str,
    output_dir: str,
    database_path: Optional[str] = None,
    total_threads: int = 16,
    max_workers: Optional[int] = None,
    threads_per_worker: Optional[int] = None,
    tree_method: str = "fasttree",
    mode_fast: bool = True,
    genetic_codes: List[int] = [0, 1, 4, 6, 11, 15, 106, 129],
    cluster_type: str = "local",
    cluster_config: Optional[Dict[str, Any]] = None,
    resume: bool = False,
):
    """
    Main GVClass pipeline with proper Prefect flow orchestration.

    This flow dynamically configures the DaskTaskRunner based on:
    - Number of queries to process
    - Available threads
    - Cluster type (local, slurm, etc.)

    Args:
        query_dir: Directory containing query sequences
        output_dir: Output directory for results
        database_path: Path to GVClass database
        total_threads: Total threads available for computation
        max_workers: Maximum number of parallel workers (None = auto)
        threads_per_worker: Threads per worker (None = auto)
        tree_method: Tree building method (fasttree or iqtree)
        mode_fast: Skip order-level marker trees if True
        genetic_codes: List of genetic codes to test
        cluster_type: Type of cluster (local, slurm, pbs, sge)
        cluster_config: Additional cluster configuration
    """
    logger = get_run_logger()

    # Step 1: Validate inputs and setup
    logger.info("Validating inputs and setting up directories")
    config = validate_and_setup_task(query_dir, output_dir, database_path)

    # Filter out completed queries if resume is enabled
    if resume:
        output_path = config["output_path"]
        original_queries = list(config["query_files"])
        filtered_queries = []
        skipped_count = 0

        for query_file in original_queries:
            query_name = query_file.stem
            summary_file = output_path / f"{query_name}.summary.tab"
            tar_file = output_path / f"{query_name}.tar.gz"

            # Skip if both summary.tab and tar.gz exist
            if summary_file.exists() and tar_file.exists():
                logger.info(f"Skipping completed query: {query_name}")
                skipped_count += 1
            else:
                filtered_queries.append(query_file)

        config["query_files"] = filtered_queries
        logger.info(
            f"Resume mode: skipped {skipped_count} completed queries, {len(filtered_queries)} remaining"
        )

    n_queries = len(config["query_files"])
    logger.info(f"Found {n_queries} queries to process")

    # Step 2: Calculate optimal worker distribution
    if threads_per_worker is None:
        n_workers, threads_per_worker = calculate_optimal_workers(
            n_queries=n_queries, total_threads=total_threads, max_workers=max_workers
        )
    else:
        n_workers = min(
            n_queries, total_threads // threads_per_worker, max_workers or n_queries
        )

    logger.info(
        f"Parallelization strategy: {n_workers} workers × {threads_per_worker} threads"
    )

    # Step 3: Process queries in parallel using ThreadPoolExecutor
    # This avoids the multiprocessing spawn issues with Dask/Prefect
    logger.info("Starting parallel query processing")

    results = []
    query_files = list(config["query_files"])
    total_queries = len(query_files)

    logger.info(
        f"Processing {total_queries} queries with {n_workers} workers ({threads_per_worker} threads each)"
    )

    # Use ThreadPoolExecutor directly for simpler, more reliable parallelization
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Submit all queries
        future_to_query = {}
        for query_file in query_files:
            future = executor.submit(
                lambda qf: process_query_task(
                    query_file=qf,
                    output_base=config["output_path"],
                    database_path=config["database_path"],
                    genetic_codes=genetic_codes,
                    tree_method=tree_method,
                    mode_fast=mode_fast,
                    threads=threads_per_worker,
                ),
                query_file,
            )
            future_to_query[future] = query_file

        # Collect results as they complete
        for future in as_completed(future_to_query):
            query_file = future_to_query[future]
            try:
                result = future.result()
                results.append(result)
                logger.info(f"Completed: {result['query']} - {result['status']}")
            except Exception as e:
                logger.error(f"Query {query_file.stem} failed: {e}")
                import traceback
                logger.error(f"Traceback: {traceback.format_exc()}")
                results.append(
                    {"query": query_file.stem, "status": "failed", "error": str(e)}
                )

    logger.info(
        f"All queries complete. Successfully processed {len([r for r in results if r['status'] == 'complete'])}/{total_queries} queries"
    )

    # Create final summary
    logger.info("Creating final summary...")
    summary_file = create_final_summary_task(results, config["output_path"])

    logger.info(f"Pipeline completed! Results in: {config['output_path']}")
    logger.info(f"Summary written to: {summary_file}")
    logger.info(
        f"Successfully processed {len([r for r in results if r['status'] == 'complete'])}/{len(results)} queries"
    )

    return summary_file


# Deployment configurations for different environments
def create_local_deployment():
    """Create deployment for local execution."""
    from prefect.deployments import Deployment

    deployment = Deployment.build_from_flow(
        flow=gvclass_flow,
        name="gvclass-local",
        parameters={"cluster_type": "local", "total_threads": 16, "mode_fast": True},
        tags=["gvclass", "local"],
        description="GVClass pipeline for local execution",
    )

    return deployment


def create_slurm_deployment():
    """Create deployment for SLURM cluster execution."""
    from prefect.deployments import Deployment

    deployment = Deployment.build_from_flow(
        flow=gvclass_flow,
        name="gvclass-slurm",
        parameters={
            "cluster_type": "slurm",
            "total_threads": 128,
            "cluster_config": {
                "queue": "normal",
                "project": "gvclass",
                "walltime": "08:00:00",
            },
        },
        tags=["gvclass", "hpc", "slurm"],
        description="GVClass pipeline for SLURM cluster execution",
    )

    return deployment
