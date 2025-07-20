from typing import List, Tuple, Dict
from collections import defaultdict
import os
import click
import pandas as pd
import pyhmmer
from pathlib import Path

from src.utils.common import validate_file_path, GVClassError
from src.utils.error_handling import error_handler, ProcessingError


def load_cutoffs_simple(cutoffs_file: str) -> Dict[str, float]:
    """
    Load cutoff values for each HMM model (score only).
    
    Args:
        cutoffs_file: Path to cutoffs file
        
    Returns:
        Dictionary mapping model names to cutoff scores
    """
    cutoffs = {}
    try:
        with open(cutoffs_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#') and not line.startswith('Model'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        model_name = parts[0]
                        cutoff_score = float(parts[1])
                        cutoffs[model_name] = cutoff_score
    except Exception as e:
        error_handler.log_error(f"Error loading cutoffs: {e}")
        
    return cutoffs


def run_pyhmmer_search_with_filtering(
    query_file: str,
    hmm_file: str,
    output_file: str,
    filtered_output: str,
    cutoffs_file: str,
    counts_file: str,
    score_file: str,
    threads: int = 4
) -> None:
    """
    Run HMM search with cutoff filtering and generate all output files.
    
    Args:
        query_file: Path to query sequences
        hmm_file: Path to HMM models
        output_file: Path to raw output
        filtered_output: Path to filtered output
        cutoffs_file: Path to cutoff values
        counts_file: Path to counts output
        score_file: Path to score output
        threads: Number of threads
    """
    # Log inputs
    error_handler.log_info("Starting HMM search with filtering:")
    error_handler.log_info(f"  Query file: {query_file}, exists: {Path(query_file).exists()}")
    error_handler.log_info(f"  HMM file: {hmm_file}, exists: {Path(hmm_file).exists()}")
    error_handler.log_info(f"  Output file: {output_file}")
    error_handler.log_info(f"  Cutoffs file: {cutoffs_file}, exists: {Path(cutoffs_file).exists()}")
    
    # First run the search
    run_pyhmmer_search(hmm_file, query_file, output_file, threads)
    
    # Check output file
    error_handler.log_info(f"HMM search output file size: {Path(output_file).stat().st_size if Path(output_file).exists() else 'NOT FOUND'}")
    
    # Load cutoffs
    cutoffs = load_cutoffs_simple(cutoffs_file)
    error_handler.log_info(f"Loaded {len(cutoffs)} cutoff values")
    
    # Parse and filter results
    model_hits = defaultdict(list)
    filtered_lines = []
    total_lines = 0
    hit_lines = 0
    
    with open(output_file, 'r') as f:
        for line in f:
            total_lines += 1
            if line.startswith('#'):
                filtered_lines.append(line)
                continue
                
            parts = line.strip().split()
            if len(parts) >= 14:
                hit_lines += 1
                query_name = parts[0]
                model_name = parts[3]
                score = float(parts[13])
                
                # Check cutoff
                if model_name in cutoffs and score >= cutoffs[model_name]:
                    filtered_lines.append(line)
                    model_hits[model_name].append((query_name, score))
    
    error_handler.log_info(f"Parsed {total_lines} lines, {hit_lines} hits, {len(filtered_lines)-total_lines+hit_lines} passed cutoffs")
                    
    # Write filtered output
    with open(filtered_output, 'w') as f:
        f.writelines(filtered_lines)
    
    error_handler.log_info(f"Wrote {len(filtered_lines)} lines to filtered output")
    error_handler.log_info(f"Found hits for {len(model_hits)} unique models")
        
    # Write counts
    with open(counts_file, 'w') as f:
        for model, hits in sorted(model_hits.items()):
            f.write(f"{model}\t{len(hits)}\n")
            
    # Write scores
    with open(score_file, 'w') as f:
        total_score = sum(max(score for _, score in hits) for hits in model_hits.values() if hits)
        f.write(f"Total score: {total_score}\n")
        for model, hits in sorted(model_hits.items()):
            if hits:
                max_score = max(score for _, score in hits)
                f.write(f"{model}\t{max_score}\n")


def run_pyhmmer_search(hmm_file: str, query_file: str, output_file: str, threads: int = 4) -> None:
    """
    Run HMM search using pyhmmer instead of external hmmsearch command.
    
    Args:
        hmm_file: Path to HMM model file
        query_file: Path to query sequence file
        output_file: Path to output file
        threads: Number of threads to use
    """
    try:
        error_handler.log_info(f"Running pyhmmer search on {query_file} with {hmm_file}")
        
        # Read HMM profiles
        with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
            profiles = list(hmm_handle)
        
        # Read query sequences
        with pyhmmer.easel.SequenceFile(query_file, digital=True) as seq_handle:
            sequences = list(seq_handle)
        
        # Run search with E-value threshold
        results = pyhmmer.hmmsearch(
            profiles, 
            sequences, 
            cpus=threads,
            E=10.0,  # E-value threshold
            domE=10.0  # Domain E-value threshold
        )
        
        # Write results in simplified domtblout format
        with open(output_file, 'w') as out_handle:
            # Write header
            out_handle.write("# target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\tbias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\tfrom\tto\tfrom\tto\tfrom\tto\tacc\n")
            
            # Iterate through results - each TopHits corresponds to one query HMM
            for i, (hmm, top_hits) in enumerate(zip(profiles, results)):
                hmm_name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                hmm_accession = hmm.accession.decode() if hmm.accession and isinstance(hmm.accession, bytes) else (hmm.accession or "-")
                hmm_length = hmm.M
                
                # Process hits for this HMM
                for hit in top_hits:
                    target_name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                    target_accession = "-"  # Simplified
                    target_length = 1000  # Placeholder - will be parsed differently downstream
                    
                    # Process domains
                    for domain_idx, domain in enumerate(hit.domains):
                        # Get alignment for coordinate information
                        if domain.alignment:
                            ali = domain.alignment
                            hmm_from = ali.hmm_from
                            hmm_to = ali.hmm_to
                            target_from = ali.target_from
                            target_to = ali.target_to
                            target_length = ali.target_length
                        else:
                            # Fallback if no alignment
                            hmm_from = 1
                            hmm_to = hmm_length
                            target_from = domain.env_from
                            target_to = domain.env_to
                            
                        # Simplified output - focus on essential fields
                        out_handle.write(f"{target_name}\t{target_accession}\t{target_length}\t")
                        out_handle.write(f"{hmm_name}\t{hmm_accession}\t{hmm_length}\t")
                        out_handle.write(f"{hit.evalue:.2e}\t{hit.score:.1f}\t{hit.bias:.1f}\t")
                        out_handle.write(f"{domain_idx+1}\t{len(hit.domains)}\t")
                        out_handle.write(f"{domain.c_evalue:.2e}\t{domain.i_evalue:.2e}\t")
                        out_handle.write(f"{domain.score:.1f}\t{domain.bias:.1f}\t")
                        # Use 1-based coordinates
                        out_handle.write(f"{hmm_from+1}\t{hmm_to}\t")
                        out_handle.write(f"{target_from+1}\t{target_to}\t")
                        out_handle.write(f"{domain.env_from+1}\t{domain.env_to}\t")
                        out_handle.write("0.90\n")  # Default accuracy
        
        error_handler.log_info(f"pyhmmer search completed. Results written to {output_file}")
        
    except Exception as e:
        raise ProcessingError(
            f"pyhmmer search failed: {str(e)}",
            step="pyhmmer_search",
            input_file=query_file
        ) from e


def get_models(modelsin: str) -> List[str]:
    """
    Generate a list of all model names from the combined HMM model file.

    Args:
        modelsin (str): Path to the combined HMM model file.

    Returns:
        List[str]: List of model names.
    """
    with open(modelsin, "r") as f:
        models = [line.split()[1].strip() for line in f if line.startswith("NAME")]
        # only get stats for general models not order level modes ("OG")
        models_general = [x for x in models if not x.startswith("OG")]
        
    return models_general


def process_hmmout(models: List[str], hmmout: str, filtered_hmmout: str, cutoffs: Dict[str, Tuple[float, float]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process HMM output to generate count and score matrices including all models.

    Args:
        models (List[str]): List of model names.
        hmmout (str): Path to the HMM output file.
        filtered_hmmout (str): Path to the filtered HMM output file.
        cutoffs (Dict[str, Tuple[float, float]]): Dictionary mapping model names to (score_cutoff, length_cutoff) tuples.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Tuple containing count and score matrices.
    """
    scores: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    protein_hits: Dict[str, Dict[str, Tuple[str, float, int, str]]] = defaultdict(dict)
    protein_hits_order: Dict[str, Dict[str, Tuple[str, float, int]]] = defaultdict(dict)
    
    with open(hmmout, 'r') as infile, open(filtered_hmmout, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                parts = line.split()
                protein_id = parts[0]
                genomeid = protein_id.split('|')[0]
                model = parts[3]
                score = float(parts[7])
                length = int(parts[2])

                if model in cutoffs:
                    score_cutoff, length_cutoff = cutoffs[model]
                    if score >= score_cutoff and length >= length_cutoff:
                        if model.startswith("OG"):
                            key = (protein_id, model)
                            if key not in protein_hits_order[genomeid]:
                                protein_hits_order[genomeid][key] = (model, score, length, line)
                                outfile.write(line)
                            elif score > protein_hits_order[genomeid][key][1]:
                                protein_hits_order[genomeid][key] = (model, score, length, line)
                                outfile.write(line)
                        else:
                            if protein_id not in protein_hits[genomeid]:
                                protein_hits[genomeid][protein_id] = (model, score, length, line)
                            elif score > protein_hits[genomeid][protein_id][1]:
                                protein_hits[genomeid][protein_id] = (model, score, length, line)

        for genomeid, protein_dict in protein_hits.items():
            for protein_id, (model, score, length, line) in protein_dict.items():
                outfile.write(f"{line}\n")
    
    for genomeid, protein_dict in protein_hits.items():
        for protein_id, (model, score, length, line) in protein_dict.items():
            counts[genomeid][model] += 1
            scores[genomeid][model] = max(scores[genomeid][model], score)
            
    for genomeid, protein_dict in protein_hits_order.items():
        for (protein_id, model), (_, score, length, line) in protein_dict.items():
            counts[genomeid][model] += 1
            scores[genomeid][model] = max(scores[genomeid][model], score)
    
    # Ensuring all models are accounted for each genome
    for genomeid in counts.keys():
        for model in models:
            counts[genomeid].setdefault(model, 0)
            scores[genomeid].setdefault(model, 0.0)

    # Converting nested dictionaries to pandas DataFrame
    count_df = pd.DataFrame(counts).T.fillna(0).astype(int)
    score_df = pd.DataFrame(scores).T.fillna(0)

    return count_df, score_df


def load_cutoffs(cutoff_file: str) -> Dict[str, Tuple[float, float]]:
    """
    Load cutoff values from a file.

    Args:
        cutoff_file (str): Path to the cutoff file.

    Returns:
        Dict[str, Tuple[float, float]]: Dictionary mapping model names to (score_cutoff, length_cutoff) tuples.
    """
    cutoffs = {}
    try:
        with open(cutoff_file, 'r') as file:
            next(file)  # Skip header
            for line in file:
                model, score_cutoff, length_cutoff = line.strip().split('\t')
                cutoffs[model] = (float(score_cutoff), float(length_cutoff))
    except Exception as e:
        print(f"Error reading {cutoff_file}: {e}")
    return cutoffs


@click.command()
@click.option('--queryfaa', '-q', type=click.Path(exists=True), required=True, help='Input query FASTA file')
@click.option('--modelscombined', '-m', type=click.Path(exists=True), required=True, help='Combined HMM model file')
@click.option('--hmmout', '-h', type=click.Path(), required=True, help='Output file for HMM hits')
@click.option('--filtered_hmmout', '-hf', type=click.Path(), required=True, help='Output file for HMM hits filtered')
@click.option('--countout', '-c', type=click.Path(), required=True, help='Output file for marker gene counts')
@click.option('--scoreout', '-s', type=click.Path(), required=True, help='Output file for highest bitscore matrix')
@click.option('--cutoff_file', '-f', type=click.Path(exists=True), help='Cutoff file for filtering HMM hits')
@click.option('--threads', '-t', default=4, help='Number of CPU threads to use')
def main(queryfaa: str, modelscombined: str, hmmout: str, filtered_hmmout: str, countout: str, scoreout: str, cutoff_file: str = None, threads: int = 4) -> None:
    """
    Main function to run HMM search and process the output.

    Args:
        queryfaa (str): Path to the input query FASTA file.
        modelscombined (str): Path to the combined HMM model file.
        hmmout (str): Path to the output file for HMM hits.
        filtered_hmmout (str): Path to the output file for filtered HMM hits.
        countout (str): Path to the output file for marker gene counts.
        scoreout (str): Path to the output file for highest bitscore matrix.
        cutoff_file (str, optional): Path to the cutoff file for filtering HMM hits. Defaults to None.
        threads (int): Number of CPU threads to use.
    """
    try:
        # Validate input files
        query_file = validate_file_path(queryfaa, must_exist=True)
        models_file = validate_file_path(modelscombined, must_exist=True)
        output_file = validate_file_path(hmmout, must_exist=False)
        
        # Validate threads parameter
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("Threads must be a positive integer")
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Check if the hmmsearch output file already exists
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            # Run hmmsearch using pyhmmer
            run_pyhmmer_search(str(models_file), str(query_file), str(output_file), threads)
        else:
            print(f"Skipping hmmsearch as the output file {output_file} already exists.")
            
    except (ValueError, FileNotFoundError, GVClassError) as e:
        print(f"Error: {e}")
        raise click.ClickException(str(e))
        
    models = get_models(modelscombined)
    
    # Load cutoff values or set default cutoffs
    if cutoff_file:
        cutoffs = load_cutoffs(cutoff_file)
    else:
        cutoffs = {model: (0.0, 0.0) for model in models}

    # Get model names
    count_matrix, score_matrix = process_hmmout(models, hmmout, filtered_hmmout, cutoffs)

    # Check if count and score matrices are empty
    if count_matrix.empty or score_matrix.empty:
        print("No hits found in the hmmsearch output. Generating empty output files.")
        # Create empty DataFrames with 'genomeid' column and model columns
        empty_count_matrix = pd.DataFrame(columns=['genomeid'] + models).set_index('genomeid')
        empty_score_matrix = pd.DataFrame(columns=['genomeid'] + models).set_index('genomeid')
        # Save empty DataFrames to output files
        empty_count_matrix.to_csv(countout, sep='\t')
        empty_score_matrix.to_csv(scoreout, sep='\t')
    else:
        # Save count and score matrices with 'genomeid' as the first column
        count_matrix.reset_index().rename(columns={'index': 'genomeid'}).to_csv(countout, sep='\t', index=False)
        score_matrix.reset_index().rename(columns={'index': 'genomeid'}).to_csv(scoreout, sep='\t', index=False)


if __name__ == '__main__':
    main()