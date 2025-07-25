import os
from collections import defaultdict
from typing import List, Dict
import click
from Bio import SeqIO
import pyswrd

from src.utils.common import validate_file_path


def run_blastp(queryfaa: str, refdb: str, blastpout: str, threads: int = 8) -> None:
    """
    Run pyswrd for sequence similarity search (replacement for diamond blastp).

    Args:
        queryfaa (str): Path to the query FASTA file.
        refdb (str): Path to the reference database (for pyswrd, this should be the .faa file).
        blastpout (str): Path to the output file for search results.
        threads (int): Number of threads to use.
    """
    # Validate input files
    query_file = validate_file_path(queryfaa, must_exist=True)
    output_file = validate_file_path(blastpout, must_exist=False)

    # Validate threads parameter
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("Threads must be a positive integer")

    # For pyswrd, we need the actual FASTA file, not a diamond database
    # Convert refdb path from .dmnd to .faa
    if refdb.endswith(".dmnd"):
        ref_faa = refdb.replace(".dmnd", ".faa")
        # Check in the faa directory
        if not os.path.exists(ref_faa):
            faa_dir = os.path.join(os.path.dirname(os.path.dirname(refdb)), "faa")
            ref_faa = os.path.join(faa_dir, os.path.basename(ref_faa))
    else:
        ref_faa = refdb

    ref_faa = validate_file_path(ref_faa, must_exist=True)

    print(f"Running pyswrd search: {query_file} vs {ref_faa}")

    try:
        # Load sequences
        query_seqs = list(SeqIO.parse(query_file, "fasta"))
        ref_seqs = list(SeqIO.parse(ref_faa, "fasta"))

        if not query_seqs:
            print(f"No sequences found in query file: {query_file}")
            open(output_file, "a").close()
            return

        if not ref_seqs:
            print(f"No sequences found in reference file: {ref_faa}")
            open(output_file, "a").close()
            return

        # Prepare sequences for pyswrd
        # Create lists of sequence strings
        queries = [str(seq.seq) for seq in query_seqs]
        targets = [str(seq.seq) for seq in ref_seqs]

        # Open output file
        with open(output_file, "w") as out:
            # Run pyswrd search with parameters similar to diamond blastp
            # Using BLOSUM62 matrix and e-value threshold
            hits_found = 0

            # Process hits from pyswrd.search
            for hit in pyswrd.search(
                queries,
                targets,
                scorer_name="BLOSUM62",
                score_threshold=0,  # We'll filter by e-value later
                threads=threads,
            ):

                # Get query and target info
                query_idx = hit.query_index
                target_idx = hit.target_index
                score = hit.score
                evalue = hit.evalue

                # Apply e-value filter (1e-5)
                if evalue > 1e-5:
                    continue

                # Get sequence IDs
                query_id = query_seqs[query_idx].id
                target_id = ref_seqs[target_idx].id

                # For now, we'll output simplified format since pyswrd doesn't provide
                # detailed alignment info like coverage and identity
                # Format: query, subject, identity, alignment_length, mismatches, gap_opens,
                #         q_start, q_end, s_start, s_end, evalue, bit_score

                # Approximate values for missing fields
                identity = 100.0  # Not available from pyswrd
                aln_len = min(len(queries[query_idx]), len(targets[target_idx]))
                mismatches = 0
                gaps = 0
                q_start = 1
                q_end = len(queries[query_idx])
                s_start = 1
                s_end = len(targets[target_idx])

                out.write(
                    f"{query_id}\t{target_id}\t{identity:.1f}\t{aln_len}\t"
                    f"{mismatches}\t{gaps}\t{q_start}\t{q_end}\t"
                    f"{s_start}\t{s_end}\t{evalue:.2e}\t{score:.1f}\n"
                )

                hits_found += 1

        if hits_found > 0:
            print(
                f"pyswrd search completed. Found {hits_found} hits. Results written to: {output_file}"
            )
        else:
            print("pyswrd search completed. No significant hits found.")

    except Exception as e:
        print(f"Error running pyswrd: {str(e)}")
        print("Creating empty output file to allow pipeline to continue.")
        open(output_file, "a").close()


def parse_blastp(blastpout: str) -> List[str]:
    """
    Parse the diamond blastp output and return the best hits.

    Args:
        blastpout (str): Path to the diamond blastp output file.

    Returns:
        List[str]: List of best hit subject IDs.
    """
    queryids_hits_dict: Dict[str, List[str]] = defaultdict(list)
    seen: List[str] = []

    try:
        # Count unique queries using proper file handling
        with open(blastpout, "r") as query_file:
            num_queries = len(set(line.split()[0] for line in query_file))

        with open(blastpout, "r") as infile:
            for line in infile:
                queryname, subjectname = line.split()[0], line.split()[1]

                # if subjectname not in seen and len(queryids_hits_dict[queryname]) <= 100 / num_queries and pident != 100:
                if (
                    subjectname not in seen
                    and len(queryids_hits_dict[queryname]) <= 100 / num_queries
                ):
                    queryids_hits_dict[queryname].append(subjectname)
                    seen.append(subjectname)
    except Exception as e:
        print(f"Error parsing BLAST output: {e}")
        return []

    besthits = [x for v in queryids_hits_dict.values() for x in v]
    return besthits


def reduce_ref(
    queryfaa: str, reffaa: str, reduced_reffaa_merged: str, refdb: str, blastpout: str
) -> None:
    """
    Reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
    """
    run_blastp(queryfaa, refdb, blastpout)
    besthits = parse_blastp(blastpout)

    seenids: List[str] = []
    outseqs: List[SeqIO.SeqRecord] = []

    for seq_record in SeqIO.parse(queryfaa, "fasta"):
        if seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

    for seq_record in SeqIO.parse(reffaa, "fasta"):
        if seq_record.description in besthits and seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

    SeqIO.write(outseqs, reduced_reffaa_merged, "fasta")


@click.command()
@click.option(
    "--queryfaa",
    "-q",
    type=click.Path(exists=True),
    required=True,
    help="Query sequences from hmmsearch for a single marker",
)
@click.option(
    "--reffaa",
    "-r",
    type=click.Path(exists=True),
    required=True,
    help="Reference sequences for a single marker",
)
@click.option(
    "--refdb",
    "-d",
    type=click.Path(exists=True),
    required=True,
    help="Diamond formatted database for a single marker",
)
@click.option(
    "--blastpout",
    "-b",
    type=click.Path(),
    required=True,
    help="Output file name for diamond blastp results",
)
@click.option(
    "--reduced_reffaa_merged",
    "-o",
    type=click.Path(),
    required=True,
    help="Output file name for reduced reference sequences",
)
def main(
    queryfaa: str, reffaa: str, refdb: str, blastpout: str, reduced_reffaa_merged: str
) -> None:
    """
    Main function to reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
    """
    if os.path.getsize(queryfaa) > 0:
        reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout)
    else:
        print(f"{queryfaa} is empty, no hits, omitting blast")


if __name__ == "__main__":
    main()
