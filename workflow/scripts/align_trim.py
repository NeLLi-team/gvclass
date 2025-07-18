
import os
import sys
import click
import pytrimal
import pyfamsa
from Bio import SeqIO

# Add the utils directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import validate_file_path, GVClassError
from error_handling import error_handler, ProcessingError


@click.command()
@click.option("--queryfaamerged", "-q", type=click.Path(exists=True), help="Input fasta file containing merged fasta sequences")
@click.option("--alnout", "-a", type=click.Path(), help="Output fasta file for aligned sequences")
@click.option("--trimout", "-t", type=click.Path(), help="Output fasta file for trimmed sequences")
@click.option("--mafftoption", "-o", default="linsi", help="MAFFT alignment option (e.g., 'auto', 'linsi', 'fftns')")
@click.option("--threads", "-j", default=4, help="Number of threads to use")
def main(queryfaamerged, alnout, trimout, mafftoption, threads):
    """
    Align and trim sequences using MAFFT and trimAl.
    
    Args:
        queryfaamerged: Input fasta file with merged sequences
        alnout: Output file for aligned sequences
        trimout: Output file for trimmed sequences
        mafftoption: MAFFT alignment option
        threads: Number of threads to use
    """
    try:
        # Validate input files
        input_file = validate_file_path(queryfaamerged, must_exist=True)
        output_aln = validate_file_path(alnout, must_exist=False)
        output_trim = validate_file_path(trimout, must_exist=False)
        
        # Validate threads parameter
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("Threads must be a positive integer")
        
        # Count the number of sequences in the input file
        with open(input_file, "r") as f:
            sequence_count = sum(1 for line in f if line.startswith(">"))

        # Proceed only if there are at least 3 sequences
        if sequence_count >= 3:
            # Validate mafftoption for pyfamsa
            valid_options = ["linsi", "ginsi", "auto"]
            if mafftoption not in valid_options:
                raise ValueError(f"Invalid alignment option: {mafftoption}. Must be one of: {valid_options}")
            
            # Use pyfamsa for alignment
            try:
                # Load sequences from file
                sequences = []
                for record in SeqIO.parse(input_file, "fasta"):
                    seq = pyfamsa.Sequence(record.id.encode(), str(record.seq).encode())
                    sequences.append(seq)
                
                # Configure aligner based on option
                if mafftoption == "linsi":
                    # Use more accurate alignment for linsi
                    aligner = pyfamsa.Aligner(threads=threads, guide_tree="upgma")
                elif mafftoption == "ginsi":
                    # Use global alignment for ginsi
                    aligner = pyfamsa.Aligner(threads=threads, guide_tree="upgma")
                else:
                    # Use default/auto alignment
                    aligner = pyfamsa.Aligner(threads=threads)
                
                # Run alignment
                alignment = aligner.align(sequences)
                
                # Save alignment to file
                with open(output_aln, 'w') as f:
                    for seq in alignment:
                        f.write(f">{seq.id.decode()}\n")
                        f.write(f"{seq.sequence.decode()}\n")
                
                error_handler.log_info(f"Alignment created successfully: {output_aln}")
                
            except Exception as e:
                raise ProcessingError(
                    f"pyfamsa alignment failed: {str(e)}",
                    step="sequence_alignment",
                    input_file=str(input_file)
                ) from e

            # Trimming of the alignment using pytrimal
            try:
                # Load alignment from file
                alignment = pytrimal.Alignment.load(str(output_aln))
                
                # Create trimmer with gap threshold 0.1 (keep columns with ≤10% gaps)
                trimmer = pytrimal.ManualTrimmer(gap_threshold=0.1)
                
                # Trim the alignment
                trimmed_alignment = trimmer.trim(alignment)
                
                # Save trimmed alignment
                trimmed_alignment.dump(str(output_trim), format="fasta")
                
                error_handler.log_info(f"Alignment trimmed successfully: {output_trim}")
                
            except Exception as e:
                raise ProcessingError(
                    f"pytrimal trimming failed: {str(e)}",
                    step="trim_alignment",
                    input_file=str(output_aln)
                ) from e
            
        else:
            print("Alignment not created as there is only one sequence in the input file.")
            
    except (ValueError, FileNotFoundError, GVClassError) as e:
        print(f"Error: {e}")
        raise click.ClickException(str(e))

if __name__ == "__main__":
    main()