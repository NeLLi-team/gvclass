import os
import sys
import shutil
import logging
import pandas as pd
from typing import List, Dict, Tuple
from Bio import SeqIO
import pyrodigal  # Using tomasbruna's fork with codes 106, 129
import pyhmmer
import click

from src.utils.common import validate_file_path
from src.utils.error_handling import ProcessingError


def run_pyhmmer_search_simple(hmm_file: str, query_file: str, output_file: str) -> None:
    """
    Simplified pyhmmer search function for opgecall.py.

    Args:
        hmm_file: Path to HMM model file
        query_file: Path to query sequence file
        output_file: Path to output file
    """
    try:
        # Read HMM profiles
        with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
            profiles = list(hmm_handle)

        # Read query sequences
        with pyhmmer.easel.SequenceFile(query_file, digital=True) as seq_handle:
            sequences = list(seq_handle)

        # Run search with default parameters
        results = pyhmmer.hmmsearch(profiles, sequences, cpus=4)

        # Write results in simplified domtblout format
        with open(output_file, "w") as out_handle:
            # Write header
            out_handle.write(
                "# target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\tbias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\tfrom\tto\tfrom\tto\tfrom\tto\tacc\n"
            )

            # Iterate through results - each TopHits corresponds to one query HMM
            for i, (hmm, top_hits) in enumerate(zip(profiles, results)):
                hmm_name = (
                    hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                )
                hmm_accession = (
                    hmm.accession.decode()
                    if hmm.accession and isinstance(hmm.accession, bytes)
                    else (hmm.accession or "-")
                )
                hmm_length = hmm.M

                # Process hits for this HMM
                for hit in top_hits:
                    target_name = (
                        hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                    )
                    target_accession = "-"  # Simplified
                    target_length = (
                        1000  # Placeholder - will be parsed differently downstream
                    )

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
                        out_handle.write(
                            f"{target_name}\t{target_accession}\t{target_length}\t"
                        )
                        out_handle.write(f"{hmm_name}\t{hmm_accession}\t{hmm_length}\t")
                        out_handle.write(
                            f"{hit.evalue:.2e}\t{hit.score:.1f}\t{hit.bias:.1f}\t"
                        )
                        out_handle.write(f"{domain_idx+1}\t{len(hit.domains)}\t")
                        out_handle.write(
                            f"{domain.c_evalue:.2e}\t{domain.i_evalue:.2e}\t"
                        )
                        out_handle.write(f"{domain.score:.1f}\t{domain.bias:.1f}\t")
                        # Convert 0-based coordinates to 1-based for output
                        # Note: pyhmmer uses 0-based inclusive coordinates
                        hmm_start_1based = hmm_from + 1
                        hmm_end_1based = hmm_to + 1  # Fix: was missing +1
                        target_start_1based = target_from + 1
                        target_end_1based = target_to + 1  # Fix: was missing +1
                        env_start_1based = domain.env_from + 1
                        env_end_1based = domain.env_to + 1  # Fix: was missing +1

                        out_handle.write(f"{hmm_start_1based}\t{hmm_end_1based}\t")
                        out_handle.write(
                            f"{target_start_1based}\t{target_end_1based}\t"
                        )
                        out_handle.write(f"{env_start_1based}\t{env_end_1based}\t")
                        out_handle.write("0.90\n")  # Default accuracy

    except Exception as e:
        raise ProcessingError(
            f"pyhmmer search failed: {str(e)}",
            step="pyhmmer_search",
            input_file=query_file,
        ) from e


def run_hmmsearch(
    faain: str, modelscombined: str, completeCutoff: float = 0.66
) -> Tuple[int, float, int, float, float, int, float, float]:
    """
    Run hmmsearch on the input FAA file and calculate relevant metrics.

    Args:
        faain (str): Path to the input FAA file.
        modelscombined (str): Path to the combined HMM models file.
        completeCutoff (float): Coverage cutoff for considering complete hits.

    Returns:
        Tuple: Metrics calculated from the hmmsearch results.
    """
    hmmout = faain.replace(".faa", ".hmmout")

    # Validate input files
    faa_file = validate_file_path(faain, must_exist=True)
    models_file = validate_file_path(modelscombined, must_exist=True)

    # Run pyhmmer search
    run_pyhmmer_search_simple(str(models_file), str(faa_file), hmmout)

    hit_count = 0
    total_score = 0
    unique_proteins = set()
    protein_scores: Dict[str, float] = {}
    protein_bestHit_scores: Dict[str, float] = {}
    unique_profiles = set()
    unique_protein2profiles = set()
    protein_bestHit_hmmCoverages: Dict[str, float] = {}
    protein_completeBestHit_scores: Dict[str, float] = {}

    with open(hmmout, "r") as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.split()
                protein_id = parts[0]
                model = parts[3]
                score = float(parts[7])
                hitScore = float(parts[13])
                profileLen = int(parts[5])
                profileFrom = int(parts[15])
                profileTo = int(parts[16])

                if protein_id not in unique_proteins:
                    unique_proteins.add(protein_id)
                    hit_count += 1

                if (
                    protein_id not in protein_scores
                    or score > protein_scores[protein_id]
                ):
                    protein_scores[protein_id] = score

                unique_profiles.add(model)
                unique_protein2profiles.add(f"{protein_id}_{model}")

                if (
                    protein_id not in protein_bestHit_scores
                    or hitScore > protein_bestHit_scores[protein_id]
                ):
                    protein_bestHit_scores[protein_id] = hitScore
                    coverage = (profileTo - profileFrom + 1) / profileLen
                    protein_bestHit_hmmCoverages[protein_id] = coverage
                    if coverage > completeCutoff:
                        protein_completeBestHit_scores[protein_id] = hitScore
                    else:
                        protein_completeBestHit_scores.pop(protein_id, None)

    if hit_count > 0:
        total_score = sum(protein_scores.values())
        avg_score = round(total_score / hit_count, 2)
        avg_bestHitScore = round(
            sum(protein_bestHit_scores.values()) / len(protein_bestHit_scores.values()),
            2,
        )
        avg_bestHitCoverage = round(
            sum(protein_bestHit_hmmCoverages.values())
            / len(protein_bestHit_hmmCoverages.values()),
            3,
        )
        if len(protein_completeBestHit_scores.values()) > 0:
            avg_completeBestHitScore = round(
                sum(protein_completeBestHit_scores.values())
                / len(protein_completeBestHit_scores.values()),
                2,
            )
        else:
            avg_completeBestHitScore = 0
        avg_proteinsPerProfile = round(
            len(unique_protein2profiles) / len(unique_profiles), 2
        )
    else:
        avg_score = 0
        avg_bestHitScore = 0
        avg_bestHitCoverage = 0
        avg_completeBestHitScore = 0
        avg_proteinsPerProfile = 0

    return (
        hit_count,
        avg_score,
        len(unique_profiles),
        avg_bestHitScore,
        avg_bestHitCoverage,
        len(protein_completeBestHit_scores.values()),
        avg_completeBestHitScore,
        avg_proteinsPerProfile,
    )


def calc_stats(fnain: str, faain: str) -> List[float]:
    """
    Calculate coding density, assembly size, and GC%.

    Args:
        fnain (str): Path to the input FNA file.
        faain (str): Path to the input FAA file.

    Returns:
        List[float]: List containing calculated statistics.
    """
    cumulative_gc = 0
    cumulative_len = 0
    cumulative_len_aa = 0
    gene_count = 0
    contigs = 0

    with open(fnain, "r") as fna_file:
        for seq_record in SeqIO.parse(fna_file, "fasta"):
            cumulative_gc += seq_record.seq.upper().count(
                "G"
            ) + seq_record.seq.upper().count("C")
            cumulative_len += len(seq_record.seq)
            contigs += 1

    with open(faain, "r") as faa_file:
        for seq_record in SeqIO.parse(faa_file, "fasta"):
            cumulative_len_aa += len(seq_record.seq)
            gene_count += 1

    return [
        contigs,
        cumulative_len,
        round(cumulative_gc / cumulative_len * 100, 2),
        gene_count,
        round(cumulative_len_aa * 3 / cumulative_len * 100, 2),
    ]


def check_filename(fnafile: str) -> None:
    """
    Check if the input filename has the correct format (.fna with no additional '.' in the filename).

    Args:
        fnafile (str): Input FNA filename.

    Raises:
        SystemExit: If the filename format is incorrect.
    """
    if len(fnafile.split(".")) > 2:
        logging.error("fnafile name contains additional '.'")
        sys.exit(1)
    if fnafile.split(".")[1] != "fna":
        logging.error("fnafile ending needs to be .fna")
        sys.exit(1)


def rename_header(bestcode: str, finalfaa: str, gffout: str) -> None:
    """
    Rename sequence header to ><filename>|<proteinid>.

    Args:
        bestcode (str): Path to the best code file.
        finalfaa (str): Path to the output FASTA file.
        gffout (str): Path to the output GFF file.
    """
    records = []
    infilebase = os.path.splitext(os.path.basename(bestcode))[0]

    with open(bestcode, "r") as best_file:
        for seq_record in SeqIO.parse(best_file, "fasta"):
            headerbase = seq_record.description.split("|")[0]
            if headerbase != infilebase:
                seq_record.id = (
                    f"{infilebase}|{str(seq_record.id).split()[0].replace(':', '__')}"
                )
                seq_record.description = ""
                records.append(seq_record)
            elif headerbase == infilebase:
                seq_record.id = str(seq_record.id).split()[0].replace(":", "__")
                seq_record.description = ""
                records.append(seq_record)

    with open(finalfaa, "w") as faa_out:
        SeqIO.write(records, faa_out, "fasta")

    # Keep only best gff
    gff_best = f"{gffout}_code{bestcode.split('_code')[-1]}"
    if os.path.isfile(gff_best):
        shutil.copy(gff_best, gffout)
    else:
        shutil.copy(f"{gffout}_codemeta", gffout)


def prepareTrainingSequence(fnafile: str) -> bytes:
    """
    Concatenate sequences in a multifasta file into one string that can be used for pyrodigal training.
    Adapted from https://github.com/althonos/pyrodigal/issues/14

    Args:
        fnafile (str): Path to the input FNA file.

    Returns:
        bytes: Concatenated sequences for pyrodigal training.
    """
    sequences = []
    for i, record in enumerate(SeqIO.parse(fnafile, "fasta")):
        if i > 0:
            # This string has all stops in all frames to prevent genes spanning records.
            sequences.append("TTAATTAATTAA")
        sequences.append(str(record.seq))
    if len(sequences) > 1:
        sequences.append("TTAATTAATTAA")
    return bytes("".join(sequences), "utf-8")


def run_genecalling_codes(fnafile: str, code: int, gffout: str) -> None:
    """
    Run gene calling with pyrodigal and defined translation table.
    tomasbruna's fork natively supports codes 106 and 129.

    Args:
        fnafile (str): Path to the input FNA file.
        code (int): Translation table code.
        gffout (str): Path to the output GFF file.
    """
    faaout = gffout.replace(".gff", ".faa")
    outSuffix = "meta" if code == 0 else code

    # tomasbruna's fork supports codes 106 and 129 natively
    metaFlag = code == 0
    geneFinder = pyrodigal.GeneFinder(meta=metaFlag)

    if code > 0:
        trainingSeq = prepareTrainingSequence(fnafile)
        geneFinder.train(trainingSeq, force_nonsd=True, translation_table=code)

    with (
        open(f"{faaout}_code{outSuffix}", "w") as proteinOut,
        open(f"{gffout}_code{outSuffix}", "w") as gffOut,
    ):
        for record in SeqIO.parse(fnafile, "fasta"):
            genes = geneFinder.find_genes(bytes(record.seq))
            genes.write_translations(proteinOut, record.id)
            genes.write_gff(gffOut, record.id)


@click.command()
@click.option(
    "--fnafile",
    "-f",
    type=click.Path(exists=True),
    required=True,
    help="Input FNA file",
)
@click.option(
    "--gffout", "-g", type=click.Path(), required=True, help="Output GFF file"
)
@click.option(
    "--genecalling_statsout",
    "-gs",
    "--genecalling-statsout",
    type=click.Path(),
    required=True,
    help="Genecalling statistics file",
)
@click.option(
    "--summary_statsout",
    "-ss",
    type=click.Path(),
    required=True,
    help="Genome stats outfile",
)
@click.option(
    "--finalfaa", "-fa", type=click.Path(), required=True, help="Output FASTA file"
)
@click.option(
    "--modelscombined",
    "-m",
    type=click.Path(exists=True),
    required=True,
    help="Combined HMM model file",
)
@click.option(
    "--finalfna", "-fn", type=click.Path(), required=True, help="Final FNA path"
)
@click.option(
    "--hmmout",
    "-o",
    type=click.Path(),
    required=True,
    help="Final hmmout from bestcode",
)
def main(
    fnafile: str,
    gffout: str,
    genecalling_statsout: str,
    summary_statsout: str,
    finalfaa: str,
    modelscombined: str,
    finalfna: str,
    hmmout: str,
) -> None:
    """
    Main function to run gene calling and calculate statistics.

    Args:
        fnafile (str): Path to the input FNA file.
        gffout (str): Path to the output GFF file.
        genecalling_statsout (str): Path to the genecalling statistics output file.
        summary_statsout (str): Path to the genome stats output file.
        finalfaa (str): Path to the output FASTA file.
        modelscombined (str): Path to the combined HMM model file.
        finalfna (str): Final FNA path.
        hmmout (str): Final hmmout from bestcode.
    """
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Check input
    check_filename(fnafile)
    os.makedirs(os.path.dirname(finalfaa), exist_ok=True)
    os.makedirs(os.path.dirname(finalfna), exist_ok=True)
    os.makedirs(os.path.dirname(summary_statsout), exist_ok=True)
    os.makedirs(os.path.dirname(hmmout), exist_ok=True)

    # Gene calling
    # tomasbruna's fork supports codes 106 and 129 for giant viruses
    codes = [0, 1, 4, 6, 11, 15, 29, 106, 129]
    for code in codes:
        run_genecalling_codes(fnafile, code, gffout)

    stats_dict: Dict[str, List[float]] = {}
    faaoutdir = os.path.dirname(gffout)
    hit_counts: Dict[str, int] = {}
    avg_scores: Dict[str, float] = {}
    profile_hit_counts: Dict[str, int] = {}
    avg_bestHitScores: Dict[str, float] = {}
    avg_bestHitCoverages: Dict[str, float] = {}
    completeBestHits_counts: Dict[str, int] = {}
    avg_completeBestHitScores: Dict[str, float] = {}
    avg_proteinsPerProfiles: Dict[str, float] = {}

    # Secure glob usage - validate the directory first
    faaoutdir_safe = validate_file_path(faaoutdir, must_exist=True)

    # Build search patterns safely
    meta_file = (
        faaoutdir_safe / f"{os.path.splitext(os.path.basename(gffout))[0]}.faa_codemeta"
    )

    faa_files = list(faaoutdir_safe.glob("*.faa_code*"))
    if meta_file.exists():
        faa_files.append(meta_file)

    for faain in faa_files:
        faain_str = str(faain)
        genome_id = os.path.basename(faain_str)
        stats_dict[genome_id] = calc_stats(fnafile, faain_str)
        (
            hit_count,
            avg_score,
            profile_hit_count,
            avg_bestHitScore,
            avg_bestHitCoverage,
            completeBestHits,
            avg_completeBestHitScore,
            avg_proteinsPerProfile,
        ) = run_hmmsearch(faain_str, modelscombined)
        hit_counts[faain_str] = hit_count
        avg_scores[faain_str] = avg_score
        profile_hit_counts[faain_str] = profile_hit_count
        avg_bestHitScores[faain_str] = avg_bestHitScore
        avg_bestHitCoverages[faain_str] = avg_bestHitCoverage
        completeBestHits_counts[faain_str] = completeBestHits
        avg_completeBestHitScores[faain_str] = avg_completeBestHitScore
        avg_proteinsPerProfiles[faain_str] = avg_proteinsPerProfile

    # Gene calling stats to select ttable that yields highest coding density
    stats_df = pd.DataFrame.from_dict(
        stats_dict,
        columns=["contigs", "LENbp", "GCperc", "genecount", "CODINGperc"],
        orient="index",
    )
    stats_df = stats_df.sort_values("CODINGperc", ascending=False)
    # Extract genetic code from filename: e.g., "file.faa_code4" -> "code4"
    stats_df["ttable"] = stats_df.index.map(lambda x: x.split("_")[-1])
    stats_df.insert(0, "query", stats_df.index.map(lambda x: os.path.splitext(x)[0]))
    coding_dens_series = stats_df.loc[stats_df["ttable"] == "codemeta", "CODINGperc"]
    coding_dens_codemeta = (
        float(coding_dens_series.iloc[0]) if not coding_dens_series.empty else 0.0
    )

    # Add hit counts and average scores to the stats dataframe
    stats_df["hits"] = stats_df.index.map(
        lambda x: hit_counts[os.path.join(faaoutdir, x)]
    )
    stats_df["avg_score"] = stats_df.index.map(
        lambda x: avg_scores[os.path.join(faaoutdir, x)]
    )
    stats_df["profile_hits"] = stats_df.index.map(
        lambda x: profile_hit_counts[os.path.join(faaoutdir, x)]
    )
    stats_df["avg_bestHitScore"] = stats_df.index.map(
        lambda x: avg_bestHitScores[os.path.join(faaoutdir, x)]
    )
    stats_df["avg_bestHitCoverage"] = stats_df.index.map(
        lambda x: avg_bestHitCoverages[os.path.join(faaoutdir, x)]
    )
    stats_df["complete_bestHits"] = stats_df.index.map(
        lambda x: completeBestHits_counts[os.path.join(faaoutdir, x)]
    )
    stats_df["avg_completeBestHitScore"] = stats_df.index.map(
        lambda x: avg_completeBestHitScores[os.path.join(faaoutdir, x)]
    )
    stats_df["avg_proteins_per_profile"] = stats_df.index.map(
        lambda x: avg_proteinsPerProfiles[os.path.join(faaoutdir, x)]
    )

    # Select bestcode based on coding density threshold, number of hits, and average score
    if coding_dens_codemeta > 0:
        bestcode = os.path.join(
            faaoutdir, f"{os.path.splitext(os.path.basename(gffout))[0]}.faa_codemeta"
        )
        max_hits = completeBestHits_counts[bestcode]
        max_avg_score = avg_bestHitScores[bestcode]
        for code in [
            "code1",
            "code4",
            "code6",
            "code11",
            "code15",
            "code29",
            "code106",
            "code129",
        ]:
            faa_code = os.path.join(
                faaoutdir, f"{os.path.splitext(os.path.basename(gffout))[0]}.faa_{code}"
            )

            # Skip if this code wasn't run (file doesn't exist or not in our dictionaries)
            if faa_code not in completeBestHits_counts:
                continue

            matching_rows = stats_df[stats_df["ttable"] == code]
            if len(matching_rows) == 0:
                continue  # Skip if this code wasn't run
            coding_dens_code = matching_rows["CODINGperc"].iloc[0]

            if (
                completeBestHits_counts[faa_code] > max_hits
                or (
                    completeBestHits_counts[faa_code] == max_hits
                    and avg_bestHitScores[faa_code] > max_avg_score
                )
                or (
                    completeBestHits_counts[faa_code] == max_hits
                    and avg_bestHitScores[faa_code] == max_avg_score
                    and float(coding_dens_code) > coding_dens_codemeta * 1.02
                )
            ):
                max_hits = completeBestHits_counts[faa_code]
                max_avg_score = avg_bestHitScores[faa_code]
                bestcode = faa_code
    else:
        bestcode = os.path.join(faaoutdir, stats_df.iloc[0].name)

    # Calculate delta and add it to the stats dataframe
    stats_df["delta"] = stats_df["CODINGperc"] - coding_dens_codemeta
    stats_df["delta"] = stats_df["delta"].round(2)

    # Write stats to files
    stats_df.to_csv(genecalling_statsout, sep="\t", index=False)

    # Write the best code stats, but always show "codemeta" in ttable column
    best_stats = stats_df[stats_df.index == os.path.basename(bestcode)].copy()
    best_stats["ttable"] = "codemeta"  # Always show codemeta for the selected code
    best_stats.to_csv(summary_statsout, sep="\t", index=False)

    # Rename headers and cleanup
    shutil.copy(fnafile, finalfna)
    rename_header(bestcode, finalfaa, gffout)
    # Secure cleanup - validate paths before removal
    gffout_safe = validate_file_path(gffout, must_exist=False)
    gffout_dir = gffout_safe.parent
    gffout_base = gffout_safe.stem

    for tempout in gffout_dir.glob(f"{gffout_base}_code*"):
        tempout_safe = validate_file_path(tempout, must_exist=False)
        if tempout_safe.exists():
            os.remove(tempout_safe)
        faaout = str(tempout_safe).replace(".gff", ".faa")
        if os.path.isfile(faaout) and faaout != finalfaa:
            os.remove(faaout)

    # Remove codemeta GFF and FAA files if they exist
    codemeta_gff = f"{gffout}_codemeta"
    if os.path.isfile(codemeta_gff) and codemeta_gff != gffout:
        os.remove(codemeta_gff)
    codemeta_faa = codemeta_gff.replace(".gff", ".faa")
    if os.path.isfile(codemeta_faa) and codemeta_faa != finalfaa:
        os.remove(codemeta_faa)

    # Copy the hmmsearch output file from the bestcode to the specified location
    bestcode_hmmout = bestcode.replace(".faa", ".hmmout")
    if os.path.isfile(bestcode_hmmout):
        shutil.copy(bestcode_hmmout, hmmout)

    # Remove remaining hmmout files - secure cleanup
    gffout_base = os.path.splitext(str(gffout_safe))[0]
    gffout_dir = gffout_safe.parent
    base_name = os.path.basename(gffout_base)

    for hmmout_file in gffout_dir.glob(f"{base_name}.hmmout_*"):
        hmmout_file_safe = validate_file_path(hmmout_file, must_exist=False)
        if str(hmmout_file_safe) != hmmout and hmmout_file_safe.exists():
            os.remove(hmmout_file_safe)


if __name__ == "__main__":
    main()
