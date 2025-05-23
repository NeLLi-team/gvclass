#!/usr/bin/env python
import os
import sys
import glob
import shutil
import logging
import pandas as pd
import subprocess
import argparse
from typing import List, Dict, Tuple
from Bio import SeqIO
import pyrodigal

def run_hmmsearch(faain: str, modelscombined: str, completeCutoff: float = 0.66) -> Tuple[int, float, int, float, float, int, float, float]:
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
    hmmsearch_cmd = f"hmmsearch --noali --domtblout {hmmout} --cpu 4 {modelscombined} {faain}"
    subprocess.run(hmmsearch_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    hit_count = 0
    total_score = 0
    unique_proteins = set()
    protein_scores: Dict[str, float] = {}
    protein_bestHit_scores: Dict[str, float] = {}
    unique_profiles = set()
    unique_protein2profiles = set()
    protein_bestHit_hmmCoverages: Dict[str, float] = {}
    protein_completeBestHit_scores: Dict[str, float] = {}

    with open(hmmout, 'r') as file:
        for line in file:
            if not line.startswith('#'):
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

                if protein_id not in protein_scores or score > protein_scores[protein_id]:
                    protein_scores[protein_id] = score

                unique_profiles.add(model)
                unique_protein2profiles.add(f'{protein_id}_{model}')

                if protein_id not in protein_bestHit_scores or hitScore > protein_bestHit_scores[protein_id]:
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
        avg_bestHitScore = round(sum(protein_bestHit_scores.values()) / len(protein_bestHit_scores.values()), 2)
        avg_bestHitCoverage = round(sum(protein_bestHit_hmmCoverages.values()) / len(protein_bestHit_hmmCoverages.values()), 3)
        if len(protein_completeBestHit_scores.values()) > 0:
            avg_completeBestHitScore = round(sum(protein_completeBestHit_scores.values()) / len(protein_completeBestHit_scores.values()), 2)
        else:
            avg_completeBestHitScore = 0
        avg_proteinsPerProfile = round(len(unique_protein2profiles) / len(unique_profiles), 2)
    else:
        avg_score = 0
        avg_bestHitScore = 0
        avg_bestHitCoverage = 0
        avg_completeBestHitScore = 0
        avg_proteinsPerProfile = 0

    return hit_count, avg_score, len(unique_profiles), avg_bestHitScore, avg_bestHitCoverage, len(protein_completeBestHit_scores.values()), avg_completeBestHitScore, avg_proteinsPerProfile

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
            cumulative_gc += seq_record.seq.upper().count('G') + seq_record.seq.upper().count('C')
            cumulative_len += len(seq_record.seq)
            contigs += 1

    with open(faain, "r") as faa_file:
        for seq_record in SeqIO.parse(faa_file, "fasta"):
            cumulative_len_aa += len(seq_record.seq)
            gene_count += 1

    return [contigs, cumulative_len, round(cumulative_gc / cumulative_len * 100, 2),
            gene_count, round(cumulative_len_aa * 3 / cumulative_len * 100, 2)]

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
                seq_record.id = f"{infilebase}|{str(seq_record.id).split()[0].replace(':', '__')}"
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
    return bytes("".join(sequences), 'utf-8')

def run_genecalling_codes(fnafile: str, code: int, gffout: str) -> None:
    """
    Run gene calling with pyrodigal and defined translation table.

    Args:
        fnafile (str): Path to the input FNA file.
        code (int): Translation table code.
        gffout (str): Path to the output GFF file.
    """
    faaout = gffout.replace(".gff", ".faa")

    metaFlag = code == 0
    outSuffix = "meta"

    # Map custom translation tables to standard ones
    # 106: TAA -> Q, ATG is the only start (as in code 6)
    # 129: TAG -> Y, ATG is the only start (as in code 29)
    code_mapping = {
        106: 6,
        129: 29
    }

    # Use the mapped code for pyrodigal
    pyrodigal_code = code_mapping.get(code, code)

    geneFinder = pyrodigal.GeneFinder(meta=metaFlag)

    if code > 0:
        trainingSeq = prepareTrainingSequence(fnafile)
        geneFinder.train(trainingSeq, force_nonsd=True, translation_table=pyrodigal_code)
        outSuffix = code

    with open(f'{faaout}_code{outSuffix}', "w") as proteinOut, open(f'{gffout}_code{outSuffix}', "w") as gffOut:
        for record in SeqIO.parse(fnafile, "fasta"):
            genes = geneFinder.find_genes(bytes(record.seq))
            genes.write_translations(proteinOut, record.id)
            genes.write_gff(gffOut, record.id)

def main():
    """
    Main function to run gene calling and calculate statistics.
    """
    parser = argparse.ArgumentParser(description='Run gene calling and calculate statistics')
    parser.add_argument('-f', '--fnafile', type=str, required=True, help='Input FNA file')
    parser.add_argument('-g', '--gffout', type=str, required=True, help='Output GFF file')
    parser.add_argument('-gs', '--genecalling_statsout', type=str, required=True, help='Genecalling statistics file')
    parser.add_argument('-ss', '--summary_statsout', type=str, required=True, help='Genome stats outfile')
    parser.add_argument('-fa', '--finalfaa', type=str, required=True, help='Output FASTA file')
    parser.add_argument('-m', '--modelscombined', type=str, required=True, help='Combined HMM model file')
    parser.add_argument('-fn', '--finalfna', type=str, required=True, help='Final FNA path')
    parser.add_argument('-o', '--hmmout', type=str, required=True, help='Final hmmout from bestcode')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Check input
    check_filename(args.fnafile)
    os.makedirs(os.path.dirname(args.finalfaa), exist_ok=True)
    os.makedirs(os.path.dirname(args.finalfna), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_statsout), exist_ok=True)
    os.makedirs(os.path.dirname(args.hmmout), exist_ok=True)

    # Gene calling
    # 106: TAA -> Q, ATG is the only start (as in code 6)
    # 129: TAG -> Y, ATG is the only start (as in code 29)
    codes = [0, 1, 4, 6, 106, 15, 11, 29, 129]  # Pass these from snakemake config instead
    for code in codes:
        run_genecalling_codes(args.fnafile, code, args.gffout)

    stats_dict: Dict[str, List[float]] = {}
    faaoutdir = os.path.dirname(args.gffout)
    hit_counts: Dict[str, int] = {}
    avg_scores: Dict[str, float] = {}
    profile_hit_counts: Dict[str, int] = {}
    avg_bestHitScores: Dict[str, float] = {}
    avg_bestHitCoverages: Dict[str, float] = {}
    completeBestHits_counts: Dict[str, int] = {}
    avg_completeBestHitScores: Dict[str, float] = {}
    avg_proteinsPerProfiles: Dict[str, float] = {}

    for faain in glob.glob(os.path.join(faaoutdir, "*.faa_code*")) + [os.path.join(faaoutdir, f"{os.path.splitext(os.path.basename(args.gffout))[0]}.faa_codemeta")]:
        genome_id = os.path.basename(faain)
        stats_dict[genome_id] = calc_stats(args.fnafile, faain)
        hit_count, avg_score, profile_hit_count, avg_bestHitScore, avg_bestHitCoverage, completeBestHits, avg_completeBestHitScore, avg_proteinsPerProfile = run_hmmsearch(faain, args.modelscombined)
        hit_counts[faain] = hit_count
        avg_scores[faain] = avg_score
        profile_hit_counts[faain] = profile_hit_count
        avg_bestHitScores[faain] = avg_bestHitScore
        avg_bestHitCoverages[faain] = avg_bestHitCoverage
        completeBestHits_counts[faain] = completeBestHits
        avg_completeBestHitScores[faain] = avg_completeBestHitScore
        avg_proteinsPerProfiles[faain] = avg_proteinsPerProfile

    # Gene calling stats to select ttable that yields highest coding density
    stats_df = pd.DataFrame.from_dict(stats_dict, columns=["contigs", "LENbp", "GCperc", "genecount", "CODINGperc"], orient="index")
    stats_df = stats_df.sort_values("CODINGperc", ascending=False)
    stats_df["ttable"] = stats_df.index.map(lambda x: x.split("_")[-1])
    stats_df.insert(0, "query", stats_df.index.map(lambda x: os.path.splitext(x)[0]))
    coding_dens_codemeta = stats_df[stats_df["ttable"] == "codemeta"]["CODINGperc"]

    # Add hit counts and average scores to the stats dataframe
    stats_df["hits"] = stats_df.index.map(lambda x: hit_counts[os.path.join(faaoutdir, x)])
    stats_df["avg_score"] = stats_df.index.map(lambda x: avg_scores[os.path.join(faaoutdir, x)])
    stats_df["profile_hits"] = stats_df.index.map(lambda x: profile_hit_counts[os.path.join(faaoutdir, x)])
    stats_df["avg_bestHitScore"] = stats_df.index.map(lambda x: avg_bestHitScores[os.path.join(faaoutdir, x)])
    stats_df["avg_bestHitCoverage"] = stats_df.index.map(lambda x: avg_bestHitCoverages[os.path.join(faaoutdir, x)])
    stats_df["complete_bestHits"] = stats_df.index.map(lambda x: completeBestHits_counts[os.path.join(faaoutdir, x)])
    stats_df["avg_completeBestHitScore"] = stats_df.index.map(lambda x: avg_completeBestHitScores[os.path.join(faaoutdir, x)])
    stats_df["avg_proteins_per_profile"] = stats_df.index.map(lambda x: avg_proteinsPerProfiles[os.path.join(faaoutdir, x)])

    # Select bestcode based on coding density threshold, number of hits, and average score
    if float(coding_dens_codemeta) > 0:
        bestcode = os.path.join(faaoutdir, f"{os.path.splitext(os.path.basename(args.gffout))[0]}.faa_codemeta")
        max_hits = completeBestHits_counts[bestcode]
        max_avg_score = avg_bestHitScores[bestcode]
        for code in ["code1", "code4", "code6", "code11", "code15", "code106", "code29", "code129"]:
            faa_code = os.path.join(faaoutdir, f"{os.path.splitext(os.path.basename(args.gffout))[0]}.faa_{code}")
            coding_dens_code = stats_df[stats_df["ttable"] == code]["CODINGperc"].iloc[0]
            if completeBestHits_counts[faa_code] > max_hits or (completeBestHits_counts[faa_code] == max_hits and avg_bestHitScores[faa_code] > max_avg_score) or (completeBestHits_counts[faa_code] == max_hits and avg_bestHitScores[faa_code] == max_avg_score and float(coding_dens_code) > float(coding_dens_codemeta) * 1.02):
                max_hits = completeBestHits_counts[faa_code]
                max_avg_score = avg_bestHitScores[faa_code]
                bestcode = faa_code
    else:
        bestcode = os.path.join(faaoutdir, stats_df.iloc[0].name)

    # Calculate delta and add it to the stats dataframe
    stats_df["delta"] = stats_df["CODINGperc"] - float(coding_dens_codemeta)
    stats_df["delta"] = stats_df["delta"].round(2)

    # Write stats to files
    stats_df.to_csv(args.genecalling_statsout, sep="\t", index=False)
    stats_df[stats_df.index == os.path.basename(bestcode)].to_csv(args.summary_statsout, sep="\t", index=False)

    # Rename headers and cleanup
    shutil.copy(args.fnafile, args.finalfna)
    rename_header(bestcode, args.finalfaa, args.gffout)
    for tempout in glob.glob(f"{args.gffout}_code*"):
        os.remove(tempout)
        faaout = tempout.replace(".gff", ".faa")
        if os.path.isfile(faaout) and faaout != args.finalfaa:
            os.remove(faaout)

    # Remove codemeta GFF and FAA files if they exist
    codemeta_gff = f"{args.gffout}_codemeta"
    if os.path.isfile(codemeta_gff) and codemeta_gff != args.gffout:
        os.remove(codemeta_gff)
    codemeta_faa = codemeta_gff.replace(".gff", ".faa")
    if os.path.isfile(codemeta_faa) and codemeta_faa != args.finalfaa:
        os.remove(codemeta_faa)

    # Copy the hmmsearch output file from the bestcode to the specified location
    bestcode_hmmout = bestcode.replace(".faa", ".hmmout")
    if os.path.isfile(bestcode_hmmout):
        shutil.copy(bestcode_hmmout, args.hmmout)

    # Remove remaining hmmout files
    for hmmout_file in glob.glob(f"{os.path.splitext(args.gffout)[0]}.hmmout_*"):
        if hmmout_file != args.hmmout:
            os.remove(hmmout_file)

if __name__ == '__main__':
    main()
