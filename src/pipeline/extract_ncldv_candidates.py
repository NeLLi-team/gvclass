#!/usr/bin/env python3
"""
Extract NCLDV candidates from GVClass results.

This script:
1. Reads gvclass_summary.tsv and filters for NCLDV:>=2
2. Extracts corresponding FAA files to ncldv_candidates_faa/
3. Creates filtered summary file ncldv_candidates_gvclass.tsv
4. Organizes output files (moves gz, removes summary.tab)
"""

import argparse
import pandas as pd
import shutil
from pathlib import Path
import re
import sys


def parse_ncldv_count(cell_value):
    """Extract NCLDV count from cell value like 'NCLDV:2'."""
    if pd.isna(cell_value) or not isinstance(cell_value, str):
        return 0

    match = re.search(r"NCLDV:(\d+)", cell_value)
    if match:
        return int(match.group(1))
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Extract NCLDV candidates from GVClass results"
    )
    parser.add_argument(
        "results_dir", help="GVClass results directory (e.g., allbins_results)"
    )
    parser.add_argument(
        "query_dir", help="Original query directory with FAA files (e.g., allbins)"
    )
    parser.add_argument(
        "--min-ncldv",
        type=int,
        default=2,
        help="Minimum NCLDV count to be considered a candidate (default: 2)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )

    args = parser.parse_args()

    results_path = Path(args.results_dir)
    query_path = Path(args.query_dir)

    # Check input paths
    if not results_path.exists():
        print(f"Error: Results directory not found: {results_path}")
        sys.exit(1)

    if not query_path.exists():
        print(f"Error: Query directory not found: {query_path}")
        sys.exit(1)

    summary_file = results_path / "gvclass_summary.tsv"
    if not summary_file.exists():
        print(f"Error: Summary file not found: {summary_file}")
        sys.exit(1)

    # Read summary file
    print(f"Reading {summary_file}...")
    df = pd.read_csv(summary_file, sep="\t")

    # Find the column containing NCLDV counts
    # Based on the data, NCLDV info is in the 'domain' column
    ncldv_col = None
    for col in ["domain", "Domain", "DOMAIN"]:
        if col in df.columns:
            ncldv_col = col
            break

    # If not found, search for any column containing NCLDV data
    if not ncldv_col:
        for col in df.columns:
            sample_values = df[col].dropna().head(20).astype(str)
            if any("NCLDV:" in val for val in sample_values):
                ncldv_col = col
                break

    if not ncldv_col:
        print("Error: Could not find column with NCLDV counts")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    print(f"Found NCLDV data in column: {ncldv_col}")

    # Extract NCLDV counts and filter
    df["ncldv_count"] = df[ncldv_col].apply(parse_ncldv_count)
    ncldv_candidates = df[df["ncldv_count"] >= args.min_ncldv].copy()

    print(
        f"\nFound {len(ncldv_candidates)} NCLDV candidates with count >= {args.min_ncldv}"
    )

    if len(ncldv_candidates) == 0:
        print("No candidates found. Exiting.")
        sys.exit(0)

    # Create output directories
    ncldv_faa_dir = results_path / "ncldv_candidates_faa"
    gvclass_dir = results_path / "gvclass"

    if not args.dry_run:
        ncldv_faa_dir.mkdir(exist_ok=True)
        gvclass_dir.mkdir(exist_ok=True)

    # Extract query names (assuming first column is query name)
    query_col = df.columns[0]
    candidate_queries = ncldv_candidates[query_col].tolist()

    print("\nNCLDV candidates:")
    for idx, row in ncldv_candidates.iterrows():
        print(f"  {row[query_col]} - NCLDV:{row['ncldv_count']}")

    # Step 1: Copy FAA files for NCLDV candidates
    print(f"\nCopying FAA files to {ncldv_faa_dir}...")
    copied = 0
    missing = []

    for query in candidate_queries:
        # Try different possible extensions
        source_files = (
            list(query_path.glob(f"{query}.faa"))
            + list(query_path.glob(f"{query}.fasta"))
            + list(query_path.glob(f"{query}.fa"))
        )

        if source_files:
            source_file = source_files[0]  # Take first match
            dest_file = ncldv_faa_dir / f"{query}.faa"

            if not args.dry_run:
                shutil.copy2(source_file, dest_file)
            print(f"  Copied: {source_file.name}")
            copied += 1
        else:
            missing.append(query)

    print(f"Copied {copied} FAA files")
    if missing:
        print(f"Warning: Could not find FAA files for {len(missing)} queries:")
        for m in missing[:10]:  # Show first 10
            print(f"  - {m}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")

    # Step 2: Write filtered summary
    output_summary = results_path / "ncldv_candidates_gvclass.tsv"
    print(f"\nWriting filtered summary to {output_summary}...")

    if not args.dry_run:
        # Drop the temporary ncldv_count column before saving
        ncldv_candidates_clean = ncldv_candidates.drop("ncldv_count", axis=1)
        ncldv_candidates_clean.to_csv(output_summary, sep="\t", index=False)

    print(f"Wrote {len(ncldv_candidates)} entries to summary file")

    # Step 3: Organize files
    print("\nOrganizing output files...")

    # Move all .tar.gz files to gvclass/
    gz_files = list(results_path.glob("*.tar.gz"))
    print(f"Moving {len(gz_files)} .tar.gz files to {gvclass_dir}...")

    for gz_file in gz_files:
        dest = gvclass_dir / gz_file.name
        if not args.dry_run:
            shutil.move(str(gz_file), str(dest))
        if not args.dry_run or len(gz_files) <= 5:
            print(f"  Moved: {gz_file.name}")

    if len(gz_files) > 5 and args.dry_run:
        print(f"  ... and {len(gz_files) - 5} more files")

    # Remove all .summary.tab files
    summary_tab_files = list(results_path.glob("*.summary.tab"))
    print(f"\nRemoving {len(summary_tab_files)} .summary.tab files...")

    if not args.dry_run:
        for tab_file in summary_tab_files:
            tab_file.unlink()

    if len(summary_tab_files) <= 5:
        for tab_file in summary_tab_files:
            print(f"  Removed: {tab_file.name}")
    else:
        print(f"  Removing {len(summary_tab_files)} files...")

    print("\nDone!")
    print("\nSummary:")
    print(f"  - NCLDV candidates found: {len(ncldv_candidates)}")
    print(f"  - FAA files copied: {copied}")
    print(f"  - Summary written to: {output_summary}")
    print(f"  - Moved {len(gz_files)} .tar.gz files to {gvclass_dir}")
    print(f"  - Removed {len(summary_tab_files)} .summary.tab files")


if __name__ == "__main__":
    main()
