#!/usr/bin/env python3
"""
Merge results from task array chunks into final summary files.
"""
import os
import sys
import glob
import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('--outdir', required=True, help='Output directory containing chunk results')
@click.option('--pattern', default='*_chunk*.txt', help='Pattern to match chunk result files')
@click.option('--summary-prefix', default='summary', help='Prefix for merged summary files')
def merge_results(outdir, pattern, summary_prefix):
    """Merge chunk results from GVClass array job."""
    outdir = Path(outdir)
    
    if not outdir.exists():
        click.echo(f"Error: Output directory {outdir} does not exist", err=True)
        sys.exit(1)
    
    # Find all chunk result files
    chunk_files = sorted(glob.glob(str(outdir / pattern)))
    
    if not chunk_files:
        click.echo(f"No chunk files found matching pattern '{pattern}' in {outdir}", err=True)
        sys.exit(1)
    
    click.echo(f"Found {len(chunk_files)} chunk files to merge")
    
    # Group files by type (removing chunk suffix)
    file_groups = {}
    for chunk_file in chunk_files:
        base_name = os.path.basename(chunk_file)
        # Remove _chunkN suffix
        file_type = base_name.split('_chunk')[0]
        if file_type not in file_groups:
            file_groups[file_type] = []
        file_groups[file_type].append(chunk_file)
    
    # Merge each group
    for file_type, files in file_groups.items():
        click.echo(f"\nMerging {len(files)} files of type: {file_type}")
        
        # Read and concatenate all files
        all_data = []
        for f in sorted(files):
            try:
                # Try to read as TSV/CSV
                df = pd.read_csv(f, sep='\t', header=0)
                all_data.append(df)
            except:
                # If not a table, concatenate as text
                with open(f, 'r') as infile:
                    all_data.append(infile.read())
        
        # Write merged file
        output_file = outdir / f"{summary_prefix}_{file_type}.txt"
        
        if all_data and isinstance(all_data[0], pd.DataFrame):
            # Merge dataframes
            merged_df = pd.concat(all_data, ignore_index=True)
            # Remove duplicates if any (based on first column)
            if not merged_df.empty:
                merged_df = merged_df.drop_duplicates(subset=merged_df.columns[0])
            merged_df.to_csv(output_file, sep='\t', index=False)
            click.echo(f"  Merged {len(all_data)} dataframes -> {output_file} ({len(merged_df)} rows)")
        else:
            # Merge text files
            with open(output_file, 'w') as outfile:
                for data in all_data:
                    outfile.write(data)
                    if not data.endswith('\n'):
                        outfile.write('\n')
            click.echo(f"  Merged {len(all_data)} text files -> {output_file}")
    
    # Clean up chunk directories if all merges successful
    click.echo("\nMerging complete!")
    
    # Create final summary report
    summary_files = list(outdir.glob(f"{summary_prefix}_*.txt"))
    if summary_files:
        click.echo(f"\nGenerated {len(summary_files)} merged summary files:")
        for f in sorted(summary_files):
            click.echo(f"  - {f.name}")

if __name__ == '__main__':
    merge_results()