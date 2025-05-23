#!/usr/bin/env python
import os
import argparse
import subprocess
from pathlib import Path

def check_and_download_resources(database_path):
    database_path = Path(database_path).resolve()
    required_files = [
        "models/combined.hmm",
        "ncldvApril24_labels.txt",
        "models_APRIL24--databaseApril24.cutoffs",
        "order_completeness.tab"
    ]

    # Create the directory if it doesn't exist
    if not database_path.exists():
        print(f"Resources directory not found at {database_path}. Creating directory...")
        database_path.mkdir(parents=True, exist_ok=True)
    else:
        print(f"Resources directory found at {database_path}.")

    # Check for missing files
    missing_files = [f for f in required_files if not (database_path / f).exists()]

    # Download resources if any files are missing
    if missing_files:
        print(f"The following required files are missing: {', '.join(missing_files)}")
        print("Downloading resources...")

        # Save current directory
        current_dir = os.getcwd()

        # Change to database directory
        os.chdir(str(database_path))

        # Download and extract resources
        url = "https://portal.nersc.gov/cfs/nelli/gvclassDB/resources.tar.gz"
        print(f"Downloading from {url}...")
        subprocess.run(["wget", "-q", url], check=True)

        print("Extracting resources...")
        subprocess.run(["tar", "-xzf", "resources.tar.gz", "--strip-components=1"], check=True)

        print("Cleaning up...")
        subprocess.run(["rm", "resources.tar.gz"], check=True)

        # Return to original directory
        os.chdir(current_dir)

        print(f"Resources downloaded and extracted to {database_path}.")
    else:
        print("All required resources are available.")

    return database_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check and download resources')
    parser.add_argument('--database_path', required=True, help='Path to database directory')
    parser.add_argument('--output', required=True, help='Output file to indicate completion')

    args = parser.parse_args()

    check_and_download_resources(args.database_path)

    # Create output file to indicate completion
    with open(args.output, 'w') as f:
        f.write("Resources checked and available\n")
