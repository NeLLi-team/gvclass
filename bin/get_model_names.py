#!/usr/bin/env python
import argparse

def get_model_names(hmm_file, mode_fast, output_file):
    # Read model names from HMM file
    modelnames = []
    with open(hmm_file, 'r') as f:
        for line in f:
            if "NAME" in line:
                model_name = line.split()[1]
                modelnames.append(model_name)
    
    # Filter models if in fast mode
    if mode_fast.lower() == 'true':
        modelnames = [x for x in modelnames if not x.startswith("OG")]
    
    # Write model names to output file
    with open(output_file, 'w') as f:
        for model in modelnames:
            f.write(f"{model}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get model names from HMM file')
    parser.add_argument('--hmm_file', required=True, help='Path to HMM file')
    parser.add_argument('--mode_fast', required=True, help='Whether to use fast mode (true/false)')
    parser.add_argument('--output', required=True, help='Output file for model names')
    
    args = parser.parse_args()
    
    get_model_names(args.hmm_file, args.mode_fast, args.output)
