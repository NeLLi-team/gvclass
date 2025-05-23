#!/usr/bin/env python
import argparse
import datetime
import json

def log_run_info(version, config, output_file):
    with open(output_file, "w") as f:
        f.write(f"GVClass Pipeline Version: {version}\n")
        f.write(f"Run Date and Time: {datetime.datetime.now()}\n\n")
        f.write("Configuration:\n")
        
        # Parse the config JSON
        config_dict = json.loads(config)
        
        # Write the config as formatted JSON
        json.dump(config_dict, f, indent=2)
        f.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Log run information')
    parser.add_argument('--version', required=True, help='Pipeline version')
    parser.add_argument('--config', required=True, help='Configuration as JSON string')
    parser.add_argument('--output', required=True, help='Output log file')
    
    args = parser.parse_args()
    
    log_run_info(args.version, args.config, args.output)
