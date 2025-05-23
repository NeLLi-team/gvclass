#!/usr/bin/env python
import sys
import subprocess
import os

# Try to import click, install if not available
try:
    import click
except ImportError:
    print("Installing click package...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "click"])
    import click

# Get the path to the original script
original_script = os.path.join(os.environ.get('WORKFLOW_PROJECTDIR', '.'), "workflow/scripts/reformat.py")

# Pass all arguments to the original script
cmd = [sys.executable, original_script] + sys.argv[1:]
subprocess.check_call(cmd)
