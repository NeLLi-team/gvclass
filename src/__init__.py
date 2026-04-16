"""
GVClass: Giant Virus Classification Pipeline

A bioinformatics tool for assigning taxonomy to giant virus contigs
or metagenome assembled genomes (GVMAGs).
"""

from src.__version__ import __version__

__author__ = "NeLLi Team"

# Core functionality imports
from . import config
from . import core
from . import utils

__all__ = ["config", "core", "utils", "__version__"]
