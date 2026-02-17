"""
GVClass: Giant Virus Classification Pipeline

A bioinformatics tool for assigning taxonomy to giant virus contigs
or metagenome assembled genomes (GVMAGs).
"""

__version__ = "1.2.2"
__author__ = "NeLLi Team"

# Core functionality imports
from . import config
from . import core
from . import utils

__all__ = ["config", "core", "utils", "__version__"]
