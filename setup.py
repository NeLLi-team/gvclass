#!/usr/bin/env python
"""Setup script for GVClass."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="gvclass",
    version="1.1.0",
    author="NeLLi Team",
    description="Giant Virus Classification Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NeLLi-team/gvclass",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    install_requires=[
        "click>=8.0.0",
        "pandas>=2.0",
        "numpy>=2.3",
        "biopython>=1.85",
        "pyhmmer>=0.10.0",
        "pytrimal>=0.8.0",
        "pyfamsa>=0.5.0",
        "pyrodigal>=3.0.0",
        "pyswrd>=0.1.0",
        "veryfasttree>=3.0.0",
    ],
    python_requires=">=3.11",
    entry_points={
        "console_scripts": [
            "gvclass=src.bin.gvclass_cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
    ],
)
