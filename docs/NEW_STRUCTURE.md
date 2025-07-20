# GVClass Repository Structure

## Overview

The repository has been restructured to follow modern Python best practices with a clear separation of concerns.

## Directory Structure

```
gvclass/
├── src/                    # All source code
│   ├── __init__.py        # Main package init
│   ├── config/            # Configuration modules
│   │   ├── __init__.py
│   │   ├── marker_sets.py # Marker set definitions
│   │   └── README.md      # Config documentation
│   ├── core/              # Core analysis modules
│   │   ├── __init__.py
│   │   ├── alignment.py   # Sequence alignment (pyfamsa/pytrimal)
│   │   ├── blast.py       # Sequence similarity search (pyswrd)
│   │   ├── gene_calling.py # Gene prediction (pyrodigal)
│   │   ├── hmm_search.py  # HMM searches (pyhmmer)
│   │   ├── phylogeny.py   # Tree building (veryfasttree/iqtree)
│   │   ├── summarize.py   # Result summarization
│   │   ├── nearest_neighbors.py # Nearest neighbor analysis
│   │   ├── extract_qhits.py    # Extract query hits
│   │   ├── reformat.py         # Input reformatting
│   │   ├── combinedout.py      # Combine outputs
│   │   └── merge_score_matrices.py # Merge score matrices
│   ├── utils/             # Utility functions
│   │   ├── __init__.py
│   │   ├── error_handling.py    # Error handling classes
│   │   ├── input_validation.py  # Input validation functions
│   │   └── common.py            # Common utilities
│   ├── bin/               # CLI entry points
│   │   ├── __init__.py
│   │   └── gvclass_cli.py      # Main CLI interface
│   └── test/              # Unit tests (to be added)
│       └── __init__.py
├── workflow/              # Snakemake workflow
│   ├── Snakefile         # Main workflow definition
│   ├── config.yml        # Workflow configuration
│   └── envs/             # Environment definitions
├── resources/            # Databases and references
├── docs/                 # Documentation
├── scripts/              # Standalone scripts
│   └── setup_database.py # Database setup script
├── example/              # Example data
├── pixi.toml            # Pixi dependency management
├── pixi.lock
├── setup.py             # Package installation
├── pyproject.toml       # Modern Python packaging
├── README.md
├── LICENSE
└── CHANGELOG.md
```

## Benefits of New Structure

1. **Cleaner imports**: Instead of complex relative imports, we now have:
   ```python
   from src.config import MIRUS_MODELS
   from src.core import alignment
   from src.utils import validate_file_path
   ```

2. **Testable**: Proper package structure makes unit testing easier

3. **Installable**: Can use `pip install -e .` for development

4. **CI/CD ready**: Standard structure for automated testing

5. **Maintainable**: Clear separation of concerns

## Migration Notes

### For Developers

- All workflow scripts are now in `src/core/`
- Utility functions are in `src/utils/`
- Configuration is in `src/config/`
- When updating imports, use absolute imports from `src`

### For Users

- The main entry points remain the same (bash scripts)
- Snakemake workflow usage is unchanged
- All functionality is preserved

## Next Steps

1. Update all imports in the moved files
2. Update Snakefile to use new paths
3. Add comprehensive unit tests in `src/test/`
4. Update CI/CD workflows
5. Consider adding more CLI commands to `src/bin/`