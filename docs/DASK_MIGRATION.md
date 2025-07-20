# Dask Migration Guide

## Overview

GVClass has been migrated from Snakemake to Dask for workflow orchestration. This provides:

- Better Python integration
- Improved parallelization and resource management
- Real-time monitoring via Dask dashboard
- Support for distributed computing (local, SLURM, etc.)

## Key Benefits

1. **Native Python**: No DSL to learn - everything is pure Python
2. **Dynamic scaling**: Workers scale up/down based on workload
3. **Memory efficiency**: Dask handles spilling to disk automatically
4. **Monitoring**: Real-time dashboard at http://localhost:8787
5. **Fault tolerance**: Failed tasks can be retried
6. **Distributed**: Easy scaling from laptop to HPC cluster

## Usage

### Basic Usage

```bash
# Using the wrapper script (recommended)
./gvclass_dask.sh example -j 8

# Using Python directly
pixi run python -m src.bin.gvclass_cli run example --n-workers 8

# With custom database location
./gvclass_dask.sh example -j 8 -d /path/to/database
```

### Advanced Options

```bash
# Use SLURM cluster
pixi run python -m src.bin.gvclass_cli run example \
    --cluster-type slurm \
    --n-workers 20 \
    --memory-per-worker 16GB

# Custom output directory
pixi run python -m src.bin.gvclass_cli run example \
    --output-dir results_custom \
    --verbose

# View Dask dashboard
pixi run python -m src.bin.gvclass_cli dashboard --scheduler localhost:8786
```

## Architecture

### Pipeline Structure

```python
# Task graph is built dynamically
pipeline = GVClassPipeline(
    query_dir="example",
    n_workers=8,
    cluster_type="local"
)

# Tasks are submitted as Dask delayed objects
reformat → gene_calling → hmm_search → extract_hits → blast → align → tree → summarize
```

### Task Dependencies

- Tasks automatically wait for their dependencies
- Parallel execution where possible
- Resource-aware scheduling

### File Management

- Input/output files are tracked automatically
- No need for explicit rule definitions
- Temporary files cleaned up based on configuration

## Monitoring

Access the Dask dashboard while the pipeline is running:
- http://localhost:8787 (default)
- Shows task progress, memory usage, worker status

## Performance Tips

1. **Workers**: Set workers to ~(CPU cores / 2) for local runs
2. **Memory**: Each worker should have 4-8GB RAM
3. **Cluster**: Use SLURM for large datasets
4. **Fast mode**: Enable `--mode-fast` to skip some analyses

## Differences from Snakemake

| Feature | Snakemake | Dask |
|---------|-----------|------|
| Configuration | YAML + Snakefile | Python code |
| Parallelization | Process-based | Thread/process hybrid |
| Monitoring | Terminal output | Web dashboard |
| Fault tolerance | Restart from checkpoint | Automatic retry |
| Resource management | Manual specification | Dynamic |

## Troubleshooting

### Memory Issues
```bash
# Increase memory per worker
--memory-per-worker 16GB
```

### Connection Issues
```bash
# Use a different scheduler port
export DASK_SCHEDULER_PORT=8788
```

### Debug Mode
```bash
# Enable verbose logging
--verbose
```

## Migration from Snakemake

If you have existing Snakemake workflows:

1. **Outputs are compatible**: Same file structure maintained
2. **Config translation**: 
   - `snakemake -j 8` → `--n-workers 8`
   - `config['treeoption']` → `--tree-method`
3. **Custom rules**: Can be added as new tasks in `src/pipeline/tasks.py`