# GVClass with Prefect + Dask

This document describes how to use the Prefect+Dask orchestration backend for GVClass, which provides robust workflow management with state tracking, caching, and failure recovery.

**Note**: This implementation uses Prefect 3.x, which has some API changes from Prefect 2.x.

## Overview

The Prefect integration provides:

- **State Management**: Tracks pipeline progress and allows resuming from failures
- **Task Caching**: Avoids re-running completed tasks (24-hour cache by default)
- **Parallel Execution**: Leverages Dask for distributed computing
- **Rich Monitoring**: Web UI for tracking pipeline runs
- **Deployment Options**: Local, HPC (SLURM/PBS/SGE), and scheduled execution
- **Automatic Retries**: Failed tasks retry automatically with exponential backoff

## Installation

The Prefect dependencies are already included in `pixi.toml`:

```bash
# Install all dependencies
pixi install

# Verify installation
python test_prefect_pipeline.py
```

## Quick Start

### 1. Run Pipeline Directly

The simplest way to run the pipeline:

```bash
# Run on example data
pixi run gvclass prefect run example/

# Run with custom output directory
pixi run gvclass prefect run my_sequences/ --output-dir results/

# Run with more workers
pixi run gvclass prefect run my_sequences/ -j 8
```

### 2. Using Deployments (Recommended for Production)

Create deployments for different scenarios:

```bash
# Create all deployments
pixi run gvclass prefect deploy

# This creates:
# - gvclass-local: For local execution
# - gvclass-hpc: For HPC clusters
# - gvclass-scheduled: For scheduled runs
```

Start the Prefect server:

```bash
# In one terminal
prefect server start

# In another terminal, start an agent
prefect agent start -q local
```

Run a deployment (Prefect 3.x uses different syntax):

```bash
# Serve deployments
pixi run gvclass prefect serve-local

# In another terminal, trigger a run
from prefect import get_client
async with get_client() as client:
    await client.create_flow_run_from_deployment(
        deployment_id="gvclass-local",
        parameters={"query_dir": "example", "output_dir": "results"}
    )
```

## Command-Line Interface

### Main Commands

```bash
# Run pipeline
pixi run gvclass prefect run <query_dir> [options]

# Create deployments
pixi run gvclass prefect deploy

# Run a deployment
pixi run gvclass prefect run-deployment -d <deployment> <query_dir>

# List recent runs
pixi run gvclass prefect list-runs

# Inspect a specific run
pixi run gvclass prefect inspect <flow_run_id>

# Open web dashboard
pixi run gvclass prefect dashboard

# Start interactive server
pixi run gvclass prefect serve-local
```

### Run Options

```bash
pixi run gvclass prefect run --help

Options:
  -o, --output-dir PATH         Output directory
  -d, --database PATH          Database path (default: resources)
  -j, --n-workers INTEGER      Number of Dask workers (default: 4)
  --threads-per-worker INTEGER Threads per worker (default: 4)
  --tree-method [fasttree|iqtree]
  --mode-fast/--no-mode-fast   Skip order-level trees (default: fast)
  --cluster-type [local|slurm|pbs|sge]
  --run-async/--run-sync       Async or sync execution
  -v, --verbose               Verbose output
```

## Using Different Cluster Types

### Local Cluster (Default)

```bash
pixi run gvclass prefect run example/ --cluster-type local -j 4
```

### SLURM Cluster

```bash
pixi run gvclass prefect run example/ \
  --cluster-type slurm \
  -j 8 \
  --threads-per-worker 8
```

The SLURM configuration can be customized in `src/pipeline/dask_clusters.py`.

### PBS/SGE Clusters

Similar to SLURM:

```bash
pixi run gvclass prefect run example/ --cluster-type pbs -j 8
```

## Monitoring and Management

### Web Dashboard

Access the Prefect UI at http://localhost:4200 after starting the server:

```bash
prefect server start
```

The dashboard shows:
- Flow run history
- Task status and logs
- Resource utilization
- Error messages

### Command-Line Monitoring

```bash
# List recent runs
pixi run gvclass prefect list-runs -n 20

# Filter by state
pixi run gvclass prefect list-runs --state FAILED

# Inspect a specific run
pixi run gvclass prefect inspect <flow_run_id>
```

## Advanced Usage

### Custom Deployments

Create custom deployments in `src/pipeline/deployment.py`:

```python
from prefect.deployments import Deployment
from src.pipeline.prefect_flow import gvclass_flow

deployment = Deployment.build_from_flow(
    flow=gvclass_flow,
    name="my-custom-deployment",
    parameters={
        "n_workers": 16,
        "tree_method": "iqtree",
        "mode_fast": False
    },
    tags=["custom", "high-memory"]
)
deployment.apply()
```

### Scheduled Runs

The `gvclass-scheduled` deployment runs daily at 2 AM UTC. Modify in `deployment.py`:

```python
from prefect.server.schemas.schedules import CronSchedule

schedule=CronSchedule(
    cron="0 */6 * * *",  # Every 6 hours
    timezone="US/Pacific"
)
```

### Task Caching

Tasks are cached for 24 hours by default. Modify in `prefect_flow.py`:

```python
from datetime import timedelta

CACHE_EXPIRATION = timedelta(hours=48)  # 48-hour cache
```

### Retry Configuration

Tasks retry automatically on failure. Customize per task:

```python
@task(
    retries=5,  # Maximum retry attempts
    retry_delay_seconds=30,  # Initial delay
    retry_exponential_backoff=True  # Exponential backoff
)
def my_task():
    pass
```

## Troubleshooting

### Common Issues

1. **Import Errors**
   ```bash
   # Reinstall in development mode
   pixi run pip install -e .
   ```

2. **Prefect Server Issues**
   ```bash
   # Reset Prefect database
   prefect server database reset -y
   ```

3. **Dask Worker Issues**
   ```bash
   # Check Dask dashboard (usually port 8787)
   # Increase worker memory if needed
   --memory-per-worker 16GB
   ```

4. **Task Failures**
   - Check logs in Prefect UI
   - Inspect specific run: `pixi run gvclass prefect inspect <id>`
   - Tasks will retry automatically (default: 2-3 attempts)

### Cleanup

Remove all Prefect artifacts:

```bash
pixi run gvclass prefect cleanup --force
```

## Best Practices

1. **Use Deployments for Production**: Direct runs are great for testing, but deployments provide better management

2. **Monitor Resource Usage**: Check Dask dashboard for memory/CPU usage

3. **Leverage Caching**: Re-runs will skip completed tasks automatically

4. **Set Appropriate Retries**: Increase retries for flaky operations (network, I/O)

5. **Use Fast Mode for Testing**: `--mode-fast` skips some analyses for quicker results

## Example Workflows

### Process Multiple Directories

```python
#!/usr/bin/env python
from src.pipeline.prefect_flow import gvclass_flow
from pathlib import Path

# Process multiple sample directories
for sample_dir in Path("samples").iterdir():
    if sample_dir.is_dir():
        gvclass_flow(
            query_dir=str(sample_dir),
            output_dir=f"results/{sample_dir.name}",
            n_workers=4
        )
```

### Custom HPC Submission

```python
from src.pipeline.dask_clusters import create_slurm_cluster
from src.pipeline.prefect_flow import gvclass_flow
from dask.distributed import Client

# Create SLURM cluster with custom settings
cluster = create_slurm_cluster(
    queue="gpu",
    project="myproject",
    walltime="08:00:00",
    cores=16,
    memory="64GB"
)

# Scale to 10 nodes
cluster.scale(10)

# Run pipeline
with Client(cluster) as client:
    result = gvclass_flow(
        query_dir="large_dataset/",
        output_dir="results/",
        n_workers=10
    )
```

## Migration from Snakemake

The Prefect pipeline provides the same functionality as the Snakemake workflow with additional benefits:

| Snakemake | Prefect |
|-----------|---------|
| `snakemake -j 8` | `gvclass prefect run -j 8` |
| `--forceall` | Delete cache or use new output dir |
| `--dryrun` | Use `--run-async` to submit without waiting |
| Rule graph | View in Prefect UI |
| Log files | Centralized in Prefect UI + local logs |

## Support

For issues specific to the Prefect integration:

1. Check this documentation
2. Run the test script: `python test_prefect_pipeline.py`
3. Check Prefect UI for detailed error messages
4. Open an issue on GitHub with the flow run ID