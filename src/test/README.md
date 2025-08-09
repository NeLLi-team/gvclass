# GVClass Performance Testing and Analysis Tools

This directory contains tools for analyzing the performance and parallelization of the GVClass pipeline.

## Tools

### 1. CPU Monitoring (`monitor_cpu_usage.py`)

Monitor CPU usage and process activity during pipeline execution.

**Usage:**
```bash
# Monitor a pipeline run
python src/test/monitor_cpu_usage.py pixi run gvclass example -o output -t 16

# With custom output directory and interval
python src/test/monitor_cpu_usage.py -o my_cpu_logs -i 0.5 pixi run gvclass example -o output -t 16
```

**Output:**
- `cpu_usage_TIMESTAMP.jsonl` - Detailed CPU samples in JSON format
- `cpu_summary_TIMESTAMP.txt` - Human-readable summary with:
  - Overall CPU utilization
  - Per-core usage statistics
  - Process activity breakdown
  - Parallelization analysis

### 2. Parallelization Analysis (`analyze_parallelization.py`)

Analyze and explain the pipeline's parallelization architecture.

**Usage:**
```bash
# Show parallelization analysis
python src/test/analyze_parallelization.py

# Run performance tests with different thread counts
python src/test/analyze_parallelization.py --test --threads 4 8 16 32

# Test with custom query directory
python src/test/analyze_parallelization.py --test --query-dir my_genomes --threads 8 16
```

**Output:**
- Detailed explanation of parallelization levels
- Performance implications for different thread counts
- Identification of bottlenecks
- Recommendations for optimization

## Understanding GVClass Parallelization

### Current Architecture

1. **Query-Level Parallelization (Dask Workers)**
   - Each genome/query runs on a separate Dask worker
   - Workers = min(n_queries, threads // threads_per_worker)
   - Limited by number of input files

2. **Within-Query Parallelization**
   - **Gene Calling**: Sequential testing of genetic codes (BOTTLENECK)
   - **HMM Search**: Single-threaded pyhmmer (BOTTLENECK)
   - **Marker Processing**: Parallel using ThreadPoolExecutor (GOOD)
     - BLAST, alignment, trimming, tree building for each marker

### Performance Expectations

For the example dataset (2 queries):

| Threads | Workers | Parallelism | Expected Speedup |
|---------|---------|-------------|------------------|
| 4       | 2×2     | 4 tasks max | Baseline         |
| 8       | 2×4     | 8 tasks max | ~1.2-1.4x        |
| 16      | 2×8     | 16 tasks max| ~1.4-1.6x        |

The limited speedup is due to:
- Only 2 queries (max 2 workers)
- Sequential genetic code testing
- Single-threaded HMM search

### Identified Bottlenecks

1. **Genetic Code Optimization** (90-135 seconds)
   - Tests 9 codes sequentially
   - Could be parallelized for 6-9x speedup

2. **HMM Search** (10-20 seconds)
   - Uses 1 thread
   - pyhmmer supports multi-threading

3. **Small Query Sets**
   - With 2 queries, max 2 workers
   - Extra threads only help marker processing

## Running Performance Analysis

### Quick Test
```bash
# Analyze current parallelization
python src/test/analyze_parallelization.py

# Monitor a single run
python src/test/monitor_cpu_usage.py pixi run gvclass example -o test_output -t 16
```

### Comprehensive Test
```bash
# Run performance tests with CPU monitoring
python src/test/analyze_parallelization.py --test --threads 4 8 16 32

# Check CPU logs
ls cpu_monitor_*threads/cpu_summary_*.txt
```

### Interpreting Results

Look for:
- **CPU Utilization**: Should approach (n_workers × 100)%
- **Active Cores**: Should match allocated threads
- **Parallel Samples**: Higher percentage = better parallelization
- **Process Distribution**: Even distribution across workers

## Recommendations

1. **For Better Performance**:
   - Process more queries at once (4+ files)
   - Use 4-8 threads per worker
   - Ensure sufficient memory (2GB per worker)

2. **For Pipeline Optimization**:
   - Parallelize genetic code testing
   - Enable multi-threaded HMM search
   - Consider task-level parallelization for small datasets