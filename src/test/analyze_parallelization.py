#!/usr/bin/env python
"""
Analyze parallelization in GVClass pipeline to understand performance bottlenecks.
"""

import sys
import time
import subprocess
import argparse
from datetime import datetime


def analyze_pipeline_parallelization():
    """Analyze the parallelization structure of the GVClass pipeline."""
    
    print("="*80)
    print("GVClass Pipeline Parallelization Analysis")
    print("="*80)
    print()
    
    # Level 1: Query-level parallelization
    print("📊 PARALLELIZATION LEVELS:")
    print()
    print("1. QUERY-LEVEL PARALLELIZATION (Dask Workers)")
    print("   - Each query (genome) is processed by a separate Dask worker")
    print("   - Number of workers = min(n_queries, threads // threads_per_worker)")
    print("   - Example: 2 queries, 16 threads → 2 workers × 8 threads each")
    print()
    
    # Level 2: Within-query parallelization
    print("2. WITHIN-QUERY PARALLELIZATION:")
    print()
    print("   A. Gene Calling (SEQUENTIAL per query)")
    print("      - Tests multiple genetic codes (0,1,4,6,11,15,29,106,129)")
    print("      - Each code tested sequentially (not parallel)")
    print("      - ⚠️ BOTTLENECK: Could parallelize genetic code testing")
    print()
    
    print("   B. HMM Search (SINGLE-THREADED per query)")
    print("      - pyhmmer search against all HMM models")
    print("      - Uses 1 thread (pyhmmer is fast but could use more threads)")
    print("      - ⚠️ BOTTLENECK: Could use multi-threaded HMM search")
    print()
    
    print("   C. Marker Processing (PARALLEL within each query)")
    print("      - Each marker (up to 99 GVOGs) processed in parallel")
    print("      - ThreadPoolExecutor with min(threads_per_worker, n_markers)")
    print("      - ✅ GOOD: This is properly parallelized")
    print("      - Each marker does:")
    print("        • BLAST search (1 thread)")
    print("        • Alignment (pyfamsa, 1 thread)")
    print("        • Trimming (pytrimal, 1 thread)")
    print("        • Tree building (fasttree/iqtree, 1 thread)")
    print()
    
    # Performance implications
    print("📈 PERFORMANCE IMPLICATIONS:")
    print()
    print("For 2 queries with different thread allocations:")
    print()
    print("- 4 threads total:")
    print("  • 2 workers × 2 threads")
    print("  • Each query can process 2 markers in parallel")
    print("  • Total parallelism: 2 queries × 2 markers = 4 parallel tasks max")
    print()
    print("- 8 threads total:")
    print("  • 2 workers × 4 threads")
    print("  • Each query can process 4 markers in parallel")
    print("  • Total parallelism: 2 queries × 4 markers = 8 parallel tasks max")
    print()
    print("- 16 threads total:")
    print("  • 2 workers × 8 threads")
    print("  • Each query can process 8 markers in parallel")
    print("  • Total parallelism: 2 queries × 8 markers = 16 parallel tasks max")
    print()
    
    # Bottlenecks
    print("🚧 IDENTIFIED BOTTLENECKS:")
    print()
    print("1. GENETIC CODE OPTIMIZATION (Sequential)")
    print("   - Currently: ~10-15 seconds per code × 9 codes = ~90-135 seconds")
    print("   - If parallelized: ~15-20 seconds total (6-9x speedup possible)")
    print()
    print("2. HMM SEARCH (Single-threaded)")
    print("   - Currently uses 1 thread")
    print("   - pyhmmer supports multi-threading")
    print("   - Could use 2-4 threads for 2-4x speedup")
    print()
    print("3. WORKER ALLOCATION")
    print("   - With only 2 queries, max 2 workers regardless of threads")
    print("   - Extra threads only help with marker parallelization")
    print("   - For better scaling: need more queries or different architecture")
    print()
    
    # Recommendations
    print("💡 RECOMMENDATIONS FOR BETTER PARALLELIZATION:")
    print()
    print("1. Parallelize genetic code testing:")
    print("   - Use ProcessPoolExecutor for genetic code optimization")
    print("   - Test all codes in parallel, select best")
    print()
    print("2. Enable multi-threaded HMM search:")
    print("   - Pass threads parameter to pyhmmer")
    print("   - Use 2-4 threads per HMM search")
    print()
    print("3. For small query sets (< 4):")
    print("   - Consider task-level parallelization instead of query-level")
    print("   - Submit each step as separate Dask task")
    print("   - Allow multiple steps to run in parallel")
    print()
    print("4. Memory considerations:")
    print("   - Each worker uses ~1-2GB RAM")
    print("   - Monitor memory to avoid swapping")
    print("   - Adjust workers based on available RAM")
    print()


def run_performance_test(query_dir, threads_list=[4, 8, 16]):
    """Run performance tests with different thread counts."""
    
    print("\n" + "="*80)
    print("Running Performance Tests")
    print("="*80)
    print()
    
    results = []
    
    for threads in threads_list:
        print(f"\nTesting with {threads} threads...")
        
        # Create unique output dir
        output_dir = f"test_output_{threads}threads_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Run with CPU monitoring
        cmd = [
            sys.executable, 
            "src/test/monitor_performance.py",
            "--",
            "pixi", "run", "gvclass",
            query_dir,
            "-o", output_dir,
            "-t", str(threads)
        ]
        
        start_time = time.time()
        try:
            subprocess.run(cmd, check=True)
            runtime = time.time() - start_time
            
            results.append({
                'threads': threads,
                'runtime': runtime,
                'output_dir': output_dir
            })
            
            print(f"✅ Completed in {runtime:.1f} seconds")
            
        except subprocess.CalledProcessError:
            print(f"❌ Failed with {threads} threads")
    
    # Summary
    if results:
        print("\n" + "="*80)
        print("Performance Summary")
        print("="*80)
        print()
        print("Threads | Runtime | Speedup")
        print("--------|---------|--------")
        
        baseline = results[0]['runtime']
        for r in results:
            speedup = baseline / r['runtime']
            print(f"{r['threads']:7d} | {r['runtime']:6.1f}s | {speedup:6.2f}x")


def main():
    parser = argparse.ArgumentParser(description='Analyze GVClass pipeline parallelization')
    parser.add_argument('--test', action='store_true', 
                       help='Run performance tests')
    parser.add_argument('--query-dir', default='example',
                       help='Query directory for tests')
    parser.add_argument('--threads', nargs='+', type=int, default=[4, 8, 16],
                       help='Thread counts to test')
    
    args = parser.parse_args()
    
    # Always show analysis
    analyze_pipeline_parallelization()
    
    # Run tests if requested
    if args.test:
        run_performance_test(args.query_dir, args.threads)
    else:
        print("\nTo run performance tests, use: --test")
        print("To monitor CPU usage during a run:")
        print("  python src/test/monitor_cpu_usage.py pixi run gvclass example -o output -t 16")


if __name__ == '__main__':
    main()