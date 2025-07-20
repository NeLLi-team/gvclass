#!/usr/bin/env python
"""
Enhanced CPU usage monitor for GVClass pipeline with better process tracking.
"""

import sys
import psutil
import time
import threading
import subprocess
from datetime import datetime
from pathlib import Path
import argparse
import json
import signal
from collections import defaultdict, deque


class EnhancedCPUMonitor:
    """Enhanced monitor with better process tracking and analysis."""
    
    def __init__(self, output_dir, interval=0.5, target_process=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.interval = interval
        self.target_process = target_process
        self.running = True
        self.start_time = datetime.now()
        
        # Output files
        timestamp = self.start_time.strftime('%Y%m%d_%H%M%S')
        self.cpu_log = self.output_dir / f"cpu_usage_{timestamp}.jsonl"
        self.summary_file = self.output_dir / f"cpu_summary_{timestamp}.txt"
        self.timeline_file = self.output_dir / f"cpu_timeline_{timestamp}.txt"
        
        # Enhanced tracking
        self.samples = []
        self.process_history = defaultdict(lambda: {
            'samples': [],
            'cmdline': '',
            'name': '',
            'parent_pid': None,
            'start_time': None,
            'end_time': None
        })
        self.phase_detection = {
            'gene_calling': deque(maxlen=20),
            'hmm_search': deque(maxlen=20),
            'marker_processing': deque(maxlen=20),
            'tree_building': deque(maxlen=20)
        }
        
        # Pre-calculate CPU percentages for all processes
        for proc in psutil.process_iter(['pid', 'name']):
            try:
                proc.cpu_percent(interval=0)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
    
    def detect_phase(self, processes):
        """Detect which phase of the pipeline is running."""
        phase_indicators = {
            'gene_calling': ['pyrodigal', 'gene_calling', 'optimize_genetic_code'],
            'hmm_search': ['pyhmmer', 'hmmsearch', 'hmm_search'],
            'marker_processing': ['blast', 'pyfamsa', 'pytrimal', 'marker'],
            'tree_building': ['VeryFastTree', 'iqtree', 'fasttree']
        }
        
        current_phases = set()
        for proc in processes:
            cmdline = proc.get('cmdline', '')
            name = proc.get('name', '')
            
            for phase, indicators in phase_indicators.items():
                if any(ind in cmdline or ind in name for ind in indicators):
                    current_phases.add(phase)
        
        return current_phases
    
    def get_process_info(self, proc):
        """Get detailed process information with CPU tracking."""
        try:
            # Get CPU with interval=0 for accurate measurement
            cpu_percent = proc.cpu_percent(interval=0)
            
            # Get memory info
            mem_info = proc.memory_info()
            
            # Get IO stats if available
            try:
                io_stats = proc.io_counters()
                io_read_mb = io_stats.read_bytes / 1024 / 1024
                io_write_mb = io_stats.write_bytes / 1024 / 1024
            except (AttributeError, psutil.AccessDenied):
                io_read_mb = 0
                io_write_mb = 0
            
            # Get thread info
            threads = proc.threads()
            
            return {
                'pid': proc.pid,
                'name': proc.name(),
                'cmdline': ' '.join(proc.cmdline())[:300],
                'cpu_percent': cpu_percent,
                'cpu_num': proc.cpu_num() if hasattr(proc, 'cpu_num') else -1,
                'memory_mb': mem_info.rss / 1024 / 1024,
                'memory_vms_mb': mem_info.vms / 1024 / 1024,
                'num_threads': proc.num_threads(),
                'thread_ids': [t.id for t in threads],
                'status': proc.status(),
                'create_time': proc.create_time(),
                'parent_pid': proc.ppid(),
                'io_read_mb': io_read_mb,
                'io_write_mb': io_write_mb
            }
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return None
    
    def monitor(self):
        """Enhanced monitoring loop."""
        print("Starting enhanced CPU monitoring...")
        
        while self.running:
            timestamp = datetime.now()
            elapsed = (timestamp - self.start_time).total_seconds()
            
            # System-wide stats
            cpu_overall = psutil.cpu_percent(interval=0)
            cpu_per_core = psutil.cpu_percent(interval=0, percpu=True)
            cpu_freq = psutil.cpu_freq(percpu=True) if hasattr(psutil, 'cpu_freq') else None
            
            # Memory and IO
            memory = psutil.virtual_memory()
            swap = psutil.swap_memory()
            io_counters = psutil.disk_io_counters()
            
            # Get all relevant processes
            processes = []
            process_tree = {}
            
            if self.target_process:
                try:
                    # Get target process and all descendants
                    parent = psutil.Process(self.target_process)
                    all_procs = [parent] + parent.children(recursive=True)
                    
                    for proc in all_procs:
                        proc_info = self.get_process_info(proc)
                        if proc_info and (proc_info['cpu_percent'] > 0 or proc_info['memory_mb'] > 10):
                            processes.append(proc_info)
                            process_tree[proc_info['pid']] = proc_info
                except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                    pass
            
            # Also monitor pipeline-related processes globally
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                try:
                    cmdline = ' '.join(proc.info.get('cmdline', []))
                    if any(keyword in cmdline for keyword in 
                          ['gvclass', 'prefect', 'dask', 'pyhmmer', 'pyrodigal',
                           'fasttree', 'VeryFastTree', 'iqtree', 'pyfamsa', 'pytrimal',
                           'blast', 'diamond', 'genetic_code', 'marker']):
                        if proc.pid not in process_tree:
                            proc_info = self.get_process_info(proc)
                            if proc_info and (proc_info['cpu_percent'] > 0 or proc_info['memory_mb'] > 20):
                                processes.append(proc_info)
                                
                                # Track process history
                                pid = proc_info['pid']
                                self.process_history[pid]['samples'].append({
                                    'time': elapsed,
                                    'cpu': proc_info['cpu_percent'],
                                    'memory': proc_info['memory_mb'],
                                    'cpu_num': proc_info['cpu_num']
                                })
                                self.process_history[pid]['cmdline'] = proc_info['cmdline']
                                self.process_history[pid]['name'] = proc_info['name']
                                self.process_history[pid]['parent_pid'] = proc_info['parent_pid']
                                if not self.process_history[pid]['start_time']:
                                    self.process_history[pid]['start_time'] = elapsed
                                self.process_history[pid]['end_time'] = elapsed
                except (psutil.NoSuchProcess, psutil.AccessDenied, AttributeError, TypeError):
                    continue
            
            # Detect current phase
            current_phases = self.detect_phase(processes)
            for phase in current_phases:
                self.phase_detection[phase].append(elapsed)
            
            # Calculate active cores (>25% usage)
            active_cores = sum(1 for cpu in cpu_per_core if cpu > 25)
            high_usage_cores = sum(1 for cpu in cpu_per_core if cpu > 75)
            
            # Calculate process totals
            total_process_cpu = sum(p['cpu_percent'] for p in processes)
            active_processes = len([p for p in processes if p['cpu_percent'] > 1])
            
            # Create sample
            sample = {
                'timestamp': timestamp.isoformat(),
                'elapsed_seconds': elapsed,
                'cpu': {
                    'overall_percent': cpu_overall,
                    'per_core': cpu_per_core,
                    'active_cores': active_cores,
                    'high_usage_cores': high_usage_cores,
                    'frequencies': [f.current for f in cpu_freq] if cpu_freq else []
                },
                'memory': {
                    'total_mb': memory.total / 1024 / 1024,
                    'used_mb': memory.used / 1024 / 1024,
                    'percent': memory.percent,
                    'available_mb': memory.available / 1024 / 1024
                },
                'swap': {
                    'used_mb': swap.used / 1024 / 1024,
                    'percent': swap.percent
                },
                'io': {
                    'read_mb': io_counters.read_bytes / 1024 / 1024 if io_counters else 0,
                    'write_mb': io_counters.write_bytes / 1024 / 1024 if io_counters else 0
                },
                'processes': processes,
                'process_summary': {
                    'total': len(processes),
                    'active': active_processes,
                    'total_cpu': total_process_cpu
                },
                'phases': list(current_phases)
            }
            
            # Write to log
            with open(self.cpu_log, 'a') as f:
                f.write(json.dumps(sample) + '\n')
            
            self.samples.append(sample)
            
            # Print live summary
            if int(elapsed) % 5 == 0:
                print(f"\r[{elapsed:6.0f}s] CPU: {cpu_overall:5.1f}% ({active_cores} cores active) | "
                      f"Processes: {active_processes} active ({total_process_cpu:5.0f}% total) | "
                      f"Memory: {memory.percent:4.1f}% | "
                      f"Phase: {', '.join(current_phases) if current_phases else 'idle':<20s}", 
                      end='', flush=True)
            
            time.sleep(self.interval)
    
    def create_enhanced_summary(self):
        """Create comprehensive summary with timeline and analysis."""
        if not self.samples:
            return
        
        print("\n\nGenerating enhanced summary...")
        
        # Basic statistics
        duration = self.samples[-1]['elapsed_seconds']
        cpu_values = [s['cpu']['overall_percent'] for s in self.samples]
        active_cores_values = [s['cpu']['active_cores'] for s in self.samples]
        memory_values = [s['memory']['percent'] for s in self.samples]
        
        with open(self.summary_file, 'w') as f:
            f.write("="*100 + "\n")
            f.write("ENHANCED CPU UTILIZATION SUMMARY\n")
            f.write("="*100 + "\n\n")
            
            f.write(f"Pipeline Runtime: {duration:.1f} seconds ({duration/60:.1f} minutes)\n")
            f.write(f"Monitoring Samples: {len(self.samples)} (interval: {self.interval}s)\n")
            f.write(f"CPU Cores: {len(self.samples[0]['cpu']['per_core'])} logical cores\n\n")
            
            # CPU Analysis
            f.write("CPU UTILIZATION ANALYSIS:\n")
            f.write("-"*50 + "\n")
            f.write(f"Overall Average: {sum(cpu_values)/len(cpu_values):.1f}%\n")
            f.write(f"Overall Peak: {max(cpu_values):.1f}%\n")
            f.write(f"Time > 50%: {sum(1 for v in cpu_values if v > 50)/len(cpu_values)*100:.1f}%\n")
            f.write(f"Time > 75%: {sum(1 for v in cpu_values if v > 75)/len(cpu_values)*100:.1f}%\n")
            f.write(f"\nAverage Active Cores: {sum(active_cores_values)/len(active_cores_values):.1f}\n")
            f.write(f"Peak Active Cores: {max(active_cores_values)}\n\n")
            
            # Per-core statistics
            n_cores = len(self.samples[0]['cpu']['per_core'])
            f.write("PER-CORE UTILIZATION:\n")
            f.write("-"*50 + "\n")
            core_stats = []
            for i in range(n_cores):
                core_values = [s['cpu']['per_core'][i] for s in self.samples]
                core_stats.append({
                    'core': i,
                    'avg': sum(core_values) / len(core_values),
                    'max': max(core_values),
                    'active_time': sum(1 for v in core_values if v > 25) / len(core_values) * 100
                })
            
            for stat in sorted(core_stats, key=lambda x: x['avg'], reverse=True):
                f.write(f"Core {stat['core']:2d}: Avg={stat['avg']:5.1f}% Max={stat['max']:5.1f}% "
                       f"Active={stat['active_time']:5.1f}% of time\n")
            
            # Phase Analysis
            f.write("\n\nPIPELINE PHASE ANALYSIS:\n")
            f.write("-"*50 + "\n")
            for phase, times in self.phase_detection.items():
                if times:
                    phase_duration = max(times) - min(times) if len(times) > 1 else 0
                    f.write(f"{phase:20s}: {phase_duration:6.1f}s "
                           f"({phase_duration/duration*100:5.1f}% of runtime)\n")
            
            # Process Analysis
            f.write("\n\nPROCESS ANALYSIS:\n")
            f.write("-"*50 + "\n")
            
            # Group processes by type
            process_types = defaultdict(lambda: {'count': 0, 'total_cpu': 0, 'total_time': 0})
            
            for pid, data in self.process_history.items():
                if data['samples']:
                    process_type = data['name']
                    # Categorize by function
                    if 'pyrodigal' in data['cmdline']:
                        process_type = 'Gene Calling'
                    elif 'pyhmmer' in data['cmdline']:
                        process_type = 'HMM Search'
                    elif 'VeryFastTree' in process_type:
                        process_type = 'Tree Building'
                    elif 'pyfamsa' in data['cmdline']:
                        process_type = 'Alignment'
                    elif 'pytrimal' in data['cmdline']:
                        process_type = 'Trimming'
                    elif 'blast' in data['cmdline'] or 'diamond' in data['cmdline']:
                        process_type = 'BLAST Search'
                    
                    process_types[process_type]['count'] += 1
                    process_types[process_type]['total_cpu'] += sum(s['cpu'] for s in data['samples'])
                    process_types[process_type]['total_time'] += data['end_time'] - data['start_time']
            
            f.write("\nProcess Type Distribution:\n")
            for ptype, stats in sorted(process_types.items(), 
                                      key=lambda x: x[1]['total_cpu'], reverse=True):
                avg_cpu = stats['total_cpu'] / max(stats['count'], 1)
                f.write(f"{ptype:20s}: {stats['count']:3d} instances, "
                       f"Avg CPU: {avg_cpu:6.1f}%, Total time: {stats['total_time']:6.1f}s\n")
            
            # Parallelization Effectiveness
            f.write("\n\nPARALLELIZATION EFFECTIVENESS:\n")
            f.write("-"*50 + "\n")
            
            # Calculate parallel efficiency
            parallel_samples = 0
            cpu_efficiency_sum = 0
            
            for sample in self.samples:
                active_procs = [p for p in sample['processes'] if p['cpu_percent'] > 5]
                if len(active_procs) > 1:
                    parallel_samples += 1
                    # Efficiency = actual CPU / (cores * 100)
                    efficiency = sample['cpu']['overall_percent'] / (sample['cpu']['active_cores'] * 100) if sample['cpu']['active_cores'] > 0 else 0
                    cpu_efficiency_sum += efficiency
            
            parallel_time_pct = parallel_samples / len(self.samples) * 100
            avg_efficiency = cpu_efficiency_sum / max(parallel_samples, 1) * 100
            
            f.write(f"Time with parallel execution: {parallel_time_pct:.1f}%\n")
            f.write(f"Average CPU efficiency during parallel phases: {avg_efficiency:.1f}%\n")
            f.write(f"Peak parallel processes: {max(s['process_summary']['active'] for s in self.samples)}\n")
            
            # Bottleneck Analysis
            f.write("\n\nBOTTLENECK ANALYSIS:\n")
            f.write("-"*50 + "\n")
            
            # Find low CPU periods
            low_cpu_periods = []
            current_period = None
            
            for i, sample in enumerate(self.samples):
                if sample['cpu']['overall_percent'] < 25:
                    if current_period is None:
                        current_period = {'start': i, 'end': i, 'min_cpu': sample['cpu']['overall_percent']}
                    else:
                        current_period['end'] = i
                        current_period['min_cpu'] = min(current_period['min_cpu'], 
                                                       sample['cpu']['overall_percent'])
                else:
                    if current_period and (current_period['end'] - current_period['start']) > 5:
                        low_cpu_periods.append(current_period)
                    current_period = None
            
            if low_cpu_periods:
                f.write(f"Detected {len(low_cpu_periods)} low CPU utilization periods:\n")
                for period in low_cpu_periods:
                    start_time = self.samples[period['start']]['elapsed_seconds']
                    end_time = self.samples[period['end']]['elapsed_seconds']
                    duration = end_time - start_time
                    f.write(f"  {start_time:6.1f}s - {end_time:6.1f}s ({duration:5.1f}s): "
                           f"Min CPU = {period['min_cpu']:.1f}%\n")
            
            # Memory Analysis
            f.write("\n\nMEMORY USAGE:\n")
            f.write("-"*50 + "\n")
            f.write(f"Average: {sum(memory_values)/len(memory_values):.1f}%\n")
            f.write(f"Peak: {max(memory_values):.1f}%\n")
            
            # Find peak memory processes
            peak_memory_procs = []
            for pid, data in self.process_history.items():
                if data['samples']:
                    peak_mem = max(s['memory'] for s in data['samples'])
                    if peak_mem > 100:  # More than 100MB
                        peak_memory_procs.append({
                            'name': data['name'],
                            'cmdline': data['cmdline'][:100],
                            'peak_mb': peak_mem
                        })
            
            if peak_memory_procs:
                f.write("\nTop Memory Consumers:\n")
                for proc in sorted(peak_memory_procs, key=lambda x: x['peak_mb'], reverse=True)[:10]:
                    f.write(f"  {proc['name']:20s}: {proc['peak_mb']:6.0f} MB\n")
                    f.write(f"    {proc['cmdline']}\n")
        
        # Create timeline visualization
        self.create_timeline()
        
        print(f"\nEnhanced summary written to: {self.summary_file}")
        print(f"Timeline visualization: {self.timeline_file}")
        print(f"Detailed log: {self.cpu_log}")
    
    def create_timeline(self):
        """Create ASCII timeline of CPU usage and phases."""
        with open(self.timeline_file, 'w') as f:
            f.write("CPU USAGE TIMELINE\n")
            f.write("="*100 + "\n\n")
            
            # Create buckets (5-second intervals)
            bucket_size = 5
            n_buckets = int(self.samples[-1]['elapsed_seconds'] / bucket_size) + 1
            
            cpu_buckets = [[] for _ in range(n_buckets)]
            phase_buckets = [set() for _ in range(n_buckets)]
            
            for sample in self.samples:
                bucket_idx = int(sample['elapsed_seconds'] / bucket_size)
                if bucket_idx < n_buckets:
                    cpu_buckets[bucket_idx].append(sample['cpu']['overall_percent'])
                    phase_buckets[bucket_idx].update(sample['phases'])
            
            # Draw timeline
            f.write("Time (s)  CPU Usage                                          Active Phases\n")
            f.write("-"*80 + "\n")
            
            for i in range(n_buckets):
                if cpu_buckets[i]:
                    avg_cpu = sum(cpu_buckets[i]) / len(cpu_buckets[i])
                    bar_length = int(avg_cpu / 2)  # Scale to 50 chars max
                    bar = '#' * bar_length + '-' * (50 - bar_length)
                    phases = ', '.join(phase_buckets[i]) if phase_buckets[i] else 'idle'
                    
                    f.write(f"{i*bucket_size:4d}-{(i+1)*bucket_size:3d}  "
                           f"{bar} {avg_cpu:5.1f}%  {phases}\n")
    
    def stop(self):
        """Stop monitoring."""
        self.running = False


def main():
    parser = argparse.ArgumentParser(description='Enhanced CPU monitoring for GVClass')
    parser.add_argument('command', nargs='+', help='Command to run and monitor')
    parser.add_argument('-o', '--output-dir', default='cpu_monitoring_v2', 
                       help='Output directory for logs')
    parser.add_argument('-i', '--interval', type=float, default=0.5,
                       help='Sampling interval in seconds')
    
    args = parser.parse_args()
    
    # Start monitoring
    monitor = EnhancedCPUMonitor(args.output_dir, args.interval)
    
    # Handle signals
    def signal_handler(sig, frame):
        print("\nStopping monitor...")
        monitor.stop()
        monitor.create_enhanced_summary()
        sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Start monitor thread
    monitor_thread = threading.Thread(target=monitor.monitor)
    monitor_thread.daemon = True
    monitor_thread.start()
    
    print("Starting enhanced CPU monitoring...")
    print(f"Output directory: {monitor.output_dir}")
    print(f"Running command: {' '.join(args.command)}")
    print("-" * 80)
    
    # Run the command
    try:
        process = subprocess.Popen(args.command)
        monitor.target_process = process.pid
        process.wait()
    except Exception as e:
        print(f"Error running command: {e}")
    finally:
        # Stop monitoring
        monitor.stop()
        monitor_thread.join()
        monitor.create_enhanced_summary()


if __name__ == '__main__':
    main()