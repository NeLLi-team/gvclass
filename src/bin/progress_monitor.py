"""Runtime progress and memory monitoring helpers for the GVClass CLI."""

from __future__ import annotations

import time
from pathlib import Path


class ResourceMonitor:
    """Monitor pipeline progress and memory usage."""

    def __init__(self, output_dir: Path, total_queries: int, process_pid: int):
        self.output_dir = output_dir
        self.total_queries = total_queries
        self.process_pid = process_pid
        self.running = True
        self.current_query = 0
        self.current_memory_mb = 0.0
        self.peak_memory_mb = 0.0
        self._last_output = ""
        self._process = self._load_process()

    def _load_process(self):
        try:
            import psutil
        except ModuleNotFoundError:
            self.psutil = None
            return None

        self.psutil = psutil
        try:
            return psutil.Process(self.process_pid)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return None

    def _update_progress(self) -> None:
        if not self.output_dir.exists():
            return
        completed = len(list(self.output_dir.glob("*.summary.tab")))
        self.current_query = min(completed, self.total_queries)

    def _update_memory(self) -> None:
        if not self._process or not self.psutil:
            return

        try:
            total_memory = self._process.memory_info().rss / 1024 / 1024
            for child in self._process.children(recursive=True):
                try:
                    total_memory += child.memory_info().rss / 1024 / 1024
                except (self.psutil.NoSuchProcess, self.psutil.AccessDenied):
                    continue
            self.current_memory_mb = total_memory
            self.peak_memory_mb = max(self.peak_memory_mb, total_memory)
        except (self.psutil.NoSuchProcess, self.psutil.AccessDenied):
            return

    def _render_progress(self) -> None:
        if self.total_queries <= 0:
            return

        percent = int((self.current_query / self.total_queries) * 100)
        output = (
            f"Progress: [{percent:3d}%] Query {self.current_query}/{self.total_queries} "
            f"| Memory: {self.current_memory_mb:.0f}MB (peak: {self.peak_memory_mb:.0f}MB)"
        )
        if output != self._last_output:
            print(f"\r{output:<80}", end="", flush=True)
            self._last_output = output

    def monitor(self) -> None:
        while self.running:
            self._update_progress()
            self._update_memory()
            self._render_progress()
            time.sleep(1)

    def stop(self) -> None:
        self.running = False
