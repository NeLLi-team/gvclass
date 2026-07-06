"""Runtime progress and memory monitoring helpers for the GVClass CLI."""

from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict

from src.pipeline.progress_events import PROGRESS_EVENT_PREFIX


class ResourceMonitor:
    """Track child-process progress events and peak memory usage."""

    def __init__(self, output_dir: Path, total_queries: int, process_pid: int):
        self.output_dir = output_dir
        self.total_queries = total_queries
        self.process_pid = process_pid
        self.running = True
        self.current_query = 0
        self.current_memory_mb = 0.0
        self.peak_memory_mb = 0.0
        self._last_output = ""
        self._last_percent = -1
        self._query_progress: Dict[str, float] = {}
        self._query_stage: Dict[str, str] = {}
        self._completed_queries: set[str] = set()
        self._failed_queries: set[str] = set()
        self._current_label = "starting"
        self._is_tty = sys.stdout.isatty()
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
        self.current_query = min(
            len(self._completed_queries) + len(self._failed_queries),
            self.total_queries,
        )

    def update_memory(self) -> None:
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

    def handle_progress_line(self, line: str) -> bool:
        if not line.startswith(PROGRESS_EVENT_PREFIX):
            return False
        try:
            event = json.loads(line[len(PROGRESS_EVENT_PREFIX) :])
        except json.JSONDecodeError:
            return True
        self.apply_event(event)
        return True

    def apply_event(self, event: Dict[str, Any]) -> None:
        event_type = event.get("event")
        query_name = str(event.get("query") or "")
        if not query_name:
            return

        progress = float(event.get("progress", 0))
        progress = max(self._query_progress.get(query_name, 0.0), progress)
        self._query_progress[query_name] = min(progress, 100.0)
        stage = str(event.get("stage") or event_type or "working")
        self._query_stage[query_name] = stage

        if event_type == "query_failed":
            self._failed_queries.add(query_name)
        elif progress >= 100:
            self._completed_queries.add(query_name)

        self._current_label = self._format_event_label(query_name, event)
        self._update_progress()
        self.render_progress()

    def render_progress(self, force: bool = False) -> None:
        if self.total_queries <= 0:
            return

        progress_sum = sum(self._query_progress.values())
        percent = int(progress_sum / max(1, self.total_queries))
        percent = max(0, min(100, percent))
        if not force and percent == self._last_percent and not self._is_tty:
            return

        bar_width = self._bar_width()
        filled = int(bar_width * percent / 100)
        bar = "#" * filled + "-" * (bar_width - filled)
        output = (
            f"Progress: [{bar}] {percent:3d}% | "
            f"{self.current_query}/{self.total_queries} queries | {self._current_label}"
        )
        if output != self._last_output:
            if self._is_tty:
                print(f"\r{output}\x1b[K", end="", flush=True)
            elif force or percent >= self._last_percent + 10 or percent == 100:
                print(output, flush=True)
            self._last_output = output
            self._last_percent = percent

    def finish(self, success: bool) -> None:
        if success:
            for query_name in self._query_progress:
                self._query_progress[query_name] = 100.0
            missing = self.total_queries - len(self._query_progress)
            for index in range(max(0, missing)):
                self._query_progress[f"__complete_{index}"] = 100.0
            self.current_query = self.total_queries
            self._current_label = "complete"
        self.render_progress(force=True)
        if self._last_output:
            print()

    def monitor(self) -> None:
        while self.running:
            self._update_progress()
            self.update_memory()
            self.render_progress()
            time.sleep(1)

    def monitor_memory(self) -> None:
        while self.running:
            self.update_memory()
            time.sleep(1)

    def stop(self) -> None:
        self.running = False

    def _bar_width(self) -> int:
        terminal_width = shutil.get_terminal_size((100, 20)).columns
        return max(20, min(40, terminal_width - 70))

    def _format_event_label(self, query_name: str, event: Dict[str, Any]) -> str:
        stage = str(event.get("stage") or event.get("event") or "working")
        stage = stage.replace("_", " ")
        marker = event.get("marker")
        completed_markers = event.get("completed_markers")
        total_markers = event.get("total_markers")
        query_label = query_name if len(query_name) <= 28 else query_name[:25] + "..."
        if completed_markers is not None and total_markers:
            return f"{query_label}: {stage} {completed_markers}/{total_markers}"
        if marker:
            return f"{query_label}: {stage} {marker}"
        return f"{query_label}: {stage}"
