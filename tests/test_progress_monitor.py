"""Progress-event rendering tests."""

from __future__ import annotations

import json
import os
from pathlib import Path

from src.pipeline.progress_events import PROGRESS_EVENT_PREFIX


def test_resource_monitor_uses_progress_events_not_output_files(
    tmp_path: Path, capsys
) -> None:
    from src.bin.progress_monitor import ResourceMonitor

    output_dir = tmp_path / "out"
    output_dir.mkdir()
    # A stale summary file should not initialize or advance progress.
    (output_dir / "stale.summary.tab").write_text("query\nstale\n")

    monitor = ResourceMonitor(output_dir, total_queries=2, process_pid=os.getpid())
    assert monitor.current_query == 0

    event = {
        "event": "query_progress",
        "query": "q1",
        "progress": 50,
        "stage": "building_marker_trees",
        "completed_markers": 2,
        "total_markers": 4,
    }
    assert monitor.handle_progress_line(PROGRESS_EVENT_PREFIX + json.dumps(event))

    captured = capsys.readouterr().out
    assert "Progress: [" in captured
    assert "25%" in captured
    assert "0/2 queries" in captured
    assert "2/4" in captured


def test_resource_monitor_counts_failed_query_as_done(tmp_path: Path, capsys) -> None:
    from src.bin.progress_monitor import ResourceMonitor

    monitor = ResourceMonitor(tmp_path, total_queries=1, process_pid=os.getpid())
    monitor.apply_event(
        {
            "event": "query_failed",
            "query": "q1",
            "progress": 100,
            "stage": "failed",
        }
    )

    captured = capsys.readouterr().out
    assert "100%" in captured
    assert "1/1 queries" in captured
