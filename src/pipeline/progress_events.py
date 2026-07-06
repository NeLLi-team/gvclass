"""Line-delimited progress events for the parent CLI renderer."""

from __future__ import annotations

import json
import os
import sys
import threading
from datetime import datetime, timezone
from typing import Any, Dict

PROGRESS_EVENT_PREFIX = "GVCLASS_PROGRESS "

_EVENT_LOCK = threading.Lock()


def progress_events_enabled() -> bool:
    return os.environ.get("GVCLASS_PROGRESS_EVENTS") == "1"


def emit_progress_event(event: Dict[str, Any]) -> None:
    """Emit one JSON progress event to stderr when parent-side rendering is on."""
    if not progress_events_enabled():
        return

    payload = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        **event,
    }
    line = PROGRESS_EVENT_PREFIX + json.dumps(payload, sort_keys=True)
    with _EVENT_LOCK:
        print(line, file=sys.stderr, flush=True)
