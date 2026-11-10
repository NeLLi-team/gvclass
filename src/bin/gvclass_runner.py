#!/usr/bin/env python
"""Runner entry point: run the pipeline from the ``run_config.json`` of the CLI.

The command line resolves every option into a ``RunConfig`` and writes it to the
output directory; this module reads that file, sets up logging and calls the flow.
It defines no option of the pipeline.
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

from src.pipeline.parallel_runner import gvclass_flow
from src.pipeline.run_config import RunConfig

LOG_FILENAME = "gvclass.log"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logger = logging.getLogger("gvclass_runner")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the command line of the runner."""
    parser = argparse.ArgumentParser(
        description="Run the GVClass pipeline from a run_config.json file."
    )
    parser.add_argument(
        "config_file", type=Path, help="Run configuration written by the command line"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Log debug messages"
    )
    return parser.parse_args(argv)


def configure_logging(output_dir: Path, *, resume: bool, verbose: bool) -> None:
    """Log to ``<output>/gvclass.log`` and to stderr.

    A resumed run appends to the log of the earlier run. The pipeline loggers
    (``gvclass_runner`` and ``src.*``) get the level; other libraries keep WARNING.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [
        logging.FileHandler(
            output_dir / LOG_FILENAME, mode="a" if resume else "w", encoding="utf-8"
        ),
        logging.StreamHandler(),
    ]
    root = logging.getLogger()
    for handler in handlers:
        handler.setFormatter(logging.Formatter(LOG_FORMAT))
        root.addHandler(handler)
    root.setLevel(logging.WARNING)
    level = logging.DEBUG if verbose else logging.INFO
    for name in ("gvclass_runner", "src"):
        logging.getLogger(name).setLevel(level)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the pipeline for the configuration in the file; return the exit status."""
    args = parse_args(argv)
    config = RunConfig.from_json(args.config_file.read_text(encoding="utf-8"))
    configure_logging(config.output_dir, resume=config.resume, verbose=args.verbose)
    try:
        summary_file = gvclass_flow(config)
    except Exception:
        logger.exception("Pipeline failed")
        return 1
    logger.info("Results written to: %s", summary_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
