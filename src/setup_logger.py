#!/usr/bin/env python

"""Logger writes to STDERR and optionally to a LOGFILE.

"""

from __future__ import annotations
from typing import Optional
import sys
from pathlib import Path
from loguru import logger


def formatter(record):
    """Custom formatter that allows for progress bar."""
    fmessage = (
        "{time:YYYY-MM-DD HH:mm:ss} | "
        "<level>{level:<8}</level> <white>|</white> "
        "<magenta>{file:<20}</magenta> <white>|</white> "
        "{message}"
    )
    return fmessage


def color_support():
    """Check for color support in stderr as a notebook or terminal/tty."""
    return sys.stderr.isatty()


def normalize_log_level(level: str) -> str:
    """Convert level input to full string if it is a substring, else return user value"""
    OPTIONS = ("TRACE", "DEBUG", "INFO", "WARNING", "ERROR")
    level = level.upper()
    for lvl in OPTIONS:
        if lvl.startswith(level):
            return lvl
    return level


def set_log_level(log_level: str = "DEBUG", log_file: Optional[Path] = None):
    """Add logger for ipyrad to stderr and optionally to file.
    """
    logger.remove()

    # always log to stderr
    logger.add(
        sink=sys.stderr,
        level=normalize_log_level(log_level),
        colorize=color_support(),
        format=formatter,
        enqueue=False,
        # traceback=True,
    )
    # optionally log to file
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(exist_ok=True)
        log_file.touch(exist_ok=True)
        logger.add(
            sink=str(log_file),
            level=normalize_log_level(log_level),
            colorize=False,
            format=formatter,
            enqueue=True,
            rotation="50 MB",
        )
    return logger


def setup_loguru_worker(log_level: str) -> None:
    """initialized on parallel Worker processes."""
    from loguru import logger
    import sys

    logger.remove()
    logger.add(
        sys.stderr,
        level=normalize_log_level(log_level),
        colorize=color_support(),
        format=formatter,
        enqueue=True,
    )


if __name__ == "__main__":
    pass
