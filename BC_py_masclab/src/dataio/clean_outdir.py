"""
Output directory cleanup utilities for MASC processing.

This module removes existing output contents before a new processing run.

Functions for:
- Clearing the campaign output directory

Last update : June 2026
Author : Baptiste Carmier
"""

import logging
import shutil

from pathlib import Path

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)


def clean_output_dir(outdir: Path) -> None:
    """
    Remove all contents of the output directory before a processing run.

    Args:
        outdir: Output directory to clean.

    Returns:
        None

    Notes:
        - Skips silently if the directory does not exist.
        - Logs a warning for files or folders that cannot be removed.
    """
    if not outdir.exists():
        return
    for child in outdir.iterdir():
        try:
            if child.is_file() or child.is_symlink():
                child.unlink()
            elif child.is_dir():
                shutil.rmtree(child)
        except Exception:
            logger.warning("Could not remove '%s' while cleaning %s", child, outdir)
