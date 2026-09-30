"""
Standalone runner for MASC daily quicklook figure generation.

This script aggregates classified ROI data into time bins and saves daily
recap PNGs under outdir/Quicklooks/.

Usage:
- Run classification first so label fields exist in ROI files
- Configure quicklooks options in config/config.yaml [quicklooks] section
- Run from BC_py_masclab/: python scripts/runner_quicklooks.py

Translated from MASC_process_classify_quicklooks.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

import sys
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root_dir))

from src.core.config import load_from_yaml
from src.pipeline.steps import run_quicklooks


def main() -> int:
    """
    Generate daily quicklook PNGs from classified ROI files.

    Args:
        None

    Returns:
        Exit code (0 on success, 1 on failure).

    Notes:
        - Requires matplotlib and existing ROI files with classification labels.
        - Output layout: Quicklooks/<yyyy.mm.dd>/yyyymmdd_quicklook#1.png (#2).
    """
    yaml_path = root_dir / "config" / "config.yaml"
    label_cfg, proc_cfg = load_from_yaml(yaml_path)
    try:
        run_quicklooks(label_cfg, proc_cfg, yaml_path=yaml_path)
        return 0
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
