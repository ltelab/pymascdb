"""
Standalone runner for MASC hydrometeor classification.

This script runs the three trained classifiers (class, riming, melting) on
processed ROI files and writes label fields back into each file.

Usage:
- Ensure ROI files already exist under label.outdir (run processing first)
- Configure classifiers in config/config.yaml [classification] section
- Run from BC_py_masclab/: python scripts/runner_classify.py

Translated from make_predictions_for_campaign.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import sys
import logging
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root_dir))

from src.core.config import load_from_yaml
from src.pipeline.steps import run_classification

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")


def main() -> int:
    """
    Run class, riming, and melting classifiers on all ROI files in outdir.

    Args:
        None

    Returns:
        Exit code (0 on success, 1 on failure).

    Notes:
        - Writes label_*, riming_*, and melting_* fields into each ROI file.
        - Requires config/config.yaml and existing processed ROI data.
    """
    yaml_path = root_dir / "config" / "config.yaml"
    label_cfg, proc_cfg = load_from_yaml(yaml_path)
    try:
        run_classification(label_cfg, proc_cfg, yaml_path=yaml_path, project_root=root_dir)
        return 0
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
