"""
Standalone runner for blowing-snow classification.

This script extracts features from raw MASC images, classifies each image
as precipitation or blowing snow, and writes CSV summary tables to outdir.

Usage:
- Configure label and [blowing_snow] sections in config/config.yaml
- Run from BC_py_masclab/: python scripts/runner_bs.py
- Outputs: blowing_snow_all.csv and blowing_snow_triplet.csv in label.outdir

Translated from run_bs.m (Christophe Praz 2018) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import sys
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root_dir))

from src.core.config import load_from_yaml
from src.pipeline.steps import run_blowing_snow


def main() -> int:
    """
    Extract blowing-snow features and classify precipitation vs blowing snow.

    Args:
        None

    Returns:
        Exit code (0 on success, 1 on failure).

    Notes:
        - Input folder and time window default to label.campaigndir / starthr_vec / endhr_vec.
        - GMM model and block parameters come from [blowing_snow] in config.yaml.
    """
    yaml_path = root_dir / "config" / "config.yaml"
    label_cfg, proc_cfg = load_from_yaml(yaml_path)
    try:
        run_blowing_snow(label_cfg, proc_cfg, yaml_path=yaml_path)
        return 0
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
