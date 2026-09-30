"""
Main entry point for the MASC processing pipeline with joblib parallelism.

This script loads config.yaml, runs batch ROI processing, and optionally
executes classification, quicklooks, and blowing-snow post-processing steps.

Optional post-processing (config [pipeline] section):
  run_classification — classify ROI files after processing
  run_quicklooks     — daily recap PNGs (requires classification labels)

Usage:
- Configure parameters in config/config.yaml
- Run: python scripts/runner_joblib.py
- Enable optional steps in the [pipeline] section (classification, quicklooks, blowing_snow)

Translated from MASC_process.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""
import sys
import logging
from pathlib import Path

root = logging.getLogger()
for handler in root.handlers[:]:
    root.removeHandler(handler)
logging.basicConfig(level=logging.ERROR, force=True)

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.core.process_joblib import masc_process
from src.core.config import load_from_yaml, load_pipeline_from_yaml, print_config_summary
from src.dataio.clean_outdir import clean_output_dir
from src.pipeline.steps import run_blowing_snow, run_classification, run_quicklooks


def main():
    """
    Run the MASC processing pipeline and optional post-processing steps.

    Args:
        None

    Returns:
        Exit code (0 on success, 1 on failure).

    Notes:
        - Reads config/config.yaml from the project root.
        - Optionally cleans outdir before processing when configured.
    """

    print("=" * 70)
    print("MASC Processing Pipeline - Joblib Parallel Runner")
    print("=" * 70)
    print()

    yaml_path = Path(__file__).resolve().parent.parent / "config" / "config.yaml"
    label_config, process_config = load_from_yaml(yaml_path)
    pipeline_config = load_pipeline_from_yaml(yaml_path)

    if process_config.saveresults and process_config.clean_outdir:
        clean_output_dir(label_config.outdir)

    print_config_summary(label_config, process_config, pipeline_config)

    print("-" * 70)
    print("Starting processing with joblib parallelism...")
    print("-" * 70)
    print()

    try:
        masc_process(label_config, process_config)

        print()
        print("-" * 70)
        print("Processing completed successfully!")
        print(f"Results saved to: {label_config.outdir}")
        if process_config.saveresults:
            print(f"Statistics file: {label_config.outdir / 'proc_stats.txt'}")
        print("-" * 70)

    except FileNotFoundError as e:
        print()
        print(f"ERROR: File not found - {e}")
        print("Please check your campaign directory path in config.yaml")
        return 1

    except NotImplementedError as e:
        print()
        print(f"ERROR: Feature not implemented - {e}")
        return 1

    except Exception as e:
        print()
        print("ERROR: Processing failed with exception:")
        print(f"  {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return 1

    if pipeline_config.run_classification:
        print()
        try:
            run_classification(
                label_config,
                process_config,
                yaml_path=yaml_path,
                project_root=yaml_path.parent.parent,
            )
        except Exception as e:
            print()
            print("ERROR: Classification failed:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            return 1

    if pipeline_config.run_quicklooks:
        print()
        try:
            run_quicklooks(label_config, process_config, yaml_path=yaml_path)
        except Exception as e:
            print()
            print("ERROR: Quicklooks failed:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            return 1

    if pipeline_config.run_blowing_snow:
        print()
        try:
            run_blowing_snow(label_config, process_config, yaml_path=yaml_path)
        except Exception as e:
            print()
            print("ERROR: Blowing snow failed:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
