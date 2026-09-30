"""
Configuration management for the MASC processing pipeline.

This module handles YAML loading, dataclass configs, and proc_params.txt export.

Functions for:
- Loading label, processing, and pipeline settings from config.yaml
- Writing proc_params.txt metadata beside processed data

Translated from create_proc_params_file.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Tuple, List, Optional

try:
    import yaml  # For YAML configuration
    _YAML_AVAILABLE = True
except Exception:
    _YAML_AVAILABLE = False

@dataclass
class ProcessingConfig:
    """Processing configuration parameters (loaded from YAML only)."""
    
    # Processing control
    parallel: bool 
    use_triplet_algo: bool 
    max_workers: int 

    generate_figs: bool 
    display_figs: bool 
    save_figs: bool 

    saveresults: bool 
    save_format: str
    clean_outdir: bool 
    
    flakebrighten: bool 

    # Thresholds
    backthresh: int 
    sizemin: int 
    min_area: int 
    minbright: float 
    max_intensthresh: float 
    min_hole_area: int 
    
    # Image margins to discard [top, bottom, left, right]
    discardmat: Tuple[int, int, int, int]
    
    # Triplet processing
    matching_tol_pix: int
    matching_tol_percent: int
    camera_order: List[int]
    triplet_required_views: int 
    triplet_allow_partial: bool 


@dataclass
class LabelConfig:
    """Labeling and data organization parameters (loaded from YAML only)."""
    # Data paths (required)
    campaigndir: Path
    starthr_vec: datetime
    endhr_vec: datetime
    outdir: Path


@dataclass
class PipelineConfig:
    """Orchestration flags for runner_joblib (post-processing steps)."""

    run_classification: bool = False
    run_quicklooks: bool = False
    run_blowing_snow: bool = False


def _project_root() -> Path:
    """Return the project root (folder containing 'config/')."""
    # src/core/config.py -> src/core -> src -> project root
    return Path(__file__).resolve().parents[2]


def load_from_yaml(yaml_path: Optional[Path] = None) -> Tuple[LabelConfig, ProcessingConfig]:
    """Load LabelConfig and ProcessingConfig from a YAML file (required).

    All configuration must come from YAML - no defaults.

    Args:
        yaml_path: Optional path to YAML file. If None, uses project_root/config/config.yaml

    Returns:
        (label_config, process_config)
    
    Raises:
        FileNotFoundError: If YAML file doesn't exist
        ImportError: If PyYAML is not installed
        KeyError: If required YAML sections or keys are missing
    """
    # Resolve default path if none provided
    if yaml_path is None:
        yaml_path = _project_root() / 'config' / 'config.yaml'

    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML config not found: {yaml_path}. Please create it and retry.")

    if not _YAML_AVAILABLE:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")

    with open(yaml_path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f) or {}

    # Extract sections - must exist
    label_section = data.get('label', {}) or {}
    proc_section = data.get('processing', {}) or {}

    # Process label section
    label_data = {}
    for key in ['campaigndir', 'starthr_vec', 'endhr_vec', 'outdir']:
        if key not in label_section:
            raise KeyError(f"Missing required label config key: '{key}'")
        
        value = label_section[key]
        
        # Handle outdir special case:
        #   "auto" → campaigndir.parent / f"{campaigndir.name}_PROCESSED"
        #     (keeps the processed folder alongside the input campaign directory)
        if key == 'outdir' and isinstance(value, str) and value.lower() == 'auto':
            campaigndir = Path(label_section.get('campaigndir', '.'))
            value = campaigndir.parent / f"{campaigndir.name}_PROCESSED"
        elif key in {'campaigndir', 'outdir'} and isinstance(value, str):
            value = Path(value)
        elif key in {'starthr_vec', 'endhr_vec'} and isinstance(value, str):
            # "all" → use min/max datetime to cover everything
            if value.lower() == 'all':
                value = datetime.min if key == 'starthr_vec' else datetime.max
            else:
                try:
                    value = datetime.fromisoformat(value)
                except ValueError:
                    try:
                        value = datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
                    except ValueError:
                        raise ValueError(f"Unrecognized datetime format for '{key}': {value}")
        label_data[key] = value

    # Process processing section
    proc_data = {}
    for key in proc_section.keys():
        proc_data[key] = proc_section[key]

    # Create config objects from YAML data
    label_config = LabelConfig(**label_data)
    process_config = ProcessingConfig(**proc_data)

    return label_config, process_config


def load_yaml_section(yaml_path: Path, section: str) -> dict:
    """
    Load a top-level YAML section as a dictionary.

    Args:
        yaml_path: Path to config.yaml.
        section: Section name (e.g. 'classification', 'quicklooks').

    Returns:
        Section contents as a dict, or an empty dict if the section is missing.

    Raises:
        ImportError: If PyYAML is not installed.
    """
    if not _YAML_AVAILABLE:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    block = data.get(section)
    return block if isinstance(block, dict) else {}


def load_pipeline_from_yaml(yaml_path: Optional[Path] = None) -> PipelineConfig:
    """
    Load post-processing orchestration flags from the [pipeline] YAML section.

    Args:
        yaml_path: Path to config.yaml; defaults to project config/config.yaml.

    Returns:
        PipelineConfig with run_classification, run_quicklooks, run_blowing_snow.

    Raises:
        FileNotFoundError: If the YAML file does not exist.
    """
    if yaml_path is None:
        yaml_path = _project_root() / "config" / "config.yaml"
    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML config not found: {yaml_path}")
    section = load_yaml_section(yaml_path, "pipeline")
    return PipelineConfig(
        run_classification=bool(section.get("run_classification", False)),
        run_quicklooks=bool(section.get("run_quicklooks", False)),
        run_blowing_snow=bool(section.get("run_blowing_snow", False)),
    )


def print_config_summary(
    label: LabelConfig,
    process: ProcessingConfig,
    pipeline: Optional[PipelineConfig] = None,
) -> None:
    """
    Print a concise summary of label, processing, and pipeline settings.

    Args:
        label: Label configuration.
        process: Processing configuration.
        pipeline: Optional pipeline orchestration flags.

    Returns:
        None
    """
    print("Processing Configuration:")
    print(f"  Mode: {'Triplet' if process.use_triplet_algo else 'Single-image'}")
    print(f"  Parallel: {process.parallel}")
    if process.parallel:
        print(f"  Max workers: {process.max_workers}")
    print(f"  Generate figures: {process.generate_figs}")
    if process.generate_figs:
        print(f"    - Display: {process.display_figs}")
        print(f"    - Save: {process.save_figs}")
    print(f"  Save results: {process.saveresults}")
    print()
    if process.use_triplet_algo:
        print("Triplet Parameters:")
        print(f"  Camera order: {process.camera_order}")
        print(f"  Required views: {process.triplet_required_views}")
        print(f"  Allow partial views: {process.triplet_allow_partial}")
        print(f"  Matching tolerance (pixels): {process.matching_tol_pix}")
        print(f"  Matching tolerance (percent): {process.matching_tol_percent}")
        print()
    print("Data Configuration:")
    print(f"  Campaign directory: {label.campaigndir}")
    print(f"  Time window: {label.starthr_vec} to {label.endhr_vec}")
    print(f"  Output directory: {label.outdir}")
    if pipeline is not None:
        print()
        print("Pipeline steps (runner_joblib):")
        print(f"  Run classification: {pipeline.run_classification}")
        print(f"  Run quicklooks:     {pipeline.run_quicklooks}")
        print(f"  Run blowing snow:   {pipeline.run_blowing_snow}")
    print()

def create_proc_params_file(pathname: Path, label: LabelConfig, process: ProcessingConfig) -> None: 
    """
    Write proc_params.txt beside processed ROI data.

    Args:
        pathname: Directory where proc_params.txt will be saved.
        label: Label configuration (campaign paths and time window).
        process: Processing configuration (thresholds and flags).

    Returns:
        None

    Translated from create_proc_params_file.m (Christophe Praz 2015) and adapted for Python
    """
    filename = Path(pathname) / 'proc_params.txt'
    filename.parent.mkdir(parents=True, exist_ok=True) # Ensure directory exists


    with open(filename, 'w') as f:
        f.write('Processing parameters associated with data in this folder\n\n\n')
        
        # Write timestamps and paths
        f.write(f'creation date          : {datetime.now().strftime("%Y-%m-%d")}\n')
        f.write(f'data processed         : {label.campaigndir}\n')
        f.write(f'starting time          : {label.starthr_vec.strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'ending time            : {label.endhr_vec.strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'processed in parallel  : {int(process.parallel)}\n\n')
        
        # Write processing parameters
        f.write(f'backtresh limit        : {process.backthresh}\n')
        f.write(f'size min               : {process.sizemin}\n')
        f.write(f'min area               : {process.min_area}\n')
        f.write(f'min brightness         : {process.minbright:.2f}\n')
        f.write(f'max intens thresh      : {process.max_intensthresh:.2f}\n')
        f.write(f'min hole area          : {process.min_hole_area}\n')
        
        # Write discard margins
        margins = process.discardmat
        f.write(f'discardmat [t,b,l,r]   : [{margins[0]} {margins[1]} {margins[2]} {margins[3]}]\n\n')
        
        # Write triplet parameters
        f.write(f'use triplet algo       : {int(process.use_triplet_algo)}\n')
        if process.use_triplet_algo:
            f.write(f'triplet required views  : {process.triplet_required_views}\n')
            f.write(f'triplet allow partial   : {int(process.triplet_allow_partial)}\n')
            f.write(f'matching tol. pix.     : {process.matching_tol_pix}\n')
            f.write(f'matching tol. percent. : {process.matching_tol_percent}\n')
            f.write(f'camera order           : {process.camera_order}\n')