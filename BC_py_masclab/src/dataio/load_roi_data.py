"""
ROI data persistence utilities for MASC processing outputs.

This module saves and loads ROI feature dictionaries in joblib, pickle,
or MATLAB-compatible .mat formats.

Functions for:
- Saving and loading ROI files
- Batch folder loading and scalar-field export

Last update : June 2026
Author : Baptiste Carmier
"""

import pickle
import numpy as np
import logging

from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

try:
    import joblib # for .joblib format
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False

try:
    import scipy.io # for .mat format
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def save_roi_data(
    roi: Dict,
    filepath: Path,
    format: str = "joblib",
    compress: Union[bool, int] = 3
) -> None:
    """
    Save a ROI dictionary to disk in the requested format.

    Args:
        roi: ROI dictionary with features and metadata.
        filepath: Output file path (extension should match format).
        format: One of 'joblib', 'pkl', or 'mat'.
        compress: joblib compression level (0–9) or bool for other formats.

    Returns:
        None

    Notes:
        - Matplotlib figure handles (keys starting with '_figure') are stripped.

    Raises:
        ImportError: If joblib or scipy is required but not installed.
        ValueError: If format is not supported.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    # Remove non-serializable items (matplotlib figures)
    roi_clean = {k: v for k, v in roi.items() if not k.startswith('_figure')}
    
    if format == "joblib":
        if not JOBLIB_AVAILABLE:
            raise ImportError("joblib not installed. Install with: pip install joblib")
        
        # joblib is optimized for numpy arrays with compression
        joblib.dump(roi_clean, filepath, compress=compress)
        
    elif format == "pkl":
        # Standard pickle
        with open(filepath, 'wb') as f:
            pickle.dump(roi_clean, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    elif format == "mat":
        if not SCIPY_AVAILABLE:
            raise ImportError("scipy not installed. Install with: pip install scipy")
        
        # Convert to MATLAB-compatible format
        mat_data = _convert_roi_to_matlab(roi_clean)
        scipy.io.savemat(filepath, mat_data)
    
    else:
        raise ValueError(f"Unknown format: {format}. Use 'joblib', 'pkl', or 'mat' in config.py")


def _mat_struct_to_dict(obj):
    """Recursively convert scipy mat_struct objects to plain Python dicts."""
    if hasattr(obj, '_fieldnames'):
        return {f: _mat_struct_to_dict(getattr(obj, f)) for f in obj._fieldnames}
    elif isinstance(obj, dict):
        return {k: _mat_struct_to_dict(v) for k, v in obj.items()}
    elif isinstance(obj, np.ndarray) and obj.dtype.names:
        # structured array → dict
        return {name: _mat_struct_to_dict(obj[name]) for name in obj.dtype.names}
    else:
        return obj


def load_roi_data(filepath: Path, format: Optional[str] = None) -> Dict:
    """
    Load a ROI dictionary from disk.

    Args:
        filepath: Path to the ROI file.
        format: File format; auto-detected from extension when None.

    Returns:
        ROI dictionary.

    Raises:
        FileNotFoundError: If the file does not exist.
        ImportError: If joblib or scipy is required but not installed.
        ValueError: If format is not supported.
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Auto-detect format from extension
    if format is None:
        ext = filepath.suffix.lower()
        if ext in ['.joblib', '.jl']:
            format = 'joblib'
        elif ext == '.pkl':
            format = 'pkl'
        elif ext == '.mat':
            format = 'mat'
        else:
            # Try joblib first (most common), fallback to pickle
            format = 'joblib'
    
    if format == "joblib":
        if not JOBLIB_AVAILABLE:
            raise ImportError("joblib not installed. Install with: pip install joblib")
        return joblib.load(filepath)
    
    elif format == "pkl":
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    
    elif format == "mat":
        if not SCIPY_AVAILABLE:
            raise ImportError("scipy not installed. Install with: pip install scipy")
        mat_data = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
        # Remove MATLAB metadata keys
        result = {k: v for k, v in mat_data.items() if not k.startswith('__')}
        # Unwrap the 'roi' envelope added by _convert_roi_to_matlab
        if list(result.keys()) == ['roi']:
            inner = result['roi']
            if hasattr(inner, '_fieldnames'):
                result = {f: getattr(inner, f) for f in inner._fieldnames}
            elif isinstance(inner, dict):
                result = inner
        # Recursively convert any remaining mat_struct objects to plain dicts
        return _mat_struct_to_dict(result)
    
    else:
        raise ValueError(f"Unknown format: {format}")

def _convert_roi_to_matlab(roi: Dict) -> Dict:
    """Convert ROI dict to MATLAB-compatible format wrapped in single 'roi' structure."""
    mat_data = {}
    
    for key, value in roi.items():
        if isinstance(value, dict):
            # Nested dict -> struct in MATLAB
            mat_data[key] = _convert_roi_to_matlab_nested(value)
        elif isinstance(value, np.ndarray):
            mat_data[key] = value
        elif isinstance(value, (list, tuple)):
            mat_data[key] = np.array(value)
        elif value is None:
            mat_data[key] = np.nan
        else:
            mat_data[key] = value
    
    # Wrap entire structure in single 'roi' key for MATLAB compatibility
    return {'roi': mat_data}


def _convert_roi_to_matlab_nested(roi: Dict) -> Dict:
    """Recursive helper for nested dicts (for sub-structures like E, E_in, E_out)."""
    mat_data = {}
    
    for key, value in roi.items():
        if isinstance(value, dict):
            mat_data[key] = _convert_roi_to_matlab_nested(value)
        elif isinstance(value, np.ndarray):
            mat_data[key] = value
        elif isinstance(value, (list, tuple)):
            mat_data[key] = np.array(value)
        elif value is None:
            mat_data[key] = np.nan
        else:
            mat_data[key] = value
    
    return mat_data


def load_all_roi_in_folder(
    folder_path: Path,
    quality: Optional[str] = None,
    format: str = "joblib"
) -> List[Dict]:
    """
    Load all ROI files from a folder, optionally filtered by quality.

    Args:
        folder_path: Path to DATA/ or DATA/GOOD or DATA/BAD.
        quality: Filter by 'GOOD' or 'BAD', or None for all files in folder_path.
        format: File format extension to load.

    Returns:
        List of ROI dictionaries.

    Raises:
        FileNotFoundError: If the target folder does not exist.
        ValueError: If format is not supported.
    """
    folder_path = Path(folder_path)
    
    if quality:
        folder_path = folder_path / quality
    
    if not folder_path.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    
    # Find all files with correct extension
    if format == "joblib":
        pattern = "*.joblib"
    elif format == "pkl":
        pattern = "*.pkl"
    elif format == "mat":
        pattern = "*.mat"
    else:
        raise ValueError(f"Unknown format: {format}")
    
    roi_files = sorted(folder_path.glob(pattern))
    
    rois = []
    for file in roi_files:
        try:
            roi = load_roi_data(file, format=format)
            rois.append(roi)
        except Exception as e:
            logger.warning("Could not load %s: %s", file, e)
    
    return rois


def export_roi_to_dict(roi: Dict, fields: Optional[List[str]] = None) -> Dict:
    """
    Export selected scalar fields from a ROI dictionary.

    Args:
        roi: ROI dictionary.
        fields: Field names to export; auto-selects scalar fields when None.

    Returns:
        Dictionary containing only the requested fields.
    """
    if fields is None:
        # Auto-select scalar fields (exclude arrays and dicts)
        fields = [k for k, v in roi.items() 
                  if not isinstance(v, (dict, np.ndarray, list)) or k == 'name']
    
    result = {}
    for field in fields:
        if field in roi:
            value = roi[field]
            # Convert numpy scalars to Python types
            if isinstance(value, np.generic):
                value = value.item()
            result[field] = value
    
    return result


def print_roi_summary(roi: Dict) -> None:
    """
    Print a human-readable summary of key ROI features.

    Args:
        roi: ROI dictionary loaded from disk.

    Returns:
        None
    """
    print("=" * 60)
    print("ROI Summary")
    print("=" * 60)
    
    # Basic info
    print(f"Filename: {roi.get('name', 'N/A')}")
    print(f"Camera: {roi.get('cam', 'N/A')}")
    print(f"Flake ID: {roi.get('id', 'N/A')}")
    print(f"Timestamp: {roi.get('tnum', 'N/A')}")
    print()
    
    # Geometric features
    print("Geometric Features:")
    for key in ['Area', 'Perimeter', 'width', 'height', 'Dmax', 'AspectRatio', 'Rectangularity', 'Nb_holes']:
        if key in roi:
            val = roi[key]
            if isinstance(val, (int, float, np.number)):
                print(f"  {key:20s}: {val:.2f}")
    print()
    
    # Intensity features
    print("Intensity Features:")
    for key in ['mean_intens', 'max_intens']:
        if key in roi:
            print(f"  {key:20s}: {roi[key]:.4f}")
    print()
    
    # Ellipse fit
    if 'E' in roi and isinstance(roi['E'], dict):
        print("Ellipse Fit:")
        print(f"  Semi-major axis (a): {roi['E'].get('a', 'N/A')}")
        print(f"  Semi-minor axis (b): {roi['E'].get('b', 'N/A')}")
        print(f"  Orientation (theta): {roi['E'].get('theta', 'N/A')}°")
        if 'EllipseAreaRatio' in roi:
            print(f"  Ellipse Area Ratio: {roi['EllipseAreaRatio']:.4f}")
        print()
    
    # Quality metrics
    print("Quality Metrics:")
    for key in ['Complexity', 'Compactness', 'Convexity', 'Solidity']:
        if key in roi:
            print(f"  {key:20s}: {roi[key]:.4f}")
    print()
    
    print("=" * 60)


if __name__ == "__main__":
    """Example usage from command line"""
    import sys
    
    if len(sys.argv) > 1:
        pkl_file = Path(sys.argv[1])
        if pkl_file.exists():
            roi = load_roi_data(pkl_file)
            print_roi_summary(roi)
        else:
            logger.error("File not found: %s", pkl_file)    
    else:
        print("Usage: python load_roi_data.py <path_to_file>")
        print("Supported formats: .joblib, .pkl, .mat")
        print("Example: python load_roi_data.py output/2015.06.20/09/DATA/GOOD/2015.06.20_09.00.01_cam0.joblib")
