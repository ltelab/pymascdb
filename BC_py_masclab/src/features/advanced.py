"""
Advanced feature extraction module for MASC analysis.

This module provides functions for extracting advanced descriptors from
snowflake ROIs, including blur indices and back-filling of missing fields.

Functions for:
- Blur index computation
- Conditional back-filling of hull, symmetry, D90, and texture features

Status: STANDALONE UTILITY, not called by the main pipeline. Like its
MATLAB counterpart process_new_descriptors.m, it is meant to be run
ad hoc on ROI files produced by older code versions, to back-fill
descriptors that did not exist at the time they were processed.

Translated from process_new_descriptors.m and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from typing import Dict

from src.features.fmeasure import fmeasure
from src.utils.descriptors_heplers import (
    compute_convex_hull,
    compute_Dmax,
    compute_D90,
    compute_symmetry_features,
    compute_blur_index
)


def process_new_descriptors(roi: Dict) -> Dict:
    """
    Add missing advanced descriptors to an existing ROI dictionary.

    Checks for missing fields and computes them when the required inputs
    are present. Useful for updating older ROI files processed with earlier
    code versions.

    Args:
        roi: ROI dictionary from basic processing or loaded from disk.

    Returns:
        Updated ROI dictionary with advanced descriptors filled in when possible.

    Notes:
        - Only computes a field if it is missing and its dependencies exist.
        - May recompute Dmax and D90 from the convex hull when D90 is absent.
    """
    is_modified = False

    # Convex hull (if missing)
    if 'hull' not in roi:
        if 'x' in roi and 'y' in roi and 'perim' in roi:
            roi['hull'] = compute_convex_hull(roi['x'], roi['y'])
            is_modified = True

    # Symmetry features (if missing)
    if 'Sym' not in roi:
        if 'bw_mask_filled' in roi and 'Dmax' in roi and 'eq_radius' in roi:
            roi['Sym'] = compute_symmetry_features(
                roi['bw_mask_filled'],
                roi['Dmax'],
                roi['eq_radius']
            )
            is_modified = True

    # Wavelet sharpness (if missing)
    if 'wavs' not in roi:
        if 'data' in roi:
            roi['wavs'] = fmeasure(roi['data'], 'WAVS', None)
            is_modified = True

    # Histogram entropy (if missing)
    if 'hist_entropy' not in roi:
        if 'data' in roi:
            roi['hist_entropy'] = fmeasure(roi['data'], 'HISE', None)
            is_modified = True

    # Area × range (if missing)
    if 'area_range' not in roi:
        if 'area' in roi and 'range_intens' in roi:
            roi['area_range'] = roi['area'] * roi['range_intens']
            is_modified = True

    # Blur index (if missing)
    if 'blur_idx' not in roi:
        if 'data' in roi:
            roi['blur_idx'] = compute_blur_index(roi['data'])
            is_modified = True

    # D90 and updated Dmax (if missing)
    if 'D90' not in roi:
        if 'hull' in roi and 'bw_mask_filled' in roi:
            # Recompute Dmax from hull
            Dmax, Dmax_theta, DmaxA, DmaxB = compute_Dmax(
                roi['hull']['xh'],
                roi['hull']['yh']
            )
            roi['Dmax'] = Dmax
            roi['Dmax_theta'] = Dmax_theta
            roi['DmaxA'] = DmaxA
            roi['DmaxB'] = DmaxB

            # Compute D90
            roi['D90'] = compute_D90(roi['bw_mask_filled'], Dmax, Dmax_theta)
            is_modified = True

    return roi
