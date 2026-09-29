"""
Basic feature extraction module for MASC analysis.

This module provides functions for extracting basic descriptors from
snowflake ROIs, including intensity, geometry, topology, and texture features.

Functions for:
- Intensity, geometric, shape, topological, and texture descriptor computation
- Optional brightening and recomputation of textural features

Translated from process_basic_descriptors.m and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import cv2
import logging 

from typing import Dict
from scipy import ndimage
from skimage.measure import regionprops, label

from src.features.fmeasure import fmeasure
from src.preprocessing.brightening import brightening
from src.preprocessing.roi import ROI
from src.core.config import ProcessingConfig
from src.utils.image import rangefilt, stdfilt
from src.utils.descriptors_heplers import (
    compute_convex_hull,
    compute_Dmax,
    compute_D90,
    fit_circle_around,
    fit_ellipse_inside,
    fit_ellipse_around,
    compute_rectangularity,
    skeleton_props,
    fractal_dim,
    haralick_props,
    compute_symmetry_features,
    compute_blur_index,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)


def _shift_coords_to_matlab_1based(roi: Dict) -> None:
    """
    Convert stored coordinates from 0-based (Python) to 1-based (MATLAB).

    All geometric descriptors are computed in 0-based crop/image indices, then
    shifted once here so exported fields match MATLAB regionprops / find /
    bwboundaries conventions. Centroid fields and E.X0/Y0 are already 1-based
    and must not be shifted again.
    """
    roi['x'] = np.asarray(roi['x']) + 1
    roi['y'] = np.asarray(roi['y']) + 1
    if len(roi['x_perim']) > 0:
        roi['x_perim'] = np.asarray(roi['x_perim']) + 1
        roi['y_perim'] = np.asarray(roi['y_perim']) + 1
    roi['x_loc'] = int(roi['x_loc']) + 1
    roi['y_loc'] = int(roi['y_loc']) + 1

    hull = roi.get('hull')
    if hull is not None:
        if len(hull.get('xh', [])) > 0:
            hull['xh'] = np.asarray(hull['xh']) + 1
            hull['yh'] = np.asarray(hull['yh']) + 1

    for key in ('DmaxA', 'DmaxB'):
        pt = roi.get(key)
        if pt is not None and len(pt) == 2:
            roi[key] = (float(pt[0]) + 1.0, float(pt[1]) + 1.0)

    for key in ('C_out', 'E_in', 'E_out'):
        fit = roi.get(key)
        if fit is not None and 'X0' in fit and 'Y0' in fit:
            fit['X0'] = float(fit['X0']) + 1.0
            fit['Y0'] = float(fit['Y0']) + 1.0

    rect = roi.get('Rect')
    if rect is not None and len(rect.get('rectx', [])) > 0:
        rect['rectx'] = np.asarray(rect['rectx']) + 1
        rect['recty'] = np.asarray(rect['recty']) + 1


def process_basic_descriptors(
    data_in: np.ndarray,
    regionprops_roi: Dict,
    process: ProcessingConfig
) -> Dict:
    """
    Process basic descriptors for a detected ROI.

    Builds the ROI dictionary with intensity, geometry, topology, texture,
    and symmetry features. Optionally brightens the cropped image and
    recomputes selected textural descriptors.

    Args:
        data_in: Input grayscale image (uint8).
        regionprops_roi: Regionprops-like dict with BoundingBox, image, coords,
            Centroid, MajorAxisLength, MinorAxisLength, and Orientation.
        process: Processing configuration (min_hole_area, flakebrighten, etc.).

    Returns:
        ROI dictionary with all computed basic descriptors.

    Notes:
        - Masks pixels outside the detection before cropping.
        - Computes hull, Dmax, D90, skeleton, fractal, Haralick, and fmeasure features.
        - When process.flakebrighten is True, stores brightened variants under roi['new'].
    """
    roi = {}

    bbox = regionprops_roi['BoundingBox']
    min_row, min_col, max_row, max_col = bbox

    # Zero pixels outside detection (masclab: mask(PixelIdxList), data_in(~mask)=0)
    data_masked = np.zeros_like(data_in)
    coords = regionprops_roi.get('coords')
    if coords is not None and len(coords) > 0:
        cr, cc = coords[:, 0], coords[:, 1]
        data_masked[cr, cc] = data_in[cr, cc]
    elif regionprops_roi.get('image') is not None:
        crop = np.asarray(data_in[min_row:max_row, min_col:max_col], dtype=np.uint8).copy()
        crop[~regionprops_roi['image'].astype(bool)] = 0
        data_masked[min_row:max_row, min_col:max_col] = crop
    else:
        data_masked = np.asarray(data_in, copy=True)

    roi_data = np.ascontiguousarray(data_masked[min_row:max_row, min_col:max_col], dtype=np.uint8)
    roi['data'] = roi_data.copy()
    local_mask = roi_data > 0
    
    # Generate masks of the snowflake + holes
    roi['bw_mask'] = local_mask
    roi['bw_mask_filled'] = ndimage.binary_fill_holes(roi['bw_mask'])
    
    # Detect holes
    bw_holes_mask = roi['bw_mask_filled'] & ~roi['bw_mask']

    if np.any(bw_holes_mask):
        # Use regionprops to detect holes
        holes_labeled = label(bw_holes_mask)
        holes_props = regionprops(holes_labeled)
        
        # Filter holes by minimum area

        valid_holes = [h for h in holes_props if h.area > process.min_hole_area]

        roi['nb_holes'] = len(valid_holes)
        roi['holes_mask'] = np.zeros_like(roi_data, dtype=bool)
            
        if roi['nb_holes'] > 0:
            all_coords = np.vstack([hole.coords for hole in valid_holes])
            roi['holes_mask'][tuple(all_coords.T)] = True 

    else:
        roi['nb_holes'] = 0
        roi['holes_mask'] = np.zeros_like(roi_data, dtype=bool)
    
    # Final mask: filled - holes
    roi['bw_mask'] = roi['bw_mask_filled'] & ~roi['holes_mask']

    # Get coordinates of filled mask
    y_idx, x_idx = np.where(roi['bw_mask_filled']>0)
    roi['y'] = y_idx.reshape(-1, 1)
    roi['x'] = x_idx.reshape(-1, 1)
    roi['area'] = float(np.sum(roi['bw_mask_filled']))
    roi['area2'] = roi['area']  # MATLAB compatibility
    roi['area_porous'] = float(np.sum(roi['bw_mask']))
    
    # Intensity (mean/max normalized; contrast later — masclab order)
    pixels_filled = roi['data'][roi['bw_mask_filled']]
    if pixels_filled.size > 0:
        roi['mean_intens'] = float(np.mean(pixels_filled, dtype=np.float64)) / 255.0
        roi['max_intens'] = float(np.max(pixels_filled)) / 255.0
        range_array = rangefilt(roi['data'], size=3)
        roi['range_intens'] = float(np.mean(range_array[roi['bw_mask_filled']], dtype=np.float64)) / 255.0
        roi['focus'] = roi['mean_intens'] * roi['range_intens']
        roi['area_focus'] = roi['area'] * roi['focus']
        roi['area_range'] = roi['area'] * roi['range_intens']
    else:
        roi['mean_intens'] = roi['max_intens'] = 0.0
        roi['range_intens'] = roi['focus'] = 0.0
        roi['area_focus'] = roi['area_range'] = 0.0
    
    # Location in original image (0-based here; +1 applied at the end like MATLAB
    # x_loc = ceil(BoundingBox(2)), y_loc = ceil(BoundingBox(1)))
    roi['x_loc'] = min_row
    roi['y_loc'] = min_col
    
    # Fitted ellipse parameters (from regionprops, using scikit-image convention)
    # scikit-image orientation: angle from 0th axis (rows/vertical), in radians, range [-pi/2, pi/2]
    # MATLAB regionprops Orientation: angle from x-axis (horizontal), in degrees, range [-90, 90]
    # Image coords (y down): matlab_deg = degrees(skimage_rad) - 90, which lands
    # in [-180, 0]; wrap by +180 to fold back into MATLAB's (-90, 90] range
    # (an ellipse orientation is defined modulo 180 degrees).
    maj_axis = regionprops_roi.get('MajorAxisLength', 0.0)
    min_axis = regionprops_roi.get('MinorAxisLength', 0.0)
    orientation_rad = regionprops_roi.get('Orientation', 0.0)
    orientation_deg = np.degrees(orientation_rad) - 90.0
    if orientation_deg < -90.0:
        orientation_deg += 180.0
    
    roi['E'] = {
        'a': maj_axis / 2.0,
        'b': min_axis / 2.0,
        'theta': orientation_deg
    }
    roi['orientation'] = orientation_deg
    roi['major_axis_length'] = maj_axis
    roi['minor_axis_length'] = min_axis
    
    # Centroid (global and local versions)
    # skimage Centroid: (row, col) 0-based in full image; MATLAB regionprops: [x, y] 1-based
    centroid_full = regionprops_roi.get('Centroid', (min_row, min_col))
    row_c = float(centroid_full[0])
    col_c = float(centroid_full[1])
    centroid_global = (col_c + 1.0, row_c + 1.0)  # MATLAB [x=col, y=row]
    centroid_local_init = (col_c - min_col + 1.0, row_c - min_row + 1.0)
    roi['centroid_global'] = centroid_global
    roi['centroid'] = centroid_global
    roi['centroid_local_init'] = centroid_local_init
    
    # Perimeter using boundary
    contours, _ = cv2.findContours(roi['bw_mask_filled'].astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    if len(contours) > 0: 
        if len(contours) > 1:
            logger.warning('Warning: more than 1 particle detected on bw_mask_filled !!!')

        contour = contours[0].squeeze()
        if contour.ndim == 1:
            contour = contour.reshape(-1, 2)

        # Drop a possible duplicate closing point before re-ordering
        if len(contour) > 1 and np.array_equal(contour[0], contour[-1]):
            contour = contour[:-1]

        # MATLAB bwboundaries starts at the leftmost, topmost boundary pixel
        # (min column / x, then min row / y), then walks the contour.
        if len(contour) > 0:
            start = int(np.lexsort((contour[:, 1], contour[:, 0]))[0])
            contour = np.vstack([contour[start:], contour[:start]])

        # MATLAB bwboundaries(..., 'noholes') returns a CLOSED contour
        # (first point repeated at the end) -> perim = n_boundary_pixels + 1
        if len(contour) > 0:
            contour = np.vstack([contour, contour[0]])

        roi['x_perim'] = contour[:, 0].reshape(-1, 1)  # Field 24
        roi['y_perim'] = contour[:, 1].reshape(-1, 1)  # Field 23
        roi['perim'] = float(len(roi['x_perim']))

        # Build a perimeter mask (bw_perim) similar to MATLAB output
        bw_perim = np.zeros_like(roi['bw_mask_filled'], dtype=bool)
        y_idx = np.clip(roi['y_perim'].astype(int).ravel(), 0, bw_perim.shape[0]-1)
        x_idx = np.clip(roi['x_perim'].astype(int).ravel(), 0, bw_perim.shape[1]-1)
        bw_perim[y_idx, x_idx] = True
        roi['bw_perim'] = bw_perim  # Field 20

    else:
        roi['x_perim'] = np.empty((0, 1))
        roi['y_perim'] = np.empty((0, 1))
        roi['perim'] = 0.0
        roi['bw_perim'] = np.zeros_like(roi['bw_mask_filled'], dtype=bool)
    
    # Convex hull — pass perim so convexity = hull.perim / roi.perim (like MATLAB)
    roi['hull'] = compute_convex_hull(roi['x'], roi['y'], roi['perim'])
    
    # Dimensions
    roi['width'] = float(roi['bw_mask'].shape[1])
    roi['height'] = float(roi['bw_mask'].shape[0])
    roi['Dmean'] = 0.5 * (roi['width'] + roi['height'])
    
    # Dmax and related
    roi['Dmax'], roi['Dmax_theta'], roi['DmaxA'], roi['DmaxB'] = compute_Dmax(roi['hull']['xh'], roi['hull']['yh'])
    
    roi['eq_radius'] = np.sqrt(roi['area'] / np.pi)
    roi['D90'] = compute_D90(roi['bw_mask_filled'], roi['Dmax'], roi['Dmax_theta'])
    
    # Complexity
    roi['complex'] = len(roi['y_perim'])/ (2.0 * np.pi * roi['eq_radius'])

    # Circumscribed circle
    roi['C_out'] = fit_circle_around(roi['x_perim'], roi['y_perim'], roi['hull']['xh'], roi['hull']['yh'])
    
    # Local centroid (in cropped coordinates) — equivalent to MATLAB:
    #   tmp = regionprops(roi.bw_mask_filled, 'Centroid'); roi.centroid_local = tmp.Centroid;
    #   roi.E.X0 = roi.centroid_local(1); roi.E.Y0 = roi.centroid_local(2);
    # skimage centroid is (row, col); MATLAB Centroid is [x=col, y=row]
    local_props = regionprops(roi['bw_mask_filled'].astype(np.uint8))
    if local_props:
        _lc = local_props[0].centroid  # (row, col) 0-indexed
        roi['centroid_local'] = (float(_lc[1]) + 1.0, float(_lc[0]) + 1.0)  # MATLAB [x, y]
    else:
        roi['centroid_local'] = (0.0, 0.0)

    roi['E']['X0'] = roi['centroid_local'][0]  # x = col
    roi['E']['Y0'] = roi['centroid_local'][1]  # y = row
    
    # Inscribed and circumscribed ellipses
    theta_rad = roi['E']['theta'] * np.pi / 180.0
    roi['E_in'] = fit_ellipse_inside(roi['x'], roi['y'], roi['x_perim'], roi['y_perim'], theta_rad)
    roi['E_out'] = fit_ellipse_around(roi['x'], roi['y'], roi['hull']['xh'], roi['hull']['yh'], theta_rad)
    
    # Rectangularity
    roi['Rect'] = compute_rectangularity(roi['x'], roi['y'], roi['perim'])
    
    # Shape descriptors
    if roi['E_out']['a'] > 0 and roi['E_out']['b'] > 0:
        roi['compactness'] = roi['area'] / (np.pi * roi['E_out']['a'] * roi['E_out']['b'])
    else:
        roi['compactness'] = 0.0
    
    if roi['C_out']['A'] > 0:
        roi['roundness'] = roi['area'] / roi['C_out']['A']
    else:
        roi['roundness'] = 0.0
    
    # Skeleton properties
    roi['skel'] = skeleton_props(roi['data'])
    
    # Fractal dimension
    roi['F'] = fractal_dim(roi['data'])
    roi['F_jac'] = 2.0 * np.log(roi['perim'] / 4.0) / np.log(roi['area']) 
    
    # Symmetry features
    roi['Sym'] = compute_symmetry_features(roi['bw_mask_filled'], roi['Dmax'], roi['eq_radius'])
    
    # Haralick texture features
    roi['H'] = haralick_props(roi['data'])

    # Blur index (imgaussfilt sigma=2/4 + std2) — process_new_descriptors.m
    roi['blur_idx'] = compute_blur_index(roi['data'])

    # More textural descriptors using fmeasure
    roi['lap'] = fmeasure(roi['data'], 'LAPM', None)
    roi['area_lap'] = roi['lap'] * roi['area']
    roi['hist_entropy'] = fmeasure(roi['data'], 'HISE', None)
    roi['wavs'] = fmeasure(roi['data'], 'WAVS', None)
    
    filled = roi['bw_mask_filled']
    pixels_mask = roi['data'][filled]
    if pixels_mask.size > 1:
        # MATLAB std2: sample std (N-1 normalization)
        roi['std'] = float(np.std(pixels_mask, dtype=np.float64, ddof=1))
    else:
        roi['std'] = 0.0

    local_std = stdfilt(roi['data'], size=3)
    roi['local_std'] = float(np.mean(local_std[filled], dtype=np.float64)) if np.any(filled) else 0.0

    local_std5 = stdfilt(roi['data'], size=5)
    roi['local_std5'] = float(np.mean(local_std5[filled], dtype=np.float64)) if np.any(filled) else 0.0

    local_std7 = stdfilt(roi['data'], size=7)
    roi['local_std7'] = float(np.mean(local_std7[filled], dtype=np.float64)) if np.any(filled) else 0.0

    roi['range_complex'] = roi['complex'] * roi['range_intens']
    # masclab: roi.contrast = double(roi.max_intens - roi.min_intens)/roi.mean_intens
    # where max_intens is a double in [0,1] but min_intens is a raw uint8.
    # In MATLAB, double-uint8 arithmetic rounds the double and SATURATES at 0,
    # so the numerator is round(max_intens) - min_intens clipped to >= 0:
    # almost always 0 (min>=1), or 1 when the mask contains a 0 pixel (holes)
    # and max_intens >= 0.5. We reproduce that quirk to match MATLAB output.
    if pixels_mask.size > 0 and roi['mean_intens'] > 0:
        roi['min_intens'] = float(np.min(pixels_mask))
        max_rounded = 1.0 if roi['max_intens'] >= 0.5 else 0.0
        roi['contrast'] = max(0.0, max_rounded - roi['min_intens']) / roi['mean_intens']
    else:
        roi['min_intens'] = 0.0
        roi['contrast'] = 0.0
    
    # Brightness adjustment if requested
    roi['new'] = {}
    if process.flakebrighten:
        roi['new']['data'] = brightening(roi['data'])
        
        # Recompute textural descriptors on brightened image
        roi['new']['lap'] = fmeasure(roi['new']['data'], 'LAPM', None)
        roi['new']['area_lap'] = roi['new']['lap'] * roi['area']
        
        range_array_new = rangefilt(roi['new']['data'], size=3)
        roi['new']['range_intens'] = float(np.mean(range_array_new[roi['bw_mask_filled']])) / 255.0 if np.any(roi['bw_mask_filled']) else 0.0
        roi['new']['range_complex'] = roi['complex'] * roi['new']['range_intens']
        
        pixels_new = roi['new']['data'][roi['bw_mask_filled']]
        if len(pixels_new) > 1:
            # MATLAB std2 on brightened crop
            roi['new']['std'] = float(np.std(pixels_new, ddof=1))
            
            local_std_new = stdfilt(roi['new']['data'], size=3)
            roi['new']['local_std'] = float(np.mean(local_std_new[roi['bw_mask_filled']])) if np.any(roi['bw_mask_filled']) else 0.0
            
            local_std5_new = stdfilt(roi['new']['data'], size=5)
            roi['new']['local_std5'] = float(np.mean(local_std5_new[roi['bw_mask_filled']])) if np.any(roi['bw_mask_filled']) else 0.0
            
            local_std7_new = stdfilt(roi['new']['data'], size=7)
            roi['new']['local_std7'] = float(np.mean(local_std7_new[roi['bw_mask_filled']])) if np.any(roi['bw_mask_filled']) else 0.0
            
            mean_new = float(np.mean(pixels_new, dtype=np.float64))
            roi['new']['contrast'] = (
                (float(np.max(pixels_new)) - float(np.min(pixels_new))) / mean_new
                if mean_new > 0 else 0.0
            )
        
        else:
            roi['new']['std'] = 0.0
            roi['new']['local_std'] = 0.0
            roi['new']['local_std5'] = 0.0
            roi['new']['local_std7'] = 0.0
            roi['new']['contrast'] = 0.0
    
    else: ## ATTENTION Pas forcément nécessaire, à vérifier le comportement
        # Brightening disabled: copy original values
        roi['new']['data'] = roi['data'].copy()
        roi['new']['lap'] = roi['lap']
        roi['new']['area_lap'] = roi['area_lap']
        roi['new']['range_intens'] = roi['range_intens']
        roi['new']['range_complex'] = roi['range_complex']
        roi['new']['std'] = roi['std']
        roi['new']['local_std'] = roi['local_std']
        roi['new']['local_std5'] = roi['local_std5']
        roi['new']['local_std7'] = roi['local_std7']
        roi['new']['contrast'] = roi['contrast']
    
    # Compute magic quality parameter (xhi)
    avg_lap = (roi['lap'] + roi['new']['lap']) / 2.0
    avg_local_std = (roi['local_std'] + roi['new']['local_std']) / 2.0
    
    # Ensure all components are positive before taking log (matches MATLAB behavior)
    if avg_lap > 0 and roi['complex'] > 0 and avg_local_std > 0 and roi['Dmean'] > 0:
        roi['xhi'] = np.log(avg_lap * roi['complex'] * avg_local_std * roi['Dmean'])
    else:
        roi['xhi'] = 0.0  # Default value when log cannot be computed

    # MATLAB uses 1-based coordinates for all exported geometry fields.
    # Apply the shift once at the end (bw_perim was built with 0-based indices).
    _shift_coords_to_matlab_1based(roi)

    return roi
