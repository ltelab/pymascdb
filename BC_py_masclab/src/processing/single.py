"""
Single image processing module for MASC analysis.

This module handles the processing of individual MASC snowflake images
through masking, ROI detection, feature extraction, and result saving.

The processing pipeline includes:
1. Image loading and validation
2. Background removal (masking)
3. ROI detection
4. Feature extraction
5. Quality assessment

Translated from MASC_picture_process.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import time
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cv2

from scipy import ndimage
from pathlib import Path
from typing import Dict, Tuple, Union
from datetime import datetime

from src.core.config import ProcessingConfig, LabelConfig, create_proc_params_file
from src.preprocessing.masking import masking
from src.preprocessing.roi import ROIDetector
from src.features.basic import process_basic_descriptors
from src.utils.identifiers import get_snowflake_id_or_fallback
from src.dataio.load_roi_data import save_roi_data

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

# Canonical quality flags used across processing/stats:
# -1 = load/processing error, 0 = no ROI or bad detection, 1 = blurry, 2 = good
FLAG_ERROR = -1
FLAG_BAD_OR_NO_ROI = 0
FLAG_BLURRY = 1
FLAG_GOOD = 2


def process_single_image(
    image_path: Union[str, Path],
    pic_info: Dict,
    label: LabelConfig,
    process: ProcessingConfig,
) -> Tuple[Dict, int, Dict]:
    """
    Process a single MASC image through the complete pipeline.

    Loads the image, applies masking, builds a binary mask, detects ROIs,
    extracts basic descriptors on the best ROI, validates quality, and
    optionally saves images, figures, and ROI data.

    Args:
        image_path: Path to the image file.
        pic_info: Metadata dict with keys filename, cam, id, fallspeed, time_num.
        label: Label configuration with output paths and campaign settings.
        process: Processing configuration object with thresholds and save options.

    Returns:
        Tuple of (roi, flag, timing):
            - roi: Dictionary with extracted features and metadata, or None on failure.
            - flag: Quality flag (-1=error, 0=no/bad ROI, 1=blurry reserved, 2=good).
            - timing: Per-step elapsed times in seconds.

    Notes:
        - Returns (None, -1, timing) if the image cannot be loaded.
        - Returns (None, 0, timing) if no ROI is detected.
        - Current single-image path returns {-1, 0, 2}; flag=1 is reserved for blurry.
    """
    timing = {}
    
    # ========================================================================
    # 1. LOAD IMAGE
    # ========================================================================
    t_loading = time.time()
    
    try:
        image_path = Path(image_path)
        
        # Check if file exists
        if not image_path.exists():
            raise FileNotFoundError(f"Image file not found: {image_path}")
        
        # Load image
        flake_data = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
        if flake_data is None:
            raise ValueError(f"Failed to load image: {image_path}")
        
        # Ensure uint8
        if flake_data.dtype != np.uint8:
            flake_data = flake_data.astype(np.uint8)
        
        # Extract metadata from pic_info
        flake_filename = pic_info.get('filename', image_path.name)
        flake_cam = pic_info.get('cam', 0)
        flake_id = pic_info.get('id', 0)
        flake_fallspeed = pic_info.get('fallspeed', np.nan)
        flake_tnum = pic_info.get('time_num', None)
        
    except Exception as e:
        logger.warning("Warning: image %s is skipped (couldn't load the image): %s", image_path.name, e)
        timing['loading'] = time.time() - t_loading
        flag = FLAG_ERROR
        return None, flag, timing
    
    timing['loading'] = time.time() - t_loading
    
    # ========================================================================
    # 2. REMOVE CLUTTER (MASKING)
    # ========================================================================
    t_clutter = time.time()
    
    data = masking(flake_data, flake_cam, process)
    
    timing['clutter'] = time.time() - t_clutter
    
    # ========================================================================
    # 3. CREATE BINARY MASK (EDGE DETECTION)
    # ========================================================================
    t_edging = time.time()
    
    # Simply create a B&W mask
    data_eroded = np.zeros_like(data, dtype=np.uint8)
    data_eroded[data > 0] = 1
    
    # Fill holes
    data_eroded = ndimage.binary_fill_holes(data_eroded).astype(np.uint8)
    
    timing['edging'] = time.time() - t_edging
    
    # ========================================================================
    # 4. DETECT ROIs
    # ========================================================================
    t_roiying = time.time()
    
    detector = ROIDetector(process)
    all_roi, idx_best, area_focus_ratio, flag_roi, status = detector.detect(data, data_eroded, flake_cam, compute_best=True)
    
    timing['roiying'] = time.time() - t_roiying
    
    # If no ROI found at all
    if len(all_roi) == 0:
        logger.warning("WARNING: no ROI candidate found in %s", flake_filename)

        timing['feature'] = 0
        timing['plotting'] = 0
        timing['saving'] = 0
        flag = FLAG_BAD_OR_NO_ROI
        return None, flag, timing 
    
    # ========================================================================
    # 5. COMPUTE FEATURES ON BEST ROI
    # ========================================================================
    t_feature = time.time()
    
    # Convert ROI object to regionprops-like dict for process_basic_descriptors
    best_roi = all_roi[idx_best]
    
    # scikit-image BoundingBox format: (min_row, min_col, max_row, max_col)
    regionprops_roi = {
        'BoundingBox': best_roi.bbox,  # (min_row, min_col, max_row, max_col)
        'image': best_roi.mask,  # Binary mask local to bbox
        'coords': best_roi.coords,  # (n_pixels, 2) array of (row, col)
        'MajorAxisLength': best_roi.major_axis_length,
        'MinorAxisLength': best_roi.minor_axis_length,
        'Orientation': best_roi.orientation,  # in radians
        'Centroid': best_roi.centroid  # (row, col)
    }
    
    # Process basic descriptors
    roi = process_basic_descriptors(data, regionprops_roi, process)
    
    # Add additional descriptors not based on the image itself
    roi['area_focus_ratio'] = area_focus_ratio
    roi['flag_roi'] = flag_roi
    roi['flag'] = flag_roi  # MATLAB compatibility
    roi['status'] = status if status else ''
    roi['n_roi'] = len(all_roi)
    roi['name'] = flake_filename
    roi['id'] = get_snowflake_id_or_fallback(flake_filename, flake_id)
    roi['cam'] = flake_cam
    roi['tnum'] = flake_tnum
    roi['fallspeed'] = flake_fallspeed
    
    timing['feature'] = time.time() - t_feature
    
    # ========================================================================
    # 6. VALIDATE FLAG
    # ========================================================================
    
    # The flag is used to count good/bad images
    if roi['flag_roi'] == 'GOOD':
        roi['status'] = 'good detection.'
        flag = FLAG_GOOD
    else:
        # Keep low-quality detections in BAD/NO-ROI bucket for now.
        # Blurry (flag=1) can be introduced later once the dedicated rule is finalized.
        flag = FLAG_BAD_OR_NO_ROI
    
    # ========================================================================
    # 7. PLOTTING (if desired)
    # ========================================================================
    t_plotting = time.time()
    
    if process.generate_figs:
        try:
            # Set matplotlib backend based on display_figs setting
            if not process.display_figs:
                matplotlib.use('Agg')  # Non-interactive backend
            else :
                matplotlib.use('TkAgg')  # Interactive backend
            
            # Create figure similar to MATLAB fig2
            fig, ax = plt.subplots(figsize=(8, 8))
            
            # Display ROI image in grayscale
            ax.imshow(roi['data'], cmap='gray', vmin=0, vmax=255)
            ax.set_aspect('equal')
            ax.set_facecolor('black')
            
            # Generate ellipse points (t from 0 to 2*pi)
            t = np.linspace(0, 2*np.pi, 100)
            
            # Convert orientation from degrees to radians (negative as per MATLAB)
            theta = -roi['E']['theta'] * np.pi / 180.0
            
            # Inner ellipse (E)
            if 'E' in roi and all(k in roi['E'] for k in ['X0', 'Y0', 'a', 'b']):
                X0, Y0, a, b = roi['E']['X0'], roi['E']['Y0'], roi['E']['a'], roi['E']['b']
                xt1 = X0 + np.cos(theta) * a * np.cos(t) - np.sin(theta) * b * np.sin(t)
                yt1 = Y0 + np.sin(theta) * a * np.cos(t) + np.cos(theta) * b * np.sin(t)
                ax.plot(xt1, yt1, 'r-', linewidth=2, label='Inner ellipse (E)')
                ax.plot(X0, Y0, 'rx', markersize=10)
            
            # Outer ellipse (E_out)
            if 'E_out' in roi and all(k in roi['E_out'] for k in ['X0', 'Y0', 'a', 'b']):
                X0, Y0, a, b = roi['E_out']['X0'], roi['E_out']['Y0'], roi['E_out']['a'], roi['E_out']['b']
                xt2 = X0 + np.cos(theta) * a * np.cos(t) - np.sin(theta) * b * np.sin(t)
                yt2 = Y0 + np.sin(theta) * a * np.cos(t) + np.cos(theta) * b * np.sin(t)
                ax.plot(xt2, yt2, 'c-', linewidth=2, label='Outer ellipse (E_out)')
                ax.plot(X0, Y0, 'co', markersize=8)
            
            # Inscribed ellipse (E_in)
            if 'E_in' in roi and all(k in roi['E_in'] for k in ['X0', 'Y0', 'a', 'b']):
                X0, Y0, a, b = roi['E_in']['X0'], roi['E_in']['Y0'], roi['E_in']['a'], roi['E_in']['b']
                xt3 = X0 + np.cos(theta) * a * np.cos(t) - np.sin(theta) * b * np.sin(t)
                yt3 = Y0 + np.sin(theta) * a * np.cos(t) + np.cos(theta) * b * np.sin(t)
                ax.plot(xt3, yt3, 'g-', linewidth=2, label='Inscribed ellipse (E_in)')
                ax.plot(X0, Y0, 'gv', markersize=8)
            
            # Circumscribed circle (C_out)
            if 'C_out' in roi and all(k in roi['C_out'] for k in ['X0', 'Y0', 'r']):
                X0, Y0, r = roi['C_out']['X0'], roi['C_out']['Y0'], roi['C_out']['r']
                xt4 = X0 + r * np.cos(t)
                yt4 = Y0 + r * np.sin(t)
                ax.plot(xt4, yt4, 'b-', linewidth=2, label='Circumscribed circle (C_out)')
            
            # Perimeter
            if 'x_perim' in roi and 'y_perim' in roi:
                ax.plot(roi['x_perim'], roi['y_perim'], 'y-', linewidth=1, label='Perimeter')
            
            # Convex hull
            if 'hull' in roi and 'xh' in roi['hull'] and 'yh' in roi['hull']:
                ax.plot(roi['hull']['xh'], roi['hull']['yh'], 'c--', linewidth=1, label='Convex hull')
            
            # Bounding rectangle
            if 'Rect' in roi and 'rectx' in roi['Rect'] and 'recty' in roi['Rect']:
                ax.plot(roi['Rect']['rectx'], roi['Rect']['recty'], 'w-', linewidth=1, label='Bounding rect')
            
            # Labels and title
            ax.set_xlabel('x axis [pixels]')
            ax.set_ylabel('y axis [pixels]')
            if 'E' in roi and 'theta' in roi['E']:
                ax.set_title(f"Orientation: {roi['E']['theta']:.2f}°")
            ax.legend(loc='upper right', fontsize=8)
            
            # Invert y-axis to match MATLAB's image coordinate system
            
            # Store figure in roi dict for later saving
            roi['_figure'] = fig
            
            if process.display_figs:
                plt.draw()
                plt.pause(0.001)  # Allow GUI to update
            
        except Exception as e:
            logger.warning("Warning: Could not generate figure: %s", e)
            
    
    timing['plotting'] = time.time() - t_plotting
    
    # ========================================================================
    # 8. SAVE RESULTS
    # ========================================================================
    t_saving = time.time()
    
    if process.saveresults:
        
        # Create metadata file if first time
        if not label.outdir.exists():
            label.outdir.mkdir(parents=True, exist_ok=True)
            create_proc_params_file(label.outdir, label, process)
            logger.info('Creation of metadata file proc_params.txt...')
        
        # Create date/hour folders from tnum
        tnum = roi.get('tnum')
        if tnum is not None and not np.isnan(tnum):
            try:
                # Check if tnum is already a datetime object
                if isinstance(tnum, datetime):
                    dt = tnum
                elif isinstance(tnum, (int, float)):
                    # Extract date from filename (format: YYYY.MM.DD_HH.MM.SS)
                    filename = roi.get('name', '')
                    if filename:
                        # Parse date from filename
                        parts = filename.split('_')
                        if len(parts) >= 2:
                            date_str = parts[0]  # YYYY.MM.DD
                            time_str = parts[1]  # HH.MM.SS
                            
                            date_parts = date_str.split('.')
                            time_parts = time_str.split('.')
                            
                            if len(date_parts) >= 3 and len(time_parts) >= 3:
                                year = int(date_parts[0])
                                month = int(date_parts[1])
                                day = int(date_parts[2])
                                hour = int(time_parts[0])
                                dt = datetime(year, month, day, hour)
                            else:
                                raise ValueError("Cannot parse date from filename")
                        else:
                            raise ValueError("Cannot parse date from filename")
                    else:
                        raise ValueError("No filename available")
                else:
                    raise ValueError(f"Unexpected tnum type: {type(tnum)}")
                
                date_folder = dt.strftime('%Y.%m.%d')
                hour_folder = dt.strftime('%H')
                
            except Exception as e:
                logger.warning("Could not parse date from tnum or filename: %s", e)
                date_folder = 'unknown_date'
                hour_folder = '00'
        else:
            # Fallback if no time info
            date_folder = 'unknown_date'
            hour_folder = '00'
        
        path2save = label.outdir / date_folder / hour_folder
        
        # Determine if GOOD or BAD based on flag from roi_detection
        # The flag was already validated in roi_detection with the same thresholds
        if flag == FLAG_GOOD:
            save_path = path2save / 'IMAGES' / 'GOOD'
        else:  # BAD/NO-ROI or BLURRY
            save_path = path2save / 'IMAGES' / 'BAD'
        
        save_path.mkdir(parents=True, exist_ok=True)
        
        # Save ROI image
        if 'data' in roi and 'name' in roi:
            img_path = save_path / roi['name']
            cv2.imwrite(str(img_path), roi['data'])
            logger.info("Saved ROI to: %s", img_path)
        
        # Save ROI data
        # Determine save path based on flag from roi_detection
        if flag == FLAG_GOOD:
            data_path = path2save / 'DATA' / 'GOOD'
        else:  # BAD/NO-ROI or BLURRY
            data_path = path2save / 'DATA' / 'BAD'
        
        data_path.mkdir(parents=True, exist_ok=True)
        
        if 'name' in roi:
            # Determine file extension based on save format
            if process.save_format == "joblib":
                ext = "joblib"
            elif process.save_format == "pkl":
                ext = "pkl"
            elif process.save_format == "mat":
                ext = "mat"
            else:
                ext = "joblib"  # Default fallback
            
            # Save ROI data using save_roi_data utility
            base_name = Path(roi['name']).stem
            data_file_path = data_path / f'{base_name}.{ext}'
            
            save_roi_data(roi, data_file_path, format=process.save_format, compress=3)
            
            logger.info("Saved ROI data to: %s", data_file_path)
        
        # Save figure if requested
        if process.save_figs and '_figure' in roi:
            # Determine figure save path based on flag from roi_detection
            if flag == FLAG_GOOD:
                fig_dir = path2save / 'FIGURES' / 'GOOD'
            else:  # BAD/NO-ROI or BLURRY
                fig_dir = path2save / 'FIGURES' / 'BAD'
            
            fig_dir.mkdir(parents=True, exist_ok=True)
            
            # Extract base filename without extension
            base_name = Path(roi['name']).stem
            fig_path = fig_dir / f'{base_name}_fig.png'
            
            roi['_figure'].savefig(fig_path, dpi=100, bbox_inches='tight')
            
            logger.info("Saved figure to: %s", fig_path)
    
    # Close figures if not displaying
    if process.generate_figs and not process.display_figs:
        plt.close('all')
    
    timing['saving'] = time.time() - t_saving
    
    return roi, flag, timing
