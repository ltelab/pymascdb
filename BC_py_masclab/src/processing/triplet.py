"""
Triplet image processing module for MASC analysis.

This module handles the processing of triplet MASC snowflake images
(3 cameras capturing the same snowflake simultaneously).

The processing pipeline includes:
1. Load all 3 (or more) camera images
2. Apply masking to remove clutter
3. Create binary masks
4. Detect ROIs on each image
5. Match ROIs across cameras based on vertical position and size
6. Select best matching triplet based on focus metric
7. Extract features for each camera view
8. Save results

Translated from MASC_triplet_process.m (Christophe Praz 2015) and adapted for Python
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
from typing import Dict, Tuple, List
from datetime import datetime

from src.core.config import ProcessingConfig, LabelConfig, create_proc_params_file
from src.preprocessing.masking import masking
from src.preprocessing.roi import ROIDetector
from src.features.basic import process_basic_descriptors
from src.utils.image import rangefilt
from src.utils.identifiers import get_snowflake_id_or_fallback
from src.dataio.load_roi_data import save_roi_data
from src.dataio.upload import lookup_fallspeed


logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)


def process_triplet_image(
    flake_id: int,
    idx_pics: List[int],
    current_dir: Path,
    pic_list,
    label: LabelConfig,
    process: ProcessingConfig
) -> Tuple[Dict, int]:
    """
    Process a multi-view snowflake event (2, 3, 4+ cameras).

    Loads all configured camera images, detects ROIs per view, matches ROIs
    across the core triplet cameras, extracts descriptors, and optionally
    saves images, figures, and ROI data.

    Args:
        flake_id: Snowflake event identifier.
        idx_pics: Indices into pic_list for this flake's camera images.
        current_dir: Directory containing the image files.
        pic_list: Campaign picture list with cam, files, time_num, fallspeed.
        label: Label configuration with output paths.
        process: Processing configuration with matching and save options.

    Returns:
        Tuple of (timing, flag):
            - timing: Per-step elapsed times in seconds (empty dict on load failure).
            - flag: Quality flag (-1=load error, 0=no ROI, 1=bad match/quality, 2=good).

    Notes:
        - Matching is performed on the first 3 cameras in camera_order.
        - Additional cameras (4th, 5th, ...) are processed independently.
        - The pivot/reference camera is the middle camera in the core triplet.


    Translated from MASC_triplet_process.m (Christophe Praz 2015) and adapted for Python
    """
    timing = {}
    camera_order = getattr(process, "camera_order", [0, 1, 2])
    required_views = max(2, int(getattr(process, "triplet_required_views", 3)))
    allow_partial = bool(getattr(process, "triplet_allow_partial", True))

    # 1) Load all configured camera images for this flake ID
    t_loading = time.time()
    camera_data = {}
    try:
        for i in idx_pics:
            cam_id = pic_list.cam[i]
            if cam_id not in camera_order:
                continue
            img_path = Path(current_dir) / pic_list.files[i]
            if not img_path.exists():
                raise FileNotFoundError(f"Image file not found: {img_path}")
            img_data = cv2.imread(str(img_path), cv2.IMREAD_GRAYSCALE)
            if img_data is None:
                raise ValueError(f"Failed to load image: {img_path}")
            if img_data.dtype != np.uint8:
                img_data = img_data.astype(np.uint8)
            camera_data[cam_id] = {
                "data": img_data,
                "filename": pic_list.files[i],
                "cam": cam_id,
                "time_num": pic_list.time_num[i] if i < len(pic_list.time_num) else None,
            }
    except Exception as exc:
        logger.warning("Couldn't load flake %s. Event discarded: %s", flake_id, exc)
        return {}, -1
    timing["loading"] = time.time() - t_loading

    available_cams = [c for c in camera_order if c in camera_data]
    min_views = required_views if not allow_partial else min(required_views, len(available_cams))
    if len(available_cams) < max(2, min_views):
        return timing, -1

    # MATLAB behavior (MASC_triplet_process.m):
    # - Matching is performed only on the first 3 cameras in camera_order (left/mid/right).
    # - The pivot/reference camera is the middle one (camera_order(2) in MATLAB, i.e. camera_order[1] here).
    # - Any additional cameras (4th, 5th, ...) are processed independently (no multiview matching).
    #
    # We keep support for "n-plets" by allowing extra cameras to be present, but we match only the core triplet.
    core_cams = available_cams[:3]
    additional_cams = available_cams[3:]

    # If user requests fewer than 3 views (e.g. 2), still do reference-based matching on the core subset.
    if len(core_cams) < max(2, min_views):
        return timing, -1

    # Pivot camera (MATLAB uses "mid" camera as reference)
    ref_cam = core_cams[1] if len(core_cams) >= 3 else core_cams[0]

    # Fallspeed metadata (indexed by flake id / dataInfo.txt, not by image index)
    flake_fallspeed = lookup_fallspeed(pic_list, flake_id)

    # 2) Masking
    t_clutter = time.time()
    for cam_id in available_cams:
        camera_data[cam_id]["data"] = masking(camera_data[cam_id]["data"], cam_id, process)
    timing["clutter"] = time.time() - t_clutter

    # 3) Binary masks
    t_edging = time.time()
    for cam_id in available_cams:
        eroded = np.zeros_like(camera_data[cam_id]["data"], dtype=np.uint8)
        eroded[camera_data[cam_id]["data"] > 0] = 1
        camera_data[cam_id]["eroded"] = ndimage.binary_fill_holes(eroded).astype(np.uint8)
    timing["edging"] = time.time() - t_edging

    # 4) ROI detection (no duplicate descriptor logic)
    t_roiying = time.time()
    detector = ROIDetector(process)
    for cam_id in available_cams:
        compute_best = cam_id in additional_cams
        all_roi, idx_best, _, _, _ = detector.detect(
            camera_data[cam_id]["data"],
            camera_data[cam_id]["eroded"],
            cam_id,
            compute_best=compute_best,
        )
        camera_data[cam_id]["all_roi"] = all_roi
        camera_data[cam_id]["idx_best"] = idx_best
        if cam_id in core_cams and len(all_roi) == 0:
            timing["roiying"] = time.time() - t_roiying
            timing["matching"] = 0
            timing["feature"] = 0
            timing["plotting"] = 0
            timing["saving"] = 0
            return timing, 0
    timing["roiying"] = time.time() - t_roiying

    # 5) Multiview matching against reference camera (MATLAB-like logic)
    t_matching = time.time()
    ref_rois = camera_data[ref_cam]["all_roi"]
    best_assignment = None
    best_score = -np.inf
    for ref_idx, ref_roi in enumerate(ref_rois):
        assignment = {ref_cam: ref_idx}
        ref_h = ref_roi.bbox[2] - ref_roi.bbox[0]
        ref_y = ref_roi.bbox[0]
        # MATLAB uses: matching_tol = max(pix, matching_tol_percent * mid_height)
        # where matching_tol_percent is a FRACTION (1 = 100% of particle height).
        # The Python config stores a PERCENT (0-100, e.g. 100 for 100%), so it
        # is always divided by 100 here. Do not copy MATLAB values (fractions)
        # into the YAML without converting them.
        tol_percent = float(getattr(process, "matching_tol_percent", 0.0))
        tol_factor = (tol_percent / 100.0)
        tol = max(float(process.matching_tol_pix), float(tol_factor) * float(ref_h))
        valid = True
        for cam_id in core_cams:
            if cam_id == ref_cam:
                continue
            match_idx = None
            best_var = float("inf")
            for cand_idx, cand_roi in enumerate(camera_data[cam_id]["all_roi"]):
                cand_h = cand_roi.bbox[2] - cand_roi.bbox[0]
                cand_y = cand_roi.bbox[0]
                if abs(cand_h - ref_h) <= tol and abs(cand_y - ref_y) <= tol:
                    var = abs(cand_h - ref_h) + abs(cand_y - ref_y)
                    if var < best_var:
                        best_var = var
                        match_idx = cand_idx
            if match_idx is None:
                valid = False
                break
            assignment[cam_id] = match_idx
        if not valid:
            continue
        area_focus = []
        for cam_id, roi_idx in assignment.items():
            roi_obj = camera_data[cam_id]["all_roi"][roi_idx]
            f = _compute_roi_focus(camera_data[cam_id]["data"], roi_obj)
            area_focus.append(f * roi_obj.area)
        score = float(np.mean(area_focus))
        if score > best_score:
            best_score = score
            best_assignment = assignment
    timing["matching"] = time.time() - t_matching
    if best_assignment is None:
        timing["feature"] = 0
        timing["plotting"] = 0
        timing["saving"] = 0
        return timing, 1

    # 6) Build ROI descriptors for selected views
    t_feature = time.time()
    selected_rois = []
    for cam_id in core_cams:
        roi_obj = camera_data[cam_id]["all_roi"][best_assignment[cam_id]]
        roi_dict = _roi_obj_to_regionprops_dict(roi_obj)
        roi = process_basic_descriptors(camera_data[cam_id]["data"], roi_dict, process)
        _attach_common_metadata(
            roi=roi,
            camera_entry=camera_data[cam_id],
            flake_id=flake_id,
            fallspeed=flake_fallspeed,
            n_roi=len(camera_data[cam_id]["all_roi"]),
        )
        selected_rois.append(roi)
    for cam_id in additional_cams:
        all_roi = camera_data[cam_id]["all_roi"]
        idx_best = camera_data[cam_id]["idx_best"]
        if not all_roi or idx_best is None:
            continue
        roi_obj = all_roi[idx_best]
        roi_dict = _roi_obj_to_regionprops_dict(roi_obj)
        roi = process_basic_descriptors(camera_data[cam_id]["data"], roi_dict, process)
        _attach_common_metadata(
            roi=roi,
            camera_entry=camera_data[cam_id],
            flake_id=flake_id,
            fallspeed=flake_fallspeed,
            n_roi=len(all_roi),
        )
        selected_rois.append(roi)
    timing["feature"] = time.time() - t_feature

    # 7) Optional figures
    t_plotting = time.time()
    if process.generate_figs:
        try:
            if not process.display_figs:
                matplotlib.use("Agg")
            else:
                matplotlib.use("TkAgg")
            for roi in selected_rois:
                roi["_figure"] = _create_individual_figure(roi)
        except Exception as exc:
            logger.warning("Could not generate triplet figures for flake %s: %s", flake_id, exc)
    timing["plotting"] = time.time() - t_plotting

    # 8) Save
    t_saving = time.time()
    is_good = _is_multiview_good(selected_rois, process)
    if process.saveresults and selected_rois:
        if not label.outdir.exists():
            label.outdir.mkdir(parents=True, exist_ok=True)
            create_proc_params_file(label.outdir, label, process)
        dt = _resolve_output_datetime(selected_rois[0])
        date_folder = dt.strftime("%Y.%m.%d") if dt else "unknown_date"
        hour_folder = dt.strftime("%H") if dt else "00"
        path2save = label.outdir / date_folder / hour_folder
        quality = "GOOD" if is_good else "BAD"
        img_dir = path2save / "IMAGES" / quality
        data_dir = path2save / "DATA" / quality
        fig_dir = path2save / "FIGURES" / quality
        img_dir.mkdir(parents=True, exist_ok=True)
        data_dir.mkdir(parents=True, exist_ok=True)
        if process.save_figs:
            fig_dir.mkdir(parents=True, exist_ok=True)
        if process.save_format == "joblib":
            ext = "joblib"
        elif process.save_format == "pkl":
            ext = "pkl"
        elif process.save_format == "mat":
            ext = "mat"
        else:
            ext = "joblib"
        for roi in selected_rois:
            _save_roi_image(roi, img_dir)
            if "name" in roi:
                base_name = Path(roi["name"]).stem
                save_roi_data(roi, data_dir / f"{base_name}.{ext}", format=process.save_format, compress=3)
                if process.save_figs and "_figure" in roi:
                    roi["_figure"].savefig(fig_dir / f"{base_name}_fig.png", dpi=100, bbox_inches="tight")

    if process.generate_figs and not process.display_figs:
        plt.close("all")
    timing["saving"] = time.time() - t_saving
    # MATLAB (MASC_triplet_process.m:555): once a triplet is matched and saved,
    # flag is always 2 — even when the quality check routed it to BAD/.
    return timing, 2


def _roi_obj_to_regionprops_dict(roi_obj) -> Dict:
    """Convert an ROI dataclass instance to a regionprops-like dictionary."""
    return {
        "BoundingBox": roi_obj.bbox,
        "image": roi_obj.mask,
        "coords": roi_obj.coords,
        "MajorAxisLength": roi_obj.major_axis_length,
        "MinorAxisLength": roi_obj.minor_axis_length,
        "Orientation": roi_obj.orientation,
        "Centroid": roi_obj.centroid,
    }


def _attach_common_metadata(
    roi: Dict,
    camera_entry: Dict,
    flake_id: int,
    fallspeed: float,
    n_roi: int,
) -> None:
    """Attach filename, camera, and flake metadata to a ROI dict (in-place)."""
    roi["n_roi"] = n_roi
    roi["name"] = camera_entry["filename"]
    roi["id"] = get_snowflake_id_or_fallback(camera_entry["filename"], flake_id)
    roi["cam"] = camera_entry["cam"]
    roi["tnum"] = camera_entry["time_num"]
    roi["fallspeed"] = fallspeed
    roi["flag_roi"] = "GOOD"
    roi["flag"] = "GOOD"
    roi["status"] = "good detection."


def _is_multiview_good(rois: List[Dict], process: ProcessingConfig) -> bool:
    """Return True if mean brightness, max intensity, and size pass thresholds."""
    if not rois:
        return False
    mean_int = float(np.mean([r.get("mean_intens", 0.0) for r in rois]))
    max_int = float(np.mean([r.get("max_intens", 0.0) for r in rois]))
    max_dim = float(np.mean([max(r.get("width", 0.0), r.get("height", 0.0)) for r in rois]))
    return (
        mean_int >= process.minbright
        and max_int >= process.max_intensthresh
        and max_dim >= process.sizemin
    )


def _resolve_output_datetime(roi: Dict):
    """Parse a datetime from ROI tnum or filename for output folder layout."""
    tnum = roi.get("tnum")
    if isinstance(tnum, datetime):
        return tnum
    filename = roi.get("name", "")
    if filename:
        try:
            parts = filename.split("_")
            date_parts = parts[0].split(".")
            time_parts = parts[1].split(".")
            return datetime(
                int(date_parts[0]),
                int(date_parts[1]),
                int(date_parts[2]),
                int(time_parts[0]),
            )
        except Exception:
            return None
    return None


def _create_individual_figure(roi: Dict) -> matplotlib.figure.Figure:
    """Create a matplotlib figure with ROI image, ellipses, and contours."""
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Display ROI image in grayscale
    ax.imshow(roi['data'], cmap='gray', vmin=0, vmax=255)
    ax.set_aspect('equal')
    ax.set_facecolor('black')
    
    # Generate ellipse points (t from 0 to 2*pi)
    t = np.linspace(0, 2*np.pi, 100)
    
    # Convert orientation from degrees to radians (negative as per MATLAB)
    if 'E' in roi and 'theta' in roi['E']:
        theta = -roi['E']['theta'] * np.pi / 180.0
    else:
        theta = 0
    
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
        ax.set_title(f"Camera {roi.get('cam', '?')} - Orientation: {roi['E']['theta']:.2f}°")
    else:
        ax.set_title(f"Camera {roi.get('cam', '?')}")
    ax.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    
    return fig


def _compute_roi_focus(img_data: np.ndarray, roi) -> float:
    """Compute focus metric as mean_intensity × range_intensity² for an ROI."""
    # bbox format is (min_row, min_col, max_row, max_col)
    min_row, min_col, max_row, max_col = roi.bbox
    roi_data = img_data[min_row:max_row, min_col:max_col]
    
    # Get pixels within ROI mask
    roi_pixels = roi_data[roi.mask > 0]
    if roi_pixels.size == 0:
        return 0.0
    
    # Compute mean intensity
    mean_intens = np.mean(roi_pixels) / 255.0
    
    # Compute range intensity (local variability)
    range_array = rangefilt(roi_data)
    range_intens = np.mean(range_array[roi.mask > 0]) / 255.0
    
    # Focus metric
    focus = mean_intens * (range_intens ** 2)
    
    return focus


def _save_roi_image(roi: Dict, save_path: Path) -> None:
    """Save the ROI grayscale image to disk if data and name are present."""
    if 'data' not in roi or 'name' not in roi:
        return
    
    img_path = save_path / roi['name']
    
    # Save image using cv2
    cv2.imwrite(str(img_path), roi['data'])
