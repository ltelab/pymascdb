"""
Region-of-interest (ROI) detection for MASC snowflake images.

This module handles the detection and validation of regions of interest
in MASC snowflake images.
Functions for:
- Connected-component labelling and ROI quality metrics
- Border filtering, best-ROI selection, and validation

Translated from roi_detection.m (Christophe Praz 2016) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import cv2

from dataclasses import dataclass
from typing import List, Tuple, Optional
from skimage.measure import label, regionprops

from src.core.config import ProcessingConfig
from src.utils.image import rangefilt


@dataclass
class ROI:
    """Detected region of interest with geometry and quality metrics."""
    bbox: Tuple[int, int, int, int]  # min_row, min_col, max_row, max_col
    mask: np.ndarray  # Binary mask local to bbox
    area: float
    perimeter: float
    centroid: Tuple[float, float]  # (row, col)
    major_axis_length: float
    minor_axis_length: float
    orientation: float  # in radians
    mean_intensity: float = 0.0
    max_intensity: float = 0.0
    range_intensity: float = 0.0
    focus: float = 0.0
    area_focus: float = 0.0
    coords: Optional[np.ndarray] = None  # (n_pixels, 2) array of (row, col)


class ROIDetector:
    """Detect, filter, and validate ROIs in MASC binary edge masks."""

    def __init__(self, process_config: ProcessingConfig):
        """
        Initialize the ROI detector with processing thresholds.

        Args:
            process_config: Processing configuration with ROI thresholds.

        Returns:
            None
        """
        self.config = process_config

    def detect(
        self,
        image: np.ndarray,
        edge_mask: np.ndarray,
        cam_id: int,
        compute_best: bool = True
    ) -> Tuple[List[ROI], Optional[int], Optional[float], Optional[str], Optional[str]]:
        """
        Detect all ROIs in an image using a binary edge mask.

        Labels connected components, filters by area and border proximity,
        computes quality metrics, selects the best ROI (area × focus), and
        validates it against configured thresholds.

        Args:
            image: Original grayscale image (uint8).
            edge_mask: Binary mask from edge detection (uint8, 0 or 1).
            cam_id: Camera ID, matched against process config camera_order
                (camera_order[0]=left LED, [1]=center, [2]=right LED).
            compute_best: Whether to select and validate the best ROI.

        Returns:
            Tuple of (all_roi, idx_best, area_focus_ratio, flag, status).

        Notes:
            - Border filtering depends on the camera position resolved
              through camera_order (left LED / center / right LED).
            - Best ROI is selected by maximum area × focus.

        Translated from roi_detection.m (Christophe Praz 2016) and adapted for Python
        """
        # Label connected components using scikit-image
        labeled_img = label(edge_mask, connectivity=2)

        all_roi = []

        # Extract properties for each labeled region (skip label 0 which is background)
        props_list = regionprops(labeled_img)

        for prop in props_list:
            # Filter by minimum area early
            if prop.area <= self.config.min_area:
                continue

            # Create ROI object with scikit-image properties
            # bbox format: (min_row, min_col, max_row, max_col)
            roi = ROI(
                bbox=prop.bbox,
                mask=prop.image,  # Binary mask local to bbox
                area=float(prop.area),
                perimeter=float(prop.perimeter),
                centroid=prop.centroid,  # (row, col) in full image
                major_axis_length=float(prop.major_axis_length),
                minor_axis_length=float(prop.minor_axis_length),
                orientation=float(prop.orientation),  # in radians
                coords=prop.coords  # (n_pixels, 2) array of (row, col)
            )

            all_roi.append(roi)

        # Filter ROIs touching the discarded borders
        all_roi = self._filter_border_rois(all_roi, edge_mask.shape, cam_id)

        # If not computing best or no ROIs, return early
        if not compute_best or len(all_roi) == 0:
            return all_roi, None, None, None, None

        # Compute quality metrics for all ROIs
        for roi in all_roi:
            self._compute_roi_metrics(roi, image)

        # Select best ROI based on area × focus
        all_area_focus = [roi.area_focus for roi in all_roi]
        idx_best = int(np.argmax(all_area_focus))

        # Compute area_focus ratio (best / second best)
        if len(all_roi) > 1:
            sorted_area_focus = sorted(all_area_focus, reverse=True)
            if sorted_area_focus[1] > 0:
                area_focus_ratio = sorted_area_focus[0] / sorted_area_focus[1]
            else:
                area_focus_ratio = float('inf')  # MATLAB x/0 -> Inf
        else:
            area_focus_ratio = float('inf')

        # Validate best ROI
        flag_roi, status = self._validate_roi(all_roi[idx_best])

        return all_roi, idx_best, area_focus_ratio, flag_roi, status

    def _filter_border_rois(
        self,
        rois: List[ROI],
        image_shape: Tuple[int, int],
        cam_id: int
    ) -> List[ROI]:
        """Drop ROIs touching discarded image borders (camera-dependent)."""
        if len(rois) == 0:
            return rois

        height, width = image_shape
        filtered_rois = []
        t, b, l, r = self.config.discardmat
        camera_order = self.config.camera_order

        for roi in rois:
            # skimage bbox: (min_row, min_col, max_row, max_col) with max exclusive.
            # MATLAB regionprops BoundingBox is [x, y, w, h] with half-pixel
            # corners, so distances are offset by ±0.5 relative to integer
            # skimage edges (roi_detection.m).
            min_row, min_col, max_row, max_col = roi.bbox
            dist_from_top = min_row + 0.5
            dist_from_bot = height - max_row - 0.5
            dist_from_left = min_col + 0.5
            dist_from_right = width - max_col - 0.5

            # camera_order[0]=left LED, [1]=center, [2]=right LED
            if cam_id == camera_order[0]:  # Left LED camera
                valid = (dist_from_top > t + 1 and
                         dist_from_bot > b + 1 and
                         dist_from_left > l + 1 and
                         dist_from_right > 1)
            elif cam_id == camera_order[2]:  # Right LED camera
                valid = (dist_from_top > t + 1 and
                         dist_from_bot > b + 1 and
                         dist_from_right > r + 1 and
                         dist_from_left > 1)
            elif cam_id == camera_order[1]:  # Center camera
                valid = (dist_from_top > t + 1 and
                         dist_from_bot > b + 1 and
                         dist_from_left > 1 and
                         dist_from_right > 1)
            else:
                # Additional cameras (e.g. 4th/5th view): check all borders
                valid = (dist_from_top > t + 1 and
                         dist_from_bot > b + 1 and
                         dist_from_left > 1 and
                         dist_from_right > 1)

            if valid:
                filtered_rois.append(roi)

        return filtered_rois

    def _compute_roi_metrics(self, roi: ROI, image: np.ndarray) -> None:
        """Compute intensity and focus metrics for an ROI (in-place)."""
        # bbox is (min_row, min_col, max_row, max_col)
        min_row, min_col, max_row, max_col = roi.bbox

        # Crop image around ROI bounding box
        cropped_image = image[min_row:max_row, min_col:max_col]

        # Get pixels within ROI mask (roi.mask is already local to the bbox)
        roi_pixels = cropped_image[roi.mask > 0]

        if len(roi_pixels) == 0:
            return

        # Maximal brightness [0, 1]
        roi.max_intensity = float(np.max(roi_pixels)) / 255.0

        # Average brightness [0, 1]
        roi.mean_intensity = float(np.mean(roi_pixels)) / 255.0

        # Local variability using rangefilt (3x3 box)
        range_array = rangefilt(cropped_image, size=3)
        roi.range_intensity = float(np.mean(range_array[roi.mask > 0])) / 255.0

        # Focus parameter for ROI selection: mean_intensity × range_intensity²
        # (roi_detection.m:86 — note the ^2, unlike the stored descriptor
        # roi.focus in process_basic_descriptors which is mean × range)
        roi.focus = roi.mean_intensity * roi.range_intensity ** 2

        # Area × focus (used for ROI selection)
        roi.area_focus = roi.focus * roi.area

    def _validate_roi(self, roi: ROI) -> Tuple[str, str]:
        """Validate an ROI against size and brightness thresholds."""
        flag = 'GOOD'
        status = ''

        # Check maximum size of ROI bounding box
        # bbox is (min_row, min_col, max_row, max_col)
        min_row, min_col, max_row, max_col = roi.bbox
        height = max_row - min_row
        width = max_col - min_col
        max_size = max(width, height)

        if max_size < self.config.sizemin:
            flag = 'BAD'
            status += ' too small.'

        # Check mean brightness
        if roi.mean_intensity < self.config.minbright:
            flag = 'BAD'
            status += ' too dark (mean).'

        # Check maximum brightness
        if roi.max_intensity < self.config.max_intensthresh:
            flag = 'BAD'
            status += ' too dark (max).'

        return flag, status
