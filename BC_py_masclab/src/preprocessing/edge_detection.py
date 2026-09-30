"""
Edge detection and morphological cleanup for MASC snowflake images.

This module provides functions for detecting and enhancing edges
in MASC snowflake images.
Functions for:
- Sobel edge detection with dilation, hole filling, and erosion

Status: LEGACY, not used by the current pipeline. The MATLAB pipeline
(MASC_picture_process.m) replaced edge_detection.m with a simple binary
mask (data > 0) followed by hole filling; single.py and triplet.py mirror
that behaviour. This module is kept as a faithful translation of
edge_detection.m for reference and possible reuse.

Translated from edge_detection.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import cv2
import matplotlib.pyplot as plt

from scipy import ndimage


def _extend_bw_mask(mask: np.ndarray, pad_size: int) -> np.ndarray:
    """Pad a binary mask on all sides to reduce border effects."""
    return np.pad(mask, pad_size, mode='edge')


def _shrink_bw_mask(mask: np.ndarray, pad_size: int) -> np.ndarray:
    """Remove symmetric padding added by ``_extend_bw_mask``."""
    if pad_size == 0:
        return mask
    return mask[pad_size:-pad_size, pad_size:-pad_size]


def detect_edges(
    image: np.ndarray,
    strel_size: int = 5,
    sobel_threshold: float = 0.008,
    illustration: bool = False
) -> np.ndarray:
    """
    Detect edges and filled areas in a MASC picture.

    Extends the image to avoid border effects, applies Sobel edge detection,
    dilates with line structuring elements, fills holes, erodes, then crops
    back to the original size.

    Args:
        image: Input grayscale image (2D uint8 array).
        strel_size: Structuring-element length for morphological operations.
        sobel_threshold: Sobel gradient threshold in [0, 1] (MATLAB convention).
        illustration: Whether to display intermediate processing steps.

    Returns:
        Binary mask with detected edges filled and eroded.

     Notes:
        - The MATLAB version uses strel('line', strel_size, angle) which creates line-shaped structuring elements
        - In OpenCV, we use cv2.getStructuringElement with appropriate shapes
        - Sobel threshold of 0.008 is scaled to 0-255 range (~ 2.04)

    Translated from edge_detection.m (Christophe Praz 2015) and adapted for Python
    """
    # Ensure input is grayscale uint8
    if image.dtype != np.uint8:
        image = image.astype(np.uint8)

    # Temporary extend mask to avoid border effects
    pad_size = 10 * strel_size
    extended = _extend_bw_mask(image, pad_size)

    # Detect edges using Sobel
    # MATLAB edge(I, 'Sobel', threshold) with threshold in [0,1]
    # OpenCV Sobel returns gradient magnitude, we need to threshold it
    sobelx = cv2.Sobel(extended, cv2.CV_64F, 1, 0, ksize=3)
    sobely = cv2.Sobel(extended, cv2.CV_64F, 0, 1, ksize=3)

    # Compute gradient magnitude
    gradient_magnitude = np.sqrt(sobelx**2 + sobely**2)

    # Normalize to [0, 1] range (MATLAB edge normalizes by max possible gradient)
    # For 8-bit images, max gradient is ~1020 (4*255 for Sobel)
    gradient_magnitude = gradient_magnitude / gradient_magnitude.max() if gradient_magnitude.max() > 0 else gradient_magnitude

    # Apply threshold to create binary edge map
    edged = (gradient_magnitude > sobel_threshold).astype(np.uint8)

    # Create line structuring elements
    # MATLAB strel('line', length, angle) creates line at specified angle
    # Angle 0 = horizontal, angle 90 = vertical
    se0 = cv2.getStructuringElement(cv2.MORPH_RECT, (strel_size, 1))   # horizontal line
    se90 = cv2.getStructuringElement(cv2.MORPH_RECT, (1, strel_size))  # vertical line

    # Dilate with both structuring elements
    # MATLAB: imdilate(edged, [se0 se90]) applies both sequentially
    dilated = cv2.dilate(edged, se0, iterations=1)
    dilated = cv2.dilate(dilated, se90, iterations=1)

    # Fill holes
    # MATLAB imfill(BW, 'holes') fills holes in binary image
    filled = ndimage.binary_fill_holes(dilated).astype(np.uint8) # From scipy

    # Erode with both structuring elements
    eroded = cv2.erode(filled, se0, iterations=1)
    eroded = cv2.erode(eroded, se90, iterations=1)

    # Shrink back to original size
    data_out = _shrink_bw_mask(eroded, pad_size)

    # Illustration if desired
    if illustration:
        fig, axes = plt.subplots(3, 2, figsize=(10, 12))

        axes[0, 0].imshow(image, cmap='gray')
        axes[0, 0].set_title('raw')
        axes[0, 0].axis('off')

        axes[0, 1].imshow(extended, cmap='gray')
        axes[0, 1].set_title('extended')
        axes[0, 1].axis('off')

        axes[1, 0].imshow(edged, cmap='gray')
        axes[1, 0].set_title('edges')
        axes[1, 0].axis('off')

        axes[1, 1].imshow(dilated, cmap='gray')
        axes[1, 1].set_title('dilated')
        axes[1, 1].axis('off')

        axes[2, 0].imshow(filled, cmap='gray')
        axes[2, 0].set_title('filled')
        axes[2, 0].axis('off')

        axes[2, 1].imshow(data_out, cmap='gray')
        axes[2, 1].set_title('final')
        axes[2, 1].axis('off')

        plt.tight_layout()
        plt.show()

    return data_out
