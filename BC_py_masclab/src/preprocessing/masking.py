"""
Dark-background masking for MASC snowflake images.

This module creates a dark background for MASC images by removing clutter
and applying thresholds based on camera position.
Functions for:
- Camera-dependent border discarding and luminosity thresholding

Translated from masking.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np

from src.core.config import ProcessingConfig


def masking(flake_data: np.ndarray, flake_cam: int, process: ProcessingConfig) -> np.ndarray:
    """
    Create a dark background for MASC images.

    Discards border strips according to camera position and zeroes pixels
    below the configured luminosity threshold.

    Args:
        flake_data: MASC picture (2D grayscale array).
        flake_cam: Camera ID, matched against process config camera_order
            (camera_order[0]=left LED, [1]=center, [2]=right LED).
        process: Processing configuration (camera_order, discardmat, backthresh).

    Returns:
        MASC picture with dark background mask applied
        
    Notes:
        - For the left LED camera (camera_order[0]): masks left side
        - For the right LED camera (camera_order[2]): masks right side
        - For the core triplet cameras (camera_order[:3]): masks top and bottom
        - Applies threshold to remove low luminosity pixels
    """
    # Make a copy to avoid modifying the input
    data = flake_data.copy()

    # Remove clutter according to the discard matrix based on camera position.
    # Camera position is resolved through camera_order, like MATLAB masking.m.
    # Index parity with MATLAB:
    # - top/left strips: 1:N (1-based) == [:N] (0-based) -> N pixels
    # - bottom/right strips: end-N:end covers N+1 pixels -> [-(N+1):] here
    #   (including the N == 0 case, where MATLAB still blanks the last line)

    t, b, l, r = process.discardmat[0], process.discardmat[1], process.discardmat[2], process.discardmat[3]
    camera_order = process.camera_order

    if flake_cam == camera_order[0] and l > 0:  # LED on left
        data[:, :l] = 0

    elif flake_cam == camera_order[2]:  # LED on right
        data[:, -(r + 1):] = 0

    if flake_cam in camera_order[:3]:
        if t > 0:
            data[:t, :] = 0  # Top
        data[-(b + 1):, :] = 0  # Bottom

    # Remove clutter below luminosity threshold
    data[data <= process.backthresh] = 0

    return data
