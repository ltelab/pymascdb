"""
Local neighborhood filters for MASC image analysis.

This module provides MATLAB-compatible range and standard-deviation filters
used in ROI quality and texture metrics.

Functions for:
- Local range filtering (max − min in a window)
- Local standard-deviation filtering

Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
from scipy import ndimage


def rangefilt(image: np.ndarray, size: int = 3) -> np.ndarray:
    """
    Compute the local range (max − min) in a square neighborhood.

    Args:
        image: Input grayscale image.
        size: Side length of the square neighborhood.

    Returns:
        Range-filtered image with the same shape as the input.

    Notes:
        - Equivalent to MATLAB rangefilt with a default 3×3 neighborhood.
    """
    def range_func(values):
        return values.max() - values.min()

    rangefilt_im = ndimage.generic_filter(image, range_func, size=size, mode='nearest')
    return rangefilt_im


def stdfilt(image: np.ndarray, size: int = 3) -> np.ndarray:
    """
    Compute the local standard deviation in a square neighborhood.

    Args:
        image: Input grayscale image.
        size: Side length of the square neighborhood.

    Returns:
        Standard-deviation-filtered image with the same shape as the input.

    Notes:
        - Equivalent to MATLAB stdfilt: sample std (N-1) in each neighborhood.
        - Border padding uses reflect (matches masclab GT; nearest underestimates
          local_std7 on small crops).
    """
    def std_func(values):
        return np.std(values, ddof=1)

    stdfilt_im = ndimage.generic_filter(
        image.astype(float), std_func, size=size, mode='reflect'
    )
    return stdfilt_im
