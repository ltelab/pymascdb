"""
Blowing-snow classification package for MASC campaigns.

This module exposes feature extraction, GMM loading, and blowing-snow
classification utilities ported from the SharedSchaerMBL MATLAB library.

Last update : June 2026
Author : Baptiste Carmier
"""

from blowing_snow.classify_bs import classify_blowing_snow
from blowing_snow.gmm_bs import BlowingSnowGMM, load_gmm

__all__ = ["classify_blowing_snow", "load_gmm", "BlowingSnowGMM"]
