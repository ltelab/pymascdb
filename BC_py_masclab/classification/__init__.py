"""
MASC snowflake classification package.

This module exposes ROI feature-vector building, classifier loading, and
campaign-level inference for hydrometeor class, riming, and melting labels.

Functions for:
- Building the canonical 96-element feature vector from ROI files
- Loading MATLAB-trained logistic classifiers
- Running campaign-level classification

Last update : June 2026
Author : Baptiste Carmier
"""

from classification.classify import classify_campaign
from classification.feature_vector import build_feature_vector
from classification.load_classifier import Classifier, load_classifier

__all__ = [
    "build_feature_vector",
    "classify_campaign",
    "load_classifier",
    "Classifier",
]
