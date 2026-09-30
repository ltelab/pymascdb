"""
Blowing-snow classification from MASC image features.

This module applies a trained GMM to per-image descriptors and produces
per-image and per-triplet precipitation vs blowing-snow labels.

Functions for:
- GMM-based per-image classification
- Triplet-level label aggregation and mixed-class flagging

Translated from:
- classify_blowingsnow.m (Christophe Praz 2018)
- compute_labelling_final.m (Mathieu Schaer 2017)
and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from blowing_snow.gmm_bs import BlowingSnowGMM, load_gmm

# Decision thresholds on the normalized angle (Schaer et al. 2020, AMT)
TRESH_PRECIP = 0.193
TRESH_BS = 0.881

# Feature columns feeding the GMM, in order [freq, porosity, Dmax, fractal]
FEATURE_COLS = ["freq", "porosity", "dmax", "fractal"]


def _flag_mixed(angle: np.ndarray) -> np.ndarray:
    """Rescaled mixed-class probability (NaN outside the mixed band)."""
    flag = np.full(angle.shape, np.nan)
    flag[(angle > TRESH_PRECIP) & (angle < TRESH_BS)] = 1.0

    p_precip = (angle - TRESH_PRECIP) / (2 * (0.5 - TRESH_PRECIP))
    lo = angle < 0.5
    flag[lo] = flag[lo] * p_precip[lo]

    p_bs = (angle - (1 - TRESH_BS)) / (2 * (TRESH_BS - 0.5))
    hi = angle >= 0.5
    flag[hi] = flag[hi] * p_bs[hi]
    return flag


def _combine_triplets(df: pd.DataFrame) -> pd.DataFrame:
    """Average normalized angle over the three camera views per flake ID and minute."""
    key = df["date_vec"].dt.floor("min")
    grp = df.assign(_min=key).groupby(["_min", "ID"], sort=True)
    out = pd.DataFrame(
        {
            "date_vec": grp["date_vec"].first().values,
            "ID": grp["ID"].first().values,
            "Normalized_Angle": grp["Normalized_Angle"].mean().values,
        }
    )
    out["Label"] = (out["Normalized_Angle"] > 0.5).astype(int)
    out["Flag_mixed"] = _flag_mixed(out["Normalized_Angle"].to_numpy())
    return out[["date_vec", "ID", "Label", "Normalized_Angle", "Flag_mixed"]]


def classify_blowing_snow(
    features: pd.DataFrame,
    gmm: BlowingSnowGMM | None = None,
    model_path: Path | str | None = None,
    show_stats: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Classify MASC image features as precipitation (0) or blowing snow (1).

    Computes the GMM normalized angle for each image, assigns binary labels,
    flags mixed cases, and aggregates labels at triplet level.

    Args:
        features: DataFrame with columns date_vec, ID, cam, freq, porosity,
            dmax, fractal (output of features_bs.extract_all).
        gmm: Pre-loaded GMM, or None to load from model_path or the default model.
        model_path: Path to gmfit_*.mat when gmm is None.
        show_stats: Whether to print classification statistics to stdout.

    Returns:
        Tuple of (all_df, triplet_df) per-image and per-triplet result tables.

    Notes:
        - Label 0 = precipitation, 1 = blowing snow (angle > 0.5).
        - Triplet table averages the angle over the three camera views per minute.

    Translated from classify_blowingsnow.m (Christophe Praz 2018) and adapted for Python
    """
    if gmm is None:
        gmm = load_gmm(model_path) if model_path else load_gmm()

    X = features[FEATURE_COLS].to_numpy(dtype=float)
    angle = gmm.normalized_angle(X)
    labels = (angle > 0.5).astype(int)

    all_df = pd.DataFrame(
        {
            "date_vec": features["date_vec"].values,
            "ID": features["ID"].values,
            "Cam": features["cam"].values,
            "Label": labels,
            "Normalized_Angle": angle,
            "Flag_mixed": _flag_mixed(angle),
        }
    )
    triplet_df = _combine_triplets(all_df)

    if show_stats:
        _print_stats(all_df, triplet_df)
    return all_df, triplet_df


def _print_stats(all_df: pd.DataFrame, triplet_df: pd.DataFrame) -> None:
    """Print blowing-snow vs precipitation statistics for image and triplet tables."""
    a = all_df["Normalized_Angle"].to_numpy()
    u = triplet_df["Normalized_Angle"].to_numpy()
    print(
        "Percentage of Blowing Snow images: %.2f%% (triplets independent), "
        "%.2f%% (triplets combined)"
        % (100 * np.mean(a > 0.5), 100 * np.mean(u > 0.5))
    )
    print(
        "Pure BS: %.2f%%, Pure Precip: %.2f%% (triplets independent)"
        % (100 * np.mean(a >= TRESH_BS), 100 * np.mean(a <= TRESH_PRECIP))
    )
    print(
        "Pure BS: %.2f%%, Pure Precip: %.2f%% (triplets combined)"
        % (100 * np.mean(u >= TRESH_BS), 100 * np.mean(u <= TRESH_PRECIP))
    )
