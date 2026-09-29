"""
Gaussian mixture model for blowing-snow classification.

This module loads and applies the 2-component GMM that separates precipitation
from blowing snow using four standardized image descriptors.

Functions for:
- Loading the trained GMM from a MATLAB model file
- Computing the normalized classification angle for feature vectors

Translated from classify_blowingsnow.m (Christophe Praz 2018) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import scipy.io as sio

DEFAULT_MODEL = Path(__file__).resolve().parent / "models" / "gmfit_random_all.mat"


@dataclass
class BlowingSnowGMM:
    """Two-class GMM (precipitation vs blowing snow) with feature standardization."""

    weights: np.ndarray        # (K,)
    means: np.ndarray          # (K, D)
    covariances: np.ndarray    # (K, D, D)
    mu_stand: np.ndarray       # (D,)
    sigma_stand: np.ndarray    # (D,)
    ind_bs: int                # 0-based component index of blowing snow
    ind_precip: int            # 0-based component index of precipitation

    def __post_init__(self) -> None:
        self.K, self.D = self.means.shape
        self._inv = np.array([np.linalg.inv(c) for c in self.covariances])
        self._logdet = np.array([np.linalg.slogdet(c)[1] for c in self.covariances])

    def _neg_log_joint(self, Xs: np.ndarray) -> np.ndarray:
        """Return -log P(x|k) P(k) for each GMM component (MATLAB marginal_prob)."""
        out = np.empty((Xs.shape[0], self.K))
        for k in range(self.K):
            diff = Xs - self.means[k]
            maha = np.einsum("ij,jk,ik->i", diff, self._inv[k], diff)
            log_n = -0.5 * (self.D * np.log(2 * np.pi) + self._logdet[k] + maha)
            out[:, k] = -(np.log(self.weights[k]) + log_n)
        return out

    def normalized_angle(self, X: np.ndarray) -> np.ndarray:
        """
        Compute the normalized classification angle in [0, 1].

        Applies log-transform to Dmax, standardizes features, then maps the
        ratio of precipitation vs blowing-snow marginal probabilities to an
        angle. Values near 0 indicate precipitation; near 1 indicate blowing snow.

        Args:
            X: Feature matrix (N × 4) with columns [freq, porosity, Dmax, fractal].

        Returns:
            Normalized angle array of shape (N,).

        Notes:
            - Follows Schaer et al. 2020 (AMT) convention.

        Translated from classify_blowingsnow.m (Christophe Praz 2018) and adapted for Python
        """
        X = np.asarray(X, dtype=float).copy()
        X[:, 2] = np.log(X[:, 2] + 1.0)                    # log-transform Dmax
        Xs = (X - self.mu_stand) / self.sigma_stand
        marg = self._neg_log_joint(Xs)
        return np.arctan(marg[:, self.ind_precip] / marg[:, self.ind_bs]) / (np.pi / 2)


def _extract_gmm_arrays(mat: dict, D: int, K: int = 2):
    """Recover weights, means, and covariances from an opaque MATLAB MCOS workspace."""
    raw = np.array(mat["__function_workspace__"]).ravel().tobytes()
    vals = np.frombuffer(raw[: len(raw) // 8 * 8], dtype="<f8")
    ok = np.isfinite(vals) & (np.abs(vals) > 1e-3) & (np.abs(vals) < 1e3)

    runs: List[np.ndarray] = []
    i = 0
    while i < len(vals):
        if ok[i]:
            j = i
            while j < len(vals) and ok[j]:
                j += 1
            runs.append(vals[i:j])
            i = j
        else:
            i += 1

    weights = next(r for r in runs if r.size == K and abs(r.sum() - 1.0) < 1e-3)
    means = next(r for r in runs if r.size == K * D).reshape(K, D, order="F")
    cov_flat = next(r for r in runs if r.size == D * D * K).reshape(D, D, K, order="F")
    covariances = np.array([cov_flat[:, :, k] for k in range(K)])
    return weights, means, covariances


def load_gmm(model_path: Path | str = DEFAULT_MODEL) -> BlowingSnowGMM:
    """
    Load the blowing-snow GMM from a MATLAB gmfit_*.mat model file.

    Args:
        model_path: Path to the trained model file (default: models/gmfit_random_all.mat).

    Returns:
        BlowingSnowGMM instance with weights, means, covariances, and standardization.

    Notes:
        - GMM parameters are recovered from the opaque MATLAB MCOS workspace.
        - mu_stand, sigma_stand, and component indices are read as plain arrays.

    Translated from classify_blowingsnow.m (Christophe Praz 2018) and adapted for Python
    """
    mat = sio.loadmat(str(model_path), squeeze_me=True, struct_as_record=False)
    mu_stand = np.asarray(mat["mu_stand_M1"], dtype=float)
    sigma_stand = np.asarray(mat["sigma_stand_M1"], dtype=float)
    ind_bs = int(mat["ind_bs_M1"]) - 1
    ind_precip = int(mat["ind_precip_M1"]) - 1
    weights, means, covariances = _extract_gmm_arrays(mat, D=mu_stand.size)
    return BlowingSnowGMM(
        weights=weights,
        means=means,
        covariances=covariances,
        mu_stand=mu_stand,
        sigma_stand=sigma_stand,
        ind_bs=ind_bs,
        ind_precip=ind_precip,
    )
