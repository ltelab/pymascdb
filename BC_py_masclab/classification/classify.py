"""
Campaign-level MASC snowflake classification.

This module loads trained logistic classifiers, builds feature matrices from
ROI files, runs inference, and writes label fields back into each ROI file.

Functions for:
- Recursive ROI discovery and feature-matrix assembly
- Skewness-driven feature transforms and logistic inference
- Campaign-wide classification with optional joblib parallelism

Translated from:
- predict_snowflakes_class.m (Christophe Praz 2015)
- predict_snowflakes_riming.m (Christophe Praz 2015)
- predict_snowflakes_melting.m (Christophe Praz 2015)
- make_predictions_for_campaign.m (Christophe Praz 2015)
and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import logging
import numpy as np

from pathlib import Path
from typing import List, Tuple

import joblib as _joblib

from classification.load_classifier import Classifier, load_classifier
from classification.feature_vector import build_feature_vector
from src.dataio.load_roi_data import load_roi_data, save_roi_data

logger = logging.getLogger(__name__)


def _find_roi_files(outdir: Path, fmt: str) -> List[Path]:
    """Recursively find ROI files with the expected extension under outdir."""
    ext_map = {'joblib': '*.joblib',
               'pkl': '*.pkl',
               'mat': '*.mat'}
    pattern = ext_map.get(fmt, '*.mat')
    all_files = sorted(outdir.rglob(pattern))
    # Keep only files whose stem starts with a digit
    roi_files = [f for f in all_files if f.stem[0].isdigit()]
    return roi_files


def _transform_features(X: np.ndarray, skew: np.ndarray) -> np.ndarray:
    """Apply per-feature pre-normalisation transforms driven by skewness."""
    X = X.copy().astype(float)
    for i, s in enumerate(skew):
        if s > 1.0:
            # MATLAB (predict_snowflakes_*.m / build_classifier.m): log(abs(X+1)),
            # NOT log(abs(X)+1) — differs for negative features (e.g. E.theta).
            X[:, i] = np.log(np.abs(X[:, i] + 1.0))
        elif s > 0.75:
            X[:, i] = np.sqrt(np.abs(X[:, i]))
        elif s < -1.0:
            X[:, i] = np.exp(X[:, i])
        elif s < -0.75:
            X[:, i] = X[:, i] ** 2
    return X


def _predict(X: np.ndarray, cl: Classifier) -> Tuple[np.ndarray, np.ndarray]:
    """Run logistic inference and return 1-based class indices and probabilities."""
    if cl.type_ != 'logistic':
        raise NotImplementedError(
            f"Only 'logistic' classifiers are supported; got '{cl.type_}'"
        )

    N = X.shape[0]
    # Prepend bias column  [1 | X]  → (N, D+1)
    tX = np.hstack([np.ones((N, 1)), X])

    if cl.type_classif == 'multiclass':
        # Numerically-stable softmax (equivalent to MATLAB exp(tX*model) / sum)
        logits = tX @ cl.model                                       # (N, K)
        logits -= logits.max(axis=1, keepdims=True)
        exp_l  = np.exp(logits)
        scores = exp_l / exp_l.sum(axis=1, keepdims=True)
        pred   = np.argmax(scores, axis=1) + 1                      # 1-based

    elif cl.type_classif == 'binary':
        # Sigmoid  (MATLAB: round(sigmoid(tX*model)))
        raw    = tX @ cl.model                                       # (N, 1) or (N,)
        scores = 1.0 / (1.0 + np.exp(-raw))
        scores = scores.reshape(N, -1)                               # (N, 1)
        # MATLAB keeps binary predictions 0-based (melting_ID is 0 or 1,
        # written as-is by predict_snowflakes_melting.m)
        pred   = np.round(scores[:, 0]).astype(int)                  # 0 or 1

    else:
        raise ValueError(f"Unknown type_classif: {cl.type_classif}")

    return pred, scores


def classify_campaign(
    outdir: Path,
    classifier_path: Path,
    save_format: str = 'mat',
    blurry_threshold: float = 0.0,
    label_prefix: str = 'label',
    n_jobs: int = 1,
) -> None:
    """
    Classify all ROI files in outdir and write predicted labels back to disk.

    Loads a MATLAB-trained logistic classifier, builds the N×96 feature matrix,
    selects and transforms features, runs inference, and updates each ROI file
    with ID, name, and probability fields.

    Args:
        outdir: Root output directory searched recursively for ROI files.
        classifier_path: Path to the MATLAB .mat classifier file.
        save_format: ROI file format ('mat', 'joblib', or 'pkl').
        blurry_threshold: Minimum xhi to classify; below this, label_ID = -9.
        label_prefix: Field prefix ('label', 'riming', or 'melting').
        n_jobs: Number of parallel workers for feature extraction (1 = sequential).

    Returns:
        None

    Notes:
        - label prefix 'label'   → label_ID, label_name, label_probs
        - label prefix 'riming'  → riming_ID, riming_name, riming_probs
        - label prefix 'melting' → melting_ID, melting_name, melting_probs

    Translated from predict_snowflakes_class.m (Christophe Praz 2015) and adapted for Python
    """
    outdir = Path(outdir)

    print(f"  Loading classifier: {classifier_path.name}")
    cl = load_classifier(classifier_path)
    print(f"  Classes: {cl.N_labels}")
    print(f"  Features used: {len(cl.feat_vec)} (out of 96)")

    # ── Find files ────────────────────────────────────────────────────────────
    roi_files = _find_roi_files(outdir, save_format)
    if not roi_files:
        print(f"  No ROI files found in {outdir} (format={save_format})")
        return
    print(f"  Found {len(roi_files)} ROI files")

    # ── Build feature matrix  X_full (N × 96) ────────────────────────────────
    print("  Building feature matrix …", end="", flush=True)

    def _fv(f):
        try:
            return build_feature_vector(load_roi_data(f, save_format))
        except Exception as e:
            logger.warning("Could not extract features from %s: %s", f, e)
            return np.zeros(96)

    if n_jobs == 1:
        rows = [_fv(f) for f in roi_files]
    else:
        rows = _joblib.Parallel(n_jobs=n_jobs)(
            _joblib.delayed(_fv)(f) for f in roi_files
        )

    X_full = np.array(rows, dtype=float)   # (N, 96)
    print(f" done  ({X_full.shape})")

    # ── Select features, transform, standardise ───────────────────────────────
    X = X_full[:, cl.feat_vec]             # (N, D)
    X = _transform_features(X, cl.skew_)

    if cl.normalization == 'standardization':
        denom = np.where(cl.std_ > 0, cl.std_, 1.0)
        X = (X - cl.mean_) / denom
    else:
        print(f"  Warning: normalization type '{cl.normalization}' not implemented; skipping.")

    # ── Predict ───────────────────────────────────────────────────────────────
    print("  Running inference …", end="", flush=True)
    pred, scores = _predict(X, cl)
    print(" done")

    # ── Write labels back into each ROI file ──────────────────────────────────
    print("  Saving labels …", end="", flush=True)
    n_fine = 0
    n_blurry = 0

    id_key    = f'{label_prefix}_ID'
    name_key  = f'{label_prefix}_name'
    probs_key = f'{label_prefix}_probs'

    for i, f in enumerate(roi_files):
        try:
            roi = load_roi_data(f, save_format)
            xhi = float(roi.get('xhi', 0.0)) if isinstance(roi, dict) else 0.0

            if xhi > blurry_threshold:
                # Multiclass predictions are 1-based (MATLAB argmax), binary
                # ones are 0-based (MATLAB round(sigmoid)); N_labels lookup
                # mirrors prediction_scheme{pred} / prediction_scheme{pred+1}.
                name_idx = pred[i] if cl.type_classif == 'binary' else pred[i] - 1
                roi[id_key]    = int(pred[i])
                roi[name_key]  = cl.N_labels[name_idx]
                roi[probs_key] = scores[i].ravel()
                n_fine += 1
            else:
                roi[id_key]    = -9
                roi[name_key]  = 'blurry'
                roi[probs_key] = scores[i].ravel()
                n_blurry += 1

            save_roi_data(roi, f, format=save_format, compress=3)
        except Exception as e:
            logger.warning("Could not write labels to %s: %s", f, e)

    print(f" done  ({n_fine} classified, {n_blurry} blurry)")
