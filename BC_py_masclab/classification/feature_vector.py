"""
Canonical 96-element feature vector builder for MASC classifiers.

This module assembles the ordered feature vector used to train and apply
masclab logistic classifiers from a processed ROI dictionary.

Feature index (1-based, matching MATLAB comments):
  1- 8  : dimensions (area, Dmax, width, height, perim, eq_radius, area_porous, porous_ratio)
  9-27  : ellipse (E.a/b/theta, E_in, E_out, areas, ratios)
 28-40  : intensity (mean, max, range, focus, lap, std, local_std, contrast)
 41-48  : brightened image (new.*)
 49-52  : skeleton (p_ratio, A_ratio, N_ends, N_junctions)
 53-56  : Haralick (Contrast, Correlation, Energy, Homogeneity)
 57-58  : fractal (F, F_jac)
 59-61  : hull (solidity, convexity, nb_angles)
 62-66  : shape (complexity, roundness, compactness, aspect_ratio, eccentricity)
 67-69  : quality + holes (xhi, has_holes, nb_holes)
 70-73  : rectangularity (A_ratio, p_ratio, aspect_ratio, eccentricity)
    74  : circumscribed circle perim ratio
 75-80  : misc (wavs, hist_entropy, chi_wrong, chi_true, local_std, complex)
 81-96  : symmetry (Sym.P0..P10, mean, std, std/mean, P6/Pmax, idx_max)


Functions for:
- Building the 96-element feature vector from ROI descriptors
Translated from load_processed_data.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np


def _g(obj, key, default=0.0):
    """Safely extract a scalar float from a dict key, returning default if missing."""
    if obj is None:
        return float(default)
    if isinstance(obj, dict):
        val = obj.get(key)
    else:
        val = getattr(obj, key, None)

    if val is None:
        return float(default)
    if isinstance(val, np.ndarray):
        if val.size == 0:
            return float(default)
        v = float(val.ravel()[0])
        return float(default) if np.isnan(v) or np.isinf(v) else v
    try:
        v = float(val)
        return float(default) if np.isnan(v) or np.isinf(v) else v
    except (TypeError, ValueError):
        return float(default)


def _len_array(arr, default=0.0):
    """Return the number of elements in an array, or default if None."""
    if arr is None:
        return float(default)
    if isinstance(arr, np.ndarray):
        return float(arr.ravel().size)
    try:
        return float(len(arr))
    except TypeError:
        return float(default)


def build_feature_vector(roi: dict) -> np.ndarray:
    """
    Build the canonical 96-element feature vector from an ROI dictionary.

    Args:
        roi: ROI dict from process_basic_descriptors() or loaded from disk.

    Returns:
        Float64 array of shape (96,) with missing fields replaced by 0.0.

    Notes:
        - Feature ordering matches load_processed_data.m (dimensions, ellipse,
          intensity, brightened, skeleton, Haralick, fractal, hull, shape,
          quality, rectangularity, misc, symmetry).
        - Indices 81–96 are FFT symmetry components (Sym.P0–P10, mean, std, …).
    """
    g = _g  # shorthand

    # Sub-structures
    E    = roi.get('E', {})     if isinstance(roi, dict) else {}
    E_in = roi.get('E_in', {}) if isinstance(roi, dict) else {}
    E_out= roi.get('E_out',{}) if isinstance(roi, dict) else {}
    hull = roi.get('hull', {}) if isinstance(roi, dict) else {}
    skel = roi.get('skel', {}) if isinstance(roi, dict) else {}
    H    = roi.get('H', {})    if isinstance(roi, dict) else {}
    Sym  = roi.get('Sym', {})  if isinstance(roi, dict) else {}
    Rect = roi.get('Rect', {}) if isinstance(roi, dict) else {}
    C_out= roi.get('C_out',{}) if isinstance(roi, dict) else {}
    new  = roi.get('new', {})  if isinstance(roi, dict) else {}

    # Frequently reused values
    area        = g(roi,  'area',        0.0)
    area_porous = g(roi,  'area_porous', 0.0)
    perim       = g(roi,  'perim',       1.0)
    complex_    = g(roi,  'complex',     1.0)
    Ea          = g(E,    'a',           1.0)
    Eb          = g(E,    'b',           1.0)
    Ein_a       = g(E_in, 'a',           1.0)
    Ein_b       = g(E_in, 'b',           1.0)
    Eout_a      = g(E_out,'a',           1.0)
    Eout_b      = g(E_out,'b',           1.0)

    # Feature 65-66: aspect ratio and eccentricity
    if Ea >= Eb and Ea > 0:
        aspect_ratio = Eb / Ea
        eccentricity = float(np.sqrt(max(0.0, 1.0 - Eb / Ea)))
    else:
        aspect_ratio = 1.0
        eccentricity = 0.0

    # Feature 61: number of hull vertices
    xh = hull.get('xh', None) if isinstance(hull, dict) else None
    nb_hull_angles = _len_array(xh)

    # Feature 51: skel N_ends (NaN → 0)
    N_ends = g(skel, 'N_ends', 0.0)
    if np.isnan(N_ends):
        N_ends = 0.0

    # Features 81-96: FFT-based symmetry (compute_symmetry_features)
    Sym_P = [g(Sym, f'P{i}', 0.0) for i in range(11)]   # P0 .. P10
    Sym_mean = g(Sym, 'mean', 0.0)
    Sym_std  = g(Sym, 'std',  0.0)
    Sym_std_over_mean = Sym_std / Sym_mean if Sym_mean != 0.0 else 0.0

    # Feature 95: P6 / max(P1..P10)
    P1_to_10 = Sym_P[1:]
    max_P = max(P1_to_10) if any(p != 0.0 for p in P1_to_10) else 1.0
    Sym_P6_over_max = Sym_P[6] / max_P if max_P != 0.0 else 0.0

    # Feature 96: argmax index of P1..P10 (matches MATLAB idx_max - 1, i.e. 0-based)
    idx_max = float(int(np.argmax(P1_to_10)))

    # ── Assemble vector (96 elements) ─────────────────────────────────────────
    xvec = [
        # --- 1-8 : dimensions ---
        area,                                                                   # 1
        g(roi, 'Dmax',        0.0),                                             # 2
        g(roi, 'width',       0.0),                                             # 3
        g(roi, 'height',      0.0),                                             # 4
        perim,                                                                  # 5
        g(roi, 'eq_radius',   0.0),                                             # 6
        area_porous,                                                            # 7
        (area - area_porous) / area if area > 0 else 0.0,                      # 8

        # --- 9-27 : ellipse ---
        Ea,                                                                     # 9
        Eb,                                                                     # 10
        g(E,    'theta',      0.0),                                             # 11
        Ein_a,                                                                  # 12
        Ein_b,                                                                  # 13
        Eout_a,                                                                 # 14
        Eout_b,                                                                 # 15
        np.pi * Ein_a  * Ein_b,                                                 # 16
        np.pi * Eout_a * Eout_b,                                                # 17
        np.pi * Ea     * Eb,                                                    # 18
        Ea    / Eout_a  if Eout_a > 0 else 0.0,                                # 19
        Eb    / Eout_b  if Eout_b > 0 else 0.0,                                # 20
        Ein_a / Eout_a  if Eout_a > 0 else 0.0,                                # 21
        Ein_b / Eout_b  if Eout_b > 0 else 0.0,                                # 22
        Ein_a / Ea      if Ea     > 0 else 0.0,                                # 23
        Ein_b / Eb      if Eb     > 0 else 0.0,                                # 24
        (Ein_a * Ein_b) / (Eout_a * Eout_b) if Eout_a * Eout_b > 0 else 0.0, # 25
        (Ein_a * Ein_b) / (Ea     * Eb)      if Ea     * Eb     > 0 else 0.0, # 26
        (Ea    * Eb)    / (Eout_a * Eout_b)  if Eout_a * Eout_b > 0 else 0.0, # 27

        # --- 28-40 : intensity ---
        g(roi, 'mean_intens',  0.0),                                            # 28
        g(roi, 'max_intens',   0.0),                                            # 29
        g(roi, 'range_intens', 0.0),                                            # 30
        g(roi, 'focus',        0.0),                                            # 31
        g(roi, 'area_focus',   0.0),                                            # 32
        g(roi, 'area_range',   0.0),                                            # 33
        g(roi, 'lap',          0.0),                                            # 34
        g(roi, 'area_lap',     0.0),                                            # 35
        g(roi, 'std',          0.0),                                            # 36
        g(roi, 'local_std',    0.0),                                            # 37
        g(roi, 'local_std5',   0.0),                                            # 38
        g(roi, 'local_std7',   0.0),                                            # 39
        g(roi, 'contrast',     0.0),                                            # 40

        # --- 41-48 : brightened image (new.*) ---
        g(new, 'range_intens', 0.0),                                            # 41
        g(new, 'lap',          0.0),                                            # 42
        g(new, 'area_lap',     0.0),                                            # 43
        g(new, 'std',          0.0),                                            # 44
        g(new, 'local_std',    0.0),                                            # 45
        g(new, 'local_std5',   0.0),                                            # 46
        g(new, 'local_std7',   0.0),                                            # 47
        g(new, 'contrast',     0.0),                                            # 48

        # --- 49-52 : skeleton ---
        g(skel, 'p_ratio',     0.0),                                            # 49
        g(skel, 'A_ratio',     0.0),                                            # 50
        N_ends,                                                                 # 51
        g(skel, 'N_junctions', 0.0),                                            # 52

        # --- 53-56 : Haralick ---
        g(H, 'Contrast',       0.0),                                            # 53
        g(H, 'Correlation',    0.0),                                            # 54
        g(H, 'Energy',         0.0),                                            # 55
        g(H, 'Homogeneity',    0.0),                                            # 56

        # --- 57-58 : fractal ---
        g(roi, 'F',            0.0),                                            # 57
        g(roi, 'F_jac',        0.0),                                            # 58

        # --- 59-61 : hull ---
        g(hull, 'solidity',    0.0),                                            # 59
        g(hull, 'convexity',   0.0),                                            # 60
        nb_hull_angles,                                                         # 61

        # --- 62-66 : shape ---
        max(1.0, complex_),                                                     # 62
        g(roi, 'roundness',    0.0),                                            # 63
        g(roi, 'compactness',  0.0),                                            # 64
        aspect_ratio,                                                           # 65
        eccentricity,                                                           # 66

        # --- 67-69 : quality + holes ---
        g(roi, 'xhi',          0.0),                                            # 67
        1.0 if g(roi, 'nb_holes', 0.0) > 0 else 0.0,                           # 68
        g(roi, 'nb_holes',     0.0),                                            # 69

        # --- 70-73 : rectangularity ---
        g(Rect, 'A_ratio',     0.0),                                            # 70
        g(Rect, 'p_ratio',     0.0),                                            # 71
        g(Rect, 'aspect_ratio',0.0),                                            # 72
        g(Rect, 'eccentricity',0.0),                                            # 73

        # --- 74 : circumscribed circle ---
        (2.0 * np.pi * g(C_out, 'r', 0.0)) / perim if perim > 0 else 0.0,     # 74

        # --- 75-80 : misc ---
        g(roi, 'wavs',         0.0),                                            # 75
        g(roi, 'hist_entropy', 0.0),                                            # 76
        max(1.0, complex_) * (1.0 + g(roi, 'range_intens', 0.0)),              # 77
        max(1.0, complex_) * (1.0 + g(roi, 'local_std7',   0.0)),              # 78
        g(roi, 'local_std',    0.0),                                            # 79
        complex_,                                                               # 80

        # --- 81-96 : symmetry (FFT-based) ---
        Sym_P[0],   Sym_P[1],  Sym_P[2],  Sym_P[3],  Sym_P[4],                # 81-85
        Sym_P[5],   Sym_P[6],  Sym_P[7],  Sym_P[8],  Sym_P[9],                # 86-90
        Sym_P[10],                                                              # 91
        Sym_mean,                                                               # 92
        Sym_std,                                                                # 93
        Sym_std_over_mean,                                                      # 94
        Sym_P6_over_max,                                                        # 95
        idx_max,                                                                # 96
    ]

    return np.array(xvec, dtype=np.float64)
