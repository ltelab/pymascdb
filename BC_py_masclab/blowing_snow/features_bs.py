"""
Feature extraction for blowing-snow classification.

This module reads raw MASC images, subtracts a per-camera median sky filter,
binarises residuals, and extracts descriptors used by the blowing-snow GMM.

Functions for:
- Campaign-wide feature extraction from PNG image trees
- Per-block median filtering and binary descriptor computation
- Photo frequency and sky-mask loading

Translated from:
- extract_all_final.m (Mathieu Schaer 2017)
- extract_features_final.m (Mathieu Schaer 2017)
- photo_frequency.m (Mathieu Schaer 2017)
- importfile_MASC.m (Mathieu Schaer 2017)
and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import ConvexHull, QhullError
from skimage import measure, morphology

from src.utils.descriptors_heplers import compute_Dmax

# Extraction parameters (extract_features_final.m)
THRES = 0.02
FACTOR_FILTER = 1.1
SE_CLOSING = 11
RATIO = 1.1
FACTOR_INCR = 0.4
THRES_INCR = 0.03
AREA_OPEN = 80000
SKY_PIXELS = 80000

_MODELS_DIR = Path(__file__).resolve().parent / "models"


def _matlab_quantile(x: np.ndarray, q: float) -> float:
    """MATLAB-compatible quantile using Hazen plotting positions."""
    if x.size == 0:
        return 0.0
    return float(np.quantile(x, q, method="hazen"))


def _matlab_strel_disk(radius: int, n: int = 4) -> np.ndarray:
    """
    Flat disk structuring element matching MATLAB strel('disk', radius, n).

    Default n=4 uses the IPT periodic-line decomposition (not a Euclidean disk).
    Verified against getnhood(strel('disk', 5)) from MATLAB.
    """
    radius = int(radius)
    n = int(n)
    if radius <= 0:
        return np.ones((1, 1), dtype=bool)
    if n == 0:
        y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
        return (x * x + y * y) <= radius * radius

    # radialSectionLength = R / (csc(pi/N) + cot(pi/(2*N)))
    radial = radius / (1.0 / np.sin(np.pi / n) + 1.0 / np.tan(np.pi / (2 * n)))
    p = int(np.floor(radial + 0.5))
    if p < 1:
        y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
        return (x * x + y * y) <= radius * radius

    canvas_r = radius + p * n + 2
    center = np.zeros((2 * canvas_r + 1, 2 * canvas_r + 1), dtype=bool)
    center[canvas_r, canvas_r] = True
    out = center.copy()
    for k in range(n):
        angle = k * np.pi / n
        vr = int(np.round(np.sin(angle)))
        vc = int(np.round(np.cos(angle)))
        g = np.gcd(abs(vr), abs(vc)) or 1
        vr //= g
        vc //= g
        pts = [(i * vr, i * vc) for i in range(-p, p + 1)]
        ys = [a for a, _ in pts]
        xs = [b for _, b in pts]
        nhood = np.zeros((max(ys) - min(ys) + 1, max(xs) - min(xs) + 1), dtype=bool)
        for y, x in pts:
            nhood[y - min(ys), x - min(xs)] = True
        out = ndimage.binary_dilation(out, structure=nhood)

    coords = np.argwhere(out)
    y0, x0 = coords.min(axis=0)
    y1, x1 = coords.max(axis=0)
    return out[y0 : y1 + 1, x0 : x1 + 1]


# Cached SE for imclose(..., strel('disk', SE_CLOSING)) — MATLAB default N=4
_CLOSING_SE = _matlab_strel_disk(SE_CLOSING, n=4)


def _matlab_regionprops_perimeter(region_mask: np.ndarray) -> float:
    """
    MATLAB regionprops(..., 'Perimeter'): Euclidean length of the outer boundary.

    Walks the CHAIN_APPROX_NONE contour (bwboundaries-like) and sums segment lengths.
    """
    m = np.ascontiguousarray(region_mask.astype(np.uint8))
    if not m.any():
        return 0.0
    contours, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return 0.0
    cnt = max(contours, key=len).reshape(-1, 2).astype(np.float64)
    if len(cnt) < 2:
        return float(len(cnt))
    pts = np.vstack([cnt, cnt[:1]])
    delta = np.diff(pts, axis=0)
    return float(np.sqrt((delta * delta).sum(axis=1)).sum())


def parse_meta(name: str) -> Tuple[pd.Timestamp, float, int]:
    """
    Parse datetime, flake ID, and camera angle from a MASC filename.

    Args:
        name: Image filename (YYYY.MM.DD_HH.MM.SS_..._camX.png).

    Returns:
        Tuple of (datetime, flake_id, cam_angle).

    Translated from importfile_MASC.m (Mathieu Schaer 2017) and adapted for Python
    """
    date = pd.to_datetime(name[:19], format="%Y.%m.%d_%H.%M.%S")
    flake_id = float(name.split("_")[3])      # token between 3rd and 4th '_'
    cam = int(name[-5])                        # char before '.png'
    return date, flake_id, cam


def photo_frequency(dates: pd.DatetimeIndex, cam_angle: np.ndarray, w: int) -> np.ndarray:
    """
    Compute the mean number of photos per minute in a moving window.

    Args:
        dates: Image timestamps.
        cam_angle: Camera angle array (0, 1, 2).
        w: Half window size in number of images.

    Returns:
        Photo frequency array (images per minute per camera).

    Translated from photo_frequency.m (Mathieu Schaer 2017) and adapted for Python
    """
    n = len(dates)
    if n < 2 * w:
        w = n // 2
    freq = np.zeros(n)
    for i in range(n):
        lo = max(0, i - w)
        hi = min(n - 1, i + w)
        dt = (dates[hi] - dates[lo]).total_seconds() / 60.0
        counts = [(cam_angle[lo:hi + 1] == c).sum() for c in (0, 1, 2)]
        freq[i] = max(counts) / dt if dt > 0 else 0.0
    return freq


def load_masks(davos: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load the three per-camera sky masks as boolean arrays.

    Args:
        davos: Use MASK_davos.mat instead of MASK.mat when True.

    Returns:
        Tuple of (MASK_0, MASK_1, MASK_2) boolean arrays.
    """
    import scipy.io as sio

    fname = "MASK_davos.mat" if davos else "MASK.mat"
    m = sio.loadmat(str(_MODELS_DIR / fname))
    return m["MASK_0"] > 0, m["MASK_1"] > 0, m["MASK_2"] > 0


def _saturating_sub(image: np.ndarray, value: np.ndarray) -> np.ndarray:
    """uint8 subtraction with rounding and 0–255 saturation (MATLAB style)."""
    out = np.rint(image.astype(np.float64) - np.asarray(value, dtype=np.float64))
    return np.clip(out, 0, 255).astype(np.uint8)


def _filter_image(image, median, filter_sum, mask, reflects):
    """Subtract the median sky filter and binarise one image."""
    factor = FACTOR_FILTER
    thres = THRES
    cond = True

    # Sun is out: image globally much brighter than the median
    if filter_sum > 0 and image.sum() / filter_sum > RATIO:
        factor = FACTOR_FILTER + FACTOR_INCR
        thres = THRES + THRES_INCR
        cond = False

    filtered = _saturating_sub(image, factor * median.astype(np.float64))

    same_size = image.shape == mask.shape
    if same_size:
        # Residual in the sky only (low sun): correct factor/thres on the sky
        if cond and int((filtered[mask] > THRES).sum()) > SKY_PIXELS:
            factor_arr = np.where(mask, FACTOR_FILTER + FACTOR_INCR, FACTOR_FILTER)
            thres = np.where(mask, THRES + THRES_INCR, THRES)
            filtered = _saturating_sub(image, factor_arr * median.astype(np.float64))
        # Remove fixed light reflects
        for r0, r1, c0, c1 in reflects:
            filtered[r0:r1, c0:c1] = 0

    binary = (filtered.astype(np.float64) / 255.0) > thres
    # imclose(Binary, strel('disk', se_closing)) — MATLAB default N=4 disk
    binary = morphology.binary_closing(binary, _CLOSING_SE)
    binary = ndimage.binary_fill_holes(binary)
    binary = morphology.remove_small_objects(binary, min_size=3, connectivity=2)
    big = morphology.remove_small_objects(binary, min_size=AREA_OPEN, connectivity=2)
    return binary & ~big


# Fixed light-reflect rectangles per cam (row0, row1, col0, col1), 0-based
_REFLECTS = {
    0: [(992, 1048, 89, 280)],
    1: [],
    2: [(329, 380, 2174, 2345), (989, 1052, 2172, 2365)],
}


def _describe(binary: np.ndarray) -> Tuple[float, float, float, np.ndarray]:
    """Return porosity, Dmax q0.7, min fractal index squared, and all Dmax values."""
    labels = measure.label(binary, connectivity=2)
    dmax_all: List[float] = []
    fractal: List[float] = []
    for r in measure.regionprops(labels):
        ys, xs = r.coords[:, 0], r.coords[:, 1]
        pts = np.column_stack([xs, ys]).astype(float)
        try:
            v = pts[ConvexHull(pts).vertices]
        except (QhullError, ValueError):
            v = pts
        dmax_all.append(compute_Dmax(v[:, 0], v[:, 1])[0])
        if r.area > 1:
            # regionprops Perimeter (MATLAB) — not skimage's perimeter attribute
            perim = _matlab_regionprops_perimeter(r.image)
            fractal.append((2 * np.log(0.25 * perim) / np.log(r.area)) ** 2)

    porosity = float(ndimage.distance_transform_edt(~binary).sum())
    dmax = _matlab_quantile(np.asarray(dmax_all), 0.7)
    frac = _matlab_quantile(np.asarray(fractal), 0.0) if fractal else 0.0
    return porosity, dmax, frac, np.asarray(dmax_all)


def extract_features_block(images, cams, masks) -> List[Tuple[float, float, float, np.ndarray]]:
    """
    Extract blowing-snow descriptors for one block of mixed-camera images.

    Args:
        images: List of grayscale images (uint8).
        cams: Camera angle for each image (0, 1, 2).
        masks: Dict mapping cam angle to sky mask boolean array.

    Returns:
        List of (porosity, dmax, fractal, dmax_all) tuples, one per image.

    Translated from extract_features_final.m (Mathieu Schaer 2017) and adapted for Python
    """
    cams = np.asarray(cams)
    medians, sums = {}, {}
    for c in (0, 1, 2):
        sel = [im for im, cc in zip(images, cams) if cc == c]
        if sel:
            med = np.median(np.stack(sel, axis=2), axis=2)
            medians[c] = med
            sums[c] = float(med.sum())

    rows = []
    for image, c in zip(images, cams):
        if c not in medians:
            rows.append((0.0, 0.0, 0.0, np.array([])))
            continue
        binary = _filter_image(image, medians[c], sums[c], masks[c], _REFLECTS.get(c, []))
        rows.append(_describe(binary) if binary.any() else (0.0, 0.0, 0.0, np.array([])))
    return rows


def extract_all(
    folder: str | Path,
    block_size: int,
    w: int,
    s: float,
    tstart,
    tstop,
    davos: bool = False,
) -> pd.DataFrame:
    """
    Extract blowing-snow features from all MASC PNG images under a folder.

    Splits images into continuous events, processes them in blocks for median
    sky filtering, and returns a DataFrame ready for classify_blowing_snow.

    Args:
        folder: Root folder containing raw MASC PNG images.
        block_size: Number of images per median-filter block.
        w: Half window size for photo frequency (images).
        s: Time gap threshold in hours to separate continuous events.
        tstart: Start of the extraction time window (inclusive).
        tstop: End of the extraction time window (inclusive).
        davos: Use Davos sky masks (MASK_davos.mat) when True.

    Returns:
        DataFrame with columns date_vec, ID, cam, freq, porosity, dmax,
        fractal, dmax_all (rows with porosity == 0 removed).

    Notes:
        - Events shorter than block_size are dropped.
        - Remaining events are truncated to a multiple of block_size.

    Raises:
        FileNotFoundError: If no PNG images are found under folder.
        ValueError: If no images fall within the requested time window.

    Translated from extract_all_final.m (Mathieu Schaer 2017) and adapted for Python
    """
    folder = Path(folder)
    files = sorted(folder.rglob("*.png"), key=lambda p: p.name)
    if not files:
        raise FileNotFoundError(f"No .png images found under {folder}")

    names = [f.name for f in files]
    meta = [parse_meta(n) for n in names]
    dates = pd.to_datetime([m[0] for m in meta])
    ids = np.array([m[1] for m in meta])
    cams = np.array([m[2] for m in meta])

    keep = (dates >= pd.Timestamp(tstart)) & (dates <= pd.Timestamp(tstop))
    files = [f for f, k in zip(files, keep) if k]
    dates, ids, cams = dates[keep], ids[keep], cams[keep]
    if len(files) == 0:
        raise ValueError("No images in the requested time window")

    # Continuous events: split where the time gap exceeds s hours
    gaps = np.diff(dates.values).astype("timedelta64[s]").astype(float) / 3600.0
    starts = np.r_[0, np.where(gaps > s)[0] + 1]
    ends = np.r_[np.where(gaps > s)[0], len(dates) - 1]
    lengths = ends - starts + 1
    long = lengths >= block_size
    starts, lengths = starts[long], lengths[long]
    ends = starts + (lengths // block_size) * block_size - 1
    lengths = ends - starts + 1

    masks = load_masks(davos)
    masks = {0: masks[0], 1: masks[1], 2: masks[2]}

    rec_date, rec_id, rec_cam, rec_freq = [], [], [], []
    rec_por, rec_dmax, rec_frac, rec_dmaxall = [], [], [], []

    total = int(lengths.sum())
    print(f"Total number of MASC images selected: {total}")
    done = 0

    for st, en in zip(starts, ends):
        sl = slice(st, en + 1)
        ev_dates = dates[sl]
        ev_ids, ev_cams = ids[sl], cams[sl]
        ev_files = files[st:en + 1]
        ev_freq = photo_frequency(ev_dates, ev_cams, w)

        for b in range(0, len(ev_files), block_size):
            blk = slice(b, b + block_size)
            imgs = [cv2.imread(str(p), cv2.IMREAD_GRAYSCALE) for p in ev_files[blk]]
            descr = extract_features_block(imgs, ev_cams[blk], masks)
            for off, (por, dmax, frac, dall) in enumerate(descr):
                i = b + off
                rec_date.append(ev_dates[i])
                rec_id.append(ev_ids[i])
                rec_cam.append(int(ev_cams[i]))
                rec_freq.append(ev_freq[i])
                rec_por.append(por)
                rec_dmax.append(dmax)
                rec_frac.append(frac)
                rec_dmaxall.append(dall)
            done += len(imgs)
            print(f"{done}/{total}")

    df = pd.DataFrame(
        {
            "date_vec": rec_date,
            "ID": rec_id,
            "cam": rec_cam,
            "freq": rec_freq,
            "porosity": rec_por,
            "dmax": rec_dmax,
            "fractal": rec_frac,
            "dmax_all": rec_dmaxall,
        }
    )
    # Drop images with no detected object (porosity == 0)
    return df[df["porosity"] > 0].reset_index(drop=True)
