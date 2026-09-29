"""
Geometric and topological descriptor helpers for MASC ROI analysis.

This module provides MATLAB-compatible functions for convex hull, Dmax/D90,
ellipse and circle fitting, skeleton, fractal dimension, Haralick texture,
symmetry, and blur metrics.

Functions for:
- Convex hull, Dmax, D90, and rectangularity
- Inscribed/circumscribed circle and ellipse fitting
- Skeleton, fractal, Haralick, symmetry, and blur descriptors

Translated from:
- compute_convex_hull.m
- compute_Dmax.m
- compute_D90.m
- fit_circle_around.m
- fit_ellipse_inside.m
- fit_ellipse_around.m
- compute_rectangularity.m
- skeleton_props.m
- boxcount.m
- fractal_dim.m
- haralick_props.m
- compute_symmetry_features.m
- compute_blur_index.m
and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import cv2

from typing import Dict, Tuple
from scipy import ndimage
from scipy.spatial import ConvexHull
from skimage.feature import graycomatrix, graycoprops
from skimage.measure import regionprops

from src.utils.bwmorph_matlab import (
    bwmorph_branchpoints,
    bwmorph_close,
    bwmorph_endpoints,
    bwmorph_open,
    bwmorph_thin,
)

_SKEL_SE = np.ones((3, 3), dtype=bool)
# 4-connectivity structuring element (matches MATLAB bwperim default conn=4)
_PERIM_SE4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)


def _points_on_hull_edge(points: np.ndarray, pa: np.ndarray, pb: np.ndarray) -> np.ndarray:
    """Indices of *points* collinear with segment pa→pb (MATLAB convhull boundary)."""
    ax, ay = float(pa[0]), float(pa[1])
    bx, by = float(pb[0]), float(pb[1])
    dx, dy = bx - ax, by - ay
    seg_len2 = dx * dx + dy * dy
    if seg_len2 == 0:
        match = np.where((points[:, 0] == ax) & (points[:, 1] == ay))[0]
        return match[:1]

    cross = (points[:, 0] - ax) * dy - (points[:, 1] - ay) * dx
    dot = (points[:, 0] - ax) * dx + (points[:, 1] - ay) * dy
    mask = (np.abs(cross) < 1e-9) & (dot >= -1e-9) & (dot <= seg_len2 + 1e-9)
    idx = np.where(mask)[0]
    return idx[np.argsort(dot[idx] / seg_len2)]


def _matlab_conv_hull_indices(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return MATLAB convhull-ordered boundary pixel indices (not just corner vertices)."""
    points = np.column_stack([x.ravel(), y.ravel()]).astype(np.float64)
    hull = ConvexHull(points)
    verts = hull.vertices
    ordered: list[int] = []
    n_verts = len(verts)
    for i in range(n_verts):
        ia = verts[i]
        ib = verts[(i + 1) % n_verts]
        edge_idx = _points_on_hull_edge(points, points[ia], points[ib])
        if ordered and edge_idx.size and edge_idx[0] == ordered[-1]:
            edge_idx = edge_idx[1:]
        ordered.extend(int(j) for j in edge_idx)
    # Close the polygon once (MATLAB convhull repeats the first vertex)
    if ordered and ordered[-1] != ordered[0]:
        ordered.append(ordered[0])
    return np.asarray(ordered, dtype=int)


def boxcount(c: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform 2-D box counting for fractal dimension estimation.

    Args:
        c: Input array (binary mask or grayscale).

    Returns:
        Tuple of (box_counts, box_sizes) arrays.

    Raises:
        ValueError: If the input array shape is invalid for 2-D box counting.

    Translated from boxcount.m and adapted for Python
    """
    arr = np.asarray(c).squeeze()
    if arr.ndim == 3 and arr.shape[2] == 3 and arr.shape[0] >= 8 and arr.shape[1] >= 8:
        arr = arr.sum(axis=2)
    bw = arr.astype(bool)
    if bw.ndim != 2:
        raise ValueError("boxcount supports 1D/2D/3D arrays; got invalid shape")

    width = max(bw.shape)
    p = int(np.ceil(np.log2(width)))
    width = 2 ** p
    if bw.shape[0] != width or bw.shape[1] != width:
        padded = np.zeros((width, width), dtype=bool)
        padded[: bw.shape[0], : bw.shape[1]] = bw
        bw = padded

    n = np.zeros(p + 1, dtype=np.float64)
    n[p] = float(bw.sum())
    work = bw.copy()
    for g in range(p - 1, -1, -1):
        siz = 2 ** (p - g)
        siz2 = round(siz / 2)
        for i in range(0, width - siz + 1, siz):
            for j in range(0, width - siz + 1, siz):
                work[i, j] = (
                    work[i, j]
                    or work[i + siz2, j]
                    or work[i, j + siz2]
                    or work[i + siz2, j + siz2]
                )
        n[g] = float(work[0:width:siz, 0:width:siz].sum())

    n = n[::-1]
    r = 2.0 ** np.arange(0, p + 1, dtype=np.float64)
    return n, r


def compute_convex_hull(x: np.ndarray, y: np.ndarray, perim: float = 0.0) -> Dict[str, np.ndarray]:
    """
    Compute the convex hull of snowflake boundary pixels.

    Args:
        x: x coordinates of boundary pixels.
        y: y coordinates of boundary pixels.
        perim: Snowflake perimeter used to compute convexity.

    Returns:
        Dictionary with keys xh, yh, solidity, perim, convexity.

    Translated from compute_convex_hull.m and adapted for Python
    """
    x = x.ravel()
    y = y.ravel()

    if len(x) < 3:
        return {'xh': x, 'yh': y, 'solidity': 0.0, 'perim': 0.0, 'convexity': 0.0}

    try:
        k = _matlab_conv_hull_indices(x, y)
        xh = np.asarray(x[k], dtype=np.float64).ravel()
        yh = np.asarray(y[k], dtype=np.float64).ravel()
    except Exception:
        return {'xh': x, 'yh': y, 'solidity': 0.0, 'perim': 0.0, 'convexity': 0.0}

    # Start at leftmost-topmost vertex (min x, then min y), single close — matches
    # the usual MATLAB convhull / bwboundaries origin convention used by the GT.
    if len(xh) > 1:
        if xh[0] == xh[-1] and yh[0] == yh[-1]:
            xh, yh = xh[:-1], yh[:-1]
        start = int(np.lexsort((yh, xh))[0])
        xh = np.concatenate([xh[start:], xh[:start], [xh[start]]])
        yh = np.concatenate([yh[start:], yh[:start], [yh[start]]])

    # solidity = n_pixels / polyarea(hull)
    hull_area = float(np.abs(np.dot(xh[:-1], yh[1:]) - np.dot(xh[1:], yh[:-1])) / 2.0)
    solidity = len(x) / hull_area if hull_area > 0 else 0.0

    # hull perimeter
    hull_perim = 0.0
    for i in range(len(xh) - 1):
        hull_perim += np.sqrt((xh[i+1] - xh[i])**2 + (yh[i+1] - yh[i])**2)

    return {
        'xh': xh,
        'yh': yh,
        'solidity': float(solidity),
        'perim': float(hull_perim),
        'convexity': float(hull_perim / perim) if perim > 0 else 0.0
    }


def compute_Dmax(xh: np.ndarray, yh: np.ndarray) -> Tuple[float, float, Tuple[float, float], Tuple[float, float]]:
    """
    Compute maximum diameter (Dmax) from convex hull points.
    
    Args:
        xh: x coordinates of convex hull vertices
        yh: y coordinates of convex hull vertices
        
    Returns:
        Tuple of (Dmax, Dmax_theta, DmaxA, DmaxB):
            - Dmax: Maximum distance between any two hull points
            - Dmax_theta: Angle of Dmax line in degrees
            - DmaxA: (x, y) coordinates of first endpoint
            - DmaxB: (x, y) coordinates of the second endpoint.

    Translated from compute_Dmax.m and adapted for Python
    """
    if len(xh) < 2:
        return 0.0, 0.0, (0.0, 0.0), (0.0, 0.0)
    
    max_dist = 0.0
    idx_a = 0
    idx_b = 0
    
    # Brute force: find maximum distance between all pairs (compute_Dmax.m)
    n_pts = len(xh)
    for i in range(n_pts - 1):
        for j in range(i + 1, n_pts):
            dist = np.sqrt((xh[i] - xh[j])**2 + (yh[i] - yh[j])**2)
            if dist > max_dist:
                max_dist = dist
                idx_a = i
                idx_b = j

    # Angle vs horizontal — acos + sign rule (compute_Dmax.m), not atan2
    ax, ay = float(xh[idx_a]), float(yh[idx_a])
    bx, by = float(xh[idx_b]), float(yh[idx_b])
    vec = np.array([bx - ax, by - ay], dtype=np.float64)
    norm_vec = float(np.linalg.norm(vec))
    if norm_vec == 0.0:
        theta_deg = 0.0
    else:
        theta = float(np.arccos(np.clip(vec[0] / norm_vec, -1.0, 1.0)))
        if theta > np.pi / 2.0:
            theta = np.pi - theta
        if (ax < bx and ay < by) or (ax > bx and ay > by):
            theta = -theta
        theta_deg = theta * 180.0 / np.pi

    DmaxA = (ax, ay)
    DmaxB = (bx, by)

    return float(max_dist), float(theta_deg), DmaxA, DmaxB


def _imrotate_loose_nearest(bw: np.ndarray, angle_deg: float) -> np.ndarray:
    """
    Emulate MATLAB imrotate(bw, angle, 'nearest', 'loose') on a binary image.

    MATLAB maps pixel centers through a rotation about the image center,
    enlarges the canvas so nothing is clipped ('loose'), and resamples with
    nearest neighbour. scipy.ndimage.rotate uses a different grid convention
    (rounding/center offsets), which shifts D0/D90 counts by 1-2 pixels.
    """
    M, N = bw.shape
    a = np.deg2rad(angle_deg)
    c, s = np.cos(a), np.sin(a)

    # Output canvas size: rotate the input corner box (MATLAB spatial coords)
    xs = np.array([0.5, N + 0.5, 0.5, N + 0.5])
    ys = np.array([0.5, 0.5, M + 0.5, M + 0.5])
    xt = xs * c - ys * s
    yt = xs * s + ys * c
    xmin, xmax = xt.min(), xt.max()
    ymin, ymax = yt.min(), yt.max()
    ncols = int(np.ceil(xmax - xmin))
    nrows = int(np.ceil(ymax - ymin))

    # Center the output box on the rotated input box
    x0 = xmin - (ncols - (xmax - xmin)) / 2.0
    y0 = ymin - (nrows - (ymax - ymin)) / 2.0

    cc, rr = np.meshgrid(np.arange(1, ncols + 1), np.arange(1, nrows + 1))
    xw = x0 + (cc - 0.5)
    yw = y0 + (rr - 0.5)

    # Inverse rotation back to input coordinates, nearest-neighbour lookup
    xi = xw * c + yw * s
    yi = -xw * s + yw * c
    col = np.floor(xi + 0.5).astype(int)
    row = np.floor(yi + 0.5).astype(int)

    out = np.zeros((nrows, ncols), dtype=bool)
    ok = (col >= 1) & (col <= N) & (row >= 1) & (row <= M)
    out[ok] = bw[row[ok] - 1, col[ok] - 1]
    return out


def compute_D90(mask: np.ndarray, Dmax: float, Dmax_theta: float) -> Dict[str, float]:
    """
    Compute D90 dimensions from the Dmax orientation.

    Rotates the filled mask by -Dmax_theta and counts occupied columns and rows.

    Args:
        mask: Filled binary mask of the snowflake.
        Dmax: Maximum diameter in pixels.
        Dmax_theta: Dmax angle in degrees.

    Returns:
        Dictionary with keys Dmax_0, Dmax_90, AR, and angle.

    Translated from compute_D90.m and adapted for Python
    """
    empty = {'Dmax_0': 0.0, 'Dmax_90': 0.0, 'AR': 0.0, 'angle': float(Dmax_theta)}
    if Dmax == 0:
        return empty

    bw = (mask > 0).astype(np.uint8)
    if not np.any(bw):
        return empty

    # MATLAB: imrotate(bw, -Dmax_angle) with defaults 'nearest' + 'loose'.
    # scipy.ndimage.rotate uses a different grid/rounding convention and is
    # off by 1-2 px on D0/D90, so we emulate imrotate exactly (see
    # _imrotate_loose_nearest). The sign convention of the emulation (y-down
    # pixel grid) is opposite to MATLAB's, hence +Dmax_theta here.
    tmp_im = _imrotate_loose_nearest(bw > 0, float(Dmax_theta))

    col_sum = np.sum(tmp_im, axis=0)
    dmax_0 = int(np.count_nonzero(col_sum > 0))

    row_sum = np.sum(tmp_im, axis=1)
    dmax_90 = int(np.count_nonzero(row_sum > 0))

    ar = float(dmax_90 / dmax_0) if dmax_0 > 0 else 0.0
    return {
        'Dmax_0': float(dmax_0),
        'Dmax_90': float(dmax_90),
        'AR': ar,
        'angle': float(Dmax_theta),
    }


def fit_circle_around(x_perim: np.ndarray, y_perim: np.ndarray, xh: np.ndarray, yh: np.ndarray) -> Dict[str, float]:
    """
    Fit the smallest circumscribed circle around the convex hull.

    Args:
        x_perim: Perimeter x coordinates.
        y_perim: Perimeter y coordinates.
        xh: Convex hull x coordinates.
        yh: Convex hull y coordinates.

    Returns:
        Dictionary with keys X0, Y0, r, and A (circle area).

    Translated from fit_circle_around.m and adapted for Python
    """
    xh = xh.ravel()
    yh = yh.ravel()

    if len(xh) < 2:
        return {'X0': 0.0, 'Y0': 0.0, 'r': 0.0, 'A': 0.0, 'status': 'ok'}

    stepsize = 0.1
    tol = 1.0

    x0_new = float(np.mean(xh))
    y0_new = float(np.mean(yh))
    r_old = np.inf
    x0_old, y0_old = x0_new, y0_new

    for _ in range(10000):  # safety cap
        xd = xh - x0_new
        yd = yh - y0_new
        r_vec = np.sqrt(xd**2 + yd**2)
        r_new = float(np.max(r_vec))

        if r_new < r_old:
            idx_support = np.where(r_vec >= r_new - tol)[0]
            mx = float(np.mean(xd[idx_support]))
            my = float(np.mean(yd[idx_support]))
            n = np.sqrt(mx**2 + my**2)
            if n == 0:
                break
            normal = stepsize * np.array([mx, my]) / n
            x0_old, y0_old = x0_new, y0_new
            x0_new += normal[0]
            y0_new += normal[1]
            r_old = r_new
        else:
            break

    r = r_old
    return {
        'X0': x0_old,
        'Y0': y0_old,
        'r': r,
        'A': float(np.pi * r**2),
        'status': 'ok',
    }


def fit_ellipse_inside(x: np.ndarray, y: np.ndarray, x_perim: np.ndarray, y_perim: np.ndarray, theta: float) -> Dict[str, float]:
    """
    Fit the largest inscribed ellipse inside the snowflake perimeter.

    Args:
        x: Snowflake pixel x coordinates.
        y: Snowflake pixel y coordinates.
        x_perim: Perimeter x coordinates.
        y_perim: Perimeter y coordinates.
        theta: Ellipse orientation angle in radians.

    Returns:
        Dictionary with keys a, b, theta, X0, Y0, and A (ellipse area).

    Notes:
        - Searches focal distances to maximize area while keeping perimeter outside.

    Translated from fit_ellipse_inside.m and adapted for Python
    """
    x = x.ravel().astype(float)
    y = y.ravel().astype(float)
    xP = x_perim.ravel().astype(float)
    yP = y_perim.ravel().astype(float)

    centroid = np.array([np.mean(x), np.mean(y)])
    xPc = xP - centroid[0]
    yPc = yP - centroid[1]

    theta_eff = -theta  # MATLAB applies theta = -theta

    # Untilt the perimeter points
    cos_t = np.cos(-theta_eff)
    sin_t = np.sin(-theta_eff)
    xPc0 = cos_t * xPc - sin_t * yPc
    yPc0 = sin_t * xPc + cos_t * yPc

    maxd = float(np.max(xPc0) - np.min(xPc0))
    if maxd <= 0:
        return {
            'a': 0.0, 'b': 0.0, 'theta': -theta_eff,
            'X0': centroid[0], 'Y0': centroid[1], 'A': 0.0, 'status': 'ok',
        }

    # MATLAB: d = 0:1:maxd  → 0, 1, ..., floor(maxd)  (not arange(0, maxd+1))
    d_arr = np.arange(0.0, np.floor(maxd) + 1.0)
    s_arr = np.zeros(len(d_arr))
    A_arr = np.zeros(len(d_arr))

    for i, d in enumerate(d_arr):
        # sum of distances from each perimeter point to the two foci
        s_vals = (np.sqrt((xPc0 - 0.5*d)**2 + yPc0**2) +
                  np.sqrt((xPc0 + 0.5*d)**2 + yPc0**2))
        s = float(np.min(s_vals))  # smallest sum = tightest constraint
        s_arr[i] = s
        s2d2 = s**2 - d**2
        A_arr[i] = np.pi * s / 4.0 * np.sqrt(s2d2) if s2d2 > 0 else 0.0

    idx = int(np.argmax(A_arr))
    s = s_arr[idx]
    d = d_arr[idx]
    a = s / 2.0
    b = float(np.sqrt(max((s/2.0)**2 - (d/2.0)**2, 0.0)))

    return {
        'a': a,
        'b': b,
        'theta': -theta_eff,
        'X0': float(centroid[0]),
        'Y0': float(centroid[1]),
        'A': float(np.pi * a * b),
        'status': 'ok',
    }


def fit_ellipse_around(x: np.ndarray, y: np.ndarray, xh: np.ndarray, yh: np.ndarray, theta: float) -> Dict[str, float]:
    """
    Fit the smallest circumscribed ellipse around the convex hull.

    Args:
        x: Snowflake pixel x coordinates.
        y: Snowflake pixel y coordinates.
        xh: Convex hull x coordinates.
        yh: Convex hull y coordinates.
        theta: Ellipse orientation angle in radians.

    Returns:
        Dictionary with keys a, b, theta, X0, Y0, and A (ellipse area).

    Notes:
        - Searches focal distances to minimize area while keeping hull points inside.

    Translated from fit_ellipse_around.m and adapted for Python
    """
    x = x.ravel().astype(float)
    y = y.ravel().astype(float)
    xh = xh.ravel().astype(float)
    yh = yh.ravel().astype(float)

    centroid = np.array([np.mean(x), np.mean(y)])
    xc = x - centroid[0]
    yc = y - centroid[1]

    theta_eff = -theta  # MATLAB applies theta = -theta

    # Untilt data and hull
    cos_t = np.cos(-theta_eff)
    sin_t = np.sin(-theta_eff)
    xc0 = cos_t * xc - sin_t * yc

    xhc = xh - centroid[0]
    yhc = yh - centroid[1]
    xhc0 = cos_t * xhc - sin_t * yhc
    yhc0 = sin_t * xhc + cos_t * yhc

    maxd = float(np.max(xc0) - np.min(xc0))
    if maxd <= 0:
        return {
            'a': 0.0, 'b': 0.0, 'theta': -theta_eff,
            'X0': centroid[0], 'Y0': centroid[1], 'A': 0.0, 'status': 'ok',
        }

    # MATLAB: d = 0:1:maxd  → 0, 1, ..., floor(maxd)  (not arange(0, maxd+1))
    d_arr = np.arange(0.0, np.floor(maxd) + 1.0)
    s_arr = np.zeros(len(d_arr))
    A_arr = np.zeros(len(d_arr))

    for i, d in enumerate(d_arr):
        # sum of distances from each hull point to the two foci
        s_vals = (np.sqrt((xhc0 - 0.5*d)**2 + yhc0**2) +
                  np.sqrt((xhc0 + 0.5*d)**2 + yhc0**2))
        s = float(np.max(s_vals))  # largest sum = tightest constraint
        s_arr[i] = s
        s2d2 = s**2 - d**2
        A_arr[i] = np.pi * s / 4.0 * np.sqrt(s2d2) if s2d2 > 0 else np.inf

    idx = int(np.argmin(A_arr))
    s = s_arr[idx]
    d = d_arr[idx]
    a = s / 2.0
    b = float(np.sqrt(max((s/2.0)**2 - (d/2.0)**2, 0.0)))

    return {
        'a': a,
        'b': b,
        'theta': -theta_eff,
        'X0': float(centroid[0]),
        'Y0': float(centroid[1]),
        'A': float(np.pi * a * b),
        'status': 'ok',
    }


def compute_rectangularity(x: np.ndarray, y: np.ndarray, perim: float) -> Dict[str, np.ndarray]:
    """
    Compute minimum-area bounding rectangle metrics.

    Args:
        x: Snowflake pixel x coordinates.
        y: Snowflake pixel y coordinates.
        perim: Snowflake perimeter in pixels.

    Returns:
        Dictionary with rectangle corners, width, length, aspect ratio, and theta.

    Translated from compute_rectangularity.m and adapted for Python
    """
    _zero = {
        'rectx': np.array([]), 'recty': np.array([]),
        'width': 0.0, 'length': 0.0, 'height': 0.0,
        'A_ratio': 0.0, 'rectangularity': 0.0,
        'p_ratio': 0.0, 'aspect_ratio': 0.0, 'eccentricity': 0.0, 'theta': 0.0,
    }
    pts = np.hstack([x, y]).astype(np.float32)
    if len(pts) < 3 or perim <= 0:
        return _zero

    rect = cv2.minAreaRect(pts)
    (rw, rh) = rect[1]
    rect_a = rw * rh
    if rect_a <= 0:
        return _zero

    box = cv2.boxPoints(rect).astype(np.float64)
    # MATLAB min-area rect polygon starts at the lowest-y corner (tie: lowest x),
    # then follows the same winding; OpenCV boxPoints starts elsewhere.
    # Closing vertex repeats the start, so sorted coord compares need this too.
    i0 = int(np.lexsort((box[:, 0], box[:, 1]))[0])
    box = np.roll(box, -i0, axis=0)

    d01 = np.linalg.norm(box[1] - box[0])
    d12 = np.linalg.norm(box[2] - box[1])
    width, length = (d01, d12) if d01 <= d12 else (d12, d01)

    ar = width / length
    a_ratio = len(pts) / rect_a
    # longer edge angle vs horizontal (masclab compute_rectangularity.m)
    i_long = 0 if d01 >= d12 else 1
    vec = box[i_long + 1] - box[i_long]
    theta = float(np.arccos(np.clip(abs(vec[0]) / np.linalg.norm(vec), 0.0, 1.0)))

    # Close the rectangle polygon (first corner repeated), like MATLAB
    rectx = np.append(box[:, 0], box[0, 0])
    recty = np.append(box[:, 1], box[0, 1])

    return {
        'rectx': rectx,
        'recty': recty,
        'width': float(width),
        'length': float(length),
        'height': float(length),
        'A_ratio': float(a_ratio),
        'rectangularity': float(a_ratio),
        'p_ratio': float(2.0 * (rw + rh) / perim),
        'aspect_ratio': float(ar),
        'eccentricity': float(np.sqrt(max(1.0 - ar, 0.0))),
        'theta': theta,
    }


def _bwperim_count(bw: np.ndarray) -> int:
    """Count perimeter pixels (MATLAB bwperim, default 4-connectivity).

    A foreground pixel is on the perimeter if at least one of its 4-connected
    neighbours is background — i.e. it is removed by a 4-connectivity erosion.
    """
    eroded = ndimage.binary_erosion(bw, structure=_PERIM_SE4, border_value=0)
    return int(np.sum(bw & ~eroded))


def skeleton_props(data: np.ndarray) -> Dict[str, float]:
    """
    Compute pseudo-skeleton topology metrics for a ROI image.

    Args:
        data: Grayscale ROI image (uint8).

    Returns:
        Dictionary with p_ratio, A_ratio, N_ends, N_junctions, length, and density.

    Translated from skeleton_props.m and adapted for Python
    """
    _zero = {
        'p_ratio': 0.0, 'A_ratio': 0.0, 'N_ends': 0, 'N_junctions': 0,
        'endpoints': 0, 'branches': 0, 'length': 0.0, 'density': 0.0,
    }
    if data.size == 0:
        return _zero

    bw = data > 0
    area = int(bw.sum())
    if area == 0:
        return _zero

    # skeleton_props.m: close → open → thin(Inf) → largest component
    bw_close_open = bwmorph_open(bwmorph_close(bw))
    bw_thin = bwmorph_thin(bw_close_open)

    # MATLAB regionprops labels with 8-connectivity, then keeps the largest blob
    labeled, n_comp = ndimage.label(bw_thin, structure=_SKEL_SE)
    if n_comp == 0:
        return _zero

    sizes = np.bincount(labeled.ravel())
    sizes[0] = 0
    sk = labeled == int(np.argmax(sizes))

    n_ends = int(bwmorph_endpoints(sk).sum())
    n_junc = int(bwmorph_branchpoints(sk).sum())
    perim_skel = int(sk.sum())
    perim = _bwperim_count(bw)

    return {
        'p_ratio': perim_skel / perim if perim else 0.0,
        'A_ratio': perim_skel / area,
        'N_ends': n_ends,
        'N_junctions': n_junc,
        'endpoints': n_ends,
        'branches': n_junc,
        'length': float(perim_skel),
        'density': perim_skel / area,
    }


def fractal_dim(data: np.ndarray, nbox_max: float = np.inf) -> float:
    """
    Estimate fractal dimension using box counting.

    Args:
        data: Input ROI image or binary mask.
        nbox_max: Maximum number of box sizes to use (np.inf for all).

    Returns:
        Fractal dimension estimate as a float.

    Translated from fractal_dim.m and adapted for Python
    """
    if data.size == 0:
        return 0.0

    try:
        nbox, rbox = boxcount(data)
    except Exception:
        return 0.0

    if len(nbox) > nbox_max:
        nbox = nbox[: int(nbox_max)]
        rbox = rbox[: int(nbox_max)]
    elif len(nbox) > 5:
        nbox = nbox[:-2]
        rbox = rbox[:-2]
    elif len(nbox) == 5:
        nbox = nbox[:-1]
        rbox = rbox[:-1]

    if len(nbox) < 2 or np.any(nbox <= 0) or np.any(rbox <= 0):
        return 0.0

    log_nbox = np.log(nbox)
    log_rbox = np.log(rbox)
    design = np.column_stack([np.ones(len(nbox)), log_rbox])
    beta, _, _, _ = np.linalg.lstsq(design, log_nbox, rcond=None)
    return float(-beta[1])


def haralick_props(data: np.ndarray) -> Dict[str, float]:
    """
    Compute Haralick GLCM texture features (mean and std over 4 directions).

    Args:
        data: Grayscale ROI image (uint8).

    Returns:
        Dictionary with Contrast, Correlation, Energy, Homogeneity and their
        ``*_std`` counterparts (sample std over the 4 angles).

    Notes:
        Faithful to ``texture/haralick_props.m`` (C. Praz):
        - ``graycomatrix(..., 'NumLevels', 256, 'Symmetric', true)`` with default
          GrayLimits ``[0 255]`` for uint8 (no min–max rescaling).
        - Offsets: ``[0 1], [-1 1], [-1 0], [-1 -1]`` (0°, 45°, 90°, 135°).
        - Background slice removed: ``glcm(2:end, 2:end)`` (``pix_value_thresh=2``).
        - Energy = ASM; Homogeneity = ``sum p/(1+|i-j|)`` (MATLAB graycoprops).
        - ``std`` is MATLAB sample std (N-1).

    Translated from haralick_props.m and adapted for Python
    """
    _keys = (
        'Contrast', 'Correlation', 'Energy', 'Homogeneity',
        'Contrast_std', 'Correlation_std', 'Energy_std', 'Homogeneity_std',
    )
    _zero = {k: 0.0 for k in _keys}

    if data.size == 0 or data.shape[0] < 2 or data.shape[1] < 2:
        return _zero

    data_u8 = np.clip(data, 0, 255).astype(np.uint8)
    if not np.any(data_u8 > 0):
        return _zero

    try:
        # MATLAB default GrayLimits for uint8 is [0 255] → SI = I (0-based bins).
        # Do NOT min–max rescale: that was the source of the ~16× Contrast blow-up.
        glcm = graycomatrix(
            data_u8,
            distances=[1],
            angles=[0.0, np.pi / 4.0, np.pi / 2.0, 3.0 * np.pi / 4.0],
            levels=256,
            symmetric=True,
            normed=False,
        )
        # pix_value_thresh = 2 → drop background co-occurrences (pixel value 0)
        glcm = glcm[1:, 1:, :, :]

        contrast = graycoprops(glcm, 'contrast').ravel()
        correlation = graycoprops(glcm, 'correlation').ravel()
        # MATLAB graycoprops 'Energy' is the Angular Second Moment (sum p^2)
        energy = graycoprops(glcm, 'ASM').ravel()

        # MATLAB Homogeneity: sum p_ij / (1+|i-j|)  (skimage uses (i-j)^2)
        P = glcm.astype(np.float64)
        sums = P.sum(axis=(0, 1), keepdims=True)
        sums[sums == 0] = 1.0
        Pn = P / sums
        n_bins = Pn.shape[0]
        ii = np.arange(n_bins)[:, None]
        jj = np.arange(n_bins)[None, :]
        w = 1.0 / (1.0 + np.abs(ii - jj))
        homogeneity = (Pn * w[:, :, None, None]).sum(axis=(0, 1)).ravel()

        def _mean_std(values: np.ndarray) -> tuple[float, float]:
            return float(np.mean(values)), float(np.std(values, ddof=1))

        c_m, c_s = _mean_std(contrast)
        r_m, r_s = _mean_std(correlation)
        e_m, e_s = _mean_std(energy)
        h_m, h_s = _mean_std(homogeneity)

        return {
            'Contrast': c_m,
            'Correlation': r_m,
            'Energy': e_m,
            'Homogeneity': h_m,
            'Contrast_std': c_s,
            'Correlation_std': r_s,
            'Energy_std': e_s,
            'Homogeneity_std': h_s,
        }
    except Exception:
        return _zero


def _sym_zero() -> Dict[str, float]:
    out = {f'P{i}': 0.0 for i in range(11)}
    out['mean'] = 0.0
    out['std'] = 0.0
    return out


def _ray_polygon_distances(xp: np.ndarray, yp: np.ndarray, dmax: float,
                           fallback: float) -> np.ndarray:
    """Max centroid-to-boundary distance per degree (masclab polyxpoly loop)."""
    n = len(xp)
    x1, y1 = xp, yp
    x2 = np.r_[xp[1:], xp[0]]
    y2 = np.r_[yp[1:], yp[0]]
    out = np.full(360, fallback, dtype=np.float64)
    rad = np.deg2rad(np.arange(360, dtype=np.float64))
    cos_a, sin_a = np.cos(rad), np.sin(rad)
    for i in range(360):
        rx, ry = dmax * cos_a[i], dmax * sin_a[i]
        dxe, dye = x2 - x1, y2 - y1
        denom = rx * dye - ry * dxe
        valid = np.abs(denom) > 1e-12
        t = np.full(n, np.nan)
        tu = np.full(n, np.nan)
        t[valid] = (x1[valid] * dye[valid] - y1[valid] * dxe[valid]) / denom[valid]
        tu[valid] = (x1[valid] * ry - y1[valid] * rx) / denom[valid]
        hit = valid & (t >= 0.0) & (t <= 1.0) & (tu >= 0.0) & (tu <= 1.0)
        if np.any(hit):
            out[i] = np.hypot(t[hit] * rx, t[hit] * ry).max()
    return out


def compute_symmetry_features(mask_filled: np.ndarray, Dmax: float,
                              eq_radius: float) -> Dict[str, float]:
    """
    Compute FFT-based radial symmetry features.

    Args:
        mask_filled: Filled binary mask of the snowflake.
        Dmax: Maximum diameter in pixels.
        eq_radius: Equivalent circular radius (fallback ray length).

    Returns:
        Dictionary with mean, std, and P0–P10 FFT power components.

    Translated from compute_symmetry_features.m and adapted for Python
    """
    if mask_filled.size == 0 or Dmax <= 0:
        return _sym_zero()

    mask = mask_filled > 0
    if not np.any(mask):
        return _sym_zero()

    props = regionprops(mask.astype(np.uint8))
    if not props:
        return _sym_zero()
    cy, cx = props[0].centroid  # row, col — MATLAB Centroid [x, y]

    mask_u8 = mask.astype(np.uint8)
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return _sym_zero()
    c = max(contours, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
    if len(c) < 3:
        return _sym_zero()

    xp = c[:, 0] - cx
    yp = c[:, 1] - cy
    distances = _ray_polygon_distances(xp, yp, float(Dmax), float(eq_radius))
    mean_d = float(distances.mean())
    # MATLAB std() uses N-1 normalization
    std_d = float(distances.std(ddof=1)) if len(distances) > 1 else 0.0
    norm = (distances - mean_d) / std_d if std_d > 0 else distances - mean_d

    f = np.fft.fft(norm)
    p1 = np.abs(f) / len(norm)
    p1 = p1[: len(norm) // 2 + 1].copy()
    if len(p1) > 2:
        p1[1:-1] *= 2.0

    out = {'mean': mean_d, 'std': std_d}
    for i in range(11):
        out[f'P{i}'] = float(p1[i]) if i < len(p1) else 0.0
    return out


def compute_blur_index(data: np.ndarray) -> Dict[str, float]:
    """
    Compute blur index from Gaussian high-pass differences.

    Args:
        data: ROI grayscale image (uint8).

    Returns:
        Dictionary with xhi2 (sigma=2) and xhi4 (sigma=4) blur metrics.

    Notes:
        Faithful to ``process_new_descriptors.m``::
            blurry = imgaussfilt(roi.data, sigma);  % keeps uint8 class
            diff = roi.data - blurry;               % uint8 saturating subtract
            xhi = std2(diff);                       % sample std (N-1)

        ``imgaussfilt`` uses FilterSize = 2*ceil(2*sigma)+1 and replicate padding.

    Translated from compute_blur_index.m / process_new_descriptors.m
    """
    if data.size == 0 or data.shape[0] < 5 or data.shape[1] < 5:
        return {'xhi2': 0.0, 'xhi4': 0.0}

    data_u8 = np.clip(data, 0, 255).astype(np.uint8)

    def _imgaussfilt_u8(img: np.ndarray, sigma: float) -> np.ndarray:
        # MATLAB: FilterSize = 2*ceil(2*sigma)+1, Padding = 'replicate'
        ksize = int(2 * np.ceil(2.0 * sigma) + 1)
        trunc = ((ksize - 1) / 2.0) / float(sigma)
        out = ndimage.gaussian_filter(
            img.astype(np.float64), sigma=sigma, truncate=trunc, mode='nearest'
        )
        return np.clip(np.round(out), 0, 255).astype(np.uint8)

    blurry2 = _imgaussfilt_u8(data_u8, 2.0)
    blurry4 = _imgaussfilt_u8(data_u8, 4.0)

    # MATLAB uint8 − uint8 saturates at 0 (no negative high-pass lobes)
    diff2 = (data_u8.astype(np.int16) - blurry2.astype(np.int16)).clip(min=0).astype(np.float64)
    diff4 = (data_u8.astype(np.int16) - blurry4.astype(np.int16)).clip(min=0).astype(np.float64)

    return {
        'xhi2': float(np.std(diff2, ddof=1)),
        'xhi4': float(np.std(diff4, ddof=1)),
    }
