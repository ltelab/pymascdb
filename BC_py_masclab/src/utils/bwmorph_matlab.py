"""
MATLAB Image Processing Toolbox compatible bwmorph operations.

Ports the morphological sub-operations used by masclab/geometry/skeleton_props.m:
  close, open, thin (Inf), endpoints, branchpoints.

`thin` reproduces MATLAB bwmorph(BW,'thin',Inf) exactly: it is the Guo & Hall
(1989) two-subiteration parallel thinning, applied through two 256-entry lookup
tables (G123 / G123') derived from the G1/G2/G3/G3' conditions documented by
MathWorks. This is the same algorithm scikit-image uses in `morphology.thin`.

Endpoints use the 8-neighbour count (== 1). Branchpoints use the skeleton
crossing number (>= 3) rather than a raw neighbour count, which avoids the
massive over-counting of junctions seen with a naive ">= 3 neighbours" rule.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

_SE3 = np.ones((3, 3), dtype=bool)

# ---------------------------------------------------------------------------
# Thinning (Guo & Hall, via MATLAB G1/G2/G3 conditions)
# ---------------------------------------------------------------------------

# Neighbourhood weight mask (center excluded). Matches scikit-image / MATLAB
# bwlookup ordering used to derive the G123 LUTs below.
_THIN_MASK = np.array([[8, 4, 2], [16, 0, 1], [32, 64, 128]], dtype=np.uint8)


def _nabe(n: int) -> np.ndarray:
    """Bits 0..8 of neighbourhood code n (bit i == neighbour i)."""
    return np.array([(n >> i) & 1 for i in range(9)], dtype=bool)


def _G1(n: int) -> bool:
    bits = _nabe(n)
    s = 0
    for i in (0, 2, 4, 6):
        if (not bits[i]) and (bits[i + 1] or bits[(i + 2) % 8]):
            s += 1
    return s == 1


def _G2(n: int) -> bool:
    bits = _nabe(n)
    n1 = n2 = 0
    for k in (1, 3, 5, 7):
        if bits[k] or bits[k - 1]:
            n1 += 1
        if bits[k] or bits[(k + 1) % 8]:
            n2 += 1
    return min(n1, n2) in (2, 3)


def _G3(n: int) -> bool:
    bits = _nabe(n)
    return not ((bits[1] or bits[2] or (not bits[7])) and bits[0])


def _G3p(n: int) -> bool:
    bits = _nabe(n)
    return not ((bits[5] or bits[6] or (not bits[3])) and bits[4])


def _build_thin_luts() -> tuple[np.ndarray, np.ndarray]:
    g1 = np.array([_G1(n) for n in range(256)], dtype=bool)
    g2 = np.array([_G2(n) for n in range(256)], dtype=bool)
    g3 = np.array([_G3(n) for n in range(256)], dtype=bool)
    g3p = np.array([_G3p(n) for n in range(256)], dtype=bool)
    g12 = g1 & g2
    return g12 & g3, g12 & g3p


_G123_LUT, _G123P_LUT = _build_thin_luts()


# ---------------------------------------------------------------------------
# Endpoints / branchpoints (crossing-number based)
# ---------------------------------------------------------------------------

# Circular neighbour ordering (clockwise) for crossing-number computation.
#   (0,0)(0,1)(0,2)
#   (1,0)     (1,2)
#   (2,0)(2,1)(2,2)
_CIRC_MASK = np.array([[1, 2, 4], [128, 0, 8], [64, 32, 16]], dtype=np.uint8)


def _crossing_number(n: int) -> int:
    """0->1 transitions around the circular 8-neighbourhood."""
    bits = [(n >> i) & 1 for i in range(8)]
    return sum(1 for k in range(8) if bits[k] == 0 and bits[(k + 1) % 8] == 1)


def _popcount(n: int) -> int:
    return int(bin(n & 0xFF).count("1"))


_CROSSING_LUT = np.array([_crossing_number(n) for n in range(256)], dtype=np.uint8)
_NEIGHBOR_COUNT_LUT = np.array([_popcount(n) for n in range(256)], dtype=np.uint8)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def bwmorph_close(bw: np.ndarray) -> np.ndarray:
    """bwmorph(bw, 'close') — imclose with 3x3 ones."""
    return ndimage.binary_closing(np.asarray(bw, dtype=bool), structure=_SE3)


def bwmorph_open(bw: np.ndarray) -> np.ndarray:
    """bwmorph(bw, 'open') — imopen with 3x3 ones."""
    return ndimage.binary_opening(np.asarray(bw, dtype=bool), structure=_SE3)


def bwmorph_thin(bw: np.ndarray, n_iter: float = np.inf) -> np.ndarray:
    """
    bwmorph(bw, 'thin', Inf) — Guo & Hall two-subiteration parallel thinning.

    Iterates the two LUT sub-passes until the image stops changing (or until
    n_iter iterations have been performed).
    """
    skel = np.asarray(bw, dtype=np.uint8).copy()
    if skel.size == 0 or not skel.any():
        return skel.astype(bool)

    remaining = n_iter
    while remaining != 0:
        before = int(skel.sum())
        for lut in (_G123_LUT, _G123P_LUT):
            codes = ndimage.correlate(skel, _THIN_MASK, mode="constant", cval=0)
            skel[np.take(lut, codes)] = 0
        if int(skel.sum()) == before:
            break
        remaining -= 1
    return skel.astype(bool)


def _neighbour_codes(bw: np.ndarray) -> np.ndarray:
    return ndimage.correlate(
        np.asarray(bw, dtype=np.uint8), _CIRC_MASK, mode="constant", cval=0
    )


def bwmorph_endpoints(bw: np.ndarray) -> np.ndarray:
    """bwmorph(bw, 'endpoints') — foreground pixels with a single neighbour."""
    bw = np.asarray(bw, dtype=bool)
    codes = _neighbour_codes(bw)
    return bw & (np.take(_NEIGHBOR_COUNT_LUT, codes) == 1)


def bwmorph_branchpoints(bw: np.ndarray) -> np.ndarray:
    """
    bwmorph(bw, 'branchpoints') — foreground pixels where >= 3 skeleton branches
    meet, detected via the crossing number (0->1 transitions around the pixel).
    """
    bw = np.asarray(bw, dtype=bool)
    codes = _neighbour_codes(bw)
    return bw & (np.take(_CROSSING_LUT, codes) >= 3)
