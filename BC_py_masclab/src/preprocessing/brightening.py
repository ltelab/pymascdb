"""
Brightening and contrast enhancement for MASC snowflake images.

Functions for:
- MATLAB-faithful adapthisteq (CLAHE) with NumTiles / ClipLimit defaults
- Small-image zero-padding before CLAHE (brightening.m)

Translated from brightening.m (Christophe Praz 2015) and adapted for Python.
adapthisteq internals follow MATLAB Image Processing Toolbox (Zuiderveld CLAHE).
Last update : September 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

import numpy as np


def _clip_histogram(hist: np.ndarray, clip_limit: int, num_bins: int) -> np.ndarray:
    """Clip a tile histogram and redistribute excess counts (MATLAB clipHistogram)."""
    hist = hist.astype(np.int64, copy=True)
    clip_limit = int(clip_limit)
    total_excess = int(np.sum(np.maximum(hist - clip_limit, 0)))
    avg_bin_incr = total_excess // num_bins
    upper_limit = clip_limit - avg_bin_incr

    for k in range(num_bins):
        if hist[k] > clip_limit:
            hist[k] = clip_limit
        elif hist[k] > upper_limit:
            total_excess -= int(clip_limit - hist[k])
            hist[k] = clip_limit
        else:
            total_excess -= avg_bin_incr
            hist[k] = hist[k] + avg_bin_incr

    k = 0
    while total_excess != 0:
        if total_excess < 0:
            break
        step_size = max(num_bins // total_excess, 1)
        for m in range(k, num_bins, step_size):
            if hist[m] < clip_limit:
                hist[m] += 1
                total_excess -= 1
                if total_excess == 0:
                    break
        k += 1
        if k > num_bins - 1:
            k = 0
    return hist


def _make_mapping_uniform(hist: np.ndarray, num_pix_in_tile: int) -> np.ndarray:
    """
    Uniform-distribution tile mapping in [0, 255], rounded to uint8 levels.

    Matches MATLAB makeMapping(..., 'uniform') followed by grayxform's
    uint8 output (round(mapping01 * 255)).
    """
    hist_sum = np.cumsum(hist.astype(np.float64))
    mapping = np.minimum(hist_sum * (255.0 / float(num_pix_in_tile)), 255.0)
    return np.floor(mapping + 0.5)


def _adapthisteq(
    image: np.ndarray,
    num_tiles: tuple[int, int] = (8, 8),
    clip_limit: float = 0.01,
    num_bins: int = 256,
) -> np.ndarray:
    """
    MATLAB adapthisteq with default parameters.

    Defaults matched: NumTiles=[8 8], ClipLimit=0.01, NBins=256,
    Distribution='uniform', Range='full'.

    The image is symmetrically padded so each tile has even dimensions and
    the grid divides the padded size exactly; padding is removed at the end.
    """
    if image.ndim != 2:
        raise ValueError("adapthisteq expects a 2-D grayscale image")

    I0 = np.asarray(image, dtype=np.uint8)
    dim_I0 = np.array(I0.shape, dtype=int)
    num_tiles_arr = np.asarray(num_tiles, dtype=int)
    if num_tiles_arr.shape != (2,) or np.any(num_tiles_arr < 2):
        raise ValueError("num_tiles must be a pair of integers >= 2")

    I = I0
    dim_I = dim_I0.copy()
    dim_tile = dim_I / num_tiles_arr

    row_div = (dim_I[0] % num_tiles_arr[0]) == 0
    col_div = (dim_I[1] % num_tiles_arr[1]) == 0
    row_even = col_even = True
    if row_div and col_div:
        row_even = (dim_tile[0] % 2) == 0
        col_even = (dim_tile[1] % 2) == 0

    no_pad_rect = None
    if not (row_div and col_div and row_even and col_even):
        pad_row = 0
        pad_col = 0
        if not row_div:
            row_tile_dim = int(dim_I[0] // num_tiles_arr[0] + 1)
            pad_row = int(row_tile_dim * num_tiles_arr[0] - dim_I[0])
        else:
            row_tile_dim = int(dim_I[0] // num_tiles_arr[0])
        if not col_div:
            col_tile_dim = int(dim_I[1] // num_tiles_arr[1] + 1)
            pad_col = int(col_tile_dim * num_tiles_arr[1] - dim_I[1])
        else:
            col_tile_dim = int(dim_I[1] // num_tiles_arr[1])

        if row_tile_dim % 2 != 0:
            pad_row += int(num_tiles_arr[0])
        if col_tile_dim % 2 != 0:
            pad_col += int(num_tiles_arr[1])

        pad_row_pre = pad_row // 2
        pad_row_post = int(np.ceil(pad_row / 2.0))
        pad_col_pre = pad_col // 2
        pad_col_post = int(np.ceil(pad_col / 2.0))

        I = np.pad(
            I,
            ((pad_row_pre, pad_row_post), (pad_col_pre, pad_col_post)),
            mode="symmetric",
        )
        no_pad_rect = (
            pad_row_pre,
            pad_col_pre,
            pad_row_pre + dim_I0[0],
            pad_col_pre + dim_I0[1],
        )

    dim_I = np.array(I.shape, dtype=int)
    dim_tile = (dim_I // num_tiles_arr).astype(int)
    num_pix_in_tile = int(np.prod(dim_tile))
    min_clip_limit = int(np.ceil(num_pix_in_tile / float(num_bins)))
    clip_limit_abs = min_clip_limit + int(
        round(float(clip_limit) * (num_pix_in_tile - min_clip_limit))
    )

    n_tile_rows = int(num_tiles_arr[0])
    n_tile_cols = int(num_tiles_arr[1])
    mappings: list[list[np.ndarray]] = [
        [None] * n_tile_cols for _ in range(n_tile_rows)  # type: ignore[list-item]
    ]

    for col in range(n_tile_cols):
        for row in range(n_tile_rows):
            r0 = row * dim_tile[0]
            c0 = col * dim_tile[1]
            tile = I[r0 : r0 + dim_tile[0], c0 : c0 + dim_tile[1]]
            hist = np.bincount(tile.ravel().astype(np.intp), minlength=num_bins)
            hist = hist[:num_bins].astype(np.int64)
            hist = _clip_histogram(hist, clip_limit_abs, num_bins)
            mappings[row][col] = _make_mapping_uniform(hist, num_pix_in_tile)

    clahe = np.zeros(I.shape, dtype=np.float64)
    img_tile_row = 0
    for k in range(1, n_tile_rows + 2):
        if k == 1:
            img_tile_num_rows = int(dim_tile[0] // 2)
            map_tile_rows = (0, 0)
        elif k == n_tile_rows + 1:
            img_tile_num_rows = int(dim_tile[0] // 2)
            map_tile_rows = (n_tile_rows - 1, n_tile_rows - 1)
        else:
            img_tile_num_rows = int(dim_tile[0])
            map_tile_rows = (k - 2, k - 1)

        img_tile_col = 0
        for l in range(1, n_tile_cols + 2):
            if l == 1:
                img_tile_num_cols = int(dim_tile[1] // 2)
                map_tile_cols = (0, 0)
            elif l == n_tile_cols + 1:
                img_tile_num_cols = int(dim_tile[1] // 2)
                map_tile_cols = (n_tile_cols - 1, n_tile_cols - 1)
            else:
                img_tile_num_cols = int(dim_tile[1])
                map_tile_cols = (l - 2, l - 1)

            ul = mappings[map_tile_rows[0]][map_tile_cols[0]]
            ur = mappings[map_tile_rows[0]][map_tile_cols[1]]
            bl = mappings[map_tile_rows[1]][map_tile_cols[0]]
            br = mappings[map_tile_rows[1]][map_tile_cols[1]]

            rslice = slice(img_tile_row, img_tile_row + img_tile_num_rows)
            cslice = slice(img_tile_col, img_tile_col + img_tile_num_cols)
            region = I[rslice, cslice]

            # grayxform(uint8, map) → rounded uint8 levels, then bilinear blend
            ulv = ul[region]
            urv = ur[region]
            blv = bl[region]
            brv = br[region]

            row_w = np.arange(img_tile_num_rows, dtype=np.float64)[:, None]
            col_w = np.arange(img_tile_num_cols, dtype=np.float64)[None, :]
            row_rev = np.arange(img_tile_num_rows, 0, -1, dtype=np.float64)[:, None]
            col_rev = np.arange(img_tile_num_cols, 0, -1, dtype=np.float64)[None, :]
            norm = float(img_tile_num_rows * img_tile_num_cols)

            clahe[rslice, cslice] = (
                row_rev * (col_rev * ulv + col_w * urv)
                + row_w * (col_rev * blv + col_w * brv)
            ) / norm

            img_tile_col += img_tile_num_cols
        img_tile_row += img_tile_num_rows

    # MATLAB uint8 assignment: round half away from zero for positive values
    out = np.clip(np.floor(clahe + 0.5), 0, 255).astype(np.uint8)
    if no_pad_rect is not None:
        r0, c0, r1, c1 = no_pad_rect
        out = out[r0:r1, c0:c1]
    return out


def brightening(image_in: np.ndarray) -> np.ndarray:
    """
    Brighten a MASC snowflake image using CLAHE histogram equalization.

    Small images are zero-padded to a minimum of 8 px in each deficient
    dimension (as in brightening.m), then passed to adapthisteq and cropped.

    Args:
        image_in: Input grayscale image (uint8).

    Returns:
        Brightened image with the same shape and dtype as the input.
    """
    tile_dim = 8
    image_in = np.asarray(image_in)
    n_lines, n_cols = image_in.shape

    if n_lines < tile_dim and n_cols < tile_dim:
        image_mod = np.zeros((tile_dim, tile_dim), dtype=np.uint8)
        image_mod[:n_lines, :n_cols] = image_in
        image_out = _adapthisteq(image_mod)
        return image_out[:n_lines, :n_cols]

    if n_lines < tile_dim:
        image_mod = np.zeros((tile_dim, n_cols), dtype=np.uint8)
        image_mod[:n_lines, :n_cols] = image_in
        image_out = _adapthisteq(image_mod)
        return image_out[:n_lines, :n_cols]

    if n_cols < tile_dim:
        image_mod = np.zeros((n_lines, tile_dim), dtype=np.uint8)
        image_mod[:n_lines, :n_cols] = image_in
        image_out = _adapthisteq(image_mod)
        return image_out[:n_lines, :n_cols]

    return _adapthisteq(image_in)
