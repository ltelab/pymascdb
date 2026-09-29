"""
Daily quicklook generation for classified MASC campaigns.

Structure mirrors MATLAB masclab (kept in one module for the Python runner):
- load / triplet merge  ≈  prediction/merge_predictions_for_campaign.m
- time-bin aggregation  ≈  first half of analysis/microstructure/make_masc_time_series.m
- figure export         ≈  Fig4 (#1) and Fig5 (#2) of make_masc_time_series.m
- fill bands            ≈  tools/fill_btw_curves.m

Translated from Christophe Praz (2015) and adapted for Python.
Last update : September 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from src.dataio.load_roi_data import load_roi_data

logger = logging.getLogger(__name__)

def _compute_riming_idx(ri_degree: float) -> float:
    """Same as masclab/prediction/compute_riming_idx.m."""
    return 0.5 * (np.sin(np.pi / 4.0 * (ri_degree - 3.0)) + 1.0)


def _triplet_stem(name: str) -> str:
    """Group cam0/cam1/cam2 files for one flake (basename without _cam_*)."""
    stem = Path(name).stem
    return re.sub(r"_cam_?\d+$", "", stem, flags=re.IGNORECASE)


def _parse_roi_datetime(roi: Dict[str, Any]) -> Optional[datetime]:
    """Resolve observation time from ``roi`` (tnum or filename)."""
    tnum = roi.get("tnum")
    if isinstance(tnum, datetime):
        dt = tnum
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.astimezone(timezone.utc)
    if isinstance(tnum, (int, float)) and not np.isnan(tnum):
        try:
            return datetime.fromordinal(int(tnum) - 366).replace(tzinfo=timezone.utc)
        except Exception:
            pass
    name = roi.get("name") or ""
    if not name:
        return None
    parts = str(name).split("_")
    if len(parts) >= 2:
        date_str, time_str = parts[0], parts[1]
        dp = date_str.split(".")
        tp = time_str.split(".")
        if len(dp) >= 3 and len(tp) >= 3:
            try:
                return datetime(
                    int(dp[0]),
                    int(dp[1]),
                    int(dp[2]),
                    int(tp[0]),
                    int(tp[1]),
                    int(tp[2]),
                    tzinfo=timezone.utc,
                )
            except ValueError:
                return None
    return None


def _as_float(x: Any, default: float = np.nan) -> float:
    try:
        if x is None:
            return default
        return float(np.asarray(x).ravel()[0])
    except Exception:
        return default


def _as_prob_vec(p: Any, n: int) -> np.ndarray:
    if p is None:
        return np.full(n, np.nan)
    a = np.asarray(p, dtype=float).ravel()
    if a.size >= n:
        return a[:n].copy()
    out = np.full(n, np.nan)
    out[: a.size] = a
    return out


def _masc_point_from_roi(roi: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Build one numeric record per ROI (per-image row)."""
    dt = _parse_roi_datetime(roi)
    if dt is None:
        return None

    e = roi.get("E") or {}
    a_e = _as_float(e.get("a"), np.nan)
    b_e = _as_float(e.get("b"), np.nan)
    ar_e = (b_e / a_e) if a_e and not np.isnan(a_e) and a_e != 0 else np.nan
    theta_e = _as_float(e.get("theta"), np.nan)

    rprobs = _as_prob_vec(roi.get("riming_probs"), 5)
    if np.all(np.isnan(rprobs)):
        ri_degree = np.nan
        ri_idx = np.nan
    else:
        ri_degree = float(np.dot(np.arange(1, 6), rprobs / (np.nansum(rprobs) + 1e-12)))
        ri_idx = _compute_riming_idx(ri_degree)

    mprob = _as_float(roi.get("melting_probs"), np.nan)

    return {
        "dt": dt,
        "name": str(roi.get("name") or ""),
        "label_id": _as_float(roi.get("label_ID"), np.nan),
        "area": _as_float(roi.get("area"), np.nan),
        "dmax": _as_float(roi.get("Dmax"), np.nan),
        "dmean": _as_float(roi.get("Dmean"), np.nan),
        "complex": _as_float(roi.get("complex"), np.nan),
        "xhi": _as_float(roi.get("xhi"), np.nan),
        "fallspeed": _as_float(roi.get("fallspeed"), np.nan),
        "ar": ar_e,
        "orientation": theta_e,
        "n_roi": _as_float(roi.get("n_roi"), np.nan),
        "roundness": _as_float(roi.get("roundness"), np.nan),
        "riming_index": ri_idx,
        "melting_id": _as_float(roi.get("melting_ID"), np.nan),
        "melting_prob": mprob,
        "label_probs": _as_prob_vec(roi.get("label_probs"), 6),
    }


def _merge_triplet_group(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Merge same-flake multi-cam rows (probabilities like merge_predictions_for_campaign.m)."""
    if len(rows) == 1:
        return rows[0].copy()

    base = rows[0].copy()
    lp = np.nansum(np.stack([r["label_probs"] for r in rows], axis=0), axis=0)
    s = np.nansum(lp)
    if s > 0:
        lp = lp / s
    lid = int(np.nanargmax(lp) + 1)
    base["label_probs"] = lp
    base["label_id"] = float(lid)

    keys_mean = (
        "area",
        "dmean",
        "fallspeed",
        "ar",
        "roundness",
        "riming_index",
        "melting_prob",
        "melting_id",
        "n_roi",
    )
    for k in keys_mean:
        base[k] = float(np.nanmean([r[k] for r in rows]))

    base["dmax"] = float(np.nanmax([r["dmax"] for r in rows]))
    base["complex"] = float(np.nanmax([r["complex"] for r in rows]))
    base["ar"] = float(np.nanmin([r["ar"] for r in rows]))
    base["orientation"] = float(np.nanmin(np.abs([r["orientation"] for r in rows])))
    base["xhi"] = float(np.nanmean([r["xhi"] for r in rows]))
    if base["melting_prob"] > 0.5:
        base["melting_id"] = 1.0
    else:
        base["melting_id"] = 0.0
    return base


def _collect_roi_paths(
    outdir: Path,
    save_format: str,
    qualities: Sequence[Literal["GOOD", "BAD"]],
) -> List[Path]:
    ext = {"joblib": ".joblib", "pkl": ".pkl", "mat": ".mat"}.get(save_format, ".joblib")
    paths: List[Path] = []
    for q in qualities:
        for p in outdir.rglob(f"DATA/{q}/*{ext}"):
            if p.is_file() and p.stem and p.stem[0].isdigit():
                paths.append(p)
    return sorted(set(paths))


def discover_time_bounds_from_outdir(
    outdir: Path, save_format: str, qualities: Sequence[Literal["GOOD", "BAD"]] = ("GOOD",)
) -> Optional[Tuple[datetime, datetime]]:
    """
    Infer min/max timestamps from ROI files under outdir.

    Args:
        outdir: Processed campaign output directory.
        save_format: ROI file format extension (joblib, pkl, mat).
        qualities: Quality subfolders to scan ('GOOD', 'BAD').

    Returns:
        Tuple of (t_min, t_max) UTC datetimes, or None if no valid files found.

    Notes:
        - Used when config starthr_vec or endhr_vec is set to 'all'.
    """
    tmin: Optional[datetime] = None
    tmax: Optional[datetime] = None
    for path in _collect_roi_paths(outdir, save_format, qualities):
        try:
            roi = load_roi_data(path, save_format)
            dt = _parse_roi_datetime(roi)
            if dt is None:
                continue
            tmin = dt if tmin is None or dt < tmin else tmin
            tmax = dt if tmax is None or dt > tmax else tmax
        except Exception as ex:
            logger.debug("skip %s: %s", path, ex)
    if tmin is None or tmax is None:
        return None
    return tmin, tmax


@dataclass
class QuicklooksOptions:
    """Configuration options for quicklook time-binning and figure export."""
    pixres_mm: float = 33.5 / 1000.0
    xi_thresh: float = 9.0
    Nmin_interval: int = 30
    Nmin_shift: int = 10
    Nclasses_masc: int = 6
    MASC_classes: Tuple[str, ...] = ("SP", "CC", "PC", "AG", "GR", "CPC")
    MASC_classes_desired: Tuple[int, ...] = (1, 2, 3, 4, 5, 6)
    N_MascSamples_min: int = 0
    use_triplet: bool = False
    OR_180: bool = True
    savefigs: bool = True
    overview_classif: bool = True
    overview_microstruct: bool = True
    disp_now: bool = True  # MATLAB disp.now
    qualities: Tuple[str, ...] = ("GOOD",)


def _hsv_cmasc() -> np.ndarray:
    """MATLAB ``cmasc``: ``hsv(6)`` then AG←magenta, CPC←green tweak."""
    # MATLAB hsv(n): hues = (0:n-1)/n, S=1, V=1
    hues = np.arange(6, dtype=float) / 6.0
    c = np.array([matplotlib.colors.hsv_to_rgb((h, 1.0, 1.0)) for h in hues])
    cmasc = c.copy()
    cmasc[3] = c[5]
    cmasc[5] = np.array([49.0, 163.0, 84.0]) / 255.0
    return cmasc


def _jet15_rgb() -> Tuple[float, float, float]:
    """MATLAB ``jet(20)(15,:)`` used for Aspect/Area Ratio bands."""
    return tuple(plt.cm.jet(np.linspace(0.0, 1.0, 20))[14, :3])  # type: ignore[return-value]


def _nanstd(a: np.ndarray) -> float:
    """Sample std (ddof=1) like MATLAB ``nanstd``."""
    a = np.asarray(a, dtype=float)
    if np.count_nonzero(~np.isnan(a)) < 2:
        return float(np.nan)
    return float(np.nanstd(a, ddof=1))


def _aggregate_bins(
    rows: List[Dict[str, Any]],
    t0: datetime,
    t1: datetime,
    opt: QuicklooksOptions,
) -> Dict[str, Any]:
    """Sliding bins like MATLAB ``tstart:minutes(shift):tstop`` with closed ``[ta, tb]``."""
    if t0.tzinfo is None:
        t0 = t0.replace(tzinfo=timezone.utc)
    if t1.tzinfo is None:
        t1 = t1.replace(tzinfo=timezone.utc)
    t0 = t0.astimezone(timezone.utc)
    t1 = t1.astimezone(timezone.utc)

    shift = timedelta(minutes=opt.Nmin_shift)
    width = timedelta(minutes=opt.Nmin_interval)
    # MATLAB includes tstop when it lands on the grid
    tgrid: List[datetime] = []
    t = t0
    while t <= t1:
        tgrid.append(t)
        t += shift
    if not tgrid:
        return {}

    nbin = len(tgrid)
    tgrid2 = [tg + width for tg in tgrid]

    classes_des = np.array(opt.MASC_classes_desired, dtype=int)
    nclass = opt.Nclasses_masc

    Nmasc_all = np.zeros(nbin)
    Nmasc_filtered = np.zeros(nbin)
    dom_class = np.full(nbin, np.nan)
    sclass = np.full((nbin, nclass), np.nan)
    melting = np.full(nbin, np.nan)
    riming = np.full(nbin, np.nan)

    dmax_m = np.full(nbin, np.nan)
    dmax_med = np.full(nbin, np.nan)
    dmax_s = np.full(nbin, np.nan)
    nroi_m = np.full(nbin, np.nan)
    nroi_med = np.full(nbin, np.nan)
    nroi_s = np.full(nbin, np.nan)
    ar_m = np.full(nbin, np.nan)
    ar_med = np.full(nbin, np.nan)
    ar_s = np.full(nbin, np.nan)
    arr_m = np.full(nbin, np.nan)
    arr_med = np.full(nbin, np.nan)
    arr_s = np.full(nbin, np.nan)
    or_m = np.full(nbin, np.nan)
    or_med = np.full(nbin, np.nan)
    or_s = np.full(nbin, np.nan)
    cplx_m = np.full(nbin, np.nan)
    cplx_med = np.full(nbin, np.nan)
    cplx_s = np.full(nbin, np.nan)
    fs_m = np.full(nbin, np.nan)
    fs_med = np.full(nbin, np.nan)
    fs_s = np.full(nbin, np.nan)

    xs = np.array([r["dmax"] for r in rows], dtype=float)
    dts = [r["dt"] for r in rows]

    for i, (ta, tb) in enumerate(zip(tgrid, tgrid2)):
        # MATLAB: Xt >= tgrid(i) & Xt <= tgrid2(i)
        idx_all = [j for j, dt in enumerate(dts) if ta <= dt <= tb]
        if not idx_all:
            continue
        idx_f = [
            j
            for j in idx_all
            if not np.isnan(rows[j]["xhi"]) and rows[j]["xhi"] >= opt.xi_thresh
        ]
        idx_f_nosp = [j for j in idx_f if rows[j]["label_id"] != 1.0]

        Nmasc_all[i] = len(idx_all)
        Nmasc_filtered[i] = len(idx_f)

        xa = xs[idx_all]

        if idx_f:
            labs = [int(rows[j]["label_id"]) for j in idx_f if not np.isnan(rows[j]["label_id"])]
            if labs:
                dom_class[i] = float(max(set(labs), key=labs.count))
            counts = np.zeros(nclass)
            for k in range(1, nclass + 1):
                counts[k - 1] = sum(1 for j in idx_f if rows[j]["label_id"] == float(k))
            sclass[i, :] = counts
            melting[i] = float(np.nanmean([rows[j]["melting_id"] for j in idx_f]))
            if idx_f_nosp:
                riming[i] = float(np.nanmean([rows[j]["riming_index"] for j in idx_f_nosp]))
            else:
                riming[i] = np.nan

        if idx_all:
            dmax_m[i] = float(np.nanmean(xa))
            dmax_med[i] = float(np.nanmedian(xa))
            dmax_s[i] = _nanstd(xa)
            nroi_m[i] = float(np.nanmean([rows[j]["n_roi"] for j in idx_all]))
            nroi_med[i] = float(np.nanmedian([rows[j]["n_roi"] for j in idx_all]))
            nroi_s[i] = _nanstd(np.array([rows[j]["n_roi"] for j in idx_all], dtype=float))
            ar_m[i] = float(np.nanmean([rows[j]["ar"] for j in idx_all]))
            ar_med[i] = float(np.nanmedian([rows[j]["ar"] for j in idx_all]))
            ar_s[i] = _nanstd(np.array([rows[j]["ar"] for j in idx_all], dtype=float))
            arr_m[i] = float(np.nanmean([rows[j]["roundness"] for j in idx_all]))
            arr_med[i] = float(np.nanmedian([rows[j]["roundness"] for j in idx_all]))
            arr_s[i] = _nanstd(np.array([rows[j]["roundness"] for j in idx_all], dtype=float))
            orv = np.array([rows[j]["orientation"] for j in idx_all], dtype=float)
            if not opt.OR_180:
                orv = np.abs(orv)
            or_m[i] = float(np.nanmean(orv))
            or_med[i] = float(np.nanmedian(orv))
            or_s[i] = _nanstd(orv)
            cplx_m[i] = float(np.nanmean([rows[j]["complex"] for j in idx_all]))
            cplx_med[i] = float(np.nanmedian([rows[j]["complex"] for j in idx_all]))
            cplx_s[i] = _nanstd(np.array([rows[j]["complex"] for j in idx_all], dtype=float))
            fs_a = np.array([rows[j]["fallspeed"] for j in idx_all], dtype=float)
            fs_a[(fs_a >= 15.0) | (fs_a < 0.1)] = np.nan
            if np.all(np.isnan(fs_a)):
                fs_m[i] = np.nan
                fs_med[i] = np.nan
                fs_s[i] = np.nan
            else:
                fs_m[i] = float(np.nanmean(fs_a))
                fs_med[i] = float(np.nanmedian(fs_a))
                fs_s[i] = _nanstd(fs_a)

    sclass_norm = np.zeros((nbin, len(classes_des)))
    for i in range(nbin):
        sub = sclass[i, classes_des - 1]
        ssum = np.nansum(sub)
        if ssum > 0:
            sclass_norm[i, :] = sub / ssum
        else:
            sclass_norm[i, :] = np.nan

    mask_low = Nmasc_filtered < opt.N_MascSamples_min
    dom_class[mask_low] = np.nan
    sclass[mask_low, :] = np.nan
    melting[mask_low] = np.nan
    riming[mask_low] = np.nan

    # Strip tz for numpy datetime64 (UTC already enforced above)
    tgrid_naive = [tg.replace(tzinfo=None) for tg in tgrid]
    return {
        "tgrid": np.array(tgrid_naive, dtype="datetime64[ns]"),
        "Nmasc_all": Nmasc_all,
        "Nmasc_filtered": Nmasc_filtered,
        "dom_class": dom_class,
        "sclass": sclass,
        "sclass_norm": sclass_norm,
        "melting": melting,
        "riming": riming,
        "Dmax": (dmax_m, dmax_med, dmax_s),
        "Nroi": (nroi_m, nroi_med, nroi_s),
        "AR": (ar_m, ar_med, ar_s),
        "ArR": (arr_m, arr_med, arr_s),
        "OR": (or_m, or_med, or_s),
        "cplx": (cplx_m, cplx_med, cplx_s),
        "fs": (fs_m, fs_med, fs_s),
    }


# MATLAB fullscreen saveas ≈ 2997×1736 px @ 150 dpi
_QL_FIGSIZE = (2997 / 150.0, 1736 / 150.0)
_QL_DPI = 150


def _fill_btw_curves(ax, x, y_lo, y_hi, y_mid, color, alpha: float = 0.35) -> None:
    """Strict port of masclab ``tools/fill_btw_curves.m``.

    Segments on finite *y1* (lower bound). For each contiguous run:
    - length ≥ 2 → ``fill`` band + mean line
    - length == 1 → vertical stem + ``ko`` (MarkerFaceColor=C)
    """
    x = np.asarray(x, dtype=float)
    y1 = np.asarray(y_lo, dtype=float).copy()
    y2 = np.asarray(y_hi, dtype=float).copy()
    y3 = np.asarray(y_mid, dtype=float).copy()

    # MATLAB nanstd of n=1 is 0; keep lo/hi finite so single samples plot
    finite_mid = np.isfinite(y3)
    y1 = np.where(finite_mid & ~np.isfinite(y1), y3, y1)
    y2 = np.where(finite_mid & ~np.isfinite(y2), y3, y2)

    valid = np.isfinite(y1)
    n = len(x)
    i = 0
    while i < n:
        if not valid[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and valid[j + 1]:
            j += 1
        xd = x[i : j + 1]
        y1d = y1[i : j + 1]
        y2d = y2[i : j + 1]
        y3d = y3[i : j + 1]
        if j != i:
            # fill([xd; flipud(xd)], [y1d; flipud(y2d)], C, 'edgecolor','none','facealpha',alpha)
            ax.fill(
                np.concatenate([xd, xd[::-1]]),
                np.concatenate([y1d, y2d[::-1]]),
                color=color,
                edgecolor="none",
                alpha=alpha,
                zorder=1,
            )
            # plot(xd, y3d, 'Color', C, 'linewidth', 2)
            ax.plot(xd, y3d, color=color, linewidth=2, zorder=2)
        else:
            # plot([xd xd],[y1d y2d],'Color',C);
            # plot(xd, y3d, 'ko', 'MarkerFaceColor', C);
            ax.plot(
                [float(xd[0]), float(xd[0])],
                [float(y1d[0]), float(y2d[0])],
                color=color,
                linewidth=1.5,
                solid_capstyle="round",
                zorder=2,
            )
            ax.plot(
                float(xd[0]),
                float(y3d[0]),
                linestyle="none",
                marker="o",
                markersize=7,
                markerfacecolor=color,
                markeredgecolor="k",
                markeredgewidth=1.0,
                zorder=3,
            )
        i = j + 1


def _panel_mean_std(mean: np.ndarray, std: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build lo/mid/hi for fill_btw_curves; NaN std → 0 (MATLAB ``nanstd`` of n=1)."""
    mid = np.asarray(mean, dtype=float)
    s = np.asarray(std, dtype=float)
    s_plot = np.where(np.isfinite(mid) & ~np.isfinite(s), 0.0, s)
    lo = mid - s_plot
    hi = mid + s_plot
    return lo, hi, mid


def _style_axes(
    ax,
    fontsize: float,
    *,
    tick_length: float = 6.0,
    grid_y: bool = False,
) -> None:
    """``box on`` + ticks out; grid x always, grid y optional (MATLAB ``grid on``)."""
    ax.set_facecolor("white")
    ax.grid(False)
    ax.xaxis.grid(True, which="major", color="0.75", linewidth=0.7, zorder=0)
    if grid_y:
        ax.yaxis.grid(True, which="major", color="0.75", linewidth=0.7, zorder=0)
    else:
        ax.yaxis.grid(False, which="both")
    ax.set_axisbelow(True)
    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        labelsize=fontsize,
        length=tick_length,
        width=1.1,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def _format_day_axis(ax, t0: datetime, t1: datetime, fontsize: float) -> None:
    """Day window with 00:00…00:00 ticks and date subtitle (MATLAB datetime look)."""
    if t0.tzinfo is not None:
        t0 = t0.replace(tzinfo=None)
    if t1.tzinfo is not None:
        t1 = t1.replace(tzinfo=None)
    ticks: List[datetime] = []
    t = t0
    while t <= t1 + timedelta(seconds=1):
        ticks.append(t)
        t = t + timedelta(hours=6)
        if len(ticks) > 48:
            break
    uniq: List[datetime] = []
    for tt in ticks:
        if t0 <= tt <= t1 and (not uniq or uniq[-1] != tt):
            uniq.append(tt)
    if not uniq or uniq[0] != t0:
        uniq = [t0] + [u for u in uniq if u != t0]
    if uniq[-1] != t1:
        uniq.append(t1)
    ax.set_xlim(mdates.date2num(t0), mdates.date2num(t1))
    ax.set_xticks([mdates.date2num(tt) for tt in uniq])
    ax.set_xticklabels([tt.strftime("%H:%M") for tt in uniq], fontsize=fontsize)
    ax.set_xlabel(t0.strftime("%b %d, %Y"), fontsize=fontsize - 1, labelpad=6)


def _maybe_now_line(ax, opt: QuicklooksOptions) -> None:
    if not opt.disp_now:
        return
    tnow = datetime.now(timezone.utc)
    y0, y1 = ax.get_ylim()
    ax.plot(
        [mdates.date2num(tnow), mdates.date2num(tnow)],
        [y0, y1],
        "r--",
        lw=2,
        label="Now",
    )


def _set_count_ylim_with_top_tick(ax, *series: np.ndarray) -> None:
    """Y-lim from 0 to a nice ceiling, with ticks that include the upper limit (MATLAB-like)."""
    from matplotlib.ticker import MaxNLocator

    vals = []
    for s in series:
        a = np.asarray(s, dtype=float)
        if a.size:
            vals.append(a)
    if not vals:
        ax.set_ylim(0, 1)
        return
    ymax = float(np.nanmax(np.concatenate(vals)))
    if not np.isfinite(ymax) or ymax <= 0:
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        return
    locator = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])
    ticks = locator.tick_values(0.0, ymax)
    ticks = ticks[ticks >= 0]
    if ticks.size == 0 or ticks[-1] < ymax:
        # ensure ceiling strictly above data so the top tick is visible
        pad = ymax * 0.05 if ymax > 0 else 1.0
        ticks = locator.tick_values(0.0, ymax + pad)
        ticks = ticks[ticks >= 0]
    top = float(ticks[-1]) if ticks.size else ymax
    if top < ymax:
        top = ymax
    ax.set_ylim(0.0, top)
    ax.set_yticks(ticks[ticks <= top + 1e-12])


def _stacked_area_with_gaps(
    ax,
    tnum: np.ndarray,
    sclass_norm: np.ndarray,
    colors: Sequence,
    labels: Sequence[str],
    legend_fontsize: float = 14,
) -> None:
    """MATLAB ``area(tgrid, sclass_norm)`` with gaps where rows are all-NaN."""
    from matplotlib.patches import Patch

    y = np.asarray(sclass_norm, dtype=float)
    valid = np.any(np.isfinite(y), axis=1)
    n = len(tnum)
    # Default dx for a single-bin segment (half a typical shift)
    if n >= 2:
        dx = float(np.median(np.diff(tnum)))
    else:
        dx = 10.0 / (24.0 * 60.0)
    i = 0
    while i < n:
        if not valid[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and valid[j + 1]:
            j += 1
        seg_t = np.asarray(tnum[i : j + 1], dtype=float)
        seg_y = np.nan_to_num(y[i : j + 1, :], nan=0.0)
        if seg_t.size == 1:
            # Degenerate area: tiny horizontal span so the column is visible
            seg_t = np.array([seg_t[0], seg_t[0] + 0.5 * dx])
            seg_y = np.vstack([seg_y, seg_y])
        ax.stackplot(
            seg_t,
            seg_y.T,
            colors=colors,
            alpha=1.0,
            edgecolor="none",
        )
        i = j + 1
    ax.legend(
        handles=[Patch(facecolor=c, edgecolor="none", label=lab) for c, lab in zip(colors, labels)],
        loc="upper right",
        fontsize=legend_fontsize,
        framealpha=1.0,
    )


def _save_overview_classif(
    data: Dict[str, Any],
    t0: datetime,
    t1: datetime,
    opt: QuicklooksOptions,
    savedir: Path,
    date_tag: str,
) -> None:
    # Fig4 in make_masc_time_series.m (fs_quicklook1 = 18)
    fs = 18
    tick_len = 12.0
    cmasc = _hsv_cmasc()
    tnum = mdates.date2num(data["tgrid"])
    cmelting = (241 / 255, 105 / 255, 19 / 255)

    date_title = f"{t0.strftime('%Y.%m.%d %H:%M')} - {t1.strftime('%Y.%m.%d %H:%M')}"
    fig = plt.figure(figsize=_QL_FIGSIZE, dpi=_QL_DPI)
    fig.patch.set_facecolor("white")

    ax1 = fig.add_subplot(4, 1, 1)
    ax1.plot(tnum, data["Nmasc_all"], "k-", lw=2, label="total")
    series_for_ylim = [data["Nmasc_all"]]
    if opt.xi_thresh > 0:
        ax1.plot(
            tnum,
            data["Nmasc_filtered"],
            "k--",
            lw=2,
            label=f"filtered xi>{opt.xi_thresh:.1f}",
        )
        series_for_ylim.append(data["Nmasc_filtered"])
    ax1.set_ylabel(
        "# MASC triplets" if opt.use_triplet else "# MASC images",
        fontsize=fs,
    )
    ax1.set_title(date_title, fontsize=fs)
    _set_count_ylim_with_top_tick(ax1, *series_for_ylim)
    _style_axes(ax1, fs, tick_length=tick_len, grid_y=False)
    _format_day_axis(ax1, t0, t1, fs)
    ax1.legend(loc="upper right", fontsize=fs - 4, framealpha=1.0)

    ax2 = fig.add_subplot(4, 1, (2, 3))
    # Keep NaNs for empty bins (MATLAB area leaves gaps; do not nan_to_num globally)
    y_stack = np.asarray(data["sclass_norm"], dtype=float)
    labs = [opt.MASC_classes[k - 1] for k in opt.MASC_classes_desired]
    colors = [cmasc[k - 1] for k in opt.MASC_classes_desired]
    _stacked_area_with_gaps(ax2, tnum, y_stack, colors, labs, legend_fontsize=fs - 4)
    ax2.set_ylabel("Proportions", fontsize=fs)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([0.0, 0.5, 1.0])
    _style_axes(ax2, fs, tick_length=tick_len, grid_y=False)
    _format_day_axis(ax2, t0, t1, fs)

    ax3 = fig.add_subplot(4, 1, 4)
    ax3.plot(tnum, data["riming"], color="b", lw=2, label="R_i")
    ax3.plot(tnum, data["melting"], color=cmelting, lw=2, label="% wet")
    ax3.set_ylabel("[-]", fontsize=fs)
    ax3.set_ylim(0, 1)
    ax3.set_yticks([0.0, 0.5, 1.0])
    _style_axes(ax3, fs, tick_length=tick_len, grid_y=False)
    ax3.legend(loc="upper right", fontsize=fs - 4, framealpha=1.0)
    _format_day_axis(ax3, t0, t1, fs)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.94, bottom=0.07, hspace=0.45)
    out = savedir / f"{date_tag}_quicklook#1.png"
    fig.savefig(out, dpi=_QL_DPI, facecolor="white")
    plt.close(fig)


def _save_overview_microstruct(
    data: Dict[str, Any],
    t0: datetime,
    t1: datetime,
    opt: QuicklooksOptions,
    savedir: Path,
    date_tag: str,
) -> None:
    """Fig5 in make_masc_time_series.m — 3×2 microstructural panels."""
    fs = 16  # fs_quicklook2
    tnum = mdates.date2num(data["tgrid"])
    fig, axes = plt.subplots(3, 2, figsize=_QL_FIGSIZE, dpi=_QL_DPI)
    fig.patch.set_facecolor("white")
    jet15 = _jet15_rgb()

    dmax_m, _, dmax_s = data["Dmax"]
    fs_m, _, fs_s = data["fs"]
    ar_m, _, ar_s = data["AR"]
    or_m, _, or_s = data["OR"]
    arr_m, _, arr_s = data["ArR"]
    cplx_m, _, cplx_s = data["cplx"]

    # MATLAB multiplies by pixres at aggregation; we store px and convert here
    pix = opt.pixres_mm
    panels = [
        (axes[0, 0], dmax_m * pix, dmax_s * pix, "Dmax", "[mm]", "r"),
        (axes[0, 1], fs_m, fs_s, "Fallspeed", "[m/s]", "b"),
        (axes[1, 0], ar_m, ar_s, "Aspect Ratio", "[-]", jet15),
        (
            axes[1, 1],
            or_m,
            or_s,
            "Orientation [-90; 90]" if opt.OR_180 else "Orientation [0; 90]",
            "[°]",
            "r",
        ),
        (axes[2, 0], arr_m, arr_s, "Area Ratio", "[-]", jet15),
        (axes[2, 1], cplx_m, cplx_s, "Complexity", "[-]", "b"),
    ]
    for ax, mean, std, title, ylab, col in panels:
        # Same as MATLAB: plot(tgrid, nan); fill_btw_curves(...); grid on; box on
        ax.plot(tnum, np.full_like(tnum, np.nan))
        lo, hi, mid = _panel_mean_std(mean, std)
        _fill_btw_curves(ax, tnum, lo, hi, mid, col, 0.35)
        ax.set_title(title, fontsize=fs)
        ax.set_ylabel(ylab, fontsize=fs)
        _style_axes(ax, fs, tick_length=8.0, grid_y=True)  # MATLAB grid on → H+V
        _format_day_axis(ax, t0, t1, fs)
        _maybe_now_line(ax, opt)
        # MATLAB: legend only on Dmax when disp.now
        if opt.disp_now and ax is axes[0, 0]:
            ax.legend(loc="upper right", fontsize=fs - 4, framealpha=1.0)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.07, hspace=0.35, wspace=0.20)
    out = savedir / f"{date_tag}_quicklook#2.png"
    fig.savefig(out, dpi=_QL_DPI, facecolor="white")
    plt.close(fig)


def load_masc_points_for_interval(
    outdir: Path,
    save_format: str,
    t_start: datetime,
    t_stop: datetime,
    use_triplet: bool,
    qualities: Sequence[str] = ("GOOD",),
) -> List[Dict[str, Any]]:
    """
    Load ROI rows with timestamps in the half-open interval [t_start, t_stop).

    Args:
        outdir: Processed campaign output directory.
        save_format: ROI file format.
        t_start: Interval start (inclusive, UTC).
        t_stop: Interval end (exclusive, UTC).
        use_triplet: Merge views sharing the same filename stem when True.
        qualities: Quality subfolders to include.

    Returns:
        List of flattened ROI point dicts ready for aggregation.
    """
    if t_start.tzinfo is None:
        t_start = t_start.replace(tzinfo=timezone.utc)
    if t_stop.tzinfo is None:
        t_stop = t_stop.replace(tzinfo=timezone.utc)
    t_start = t_start.astimezone(timezone.utc)
    t_stop = t_stop.astimezone(timezone.utc)

    qtuple = tuple(q for q in qualities if q in ("GOOD", "BAD"))
    paths = _collect_roi_paths(outdir, save_format, qtuple)  # type: ignore[arg-type]
    raw: List[Dict[str, Any]] = []
    for path in paths:
        try:
            roi = load_roi_data(path, save_format)
            row = _masc_point_from_roi(roi)
            if row is None:
                continue
            if not (t_start <= row["dt"] < t_stop):
                continue
            raw.append(row)
        except Exception as ex:
            logger.warning("Could not load %s: %s", path, ex)

    if not use_triplet:
        return raw

    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in raw:
        stem = _triplet_stem(str(row["name"]))
        groups.setdefault(stem, []).append(row)
    merged: List[Dict[str, Any]] = []
    for _stem, g in groups.items():
        merged.append(_merge_triplet_group(g))
    return merged


def run_quicklooks_for_range(
    outdir: Path,
    save_format: str,
    day_start: datetime,
    day_end: datetime,
    opt: QuicklooksOptions,
    save_root: Path,
) -> List[Path]:
    """
    Build quicklook recap figures for one calendar-day window [day_start, day_end).

    Args:
        outdir: Processed campaign output directory.
        save_format: ROI file format.
        day_start: Window start (inclusive, UTC).
        day_end: Window end (exclusive, UTC).
        opt: QuicklooksOptions controlling binning and figure types.
        save_root: Root directory for yyyy.mm.dd output subfolders.

    Returns:
        List of paths to PNG files written (may be empty if no data).

    Notes:
        - Filenames match MATLAB: yyyymmdd_quicklook#1.png and #2.png.

    Translated from make_masc_time_series.m (Christophe Praz 2015) and adapted for Python
    """
    if day_start.tzinfo is None:
        day_start = day_start.replace(tzinfo=timezone.utc)
    if day_end.tzinfo is None:
        day_end = day_end.replace(tzinfo=timezone.utc)
    day_start = day_start.astimezone(timezone.utc)
    day_end = day_end.astimezone(timezone.utc)

    qualities = tuple(q for q in opt.qualities if q in ("GOOD", "BAD"))
    if not qualities:
        qualities = ("GOOD",)

    rows = load_masc_points_for_interval(
        outdir,
        save_format,
        day_start,
        day_end,
        use_triplet=opt.use_triplet,
        qualities=qualities,
    )
    created: List[Path] = []
    if not rows:
        logger.warning(
            "No ROI rows in [%s, %s] under %s — skipping quicklooks.",
            day_start,
            day_end,
            outdir,
        )
        return created

    data = _aggregate_bins(rows, day_start, day_end, opt)
    if not data:
        return created

    if float(np.nansum(data["Nmasc_all"])) == 0:
        return created

    subdir = save_root / day_start.strftime("%Y.%m.%d")
    subdir.mkdir(parents=True, exist_ok=True)
    date_tag = day_start.strftime("%Y%m%d")

    if opt.savefigs:
        if opt.overview_classif:
            _save_overview_classif(data, day_start, day_end, opt, subdir, date_tag)
            created.append(subdir / f"{date_tag}_quicklook#1.png")
        if opt.overview_microstruct:
            _save_overview_microstruct(data, day_start, day_end, opt, subdir, date_tag)
            created.append(subdir / f"{date_tag}_quicklook#2.png")

    return created


def run_quicklooks_campaign(
    outdir: Path,
    save_format: str,
    t_min: datetime,
    t_max: datetime,
    opt: QuicklooksOptions,
    save_root: Path,
) -> List[Path]:
    """
    Run quicklooks for each calendar day overlapping a time window.

    Args:
        outdir: Processed campaign output directory.
        save_format: ROI file format.
        t_min: Campaign window start (UTC).
        t_max: Campaign window end (UTC).
        opt: QuicklooksOptions.
        save_root: Root directory for daily quicklook subfolders.

    Returns:
        List of all PNG paths written across the window.
    """
    if t_min.tzinfo is None:
        t_min = t_min.replace(tzinfo=timezone.utc)
    if t_max.tzinfo is None:
        t_max = t_max.replace(tzinfo=timezone.utc)
    t_min = t_min.astimezone(timezone.utc)
    t_max = t_max.astimezone(timezone.utc)

    if t_max <= t_min:
        return []

    d0 = t_min.date()
    d1 = t_max.date()
    all_created: List[Path] = []
    d = d0
    while d <= d1:
        day_start = datetime(d.year, d.month, d.day, tzinfo=timezone.utc)
        day_end = day_start + timedelta(days=1)
        win0 = max(day_start, t_min)
        win1 = min(day_end, t_max)
        if win1 > win0:
            all_created.extend(
                run_quicklooks_for_range(
                    outdir, save_format, day_start, day_end, opt, save_root
                )
            )
        d = d + timedelta(days=1)
    return all_created
