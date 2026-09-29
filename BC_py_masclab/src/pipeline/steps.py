"""
Shared post-processing steps for the MASC pipeline.

This module orchestrates classification, daily quicklooks, and blowing-snow
classification after ROI processing.

Functions for:
- Running the three MASC classifiers on processed ROI files
- Generating daily quicklook recap figures
- Classifying blowing-snow events from raw images

Translated from:
- make_predictions_for_campaign.m (Christophe Praz 2015)
- merge_predictions_for_campaign.m (Christophe Praz 2015)
- run_bs.m (Christophe Praz 2018)
and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

from __future__ import annotations

from datetime import timezone
from pathlib import Path
from typing import List, Literal, Optional

from src.core.config import LabelConfig, ProcessingConfig, _project_root, load_yaml_section

try:
    import yaml
    _YAML_OK = True
except ImportError:
    yaml = None  # type: ignore
    _YAML_OK = False


def _resolve_yaml_path(yaml_path: Optional[Path]) -> Path:
    """Return config.yaml path, defaulting to project config/config.yaml."""
    if yaml_path is None:
        return _project_root() / "config" / "config.yaml"
    return Path(yaml_path)


def _load_classification_section(yaml_path: Path) -> dict:
    """Load and validate the [classification] section from config.yaml."""
    if not _YAML_OK:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")
    section = load_yaml_section(yaml_path, "classification")
    if not section:
        raise KeyError(
            "Missing [classification] section in config.yaml. "
            "Please add classifiers_dir, class, riming, melting keys."
        )
    return section


def _load_quicklooks_section(yaml_path: Path) -> dict:
    """Load the [quicklooks] section from config.yaml (empty dict if missing)."""
    if not _YAML_OK:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")
    return load_yaml_section(yaml_path, "quicklooks") or {}


def run_classification(
    label_cfg: LabelConfig,
    proc_cfg: ProcessingConfig,
    yaml_path: Optional[Path] = None,
    project_root: Optional[Path] = None,
) -> None:
    """
    Run the three MASC classifiers and write labels into ROI files under outdir.

    Args:
        label_cfg: Label configuration with outdir path.
        proc_cfg: Processing configuration (save_format, etc.).
        yaml_path: Optional path to config.yaml.
        project_root: Optional project root for relative classifier paths.

    Returns:
        None

    Raises:
        FileNotFoundError: If outdir or classifiers are missing.
        ImportError: If PyYAML is not installed.
        KeyError: If the [classification] section is incomplete.

    Translated from make_predictions_for_campaign.m (Christophe Praz 2015) and adapted for Python
    """
    yaml_path = _resolve_yaml_path(yaml_path)
    root_dir = Path(project_root) if project_root else _project_root()
    cl_cfg = _load_classification_section(yaml_path)

    classifiers_dir = Path(cl_cfg["classifiers_dir"])
    if not classifiers_dir.is_absolute():
        classifiers_dir = (root_dir / classifiers_dir).resolve()

    blurry_threshold = float(cl_cfg.get("blurry_threshold", 0.0))
    n_jobs = int(cl_cfg.get("n_jobs", 1))
    save_format = proc_cfg.save_format
    outdir = label_cfg.outdir

    print("=" * 70)
    print("MASC Classification Pipeline")
    print("=" * 70)
    print(f"  Output directory : {outdir}")
    print(f"  Classifiers dir  : {classifiers_dir}")
    print(f"  Save format      : {save_format}")
    print(f"  Blurry threshold : {blurry_threshold}")
    print(f"  Workers          : {n_jobs}")
    print()

    if not outdir.exists():
        raise FileNotFoundError(
            f"outdir does not exist: {outdir}. Run processing first (runner_joblib.py)."
        )

    from classification.classify import classify_campaign

    cl_class_path = classifiers_dir / cl_cfg["class"]
    print(f"[1/3] Hydrometeor class classifier: {cl_cfg['class']}")
    classify_campaign(
        outdir=outdir,
        classifier_path=cl_class_path,
        save_format=save_format,
        blurry_threshold=blurry_threshold,
        label_prefix="label",
        n_jobs=n_jobs,
    )
    print()

    cl_riming_path = classifiers_dir / cl_cfg["riming"]
    print(f"[2/3] Riming classifier: {cl_cfg['riming']}")
    classify_campaign(
        outdir=outdir,
        classifier_path=cl_riming_path,
        save_format=save_format,
        blurry_threshold=blurry_threshold,
        label_prefix="riming",
        n_jobs=n_jobs,
    )
    print()

    cl_melting_path = classifiers_dir / cl_cfg["melting"]
    print(f"[3/3] Melting classifier: {cl_cfg['melting']}")
    classify_campaign(
        outdir=outdir,
        classifier_path=cl_melting_path,
        save_format=save_format,
        blurry_threshold=blurry_threshold,
        label_prefix="melting",
        n_jobs=n_jobs,
    )
    print()

    print("=" * 70)
    print("Classification completed.")
    print(f"Labels written into ROI files in: {outdir}")
    print("=" * 70)


def _quicklooks_options_from_dict(d: dict, process_use_triplet: bool):
    """Build QuicklooksOptions from a YAML [quicklooks] section dict."""
    from src.analysis.quicklooks import QuicklooksOptions

    return QuicklooksOptions(
        pixres_mm=float(d.get("pixres_mm", 33.5 / 1000.0)),
        xi_thresh=float(d.get("xi_thresh", 9.0)),
        Nmin_interval=int(d.get("Nmin_interval", 30)),
        Nmin_shift=int(d.get("Nmin_shift", 10)),
        Nclasses_masc=int(d.get("Nclasses_masc", 6)),
        MASC_classes=tuple(d.get("MASC_classes", ("SP", "CC", "PC", "AG", "GR", "CPC"))),
        MASC_classes_desired=tuple(d.get("MASC_classes_desired", (1, 2, 3, 4, 5, 6))),
        N_MascSamples_min=int(d.get("N_MascSamples_min", 0)),
        use_triplet=(
            process_use_triplet
            if d.get("use_triplet") is None
            else bool(d.get("use_triplet"))
        ),
        OR_180=bool(d.get("OR_180", True)),
        savefigs=bool(d.get("savefigs", True)),
        overview_classif=bool(d.get("overview_classif", True)),
        overview_microstruct=bool(d.get("overview_microstruct", True)),
        disp_now=bool(d.get("disp_now", True)),
        qualities=tuple(d.get("qualities", ("GOOD",))),
    )


def run_quicklooks(
    label_cfg: LabelConfig,
    proc_cfg: ProcessingConfig,
    yaml_path: Optional[Path] = None,
) -> List[Path]:
    """
    Generate daily quicklook PNGs from classified ROI files under outdir.

    Args:
        label_cfg: Label configuration with outdir and time window.
        proc_cfg: Processing configuration (save_format, use_triplet_algo).
        yaml_path: Optional path to config.yaml.

    Returns:
        List of paths to written figure files.

    Raises:
        FileNotFoundError: If outdir is empty and time bounds cannot be discovered.

    Translated from merge_predictions_for_campaign.m (Christophe Praz 2015) and adapted for Python
    """
    from src.analysis.quicklooks import (
        discover_time_bounds_from_outdir,
        run_quicklooks_campaign,
    )

    yaml_path = _resolve_yaml_path(yaml_path)
    ql_raw = _load_quicklooks_section(yaml_path)

    outdir = Path(label_cfg.outdir)
    save_fmt = str(proc_cfg.save_format)

    savepath_raw = ql_raw.get("savepath", "auto")
    if isinstance(savepath_raw, str) and savepath_raw.lower() == "auto":
        save_root = outdir / "Quicklooks"
    else:
        save_root = Path(savepath_raw)

    opt = _quicklooks_options_from_dict(ql_raw, proc_cfg.use_triplet_algo)

    t_min = label_cfg.starthr_vec
    t_max = label_cfg.endhr_vec
    if t_min.tzinfo is None:
        t_min = t_min.replace(tzinfo=timezone.utc)
    if t_max.tzinfo is None:
        t_max = t_max.replace(tzinfo=timezone.utc)
    t_min = t_min.astimezone(timezone.utc)
    t_max = t_max.astimezone(timezone.utc)

    if t_min.year < 1900 or t_max.year > 2100:
        bounds = discover_time_bounds_from_outdir(
            outdir,
            save_fmt,
            tuple[Literal['GOOD', 'BAD'], ...](q for q in opt.qualities if q in ("GOOD", "BAD")),  # type: ignore[arg-type]
        )
        if bounds is None:
            raise FileNotFoundError(
                "No ROI files found under outdir (needed for open-ended time window)."
            )
        dmin, dmax = bounds
        dmin = dmin.astimezone(timezone.utc)
        dmax = dmax.astimezone(timezone.utc)
        if t_min.year < 1900:
            t_min = dmin
        if t_max.year > 2100:
            t_max = dmax
        print(f"Time bounds (after discovery): {t_min} … {t_max}")

    print("=" * 70)
    print("MASC quicklooks (aggregation + recap figures)")
    print("=" * 70)
    print(f"  outdir     : {outdir}")
    print(f"  save_root  : {save_root}")
    print(f"  use_triplet: {opt.use_triplet}")
    print(f"  qualities  : {opt.qualities}")
    print()

    paths = run_quicklooks_campaign(outdir, save_fmt, t_min, t_max, opt, save_root)
    print(f"Wrote {len(paths)} figure file(s).")
    for p in paths[:20]:
        print(f"  {p}")
    if len(paths) > 20:
        print(f"  … ({len(paths) - 20} more)")

    return paths


def _clamp_to_pandas_bounds(dt):
    """Keep a datetime safely within the pandas Timestamp range (1677-2262).

    Bounds are kept well inside the limits so the later pd.Timestamp() never
    overflows; MASC data is well within this window anyway.
    """
    from datetime import datetime as _dt

    return min(max(dt, _dt(1678, 1, 1)), _dt(2261, 12, 31))


def run_blowing_snow(
    label_cfg: LabelConfig,
    proc_cfg: ProcessingConfig,
    yaml_path: Optional[Path] = None,
) -> List[Path]:
    """
    Classify raw MASC images as precipitation or blowing snow.

    Args:
        label_cfg: Label configuration (campaigndir, time window, outdir).
        proc_cfg: Processing configuration (unused except for consistency).
        yaml_path: Optional path to config.yaml [blowing_snow] section.

    Returns:
        List of paths to the written CSV files (all and triplet summaries).

    Translated from run_bs.m (Christophe Praz 2018) and adapted for Python
    """
    import joblib

    from blowing_snow.classify_bs import classify_blowing_snow
    from blowing_snow.features_bs import extract_all
    from blowing_snow.gmm_bs import DEFAULT_MODEL, load_gmm

    yaml_path = _resolve_yaml_path(yaml_path)
    bs = load_yaml_section(yaml_path, "blowing_snow")

    in_path = label_cfg.campaigndir
    if str(bs.get("in_path", "auto")).lower() != "auto":
        in_path = Path(bs["in_path"])

    model_path = DEFAULT_MODEL
    if str(bs.get("model_path", "auto")).lower() != "auto":
        model_path = Path(bs["model_path"])

    davos = bs.get("davos", "auto")
    davos = "davos" in str(in_path).lower() if str(davos).lower() == "auto" else bool(davos)

    block_size = int(bs.get("block_size", 15))
    s = float(bs.get("s", 2))
    w = int(bs.get("w", 10))
    reprocess = bool(bs.get("reprocess", True))

    outdir = label_cfg.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    feat_path = outdir / "features_bs.joblib"
    all_out = outdir / "blowing_snow_all.csv"
    triplet_out = outdir / "blowing_snow_triplet.csv"

    print("=" * 70)
    print("MASC Blowing-snow classification")
    print("=" * 70)
    print(f"  Input images : {in_path}")
    print(f"  Sky mask     : {'davos' if davos else 'default'}")
    print(f"  GMM model    : {model_path.name}")
    print(f"  block/s/w    : {block_size} / {s} / {w}")
    print(f"  Output dir   : {outdir}")
    print()

    if reprocess or not feat_path.exists():
        t_start = _clamp_to_pandas_bounds(label_cfg.starthr_vec)
        t_stop = _clamp_to_pandas_bounds(label_cfg.endhr_vec)
        features = extract_all(in_path, block_size, w, s, t_start, t_stop, davos=davos)
        joblib.dump(features, feat_path, compress=3)
    else:
        print(f"Reusing existing features: {feat_path}")
        features = joblib.load(feat_path)

    gmm = load_gmm(model_path)
    all_df, triplet_df = classify_blowing_snow(features, gmm=gmm)

    all_df.to_csv(all_out, index=False, date_format="%d/%m/%Y %H:%M:%S")
    triplet_df.to_csv(triplet_out, index=False, date_format="%d/%m/%Y %H:%M:%S")
    print(f"Wrote {all_out}")
    print(f"Wrote {triplet_out}")
    return [all_out, triplet_out]
