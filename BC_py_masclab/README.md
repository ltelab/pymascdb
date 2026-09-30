# BC_py_masclab — MASC snowflake processing pipeline (Python)

`BC_py_masclab` is a Python port of the MATLAB `masclab` toolbox (C. Praz, EPFL-LTE) used to
process images from a **Multi-Angle Snowflake Camera (MASC)**. Starting from the raw PNG
frames recorded by the three MASC cameras, the pipeline:

1. discovers the campaign folders to process for a given time window,
2. isolates the snowflake in each image (masking + region-of-interest detection),
3. computes ~90 morphological, textural and quality descriptors per snowflake view,
4. optionally matches the three camera views of the same particle (triplet mode),
5. applies the pre-trained MATLAB logistic classifiers (hydrometeor type, riming degree,
   melting) to every processed view,
6. produces daily "quicklook" recap figures,
7. optionally classifies the raw image stream as precipitation vs. blowing snow (GMM).

Everything is driven by a single YAML file (`config/config.yaml`); there are no command-line
arguments.

Documentation is split in three files:

| Document                                                   | Content                                                   |
|------------------------------------------------------------|-----------------------------------------------------------|
| this `README.md`                                           | layout, installation, input data, usage, outputs           |
| [`config/README.md`](config/README.md)                     | reference of every `config.yaml` key                       |
| [`src/pipeline/README.md`](src/pipeline/README.md)         | how each processing stage works, descriptors, ROI dictionary |

---

## Table of contents

- [Repository layout](#repository-layout)
- [Requirements and installation](#requirements-and-installation)
- [Input data format](#input-data-format)
- [Usage](#usage)
- [Outputs](#outputs)
- [Lineage and credits](#lineage-and-credits)

---

## Repository layout

```
BC_py_masclab/
├── config/
│   ├── config.yaml              # single source of configuration for every runner
│   └── README.md                # configuration reference
├── environment.yaml             # conda / micromamba environment (pinned dependencies)
├── scripts/                     # entry points (run these)
│   ├── runner_joblib.py         # image processing (+ optional post-processing steps)
│   ├── runner_classify.py       # classification only
│   ├── runner_quicklooks.py     # daily quicklook figures only
│   └── runner_bs.py             # blowing-snow classification only
├── src/
│   ├── core/
│   │   ├── config.py            # dataclasses + YAML loading, proc_params.txt writer
│   │   └── process_joblib.py    # batch orchestration, joblib parallelism, statistics
│   ├── dataio/
│   │   ├── uploaddirs.py        # discover yyyy.mm.dd/HH folders in the time window
│   │   ├── upload.py            # parse imgInfo.txt / dataInfo.txt -> PictureList
│   │   ├── load_roi_data.py     # save/load ROI dicts (.joblib / .pkl / .mat)
│   │   └── clean_outdir.py      # wipe the output folder
│   ├── preprocessing/
│   │   ├── masking.py           # border strips + background threshold
│   │   ├── roi.py               # ROIDetector: connected components, filtering, scoring
│   │   ├── brightening.py       # CLAHE brightening (MATLAB adapthisteq equivalent)
│   │   └── edge_detection.py    # legacy Sobel edge detector (not used by the pipeline)
│   ├── processing/
│   │   ├── single.py            # process one image end to end
│   │   └── triplet.py           # process the 3 views of one particle end to end
│   ├── features/
│   │   ├── basic.py             # process_basic_descriptors: the full descriptor set
│   │   ├── fmeasure.py          # focus measures (LAPM, HISE, WAVS, ...)
│   │   └── advanced.py          # standalone back-fill of "new" descriptors on old ROI files
│   ├── analysis/
│   │   └── quicklooks.py        # time aggregation + recap figures
│   ├── pipeline/
│   │   ├── steps.py             # run_classification / run_quicklooks / run_blowing_snow
│   │   └── README.md            # how the pipeline works
│   └── utils/
│       ├── descriptors_heplers.py  # convex hull, Dmax/D90, ellipse & circle fits, skeleton,
│       │                           # fractal dimension, Haralick, symmetry, blur index
│       ├── bwmorph_matlab.py       # MATLAB-compatible bwmorph (thin, endpoints, branchpoints…)
│       ├── image.py                # rangefilt / stdfilt
│       ├── identifiers.py          # cam_id / flake_id parsing from filenames
│       └── upload_vanderbilt.py    # metadata parser for Vanderbilt-style FLAKE_*.png datasets
├── classification/
│   ├── logit_trained_models/    # MATLAB-trained logistic classifiers (*.mat)
│   ├── load_classifier.py       # read MATLAB .mat logistic classifiers
│   ├── feature_vector.py        # build the canonical 96-element feature vector from an ROI
│   └── classify.py              # skew transforms, standardisation, inference, write-back
└── blowing_snow/
    ├── features_bs.py           # sky-median filtering + per-image descriptors
    ├── gmm_bs.py                # load the 2-component GMM from .mat, classification angle
    ├── classify_bs.py           # per-image and per-triplet labelling
    └── models/                  # gmfit_random_all.mat, MASK.mat, MASK_davos.mat
```

All modules import each other with absolute package paths (`src.…`, `classification.…`,
`blowing_snow.…`); the runners insert the `BC_py_masclab/` root into `sys.path`, so run them
from anywhere but keep the folder structure intact.

---

## Requirements and installation

### Software

- **Python 3.11** (pinned in `environment.yaml`).
- A conda-compatible package manager: [Miniconda](https://docs.anaconda.com/miniconda/) or
  [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).
- The MATLAB-trained classifiers (`*.mat`) from the original `masclab`, shipped in
  `classification/logit_trained_models/` (the folder `config.yaml` points to by default).
  MATLAB itself is **not** required. There are 3 of them:
   - `logit_v.1.1_FINAL3_6classes_weighted_scheme.mat`
   - `logit_v.1.1_FINAL3_riming_5classes_weighted_scheme.mat`
   - `logit_v.1.1_FINAL3_melting_2classes_noweight.mat`

### Python dependencies

All dependencies are pinned in `environment.yaml` (environment name `py_masclab`). The
pipeline itself relies on `numpy`, `scipy`, `opencv`, `scikit-image`, `pywavelets`,
`matplotlib`, `joblib`, `tqdm`, `pyyaml` and `pandas`; the remaining entries (`pytest`,
`ipykernel`, `bokeh`, `shapely`, `pillow`) support testing and interactive analysis.

### Install

From the `BC_py_masclab/` folder, with Miniconda:

```bash
conda env create -f environment.yaml
conda activate py_masclab
```

or with micromamba:

```bash
micromamba create -f environment.yaml
micromamba activate py_masclab
```

Sanity check (should print nothing):

```bash
python -c "import numpy, scipy, cv2, skimage, pywt, matplotlib, joblib, tqdm, yaml, pandas"
```

To update an existing environment after `environment.yaml` changes:
`conda env update -f environment.yaml --prune` (or `micromamba install -f environment.yaml`).

---

## Input data format

The pipeline expects the standard MASC folder hierarchy produced by the acquisition software
(or by `masclab`'s `reorganize_MASC_folders.m`):

```
<campaigndir>/
└── yyyy.mm.dd/                       # one folder per day
    └── HH/                           # one folder per hour
        ├── imgInfo.txt
        ├── dataInfo.txt
        └── yyyy.mm.dd_HH.MM.SS_flake_<ID>_cam_<C>.png   # one PNG per view
```

- **Images**: 8-bit grayscale PNGs. Filenames must contain `flake_<ID>` and `cam_<C>` (with or
  without the underscore before the camera index); both are parsed by `src/utils/identifiers.py`.
- **`imgInfo.txt`** — one line per image, whitespace-separated:

  ```
  <flake_id> <cam> <MM.DD.YYYY> <HH:MM:SS.ffffff> <filename> <extra>
  45924      0     06.20.2015   09:00:00.718243   2015.06.20_09.00.00_flake_45924_cam_0.png 0
  ```

- **`dataInfo.txt`** — one line per particle (fall speed in m/s):

  ```
  <flake_id> <MM.DD.YYYY> <HH:MM:SS.ffffff> <fallspeed>
  45924      06.20.2015   09:00:00.843043   1.9198
  ```

  If `dataInfo.txt` is missing, fall speeds are set to `NaN` and processing continues.

> **NOTE — the image list comes from `imgInfo.txt`, not from a scan of the folder.**
> The standard pipeline (`runner_joblib.py` → `src/core/process_joblib.py`) never globs the
> folder for `*.png`. `src/dataio/uploaddirs.py` selects the `yyyy.mm.dd/HH` folders **by their
> name/date only** (it does not look at their contents), and `src/dataio/upload.py` builds the
> list of images to process **exclusively from `imgInfo.txt`**. Consequences:
> - If `imgInfo.txt` is missing or unreadable, `upload()` logs an error and returns an **empty**
>   `PictureList`; the folder is then **skipped silently** (no images processed, no "no PNG"
>   warning).
> - PNG files physically present in the folder but **not listed** in `imgInfo.txt` are ignored.
> - The only file-level check is **per image**: `single.py` / `triplet.py` raise
>   `FileNotFoundError` if a filename listed in `imgInfo.txt` does not exist on disk (or
>   `ValueError` if it cannot be decoded).
>
> In other words, there is **no "does this folder contain PNGs?" pre-check** before processing
> starts — the pipeline trusts `imgInfo.txt`. (This differs from the blowing-snow path in
> `blowing_snow/features_bs.py` and the Vanderbilt uploader below, which do scan for `*.png`.)

Vanderbilt-style datasets (`FLAKE_*.png` without `imgInfo.txt`) can be read with
`src/utils/upload_vanderbilt.py`, which builds the same `PictureList` structure.

---

## Usage

Edit `config/config.yaml` — at least `label.campaigndir`, and the `pipeline` flags for the
post-processing steps you want (see [`config/README.md`](config/README.md)) — then run one of
the scripts. All of them can be launched from any working directory.

### Full run (processing + selected post-processing steps)

```bash
python scripts/runner_joblib.py
```

`runner_joblib.py`:

1. loads and prints the configuration summary,
2. wipes `outdir` if `processing.clean_outdir` is `true` (**the whole folder is deleted**),
3. processes every image/triplet of the time window in parallel and writes `proc_stats.txt`,
4. runs, in this order and only if enabled in the `pipeline` section: classification →
   quicklooks → blowing snow.

### Standalone steps

Each post-processing step can be re-run on an existing `outdir` without reprocessing images:

```bash
python scripts/runner_classify.py     # (re)label existing ROI files
python scripts/runner_quicklooks.py   # figures from labelled ROI files
python scripts/runner_bs.py           # blowing-snow labels from the raw images
```

---

## Outputs

With `outdir: auto` and `campaigndir: .../RAW`, everything is written to `.../RAW_PROCESSED/`:

```
RAW_PROCESSED/
├── proc_params.txt                 # processing parameters used (human readable)
├── proc_stats.txt                  # counts and timing breakdown of the run
├── features_bs.joblib              # blowing-snow features cache (if run)
├── blowing_snow_all.csv            # blowing-snow labels per image (if run)
├── blowing_snow_triplet.csv        # blowing-snow labels per particle (if run)
├── Quicklooks/
│   └── yyyy.mm.dd/
│       ├── yyyymmdd_quicklook#1.png
│       └── yyyymmdd_quicklook#2.png
└── yyyy.mm.dd/
    └── HH/
        ├── DATA/    {GOOD,BAD}/  <image stem>.joblib      # ROI dictionary (or .pkl / .mat)
        ├── IMAGES/  {GOOD,BAD}/  <image stem>.png         # cropped ROI image
        └── FIGURES/ {GOOD,BAD}/  <image stem>_fig.png     # diagnostic figure (if save_figs)
```

`GOOD` holds views with a valid detection (`flag_roi == 'GOOD'`); `BAD` holds views with no
ROI or a rejected one. Classification and quicklooks work on the `DATA/` files.

Each ROI file is a plain Python `dict` holding the descriptors, the bookkeeping fields and,
after classification, the `label_*`, `riming_*`, `melting_*` fields; its content is detailed
in [`src/pipeline/README.md`](src/pipeline/README.md#9-the-roi-dictionary). With
`save_format: mat` the same dictionary is stored as a MATLAB struct named `roi`, loadable
from `masclab`.

---

## Lineage and credits

