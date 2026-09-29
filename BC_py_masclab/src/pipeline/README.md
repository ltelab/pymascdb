# How the pipeline works

This document describes what each stage of `BC_py_masclab` does, module by module. Parameter
names in backticks refer to keys of `config/config.yaml`; see
[`config/README.md`](../../config/README.md) for their meaning and defaults.

`src/pipeline/steps.py` holds the three post-processing entry points
(`run_classification`, `run_quicklooks`, `run_blowing_snow`) called by `scripts/runner_joblib.py`
and by the standalone runners. Image processing itself is orchestrated by
`src/core/process_joblib.py`.

## Contents

1. [Folder and metadata discovery](#1-folder-and-metadata-discovery)
2. [Single-image processing](#2-single-image-processing)
3. [Triplet processing](#3-triplet-processing)
4. [ROI detection](#4-roi-detection)
5. [Descriptor extraction](#5-descriptor-extraction)
6. [Classification](#6-classification)
7. [Quicklooks](#7-quicklooks)
8. [Blowing-snow classification](#8-blowing-snow-classification)
9. [The ROI dictionary](#9-the-roi-dictionary)

---

## 1. Folder and metadata discovery

`src/dataio/uploaddirs.py` walks `campaigndir`, keeps the `yyyy.mm.dd/HH` folders whose date
falls inside `[starthr_vec, endhr_vec]`, and returns them sorted. For each folder,
`src/dataio/upload.py` parses `imgInfo.txt` and `dataInfo.txt` into a `PictureList`
dataclass (`files`, `id`, `cam`, `time_vec`, `time_str`, `time_num`, `fallspeed`, `fallid`,
`id_unique`).

`src/core/process_joblib.py` (`masc_process`) then builds the list of work items — one per
image in single mode, one per unique flake ID in triplet mode — and dispatches them with
`joblib.Parallel` (loky backend). Per-item timings (loading, masking, ROI detection,
descriptors, matching, plotting, saving) are accumulated and summarised in
`<outdir>/proc_stats.txt`, together with the number of triplets processed, complete, matched,
unmatched and processed independently.

## 2. Single-image processing

`src/processing/single.py` (`process_single_image`) runs for one PNG:

1. **Load** the image as an 8-bit grayscale array (`cv2.imread`); intensity descriptors are
   normalised to [0, 1] later by dividing by 255.
2. **Mask** (`src/preprocessing/masking.py`): blank the border strips defined by `discardmat`
   according to the camera position in `camera_order` — top and bottom strips for the three
   core cameras, plus the left strip for the left camera (`camera_order[0]`) and the right
   strip for the right camera (`camera_order[2]`) — then set every pixel `<= backthresh` to 0.
3. **Binary mask**: `data > 0`, holes filled with `scipy.ndimage.binary_fill_holes`.
4. **ROI detection** (see §4) → the best ROI plus a `GOOD` / `BAD` flag and a status string.
5. **Descriptor extraction** (see §5) on the best ROI.
6. **Bookkeeping**: `name`, `cam`, `id`, `tnum` (UNIX timestamp), `fallspeed`, `flag_roi`,
   `status`, `n_roi`, global/local centroids.
7. **Save** (if `saveresults`): the cropped ROI image, the ROI dictionary, and optionally the
   diagnostic figure, under `GOOD/` or `BAD/` depending on the flag. The first image also
   writes `<outdir>/proc_params.txt`.

An image with no valid ROI still produces a minimal ROI file (in `BAD/`) so that every input
is accounted for.

## 3. Triplet processing

`src/processing/triplet.py` (`process_triplet_image`) handles all views of one flake ID at
once:

1. Load and mask each available view (steps 1–4 above run per camera, all candidate ROIs
   kept).
2. **Match ROIs across cameras**: because the cameras are horizontally aligned, the same
   particle must appear at a similar vertical position with a similar height in every view.
   For each candidate combination, the vertical centroid offset must be within
   `matching_tol_pix` and the height difference within `matching_tol_percent` % of the
   reference height.
3. **Choose the best combination**: among matching combinations, the one with the highest
   mean `area × focus` wins (the sharpest, largest consistent particle).
4. Extract descriptors for each selected view, store `flag = 2` (matched) and save one ROI
   file per camera, plus one figure per view.
5. **Fallback**: if fewer than `triplet_required_views` views exist and
   `triplet_allow_partial` is `false`, or if no combination matches, each view is processed
   independently as in §2 (counted as "processed indep." in `proc_stats.txt`).

## 4. ROI detection

`src/preprocessing/roi.py` (`ROIDetector.detect`):

1. Label connected components of the binary mask (`skimage.measure.label`, 8-connectivity)
   and drop those smaller than `min_area`.
2. Drop components that touch the blanked border strips (they are probably cut particles).
3. For each remaining component compute: area, bounding box, centroid, mean/max intensity,
   local range intensity (`rangefilt`), focus (Laplacian energy) and `area × focus`.
4. The **best ROI** is the one maximising `area × focus`. `area_focus_ratio` (best / second
   best) is stored as an ambiguity indicator.
5. **Validation** (`_validate_roi`): the best ROI is `GOOD` unless its bounding box is
   smaller than `sizemin`, its mean intensity is below `minbright`, or its maximum intensity
   is below `max_intensthresh`; otherwise it is `BAD` with an explanatory `status`.

## 5. Descriptor extraction

`src/features/basic.py` (`process_basic_descriptors`) fills the ROI dictionary. Grouped by
family (helper functions live in `src/utils/descriptors_heplers.py` and
`src/features/fmeasure.py`):

| Family                  | Fields (main)                                                                                                     |
|-------------------------|-------------------------------------------------------------------------------------------------------------------|
| Size                    | `area`, `area_porous`, `perim`, `width`, `height`, `eq_radius`, `Dmean`, `Dmax`, `Dmax_theta`, `D90` (`Dmax_0`, `Dmax_90`, `AR`) |
| Shape                   | `E` (fitted ellipse `a`, `b`, `theta`, centre), `E_in` (largest inscribed), `E_out` (smallest circumscribed), `C_out` (circumscribed circle), `Rect` (min-area rectangle, `rectangularity`, `width`, `height`), `roundness`, `compactness`, `complex` |
| Convex hull / holes     | `hull` (`solidity`, `convexity`, `perim`, vertices), `nb_holes`, `holes_mask`                                     |
| Topology                | `skel` (`endpoints`, `branches`, `length`, `density`), `F` (box-counting fractal dimension), `F_jac`               |
| Texture / intensity     | `mean_intens`, `max_intens`, `min_intens`, `std`, `contrast`, `range_intens`, `local_std`, `local_std5`, `local_std7`, `hist_entropy`, `H` (Haralick `Contrast`, `Correlation`, `Energy`, `Homogeneity`) |
| Focus                   | `lap` (LAPM), `wavs` (WAVS), `focus`, `area_focus`, `area_lap`, `area_range`, `range_complex`                    |
| Symmetry                | `Sym` (radial-profile symmetry descriptors)                                                                       |
| Quality                 | `xhi` — the "magic" image-quality index combining size, focus and contrast; used by classification (`blurry_threshold`) and quicklooks (`xi_thresh`) |
| Brightened recompute    | `new` — the texture/focus subset recomputed on a CLAHE-brightened crop when `flakebrighten` is `true`             |
| Raw data                | `data` (cropped grey image), `bw_mask`, `bw_mask_filled`, `bw_perim`, `x_perim`/`y_perim`, `x`/`y`/`x_loc`/`y_loc` |

Small holes (`< min_hole_area`) are filled before the geometric descriptors are computed.
Focus measures follow Pertuz et al. (2013); Haralick features are computed on a symmetric
256-level GLCM over the four standard offsets.

## 6. Classification

`classification/` applies the MATLAB-trained multinomial logistic regressions to every ROI
file under `outdir` (`classify_campaign`, called three times by `run_classification`):

1. `load_classifier.py` reads the `.mat` file and exposes weights, the indices of the
   features used (converted to 0-based), per-feature normalisation statistics, the skewness
   transform flags and the class labels.
2. `feature_vector.py` builds the canonical **96-element feature vector** from the ROI dict
   (areas, ellipse ratios, intensities, Haralick, hull, skeleton, symmetry, brightened `new`
   descriptors, etc.). Missing fields default to `0.0`.
3. `classify.py` selects the classifier's features, applies the skewness-driven transforms
   (`log`, `sqrt`, `exp`, square), standardises with the stored mean/std, then computes the
   softmax (multiclass) or sigmoid (binary) probabilities.
4. Results are written back into each ROI file:

   | Prefix    | Fields written                                         | Values                                         |
   |-----------|--------------------------------------------------------|------------------------------------------------|
   | `label`   | `label_ID`, `label_name`, `label_probs`                | 1–6: SP, CC, PC, AG, GR, CPC                   |
   | `riming`  | `riming_ID`, `riming_name`, `riming_probs`             | 1–5 riming degree                              |
   | `melting` | `melting_ID`, `melting_name`, `melting_probs`          | 0 = dry, 1 = melting                           |

   Views with `xhi <= blurry_threshold` receive `*_ID = -9` and `*_name = "blurry"`
   (probabilities are still stored).

## 7. Quicklooks

`src/analysis/quicklooks.py` (`run_quicklooks_campaign`):

1. Load every ROI file of the requested `qualities` in the time window (if the window is
   `"all"`, bounds are discovered from the files themselves).
2. Keep views with `xhi >= xi_thresh`; in triplet mode merge the three views of a particle
   (probabilities summed and renormalised, Dmax = max over views, aspect ratio = min,
   orientation = min |angle|, melting = mean probability, riming index = mean).
3. Slide `Nmin_interval`-minute bins every `Nmin_shift` minutes and compute, per bin:
   particle count, class proportions and dominant class, melting fraction, mean riming
   index, Dmax / aspect-ratio / orientation statistics and fall speed.
4. Export per day, under `<savepath>/yyyy.mm.dd/`:
   - `yyyymmdd_quicklook#1.png` — classification overview (counts, class proportions, dominant
     class, melting and riming time series);
   - `yyyymmdd_quicklook#2.png` — microstructure overview (Dmax, aspect ratio, orientation,
     fall speed).

## 8. Blowing-snow classification

`blowing_snow/` works on the **raw PNGs**, independently of ROI processing:

1. `features_bs.py` sorts the images chronologically, splits them into events separated by
   more than `s` hours, and for each block of `block_size` images per camera estimates the
   sky background as the per-pixel median. Each image is background-subtracted, binarised and
   described by porosity, Dmax (0.7 quantile of particle sizes), box-counting fractal
   dimension and the local photo frequency (`w`). Features are cached in
   `<outdir>/features_bs.joblib`.
2. `gmm_bs.py` loads the 2-component Gaussian mixture trained in MATLAB
   (`gmfit_random_all.mat`), applies the same transforms (`log(Dmax+1)`, standardisation) and
   converts the posterior into a normalised **classification angle** in [0, 1].
3. `classify_bs.py` labels each image `0` (precipitation) if the angle is `< 0.5`, `1`
   (blowing snow) otherwise, flags "mixed" cases near the boundary, and aggregates labels per
   triplet (same minute and flake ID).
4. Outputs: `<outdir>/blowing_snow_all.csv` (per image) and
   `<outdir>/blowing_snow_triplet.csv` (per particle).

## 9. The ROI dictionary

Each processed view is stored as a plain Python `dict` (`.joblib` by default; `.pkl` or a
MATLAB struct named `roi` with `save_format: mat`). It contains:

- the descriptors of §5,
- bookkeeping: `name`, `cam`, `id`, `tnum`, `fallspeed`, `flag`, `flag_roi`, `status`,
  `n_roi`, `centroid`, `centroid_global`, `centroid_local`,
- after classification: `label_*`, `riming_*`, `melting_*` (see §6).

Read them with `src/dataio/load_roi_data.py` (`load_roi_data` for one file,
`load_all_roi_in_folder` for a `DATA/` folder, `print_roi_summary` for a quick look).
