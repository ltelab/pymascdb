# Configuration reference — `config/config.yaml`

Every runner in `scripts/` reads `BC_py_masclab/config/config.yaml`; there are no
command-line arguments. The file has six sections, documented below. The configuration is
parsed by `src/core/config.py` (`load_from_yaml`, `load_pipeline_from_yaml`,
`load_yaml_section`) into the `LabelConfig`, `ProcessingConfig` and `PipelineConfig`
dataclasses; the `classification`, `quicklooks` and `blowing_snow` sections are read as
plain dictionaries by `src/pipeline/steps.py`.

All keys of the `label` and `processing` sections are **required** (no defaults in code).
Keys of the other sections fall back to the defaults given below when omitted.

For what each step does with these parameters, see [`src/pipeline/README.md`](../src/pipeline/README.md).

---

## `pipeline` — what `runner_joblib.py` does after image processing

| Key                  | Default | Description                                                    |
|----------------------|---------|----------------------------------------------------------------|
| `run_classification` | `false` | run the three classifiers on the freshly processed ROI files   |
| `run_quicklooks`     | `false` | generate daily quicklook figures                               |
| `run_blowing_snow`   | `false` | run the blowing-snow GMM classification on the raw images      |

The steps run in this order: classification → quicklooks → blowing snow. The standalone
runners (`runner_classify.py`, `runner_quicklooks.py`, `runner_bs.py`) ignore this section.

## `label` — where the data is and which period to process

| Key           | Description                                                                                         |
|---------------|-----------------------------------------------------------------------------------------------------|
| `campaigndir` | root folder containing the `yyyy.mm.dd/` day folders                                                |
| `outdir`      | output folder. `"auto"` → `<parent of campaigndir>/<campaigndir name>_PROCESSED`                     |
| `starthr_vec` | start of the time window, `"YYYY-MM-DD HH:MM:SS"` (ISO also accepted) or `"all"`                     |
| `endhr_vec`   | end of the time window, same formats. `"all"` on both keys processes every folder found              |

On Windows, write paths with forward slashes (`C:/data/RAW`) or quote them.

## `processing` — image processing parameters

| Key                        | Default | Description                                                                                      |
|----------------------------|---------|--------------------------------------------------------------------------------------------------|
| `use_triplet_algo`         | `true`  | triplet mode (match the 3 views of a particle) vs. independent single-image mode                  |
| `parallel`, `max_workers`  | `true`, `6` | joblib parallelism; `max_workers` also caps memory usage (each worker holds one full image) |
| `saveresults`              | `true`  | write ROI files, cropped images and (optionally) figures                                          |
| `save_format`              | `joblib`| `joblib` \| `pkl` \| `mat` (MATLAB-compatible struct)                                            |
| `clean_outdir`             | `true`  | **wipe `outdir` at start-up** — disable when several campaigns share an output folder            |
| `generate_figs`, `display_figs`, `save_figs` | `true`, `false`, `true` | per-ROI diagnostic figures (image + ellipse fits + skeleton). `display_figs` needs a GUI backend; disabling `save_figs` speeds up the run noticeably |
| `flakebrighten`            | `true`  | recompute texture/focus descriptors on a CLAHE-brightened crop (stored under `roi['new']`)        |
| `backthresh`               | `12`    | grey level at or under which a pixel is considered background                                    |
| `min_area`                 | `20`    | minimum connected-component area (px) to be considered an ROI                                    |
| `sizemin`                  | `0`     | minimum bounding-box size (px) for a valid ROI                                                   |
| `minbright`, `max_intensthresh` | `0`, `0` | minimum mean / maximum grey level for a valid ROI (0 disables)                              |
| `min_hole_area`            | `10`    | holes smaller than this are filled before descriptor extraction                                  |
| `discardmat`               | `[400, 400, 410, 410]` | image border strips to blank out `[top, bottom, left, right]` (px)                    |
| `matching_tol_pix`         | `50`    | triplet matching: max vertical offset between views (px)                                         |
| `matching_tol_percent`     | `100`   | triplet matching: max relative height difference between views, in **percent** (0–100)           |
| `camera_order`             | `[0, 1, 2]` | camera indices `[left, middle, right]`; drives which border strips apply to which view        |
| `triplet_required_views`   | `3`     | number of views needed for a "complete" triplet                                                  |
| `triplet_allow_partial`    | `true`  | if fewer views are available, still process the ones present                                     |

## `classification`

| Key                          | Default | Description                                                                                  |
|------------------------------|---------|----------------------------------------------------------------------------------------------|
| `classifiers_dir`            | —       | folder containing the MATLAB `.mat` classifiers; relative paths are resolved from `BC_py_masclab/` (`classification/logit_trained_models`) |
| `class`, `riming`, `melting` | —       | file names of the three classifiers inside `classifiers_dir`                                 |
| `blurry_threshold`           | `0`     | ROIs with `xhi <= blurry_threshold` get `*_ID = -9` and `*_name = "blurry"`                  |
| `n_jobs`                     | `1`     | parallel workers for feature-vector extraction (`1` = sequential)                            |

## `quicklooks`

| Key                              | Default | Description                                                                             |
|----------------------------------|---------|-----------------------------------------------------------------------------------------|
| `savepath`                       | `"auto"` | `"auto"` → `<outdir>/Quicklooks/yyyy.mm.dd/`                                           |
| `pixres_mm`                      | `0.0335` | pixel size in mm (standard MASC) used to convert Dmax to mm                            |
| `xi_thresh`                      | `9`     | minimum `xhi` (image quality) for a view to be included                                  |
| `Nmin_interval`, `Nmin_shift`    | `30`, `10` | width and step of the sliding time bins, in minutes                                   |
| `Nclasses_masc`, `MASC_classes`, `MASC_classes_desired` | `6`, `[SP, CC, PC, AG, GR, CPC]`, `[1..6]` | class scheme                              |
| `N_MascSamples_min`              | `0`     | minimum number of particles for a bin to be drawn                                        |
| `use_triplet`                    | `false` | merge the views of each particle before aggregating (`null` → `processing.use_triplet_algo`; `false` matches MATLAB GT quicklooks) |
| `OR_180`                         | `true`  | fold orientations into [-90°, 90°]                                                        |
| `savefigs`, `overview_classif`, `overview_microstruct`, `disp_now` | `true`, `true`, `true`, `false` | which figures to produce / save / show |
| `qualities`                      | `["GOOD"]` | which ROI buckets to read (`["GOOD"]` or `["GOOD", "BAD"]`)                           |

## `blowing_snow`

| Key          | Default  | Description                                                                                |
|--------------|----------|--------------------------------------------------------------------------------------------|
| `in_path`    | `"auto"` | raw image folder; `"auto"` → `label.campaigndir`                                            |
| `model_path` | `"auto"` | GMM `.mat`; `"auto"` → `blowing_snow/models/gmfit_random_all.mat`                           |
| `davos`      | `"auto"` | sky mask selection; `"auto"` → Davos mask if `"davos"` appears in `in_path`, else default   |
| `block_size` | `15`     | number of consecutive images used for the median sky estimate                              |
| `s`          | `2`      | gap (hours) above which the image sequence is split into separate events                   |
| `w`          | `10`     | half-window (images) for the photo-frequency descriptor                                    |
| `reprocess`  | `true`   | `false` → reuse `<outdir>/features_bs.joblib` if present instead of re-extracting features |

---

## Typical adjustments

- **Only produce ROI files** (no post-processing): set every key of `pipeline` to `false`.
- **Re-label an existing run** with other classifiers or another `blurry_threshold`: edit
  `classification` and run `scripts/runner_classify.py`; no reprocessing needed.
- **Change quicklook binning or plots**: edit `quicklooks` and run
  `scripts/runner_quicklooks.py`.
- **Several campaigns into one `outdir`**: set `clean_outdir: false` and give `outdir` an
  explicit path.
- **Headless server**: keep `display_figs` and `quicklooks.disp_now` at `false` (figures are
  rendered with the Agg backend).
