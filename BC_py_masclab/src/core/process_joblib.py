"""
Core batch processing for MASC data with joblib parallelism.

This module discovers campaign folders, runs single-image or triplet processing,
and writes timing and quality statistics.

Functions for:
- Parallel single-image and triplet batch processing
- Processing statistics and timing accumulation

Translated from MASC_process.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import os
import time
import logging

from pathlib import Path
from typing import List, Tuple, Optional, Dict
from joblib import Parallel, delayed
from tqdm import tqdm

from src.core.config import ProcessingConfig, LabelConfig
from src.dataio.uploaddirs import uploaddirs
from src.dataio.upload import upload, lookup_fallspeed
from src.processing.single import process_single_image
from src.processing.triplet import process_triplet_image

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)  

class ProcessingStats:
    """Statistics accumulated during batch processing."""
    def __init__(self, use_triplet_algo):
        self.N_tot: int = 0

        if use_triplet_algo:
            # Triplet-specific stats
            self.N_triplet_full: int = 0
            self.N_triplet_miss: int = 0
            self.N_matched: int = 0
            self.N_no_match: int = 0
            self.N_processed_indep: int = 0

        else:
            self.N_good: int = 0
            self.N_bad: int = 0
            self.N_blurry: int = 0


class Timings:
    """Timing statistics for various processing steps."""
    def __init__(self, use_triplet_algo):

        self.tt_uploading: float = 0.0
        self.tt_loading: float = 0.0
        self.tt_clutter: float = 0.0
        self.tt_edging: float = 0.0
        self.tt_roiying: float = 0.0
        self.tt_feature: float = 0.0
        self.tt_plotting: float = 0.0
        self.tt_saving: float = 0.0
        self.tt_program: float = 0.0

        if use_triplet_algo:
            self.tt_matching: float = 0.0


def process_single_image_wrapper(
    img_path: Path,
    pic_info: Dict,
    label: LabelConfig,
    process: ProcessingConfig,
) -> Tuple[Optional[Dict], int, Dict]:
    """
    Joblib wrapper around process_single_image with exception handling.

    Args:
        img_path: Path to the image file.
        pic_info: Metadata dict (filename, cam, id, fallspeed, time_num).
        label: Label configuration with output paths.
        process: Processing configuration.

    Returns:
        Tuple of (roi, flag, timing) passed through from process_single_image,
        or (None, -1, {}) on error.
    """
    try:
        roi, flag, timing = process_single_image(img_path, pic_info, label, process)
        return roi, flag, timing
    except Exception as exc:
        logger.warning("Error processing %s: %s", img_path, exc)
        return None, -1, {}


def process_triplet_image_wrapper(
    flake_id: int,
    idx_pics: List[int],
    current_dir: Path,
    pic_list,
    label: LabelConfig,
    process: ProcessingConfig
) -> Tuple[Dict, int, int]:
    """
    Joblib wrapper around process_triplet_image with exception handling and
    MATLAB-like fallback: if no triplet could be matched (flag < 2), each
    image of the event is reprocessed independently through
    process_single_image (MASC_process.m:179-201).

    Args:
        flake_id: Snowflake event identifier.
        idx_pics: Indices into pic_list for this flake's images.
        current_dir: Directory containing the image files.
        pic_list: Campaign picture list with metadata.
        label: Label configuration with output paths.
        process: Processing configuration.

    Returns:
        Tuple of (timing, flag, n_indep_good):
            - timing: accumulated step timings (triplet + fallback singles).
            - flag: triplet flag (-1=error, 0=no ROI, 1=no match, 2=matched).
            - n_indep_good: number of fallback images processed with flag 2.
    """
    try:
        timing, flag = process_triplet_image(flake_id, idx_pics, current_dir, pic_list, label, process)
    except Exception as exc:
        logger.warning("Error processing triplet %s: %s", flake_id, exc)
        timing, flag = {}, -1

    n_indep_good = 0
    if flag < 2:
        for idx in idx_pics:
            pic_info = {
                "filename": pic_list.files[idx] if pic_list.files else f"image_{idx}.png",
                "cam": pic_list.cam[idx],
                "id": pic_list.id[idx],
                "fallspeed": lookup_fallspeed(pic_list, pic_list.id[idx]),
                "time_num": pic_list.time_num[idx] if idx < len(pic_list.time_num) else None,
            }
            img_path = Path(current_dir) / pic_info["filename"]
            try:
                _, s_flag, s_timing = process_single_image(img_path, pic_info, label, process)
            except Exception as exc:
                logger.warning("Error processing fallback image %s: %s", pic_info["filename"], exc)
                s_flag, s_timing = -1, {}
            if s_flag == 2:
                n_indep_good += 1
            for key, value in (s_timing or {}).items():
                timing[key] = timing.get(key, 0.0) + value

    return timing, flag, n_indep_good


def masc_process(label: LabelConfig, process: ProcessingConfig) -> None:
    """
    Run the main MASC batch processing pipeline with optional joblib parallelism.

    Args:
        label: Label configuration (campaign paths, date range, output directory).
        process: Processing configuration (thresholds, algorithms, save options).

    Returns:
        None

    Notes:
        - Camera parameters are integrated into ProcessingConfig in Python.
        - Parallelism is controlled via process.parallel and process.max_workers.
        - Triplet matching is controlled via process.use_triplet_algo.

    Translated from MASC_process.m (Christophe Praz 2015) and adapted for Python
    """
    
    # Find all relevant snowflake directories
    try:
        dir_list = uploaddirs(label.campaigndir, label.starthr_vec, label.endhr_vec)
    except NotImplementedError:
        logger.error("You need to implement directory discovery based on campaign structure")
        return
    
    # Initialize statistics and timings
    stats = ProcessingStats(process.use_triplet_algo)
    timings = Timings(process.use_triplet_algo)
    
    # Start program timer
    t_startprogram = time.time()
    
    # ===== Setup joblib parallelism if enabled =====
    if process.parallel:
        # Set environment variables to avoid nested parallelism
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["BLAS_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        
        # Determine number of workers
        n_workers = min(process.max_workers, os.cpu_count() - 2 if os.cpu_count() > 2 else 1)
        print(f"Joblib parallelism enabled with {n_workers} workers")
    
    # Loop over all relevant directories
    for i_dir in range(len(dir_list)):
        current_dir = dir_list[i_dir]
        print(f"\n***** Processing directory {i_dir+1}/{len(dir_list)}: {current_dir}")
        
        # Start upload timer for this directory
        t_uploading = time.time()
        
        # Retrieve snowflakes info
        pic_list = upload(current_dir)
        # id_unique is now computed inside upload()

        # Number of new flakes added
        n_id = len(pic_list.id_unique)
        if not process.parallel:
            stats.N_tot += n_id
        
        timings.tt_uploading += time.time() - t_uploading
        
        # ===== Prepare processing tasks =====
        # Collect all tasks and metadata
        all_tasks = []
        task_metadata = []
        
        if not process.use_triplet_algo:
            # ===== Single-image mode (no triplet matching) =====
            
            for j in range(len(pic_list.id)):
                # Construct pic_info dict for this image
                pic_info = {
                    "filename": pic_list.files[j] if pic_list.files else f"image_{j}.png",
                    "cam": pic_list.cam[j],
                    "id": pic_list.id[j],
                    "fallspeed": lookup_fallspeed(pic_list, pic_list.id[j]),
                    "time_num": pic_list.time_num[j] if j < len(pic_list.time_num) else None,
                }
                
                img_path = Path(current_dir) / pic_list.files[j]
                
                if process.parallel:
                    # Create delayed task for joblib
                    task = delayed(process_single_image_wrapper)(img_path, pic_info, label, process)
                    all_tasks.append(task)
                    task_metadata.append({"type": "single", "filename": pic_info["filename"]})
                else:
                    # Sequential processing
                    if pic_list.files:
                        logger.info("Processing %s...", pic_list.files[j])
                    
                    try:
                        roi, flag, timing = process_single_image(img_path, pic_info, label, process)
                    except Exception as exc:
                        logger.error("Unknown error while processing image %s: %s", pic_list.files[j], exc)
                        flag = -1
                        timing = {}
                    
                    # Accumulate statistics
                    _accumulate_single_stats(stats, flag)
                    _accumulate_timings(timings, timing, process.use_triplet_algo)
        
        else:
            # ===== Triplet matching mode =====
            
            if not process.parallel:
                logger.info("Processing %d unique snowflake IDs with triplet matching...", len(pic_list.id_unique))
            
            # Loop over unique picture IDs (triplets or incomplete triplets)
            for j, unique_id in enumerate(pic_list.id_unique):
                # Find all images with this ID
                idx_pics = [i for i, pic_id in enumerate(pic_list.id) if pic_id == unique_id]
                min_views = max(2, getattr(process, "triplet_required_views", 3))
                allow_partial = getattr(process, "triplet_allow_partial", True)
                should_run_triplet = len(idx_pics) >= min_views or (allow_partial and len(idx_pics) >= 2)
                
                # Run multiview matching if enough views are available
                if should_run_triplet:
                    if not process.parallel:
                        stats.N_triplet_full += 1
                    
                    if process.parallel:
                        # Create delayed task for triplet
                        task = delayed(process_triplet_image_wrapper)(
                            unique_id, idx_pics, current_dir, pic_list, label, process
                        )
                        all_tasks.append(task)
                        task_metadata.append({"type": "triplet", "id": unique_id, "n_images": len(idx_pics)})
                    else:
                        # Sequential processing (wrapper includes the MATLAB
                        # fallback to independent processing when flag < 2)
                        timing, flag, n_indep_good = process_triplet_image_wrapper(
                            unique_id, idx_pics, current_dir, pic_list, label, process
                        )
                        
                        # Accumulate statistics
                        _accumulate_triplet_stats(stats, flag, n_indep_good)
                        _accumulate_timings(timings, timing, process.use_triplet_algo)
                
                else:
                    # Not enough views for multiview matching - process independently
                    if not process.parallel:
                        stats.N_triplet_miss += 1
                        stats.N_processed_indep += len(idx_pics)
                    
                    for idx in idx_pics:
                        pic_info = {
                            "filename": pic_list.files[idx] if pic_list.files else f"image_{idx}.png",
                            "cam": pic_list.cam[idx],
                            "id": pic_list.id[idx],
                            "fallspeed": lookup_fallspeed(pic_list, pic_list.id[idx]),
                            "time_num": pic_list.time_num[idx] if idx < len(pic_list.time_num) else None,
                        }
                        
                        img_path = Path(current_dir) / pic_list.files[idx]
                        
                        if process.parallel:
                            task = delayed(process_single_image_wrapper)(img_path, pic_info, label, process)
                            all_tasks.append(task)
                            task_metadata.append({"type": "incomplete_triplet", "filename": pic_info["filename"]})
                        else:
                            try:
                                roi, flag, timing = process_single_image(img_path, pic_info, label, process)
                            except Exception as exc:
                                logger.error("Unknown error while processing image %s: %s", pic_list.files[idx], exc)
                                flag = -1
                                timing = {}
                            
                            # Accumulate statistics
                            _accumulate_single_stats(stats, flag)
                            _accumulate_timings(timings, timing, process.use_triplet_algo)
        
        # ===== Execute parallel tasks if enabled =====
        
        if process.parallel and all_tasks:
            print(f"\nProcessing {len(all_tasks)} tasks with joblib ({n_workers} workers)...")
            
            # Execute all tasks in parallel with joblib with progress bar
            # backend='loky' is the default and uses process-based parallelism
            results = Parallel(n_jobs=n_workers, backend='loky', verbose=0)(
                tqdm(all_tasks, desc="Processing images", unit="image")
            )
            
            # Accumulate statistics from results
            for result, metadata in zip(results, task_metadata):
                if metadata["type"] == "single":
                    roi, flag, timing = result
                    _accumulate_single_stats(stats, flag)
                    _accumulate_timings(timings, timing, process.use_triplet_algo)
                elif metadata["type"] == "triplet":
                    timing, flag, n_indep_good = result
                    stats.N_triplet_full += 1
                    _accumulate_triplet_stats(stats, flag, n_indep_good)
                    _accumulate_timings(timings, timing, process.use_triplet_algo)
                elif metadata["type"] == "incomplete_triplet":
                    roi, flag, timing = result
                    stats.N_triplet_miss += 1
                    stats.N_processed_indep += 1
                    # Note: Don't call _accumulate_single_stats here since triplet mode
                    # ProcessingStats doesn't have N_good/N_bad/N_blurry attributes
                    _accumulate_timings(timings, timing, process.use_triplet_algo)
            
            print(f"Completed {len(all_tasks)} tasks")
    
    # ===== End of processing loop =====
    
    # Compute N_tot for parallel mode
    if process.parallel:
        if process.use_triplet_algo:
            stats.N_tot = stats.N_triplet_full + stats.N_triplet_miss
        else:
            stats.N_tot = stats.N_good + stats.N_bad + stats.N_blurry
    
    timings.tt_program = time.time() - t_startprogram
    
    # ===== Print summary statistics =====
    
    print("*********************************************************")
    print(f"***** Task finished ! Total Time spent    : {timings.tt_program:.2f} seconds")
    
    if process.use_triplet_algo:
        print(f"***** Number of triplets processed        : {stats.N_tot:d}")
        print(f"***** Number of triplet with 3+ views     : {stats.N_triplet_full:d} {stats.N_triplet_full/stats.N_tot*100:.1f} %")
        print(f"***** Number of missing triplets          : {stats.N_triplet_miss:d} {stats.N_triplet_miss/stats.N_tot*100:.1f} %")
        if stats.N_triplet_full > 0:
            print(f"***** Number of matched triplets          : {stats.N_matched:d} {stats.N_matched/stats.N_triplet_full*100:.1f} %")
            print(f"***** Number of triplets without match    : {stats.N_no_match:d} {stats.N_no_match/stats.N_triplet_full*100:.1f} %")
        print(f"***** Number of imgs processed indep.     : {stats.N_processed_indep:d} {stats.N_processed_indep/(stats.N_processed_indep+3*stats.N_triplet_full)*100:.1f} %")
    else:
        N_tot_im = stats.N_good + stats.N_bad + stats.N_blurry
        print(f"***** Number of flakes found              : {stats.N_tot:d}")
        print(f"***** Number of pictures processed        : {N_tot_im:d}")
        if N_tot_im > 0:
            print(f"***** Number of good snowflakes           : {stats.N_good:d} {stats.N_good/N_tot_im*100:.1f} %")
            print(f"***** Number of blurry snowflakes         : {stats.N_blurry:d} {stats.N_blurry/N_tot_im*100:.1f} %")
            print(f"***** Number of no/bad detections         : {stats.N_bad:d} {stats.N_bad/N_tot_im*100:.1f} %")
    
    tt_all = (timings.tt_uploading + timings.tt_loading + timings.tt_clutter + 
              timings.tt_edging + timings.tt_roiying + timings.tt_feature + 
              timings.tt_plotting + timings.tt_saving)
    if process.use_triplet_algo:
        tt_all += timings.tt_matching
    
    if tt_all > 0:
        print(f"***** Total time uploading directories    : {timings.tt_uploading/tt_all*100:.1f} %")
        print(f"***** Total time loading pictures         : {timings.tt_loading/tt_all*100:.1f} %")
        print(f"***** Total time removing clutter         : {timings.tt_clutter/tt_all*100:.1f} %")
        print(f"***** Total time edging pictures          : {timings.tt_edging/tt_all*100:.1f} %")
        print(f"***** Total time computing feature        : {timings.tt_feature/tt_all*100:.1f} %")
        print(f"***** Total time selecting flake ROI      : {timings.tt_roiying/tt_all*100:.1f} %")
        if process.use_triplet_algo:
            print(f"***** Total time matching particules      : {timings.tt_matching/tt_all*100:.1f} %")

        print(f"***** Total time plotting                 : {timings.tt_plotting/tt_all*100:.1f} %")
        print(f"***** Total time saving processed data    : {timings.tt_saving/tt_all*100:.1f} %")
    
    if process.use_triplet_algo and stats.N_tot > 0:
        print(f"***** Net average time per triplet         : {timings.tt_program/stats.N_tot:.2f}")
    elif not process.use_triplet_algo and N_tot_im > 0:
        print(f"***** Net average time per flake           : {timings.tt_program/N_tot_im:.2f}")
    
    print("*********************************************************\n")
    
    # ===== Save statistics to file =====
    
    outdir = Path(label.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stats_file = outdir / "proc_stats.txt"
    with open(stats_file, "w") as f:
        f.write("*********************************************************\n")
        f.write(f"***** Task finished ! Total Time spent    : {timings.tt_program:.2f} seconds\n")
        
        if process.use_triplet_algo:
            f.write(f"***** Number of triplets processed        : {stats.N_tot:d}\n")
            f.write(f"***** Number of triplet with 3+ views     : {stats.N_triplet_full:d} {stats.N_triplet_full/stats.N_tot*100:.1f} %\n")
            f.write(f"***** Number of missing triplets          : {stats.N_triplet_miss:d} {stats.N_triplet_miss/stats.N_tot*100:.1f} %\n")
            if stats.N_triplet_full > 0:
                f.write(f"***** Number of matched triplets          : {stats.N_matched:d} {stats.N_matched/stats.N_triplet_full*100:.1f} %\n")
                f.write(f"***** Number of triplets without match    : {stats.N_no_match:d} {stats.N_no_match/stats.N_triplet_full*100:.1f} %\n")
            f.write(f"***** Number of imgs processed indep.     : {stats.N_processed_indep:d} {stats.N_processed_indep/(stats.N_processed_indep+3*stats.N_triplet_full)*100:.1f} %\n")
        else:
            N_tot_im = stats.N_good + stats.N_bad + stats.N_blurry
            f.write(f"***** Number of flakes found              : {stats.N_tot:d}\n")
            f.write(f"***** Number of pictures processed        : {N_tot_im:d}\n")
            if N_tot_im > 0:
                f.write(f"***** Number of good snowflakes           : {stats.N_good:d} {stats.N_good/N_tot_im*100:.1f} %\n")
                f.write(f"***** Number of blurry snowflakes         : {stats.N_blurry:d} {stats.N_blurry/N_tot_im*100:.1f} %\n")
                f.write(f"***** Number of no/bad detections         : {stats.N_bad:d} {stats.N_bad/N_tot_im*100:.1f} %\n")
        
        if tt_all > 0:
            f.write(f"***** Total time uploading directories    : {timings.tt_uploading/tt_all*100:.1f} %\n")
            f.write(f"***** Total time loading pictures         : {timings.tt_loading/tt_all*100:.1f} %\n")
            f.write(f"***** Total time removing clutter         : {timings.tt_clutter/tt_all*100:.1f} %\n")
            f.write(f"***** Total time edging pictures          : {timings.tt_edging/tt_all*100:.1f} %\n")
            f.write(f"***** Total time computing feature        : {timings.tt_feature/tt_all*100:.1f} %\n")
            f.write(f"***** Total time selecting flake ROI      : {timings.tt_roiying/tt_all*100:.1f} %\n")
            if process.use_triplet_algo:
                f.write(f"***** Total time matching particules      : {timings.tt_matching/tt_all*100:.1f} %\n")
            f.write(f"***** Total time plotting                 : {timings.tt_plotting/tt_all*100:.1f} %\n")
            f.write(f"***** Total time saving processed data    : {timings.tt_saving/tt_all*100:.1f} %\n")
        
        if process.use_triplet_algo and stats.N_tot > 0:
            f.write(f"***** Net average time per triplet         : {timings.tt_program/stats.N_tot:.2f}\n")
        elif not process.use_triplet_algo and N_tot_im > 0:
            f.write(f"***** Net average time per flake           : {timings.tt_program/N_tot_im:.2f}\n")
        
        f.write("*********************************************************\n")


# ===== Helper functions for statistics accumulation =====

def _accumulate_single_stats(stats: ProcessingStats, flag: int) -> None:
    """Accumulate statistics for single image processing."""
    # Canonical quality flags:
    # -1=error, 0=no ROI or bad detection, 1=blurry, 2=good
    # Note: single.py currently emits {-1, 0, 2}; blurry (1) is reserved.
    if flag == 0:
        stats.N_bad += 1
    elif flag == 1:
        stats.N_blurry += 1
    elif flag == 2:
        stats.N_good += 1


def _accumulate_triplet_stats(stats: ProcessingStats, flag: int, n_indep_good: int) -> None:
    """Accumulate statistics for triplet processing.

    Mirrors MASC_process.m:178-203: flag == 2 counts as a matched triplet;
    any flag < 2 (error, no ROI, no match) counts as no-match, and
    n_indep_good is the number of fallback images successfully processed
    independently (MATLAB increments N_processed_indep only on flag == 2).
    """
    if flag == 2:  # GOOD - matched triplet
        stats.N_matched += 1
    else:  # error / no ROI / no match — fallback path
        stats.N_no_match += 1
        stats.N_processed_indep += n_indep_good


def _accumulate_timings(timings: Timings, timing: Dict, use_triplet_algo: bool) -> None:
    """Accumulate timing statistics."""
    if timing:
        timings.tt_loading += timing.get("loading", 0.0)
        timings.tt_clutter += timing.get("clutter", 0.0)
        timings.tt_edging += timing.get("edging", 0.0)
        timings.tt_roiying += timing.get("roiying", 0.0)
        timings.tt_feature += timing.get("feature", 0.0)
        timings.tt_plotting += timing.get("plotting", 0.0)
        timings.tt_saving += timing.get("saving", 0.0)
        if use_triplet_algo:
            timings.tt_matching += timing.get("matching", 0.0)
