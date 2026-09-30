"""
Vanderbilt-format MASC picture list loader.

This module parses FLAKE_* filenames when imgInfo.txt and dataInfo.txt are
not available (Vanderbilt triplet campaigns).

Functions for:
- Scanning FLAKE_*.png files and building a PictureList

Translated from MASC_VanderbiltTriplets_process.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import logging

from pathlib import Path
from datetime import datetime

from src.utils.dataio.upload import PictureList

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

def upload_vanderbilt(dirname: Path) -> PictureList:
    """
    Build a PictureList by parsing FLAKE_* filenames in a directory.

    Args:
        dirname: Directory containing FLAKE_*.png files.

    Returns:
        PictureList with files, id, cam, time fields; fallspeed filled with NaN.

    Notes:
        - Images sharing the same timestamp receive the same flake ID (triplets).
        - Filename format: FLAKE_YYYYMMDD_HHMMSS_camX.png.

    Translated from MASC_VanderbiltTriplets_process.m (Christophe Praz 2015) and adapted for Python
    """
    
    dirname = Path(dirname)
    pic = PictureList()
    
    # ===== Find all FLAKE_* files =====
    
    try:
        # Get all FLAKE_*.png files and sort them
        flake_files = sorted(dirname.glob("FLAKE*.png"))
        
        if not flake_files:
            logger.warning("No FLAKE*.png files found in %s", dirname)
            return pic
        
        print(f"Found {len(flake_files)} Vanderbilt FLAKE images in {dirname}")
        
    except Exception as e:
        logger.error("ERROR: Unable to read directory %s", dirname)
        logger.error("Details: %s", e)
        return pic
    
    # ===== Parse each filename to extract metadata =====
    
    current_id = 0
    previous_timestamp = None
    
    for file_path in flake_files:
        filename = file_path.name
        
        try:
            # Parse filename: FLAKE_YYYYMMDD_HHMMSS_camX.png
            #                 012345678901234567890123456789
            #                       6  13 15  20    -4
            
            # Extract date string (positions 6-13): YYYYMMDD
            date_str = filename[6:14]  # 8 characters
            
            # Extract time string (positions 15-21): HHMMSS
            time_str = filename[15:21]  # 6 characters
            
            # Extract camera ID (position end-5): X in camX.png
            # filename ends with "camX.png" so cam ID is at position -5
            cam_id = int(filename[-5])
            
            # Combine date and time for timestamp comparison
            timestamp = date_str + time_str
            
            # Parse into datetime components
            year = int(date_str[0:4])
            month = int(date_str[4:6])
            day = int(date_str[6:8])
            hour = int(time_str[0:2])
            minute = int(time_str[2:4])
            second = int(time_str[4:6])
            
            # Create datetime object
            dt = datetime(year, month, day, hour, minute, second)
            
            # Assign picture ID based on timestamp
            # Images with identical timestamps belong to the same snowflake (triplet)
            if previous_timestamp is None or timestamp != previous_timestamp:
                # New timestamp = new flake ID
                current_id += 1
            
            # Store data
            pic.files.append(filename)
            pic.cam.append(cam_id)
            pic.id.append(current_id)
            pic.time_vec.append([year, month, day, hour, minute, second])
            pic.time_str.append(timestamp)
            pic.time_num.append(dt.timestamp())
            
            # Update previous timestamp for next iteration
            previous_timestamp = timestamp
            
        except (ValueError, IndexError) as e:
            logger.warning("Could not parse filename '%s': %s", filename, e)
            continue
    
    # ===== Compute unique picture IDs =====
    pic.id_unique = sorted(set(pic.id))
    
    # ===== Fill fallspeed data with NaN (not available for Vanderbilt) =====
    pic.fallid = pic.id_unique.copy()
    pic.fallspeed = [np.nan] * len(pic.fallid)
    
    print(f"Parsed {len(pic.files)} images into {len(pic.id_unique)} unique snowflake IDs")
    
    return pic


def upload_vanderbilt_legacy(dirname: Path) -> PictureList:
    """Legacy alias for upload_vanderbilt (MATLAB naming compatibility)."""
    return upload_vanderbilt(dirname)
