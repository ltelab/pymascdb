"""
Campaign picture list loader for MASC batch processing.

This module reads imgInfo.txt and dataInfo.txt to build a PictureList with
filenames, timestamps, camera IDs, and fallspeed metadata.

Functions for:
- Parsing campaign metadata files into a PictureList container

Translated from upload.m (Christophe Praz 2015) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""


import numpy as np
import logging

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from typing import List

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

@dataclass
class PictureList:
    """Container for picture metadata extracted from imgInfo.txt and dataInfo.txt."""
    
    time_vec: List[List[int]] = field(default_factory=list)  # [[year, month, day, hour, min, sec], ...]
    id: List[int] = field(default_factory=list)
    cam: List[int] = field(default_factory=list)
    files: List[str] = field(default_factory=list)
    time_str: List[str] = field(default_factory=list)
    time_num: List[float] = field(default_factory=list)
    fallspeed: List[float] = field(default_factory=list)
    fallid: List[int] = field(default_factory=list)
    id_unique: List[int] = field(default_factory=list)  # Unique picture IDs (for triplet processing)


def _is_int(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def upload(dirname: Path) -> PictureList:
    """
    Load MASC image metadata from imgInfo.txt and dataInfo.txt.

    Args:
        dirname: Directory containing imgInfo.txt and dataInfo.txt.

    Returns:
        PictureList with files, id, cam, time_vec, time_str, time_num,
        fallspeed, fallid, and id_unique fields.

    Notes:
        - Fallspeed values are NaN when dataInfo.txt is missing or incomplete.
        - imgInfo.txt format: ID Cam Date Time Filename (whitespace-separated).
    """
    
    dirname = Path(dirname)
    pic = PictureList()
    
    # ===== Import imagelist and time from imgInfo.txt =====
    
    imginfo_path = dirname / "imgInfo.txt"
    
    try:
        with open(imginfo_path, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Split by whitespace (tabs or spaces)
            parts = line.split()
            
            # Skip the optional header line ("flake ID  camera ID  date ...")
            # and any malformed line: a data line starts with an integer ID.
            if len(parts) < 5 or not _is_int(parts[0]):
                continue
            
            # Parse fields: ID  Cam  Date  Time  Filename
            pic_id = int(parts[0])
            pic_cam = int(parts[1])
            day = parts[2]        # e.g., "06.20.2015"
            timestamp = parts[3]  # e.g., "09:00:01"
            filename = parts[4]   # e.g., "image_001_0.png"
            # il y'a un 0 en plus mais aucune idée de à quoi il correspond
            
            # Parse date: MM.DD.YYYY
            d_parts = day.split('.')
            picmo = int(d_parts[0])
            picdd = int(d_parts[1])
            picyr = int(d_parts[2])
            
            # Parse time: HH:MM:SS
            t_parts = timestamp.split(':')
            pichh = int(t_parts[0])
            picmm = int(t_parts[1])
            picss = int(float(t_parts[2]))  # Handle fractional seconds
            
            # Store data
            pic.id.append(pic_id)
            pic.cam.append(pic_cam)
            pic.files.append(filename)
            pic.time_vec.append([picyr, picmo, picdd, pichh, picmm, picss])
        
    except FileNotFoundError:
        logger.error("ERROR : imgInfo.txt not found in %s", dirname)
        return pic
    except Exception as e:
        logger.error("ERROR : unable to parse imgInfo.txt in %s", dirname)
        logger.error("Details: %s", e)
        return pic
    
    # Convert time_vec to datetime strings and numerical timestamps
    pic.time_str = []
    pic.time_num = []
    
    for tv in pic.time_vec:
        dt = datetime(tv[0], tv[1], tv[2], tv[3], tv[4], tv[5])
        pic.time_str.append(dt.strftime('%Y-%m-%d %H:%M:%S'))
        # Convert to timestamp (seconds since epoch) for MATLAB datenum compatibility
        pic.time_num.append(dt.timestamp())
    
    # ===== Import fallspeed data from dataInfo.txt =====
    
    datainfo_path = dirname / "dataInfo.txt"
    
    try:
        with open(datainfo_path, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Split by whitespace
            parts = line.split()
            
            # Skip the optional header line ("snowflake id  date  time  fall speed")
            if len(parts) < 4 or not _is_int(parts[0]):
                continue
            
            # Parse fields: ID  Date  Time  Fallspeed
            fallid = int(parts[0])
            # parts[1] is date, parts[2] is time - we can ignore these
            fallspeed = float(parts[3])
            
            pic.fallid.append(fallid)
            pic.fallspeed.append(fallspeed)
    
    except Exception as e:
        # If dataInfo.txt is missing or unreadable, fill with NaN
        pic.fallid = sorted(set(pic.id))  # unique IDs
        pic.fallspeed = [np.nan] * len(pic.fallid)
        logger.error("ERROR : unable to open/parse dataInfo.txt in %s", dirname)
        logger.error("Details: %s", e)
    
    # ===== Compute unique picture IDs (for triplet processing) =====
    pic.id_unique = sorted(set(pic.id))
    
    return pic


def lookup_fallspeed(pic_list: PictureList, flake_id: int) -> float:
    """
    Return fallspeed for a snowflake id from dataInfo.txt (fallid/fallspeed).

    Args:
        pic_list: PictureList from upload(); fallspeed is parallel to fallid,
            not to the per-image id/files lists.
        flake_id: Snowflake id (same as imgInfo / dataInfo id).

    Returns:
        Fallspeed in m/s, or NaN if the id is missing.
    """
    for fid, speed in zip(pic_list.fallid, pic_list.fallspeed):
        if fid == flake_id:
            return float(speed)
    return float(np.nan)
