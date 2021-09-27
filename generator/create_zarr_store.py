#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 09:49:54 2021

@author: ghiggi
"""
#-----------------------------------------------------------------------------.
#############################################
### Zarr store metadata/attrs generation  ###
#############################################
##--------------------------------------------.
import numpy as np 
import xarray as xr 
import os
import zarr
import dask
import pandas as pd

dir_path = "/ltenas3/MASC_DB"
zarr_array_fpath = os.path.join(dir_path,"Array.zarr")
cam0_fpath = os.path.join(dir_path, "MASCdb_cam0.parquet")  

# Read cam0
cam0_df = pd.read_parquet(cam0_fpath).convert_dtypes()
flake_ids = cam0_df['flake_id'].values 

# Open zarr array and define dask chunk size 
chunks = (512,1024,1024,3) # 1.5 MB 
dask_arr = dask.array.from_zarr(zarr_array_fpath, chunks=chunks)

# Define global attributes 
global_attr = {"title": "MASCDB" ,
               "version": "MASCDB v1.0",
               "authors": "Jacopo Grazioli, Gionata Ghiggi",
               "contacts": "jacopo.grazioli@epfl.ch; gionata.ghiggi@epfl.ch",
		       "laboratory": "Environmental Remote Sensing Laboratory (LTE)",
               "institution": "EPFL",
		       "references": "Grazioli et al.,2022. MASCDB - XXXXX, Nat. Sci. Data, doi: https://doi.org/XXXX"
               }

# Create a xr.Dataset 
ds = xr.Dataset(
    data_vars=dict(
        data=(["flake_id", "x", "y", "cam_id"], dask_arr),
    ),
    coords=dict(
        cam_id = np.array([0,1,2]),
        flake_id = flake_ids,
    ),
    attrs=global_attr,
)
ds

# ds = ds.isel(flake_id=slice(0,2048))

# Save to zarr 
p = ds.to_zarr("/ltenas3/MASC_DB/new_zarr.zarr", compute=False)

from dask.diagnostics import ProgressBar
with ProgressBar():                        
    _ = p.compute()
 
 