#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 15:41:45 2021

@author: ghiggi
"""
import numpy as np 
import xarray as xr 
import os
import zarr
import pandas as pd

dir_path = "/data/MASC_DB"

#-----------------------------------------------------------------------------.
### Plot images / triplets 
fpath = os.path.join(dir_path,"MASCdb.zarr")
ds = xr.open_zarr(fpath)

# Stack CAM_ID and Triplet_ID into Image_ID
ds_stacked = ds.stack(image_ID = ('flake_id','cam_id'))
ds_stacked.isel(image_ID = slice(0,10)) # subset image_ID (by position only !)
ds_stacked.sel(cam_id=0)                # subset CAM_ID   (by label only !)  (This loose CAM_ID info!)
ds_stacked.sel(flake_id=slice(0,10))   # subset TripletID (by label only !) (This loose TripletID info!)

# Unstack Image_ID
ds_unstacked = ds_stacked.unstack('image_ID') 
ds_unstacked['data'].dims 
ds['data'].dims  # ! Dimension order different than original ds ! 

#-----------------------------------------------------------------------------.
### CAM DB 
feature_cam0_fpath = os.path.join(dir_path, "MASCdb_cam0.parquet")
feature_cam1_fpath = os.path.join(dir_path, "MASCdb_cam1.parquet")
feature_cam2_fpath = os.path.join(dir_path, "MASCdb_cam2.parquet")

cam0 = pd.read_parquet(feature_cam0_fpath)
cam1 = pd.read_parquet(feature_cam1_fpath)
cam2 = pd.read_parquet(feature_cam2_fpath)

# First row: flake_id start at 10 ... then 11,12,13,14 ... then 1 !
print(cam0)
cam0['flake_id'].values[0:40]
np.unique(cam0['flake_id'])

# Select column 
cam0['label_id']
cam0['label_name'] 

cam0.columns.values.tolist()  # TODO: get_columns_list()

np.unique(cam0['label_id'])
np.unique(cam0['label_name'])  

np.unique(cam0['riming_id'])        # predicted?  
np.unique(cam0['riming_deg_level']) # what is that?   

# Features in xarray Dataset
cam0_ds = xr.Dataset.from_dataframe(cam0)
cam0_ds = cam0_ds.set_coords('flake_id').swap_dims({'index': 'flake_id'}) # 

cam1_ds = xr.Dataset.from_dataframe(cam1)
cam1_ds = cam1_ds.set_coords('flake_id').swap_dims({'index': 'flake_id'})

cam2_ds = xr.Dataset.from_dataframe(cam2)
cam2_ds = cam2_ds.set_coords('flake_id').swap_dims({'index': 'flake_id'})

#-----------------------------------------------------------------------------.
# Triplet DB
triplet_fpath = os.path.join(dir_path, "MASCdb_triplet.parquet")
db = pd.read_parquet(triplet_fpath)
db.shape
db.columns.values.tolist()  

#-----------------------------------------------------------------------------.
#################
## Filtering ####
#################  
# filter(marketing.AmountSpent > 2000) & (marketing.History == 'High') # now 
# filter(AmountSpent > 2000) & (History == 'High')) # ideal 
    
# https://github.com/kieferk/dfply
# https://github.com/dodger487/dplython
# https://pythonhosted.org/dplython/
# https://github.com/coursera/pandas-ply
# https://github.com/koaning/kadro
# dplython:sift --> filter        #  pandas-ply.ply_where(X.arr > 30, X.dep > 30)) #

#-----------------------------------------------------------------------------.
#############################################
### Zarr store metadata/attrs generation  ###
#############################################
##--------------------------------------------.
import matplotlib.pyplot as plt
fpath='/data/MASC_DB/MASCdb.zarr/data/'
#-------------------------------------------------------------------------------.
### Create metadata for custom zarr store 
fpath = os.path.join(dir_path, "MASCdb_orig","MASCdb.zarr")
arr = zarr.open_array(fpath, group=None)  # 1.69 TiB -->  (589416, 1024, 1024, 3) 
arr.chunks
arr.compressor

arr1 = arr[0:10,:,:,:]

da = xr.DataArray(data=arr1,
                  dims=["flake_id","y","x","cam_id"],
                  coords=dict(cam_id=[0,1,2]))

ds = da.to_dataset(name='data')
ds.to_zarr("/data/MASC_DB/zarr_metadata.zarr")
    
    
ds.to_zarr(os.path.join("/media/ghiggi/New Volume/Data/MASCDB","Campaign.zarr"))



l_stores = glob.glob("*.zarr")
l_ds = [xr.open_zarr(fpath) for fpath in l_stores]
ds =xr.concat(l_ds, dim="triplet_id")
ds.to_zarr()
## And the modify .attr, .metadata according to original zarr 

 