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

dir_path = "/media/ghiggi/New Volume/Data/MASCDB"

#-----------------------------------------------------------------------------.
### Plot images / triplets 
fpath = os.path.join(dir_path,"MASCdb.zarr")
ds = xr.open_zarr(fpath)

# Stack CAM_ID and Triplet_ID into Image_ID
ds_stacked = ds.stack(image_ID = ('TripletID','CAM_ID'))
ds_stacked.isel(image_ID = slice(0,10)) # subset image_ID (by position only !)
ds_stacked.sel(CAM_ID=0)                # subset CAM_ID   (by label only !)  (This loose CAM_ID info!)
ds_stacked.sel(TripletID=slice(0,10))   # subset TripletID (by label only !) (This loose TripletID info!)

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

# TODO: get_feature_list (without predicted classes, datetime, flake_id) 
# TODO: get_classes_list (riming_*, label_*, melting_*) 
# TODO: which rows have been used for manual classification? 

# TODO:
# get_label_name_dict(name: id) .. and viceversa --> Done in aux.py 
# get_label_id_dict(id: name)
# get_riming_id_dict(id: name)
# get_riming_name_dict(name: id) TODO: riming_name ? not provided !


np.unique(cam0['label_id'])
np.unique(cam0['label_name'])  

np.unique(cam0['riming_id'])        # predicted? 
np.unique(cam0['riming_deg_level']) # what is that?    TODO: riming_name column !!! 

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

np.unique(db['bs_precip_type']) # ! Do not match with https://github.com/jacgraz/pymascdb/blob/master/reader/database_reader.py#L57 
                                # What ['', means? --> JGR answer: indeed they do not match (the function in reader is older)
                                # '' means that for any reason the blowing snow estimation is not available. Not a good choice. Better None?
## Triplet_ID 
'flake_id'

### Predicted variables (average of the three cam?)
'riming_deg_level', # what is that? 
'riming_id',        # create riming_name ? or riming_degree?
'melting_id',
'melting_prob',  
'label_name',      # --> snowflake_label ?
'label_id',
'label_id_prob',

# Other descriptors (single estimate based on 3 images?)
'fallspeed' # not in cam DB 

# TODO: 
# get_3d_gan_variables() 
'3dgan_mass',
'3dgan_vol_ch',
'3dgan_r_g',

# get_env_variables() 
'env_T',
'env_P',
'env_DD',
'env_FF',
'env_RH'

# get_blowing_snow_variables()  
'bs_nor_angle',
'bs_mix_ind',
'bs_precip_type',

# get_data_info 
'datetime',
'campaign',
'latitude',
'longitude',
'altitude',
 
# This also present in cam DB  ... same? average? else? 
'pix_size',
'Xhi',
'n_roi', 
'Dmax'    

### TODO: check images are not all 0 !!!!

# Report somewhere 
# 0 : no values 
# > 0 : some snowflake 
# max 255 

####################
### TODO CHECKS ####
####################
### 1
# Retrieve the event leading to largest Dmax 
timedelta_thr = np.timedelta64(2, 'h')
mascdb.define_event_id(timedelta_thr=timedelta_thr)
mascdb.arrange('cam0.Dmax', decreasing=True)        # TODO CHECK WHY THIS RAISE WARNING

event_id = mascdb.arrange('cam0.Dmax', decreasing=True).cam0['event_id']
idx_event = mascdb.cam0['event_id'] == event_id[0]
mascdb_event = mascdb.isel(idx_event).arrange('cam0.datetime', decreasing=False)
print(mascdb_event)

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

#-------------------------------------------------------------------------------.
### Create metadata for custom zarr store 
fpath = os.path.join(dir_path, "MASCdb_orig","MASCdb.zarr")
arr = zarr.open_array(fpath, group=None)  # 1.69 TiB -->  (589416, 1024, 1024, 3) 
arr.chunks
arr.compressor

arr1 = arr[0:10,:,:,:]

da = xr.DataArray(data=arr1,
                  dims=["TripletID","y","x","CAM_ID"],
                  coords=dict(CAM_ID=[0,1,2], TripletID=np.arange(start, end)))


ds = da.to_dataset(name='campaign_name')
ds.to_zarr(os.path.join("/media/ghiggi/New Volume/Data/MASCDB","Campaign.zarr"))

l_stores = glob.glob("*.zarr")
l_ds = [xr.open_zarr(fpath) for fpath in l_stores]
ds =xr.concat(l_ds, dim="triplet_id")
ds.to_zarr()
## And the modify .attr, .metadata according to original zarr 

 