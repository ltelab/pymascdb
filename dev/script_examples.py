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
# get_label_name_dict(name: id) .. and viceversa
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
                                # What ['', means? 
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

#-----------------------------------------------------------------------------.
#################
## Filtering ####
################# get always the idx ... to subset all db and ds 

### Always return a copy of the object 
# arrange (by size i.e. , pixel number, diameter)  (on what cam db?)
# select_first_n 
# select_last_n 
# sample_n 

# select (column) 
# filter (row) 
 
# filter(marketing.AmountSpent > 2000) & (marketing.History == 'High') # now 
# filter(AmountSpent > 2000) & (History == 'High')) # ideal 
    
# https://github.com/kieferk/dfply
# https://github.com/dodger487/dplython
# https://pythonhosted.org/dplython/
# https://github.com/coursera/pandas-ply
# https://github.com/koaning/kadro

# select : db columns 
# dplython:sift --> filter        #  pandas-ply.ply_where(X.arr > 30, X.dep > 30)) #
# sample_n : sample 
# arrange : reorder 

# capture expression ... check cam0, cam1, cam2, triplet.env (env.T > ), 

#-----------------------------------------------------------------------------.
###########################
### Already implemented ###
###########################
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
                  coords=dict(CAM_ID=[0,1,2]))

ds = da.to_dataset(name='data')
ds.to_zarr(os.path.join("/media/ghiggi/New Volume/Data/MASCDB","my_masc.zarr"))

## And the modify .attr, .metadata according to original zarr 

#-----------------------------------------------------------------------------.
### Plot images / triplets 
fpath = os.path.join(dir_path,"MASCdb.zarr")
ds = xr.open_zarr(fpath)

# Select CAM_ID 0 
ds.sel(CAM_ID = 0)   # by label
ds.isel(CAM_ID = 0)  # by index 

# Display 9 images of CAM_ID = 0 [plot_flakes()]
p = ds['data'].sel(CAM_ID = 0, TripletID = slice(0,9)).plot(x='x',y='y', row="TripletID", col_wrap=3, 
                                                            aspect=1, 
                                                            cmap='gray', add_colorbar=False, 
                                                            vmin=0, vmax=255)
for i, ax in enumerate(p.axes.flat):
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_axis_off() 
p.fig.subplots_adjust(wspace=0.01, hspace=0.1)

# add option to report mm ... 
# pix_size = self.triplet.loc[index].pix_size*1e3  # mm
# xsize = out.shape[0]*pix_size # mm
# ysize = out.shape[1]*pix_size # mm


# Display 3 images of all CAM [plot_triplets()]
p = ds['data'].isel(TripletID = slice(0,3)).plot(x='x',y='y',
                                                 col="CAM_ID",
                                                 row="TripletID",
                                                 aspect=1,
                                                 cmap='gray', add_colorbar=False,
                                                 vmin=0, vmax=255)
 
for i, ax in enumerate(p.axes.flat):
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_axis_off()
p.fig.subplots_adjust(wspace=0.01, hspace=0.01)


### zoom function 
# --> should we respect aspect=1 and return square images 
img_triplet = ds['data'].isel(TripletID=0).values
plt.imshow(img_triplet[:,:,0], vmin=0, vmax=255, cmap="gray") 
plt.imshow(img_triplet[:,:,2], vmin=0, vmax=255, cmap="gray") 
plt.colorbar()

img = img_triplet[:,:,0]

def _internal_bbox(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return rmin, rmax, cmin, cmax

def _get_zoomed_image(img):
    rmin, rmax, cmin, cmax = _internal_bbox(img)
    zoom_img = img[rmin:rmax+1, cmin:cmax+1]
    return zoom_img

def _center_image(img, nrow, ncols): 
    r, c = img.shape 
    col_incr = int((ncol - c)/2)
    row_incr = int((nrow - r)/2)
    arr = np.zeros((nrow, ncol))
    arr[slice(row_incr,row_incr+r), slice(col_incr,col_incr+c)] = img 
    return arr 


img = img_triplet[:,:,0]
plt.imshow(_get_zoomed_image(img), vmin=0, vmax=255, cmap="gray") 
l_imgs = [img_triplet[:,:,i] for i in range(3)]

# Zoom a list of image 
l_zoomed = [_get_zoomed_image(img) for img in l_imgs]
# Ensure same shape across all zoomed images 
l_shapes = [img.shape for img in l_zoomed]
r_max, c_max = (max(n) for n in zip(*l_shapes)) # Get max number of row an columns 
l_zoomed = [_center_image(img, nrow=r_max, ncols=c_max) for img in l_zoomed]

for img in l_zoomed:
    plt.imshow(img, vmin=0, vmax=255, cmap="gray")
    plt.show()