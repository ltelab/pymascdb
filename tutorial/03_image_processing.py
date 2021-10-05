#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 11:22:25 2021

@author: ghiggi
"""
##----------------------------------------------------------------------------.
##################################
### Image processing Tutorial ####
##################################
import os
os.chdir("/home/ghiggi/Projects/pymascdb")
#os.chdir("/home/grazioli/CODES/python/pymascdb")
import math
import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt

from dask.distributed import Client
from skimage import filters, measure, morphology
import mascdb.api
from mascdb.api import MASC_DB
from mascdb.utils_img import _compute_2Dimage_descriptors 

dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
#dir_path = "/data/MASC_DB/"

# Initialize dask distributed client (for parallel processing)
client = Client(processes=False)

##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Select the largest        
mascdb = mascdb.arrange('cam0.Dmax', decreasing=True)
img = mascdb.da.isel(cam_id=0, flake_id=0).values

# 0 : no values 
# > 0 : some snowflake 
# max 255 

##----------------------------------------
# Get contours 
plt.imshow(img)
for contour in measure.find_contours(img, 0):
    # coords = contour
    coords = measure.approximate_polygon(contour, tolerance=2.5)
    plt.plot(coords[:, 1], coords[:, 0], '-r', linewidth=2)
plt.show()

##----------------------------------------
### Retrieve some standard image properties 
# Obtain a binary image
threshold = filters.threshold_otsu(img)
binary_mask = img > threshold
plt.imshow(binary_mask, interpolation='None')
plt.show()
# Remove inner holes and other objects 
binary_mask = morphology.remove_small_objects(binary_mask, 50)
plt.imshow(binary_mask, interpolation='None')
plt.show()
binary_mask = morphology.remove_small_holes(binary_mask, 50)
plt.imshow(binary_mask, interpolation='None')
plt.show()
# Labels regions 
# --> Background is label 0 
labels = measure.label(binary_mask)
plt.imshow(labels, interpolation='None')
plt.show()
# Remove spurious stuffs 
labels[labels != 1] = 0  # TODO: check if 1 is always the largest
plt.imshow(labels, interpolation='None')
plt.show()
# Compute properties 
properties = ['area','filled_area', 'perimeter', 'eccentricity', 'equivalent_diameter',
              'euler_number',  'solidity', 'feret_diameter_max',
              'orientation',
              'minor_axis_length', 'major_axis_length'] # 'axis_major_length','axis_minor_length'] 
# 'inertia_tensor',
# 'inertia_tensor_eigvals',
# 'moments',
# 'moments_central',
# 'moments_hu',
# 'moments_normalized',
    
region_props = measure.regionprops(labels)
props_dict = measure.regionprops_table(labels, properties = properties) 
props_dict = {k: v[0] for k,v in props_dict.items()}

props_list = [props_dict[prop] for prop in properties]
# region_props[0].convex_image            
# region_props[0].filled_image

# Plot minor and major axis of the snowflake 
fig, ax = plt.subplots()
ax.imshow(labels, cmap=plt.cm.gray)

for props in region_props:
    y0, x0 = props.centroid
    orientation = props.orientation
    x1 = x0 + math.cos(orientation) * 0.5 * props.minor_axis_length
    y1 = y0 - math.sin(orientation) * 0.5 * props.minor_axis_length
    x2 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
    y2 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length

    ax.plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
    ax.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
    ax.plot(x0, y0, '.g', markersize=15)

    minr, minc, maxr, maxc = props.bbox
    bx = (minc, maxc, maxc, minc, minc)
    by = (minr, minr, maxr, maxr, minr)
    ax.plot(bx, by, '-b', linewidth=2.5)

 
plt.show()

#-----------------------------------------------------------------------------.
##################################################
### Compute descriptors for all mascdb images ####
##################################################
### Define function to compute descriptors 
properties = ['area','filled_area', 'perimeter', 'eccentricity', 'equivalent_diameter',
              'euler_number',  'solidity', 'feret_diameter_max',
              'orientation',
              'minor_axis_length', 'major_axis_length'] # 'axis_major_length','axis_minor_length'] 

def _descriptor_fun(img, properties):
    # Obtain a binary image
    threshold = filters.threshold_otsu(img)
    binary_mask = img > threshold
    # Remove inner holes and other objects 
    binary_mask = morphology.remove_small_objects(binary_mask, 50)
    binary_mask = morphology.remove_small_holes(binary_mask, 50)
    # Labels regions 
    # --> Background is label 0 
    labels = measure.label(binary_mask)
    # Remove spurious stuffs 
    labels[labels != 1] = 0  # TODO: check if 1 is always the largest
    # Compute properties 
    region_props = measure.regionprops(labels)
    props_dict = measure.regionprops_table(labels, properties = properties) 
    props_dict = {k: v[0] for k,v in props_dict.items()}
    props_arr = np.array([props_dict[prop] for prop in properties])
    return props_arr

# Check it works properly for many images (do it manually)
img = mascdb.da.isel(CAM_ID=0, TripletID=0).values
props_arr = _descriptor_fun(img=img, properties = properties)

# Check that it works for a subset of mascdb (first 100, then 1000) to 
#  ensure no memory leak in the code 

# http://image.dask.org/en/latest/coverage.html for dask-optimized image processing functions

# properties = ["area"]
da = mascdb.da.isel(TripletID = slice(0,100))
x = "x"
y = "y"
fun_kwargs = {'properties': properties}
fun = _descriptor_fun
labels = properties 

da_descriptors = _compute_2Dimage_descriptors(da = mascdb.da.isel(TripletID = slice(0,100)), 
                                              fun = _descriptor_fun,
                                              labels = labels, 
                                              fun_kwargs = fun_kwargs)

da_descriptors = _compute_2Dimage_descriptors(da = mascdb.da.isel(TripletID = slice(0,500)), 
                                              fun = _descriptor_fun,
                                              labels = labels, 
                                              fun_kwargs = fun_kwargs)

da_descriptors = _compute_2Dimage_descriptors(da = mascdb.da.isel(TripletID = slice(0,1000)), 
                                              fun = _descriptor_fun,
                                              labels = labels, 
                                              fun_kwargs = fun_kwargs)

#-----------------------------------------------------------------------------.
## Run descriptors computation  
# - If the dask.distributed.Client() is called above, navigate to
#   http://localhost:8787/status to see the diagnostic dashboard
#   and control the computation progress

# - If the dask.distributed.client() is not called, it will use by default the 
#   dask.scheduler('threaded') and a progressbar will be displayed in the terminal 

da_descriptors = _compute_2Dimage_descriptors(da = da, 
                                              fun = _descriptor_fun,
                                              labels = labels, 
                                              fun_kwargs = fun_kwargs)

#-----------------------------------------------------------------------------.
# Apply descriptors manually to mascdb 
cam0 = da_descriptors.isel(CAM_ID = 0).to_dataset('descriptor').to_pandas().drop(columns='CAM_ID')
cam1 = da_descriptors.isel(CAM_ID = 1).to_dataset('descriptor').to_pandas().drop(columns='CAM_ID')
cam2 = da_descriptors.isel(CAM_ID = 2).to_dataset('descriptor').to_pandas().drop(columns='CAM_ID')
new_mascdb = mascdb.add_cam_columns(cam0=cam0, cam1=cam1, cam2=cam2, force=False, complete=True)

#-----------------------------------------------------------------------------.
# Apply directly to the full mascdb (when you are sure no memory leaks !)
new_mascdb = mascdb.compute_2Dimage_descriptors(fun = _descriptor_fun,
                                                labels = labels, 
                                                fun_kwargs = fun_kwargs)

##----------------------------------------------------------------------------.    
 
############################ 
#### Remove condensation ###
############################ 
# --> FFT of first image 
# --> FFT of second image 
# --> FFT of third image 
# --> Substract median to each
# --> Difference and reconstruct 

#-----------------------------------------------------------------------------.
############################
#### Image enhancements ####
############################
# Plot largest        
mascdb = mascdb.arrange('cam0.Dmax', decreasing=True) 
mascdb.plot_triplets(indices=list(range(0,8)), zoom=False)   # same images? 
mascdb.plot_triplets(indices=list(range(0,10)), zoom=True, enhancement=None)   

# Plot smallest
mascdb = mascdb.arrange('cam0.Dmax', decreasing=False) 
mascdb.plot_triplets(indices=list(range(0,10)), zoom=True)     
mascdb.plot_triplets(indices=list(range(0,10)), zoom=True, enhancement=None)   

#-----------------------------------------------------------------------------.
from mascdb.utils_img import xri_contrast_stretching 
from mascdb.utils_img import xri_hist_equalization 
from mascdb.utils_img import xri_local_hist_equalization 
  
mascdb = mascdb.arrange('cam0.Dmax', decreasing=True) 
da = mascdb.da.isel(CAM_ID=0, TripletID=0)

### Comparison 
img = da.values
img_rescale = xri_contrast_stretching(da, pmin=2, pmax=98).values
img_global_eq = xri_hist_equalization(da, adaptive=False).values
img_adapt_eq = xri_hist_equalization(da, adaptive=True).values
img_local_eq = xri_local_hist_equalization(da).values
 
# Plot 
l_imgs = [img, img_rescale, img_global_eq, img_adapt_eq, img_local_eq]
l_titles = ["Original", "Contrast stretching", "Global Histogram Equalization", 
            "Adaptive Histogram Equalization", "Local Histogram Equalization"]
fig, axs = plt.subplots(3,2,figsize=(7,10))
for i, ax in enumerate(axs.flatten()):
    if i <= 5:
        ax.imshow(l_imgs[i], cmap="gray", vmin=0, vmax=255)
        ax.set_title(l_titles[i])
        ax.set_axis_off()
plt.show()

# Random plots 
x = 'x'
y = 'y'
da = mascdb.da.isel(CAM_ID=0, TripletID=1000)
da.plot.imshow(x=x, y=y, row="TripletID", col="CAM_ID", vmin=0, vmax=255, cmap="gray")
xri_contrast_stretching(da, pmin=2, pmax=98).plot.imshow(x=x, y=y, row="TripletID", col="CAM_ID", vmin=0, vmax=255, cmap="gray")
xri_hist_equalization(da, adaptive=False).plot.imshow(x=x, y=y, row="TripletID", col="CAM_ID", vmin=0, vmax=255, cmap="gray")
xri_hist_equalization(da, adaptive=True).plot.imshow(x=x, y=y, row="TripletID", col="CAM_ID", vmin=0, vmax=255, cmap="gray")
# Slow !!!
xri_local_hist_equalization(da).plot.imshow(x=x, y=y, row="TripletID", col="CAM_ID", vmin=0, vmax=255, cmap="gray")

#-----------------------------------------------------------------------------.
