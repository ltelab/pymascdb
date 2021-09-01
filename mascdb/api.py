#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:56:51 2021

@author: ghiggi
"""
import os
import xarray as xr
import pandas as pd
import numpy as np
import mascdb.pd_sns_accessor 

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

def _center_image(img, nrow, ncol): 
    r, c = img.shape 
    col_incr = int((ncol - c)/2)
    row_incr = int((nrow - r)/2)
    arr = np.zeros((nrow, ncol))
    arr[slice(row_incr,row_incr+r), slice(col_incr,col_incr+c)] = img 
    return arr 

#-----------------------------------------------------------------------------.
class MASC_DB:
    """
    Masc database class to read and manipulate the 4 databases of 
    descriptors (one for each cam) as well TODO

    """
    def __init__(self, dir_path):
        """
        TO DO
        
        """
        zarr_store_fpath = os.path.join(dir_path,"MASCdb.zarr")
        cam0_fpath = os.path.join(dir_path, "MASCdb_cam0.parquet")
        cam1_fpath = os.path.join(dir_path, "MASCdb_cam1.parquet")
        cam2_fpath = os.path.join(dir_path, "MASCdb_cam2.parquet")
        triplet_fpath = os.path.join(dir_path, "MASCdb_triplet.parquet")
        
        self.da = xr.open_zarr(zarr_store_fpath)['data']

        # Read data into dataframes
        self.cam0    = pd.read_parquet(cam0_fpath)
        self.cam1    = pd.read_parquet(cam1_fpath)
        self.cam2    = pd.read_parquet(cam2_fpath)
        self.triplet = pd.read_parquet(triplet_fpath)


    def plot_triplets(self, indices=None, random = True, n_triplets = 1, zoom=True, **kwargs):
        # Check args
        if indices is not None: 
            # Check valid indices 
            # TODO: 
            random = False
        if random: # (use df.sample)
            indices = list(pd.Index(np.random.choice(self.triplet.index, n_triplets)))
        # TODO: check max n_triplets and len(indices) !!!
        # TODO: add stuff of mm @Jacobo did 
        # add option to report mm ... 
        # pix_size = self.triplet.loc[index].pix_size*1e3  # mm
        # xsize = out.shape[0]*pix_size # mm
        # ysize = out.shape[1]*pix_size # mm
        #--------------------------------------------------.
        # Subset triplet(s) images 
        da_subset = self.da.isel(TripletID = indices).transpose(...,'CAM_ID','TripletID')
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            # Get list of images 
            l_imgs = [da_subset.isel(TripletID=i,CAM_ID=j).values for i in range(len(indices)) for j in range(3)]
            # Zoom a list of image 
            l_zoomed = [_get_zoomed_image(img) for img in l_imgs]
            # Ensure same shape across all zoomed images 
            l_shapes = [img.shape for img in l_zoomed]
            r_max, c_max = (max(n) for n in zip(*l_shapes)) # Get max number of row an columns 
            l_zoomed = [_center_image(img, nrow=r_max, ncol=c_max) for img in l_zoomed]
            # Reassign to da_subset 
            da_subset = da_subset.isel(x=slice(0, c_max), y=slice(0,r_max))
            l_stack = [np.stack(l_zoomed[i*3:3*(i+1)], axis=-1) for i in range(len(indices))]
            new_arr = np.stack(l_stack, axis=-1) # TODO: check order is correct 
            da_subset.values = new_arr
        #--------------------------------------------------.
        # Plot triplet(s)
        row = "TripletID" if len(indices) > 1 else None
        p = da_subset.plot(x='x',y='y',
                           col="CAM_ID",
                           row=row,
                           aspect=1,
                           cmap='gray', add_colorbar=False,
                           vmin=0, vmax=255, 
                           **kwargs)
        # Nice layout 
        for i, ax in enumerate(p.axes.flat):
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_axis_off()
        p.fig.subplots_adjust(wspace=0.01, hspace=0.01)   
        #--------------------------------------------------. 
        return p       
            
    def plot_flake(self, CAM_ID=None, index=None, random = True, zoom=True, **kwargs):
        # Check args
        if index is not None: 
            if not isinstance(index, int):
                raise TypeError("'plot_flake' expects 'index' as an integer.")
        if random: # (use df.sample)
           if index is None:
               index = list(pd.Index(np.random.choice(self.triplet.index, 1)))
           if CAM_ID is None: 
                CAM_ID = np.random.choice([0,1,2], 1)
             
        # TODO: check max n_triplets and len(indices) !!!
        # TODO: check valid CAM_ID if not None (int value, list... )

        # TODO: add stuff of mm @Jacobo did 
        # add option to report mm ... 
        # pix_size = self.triplet.loc[index].pix_size*1e3  # mm
        # xsize = out.shape[0]*pix_size # mm
        # ysize = out.shape[1]*pix_size # mm
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If CAM_ID is an integer (instead of list length 1 ... the CAM_ID dimension is dropped)
        da_img = self.da.isel(TripletID = index, CAM_ID = CAM_ID) 
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            # Zoom a list of image 
            zoomed_img = _get_zoomed_image(da_img.values)  
            r_max, c_max = zoomed_img.shape
            # Reassign to da_subset 
            da_img = da_img.isel(x=slice(0, c_max), y=slice(0,r_max))
            da_img.values = zoomed_img
        #--------------------------------------------------.
        # Plot triplet(s)
        p = da_img.plot(x='x',y='y', # TODO no aspect without size 
                        cmap='gray', add_colorbar=False, 
                        vmin=0, vmax=255,
                        **kwargs)
        #--------------------------------------------------. 
        return p  

    def plot_flakes(self, CAM_ID=None, indices=None, random = True, n_images = 9, zoom=True, **kwargs):
        # Check args
        if indices is not None: 
            if isinstance(indices, int):
                indices = list(indices)
            # Check valid indices 
            # TODO: 
            random = False
        if random: # (use df.sample)
           if indices is None:
               indices = list(pd.Index(np.random.choice(self.triplet.index, n_images)))
           if CAM_ID is None: 
                CAM_ID = np.random.choice([0,1,2], 1)
        print(len(indices))
        if len(indices) == 1: 
            print("You should use 'plot_flake' to plot a single image.")
            return self.plot_flake(index=indices[0], CAM_ID=CAM_ID, random=False, zoom=zoom, *kwargs)
            
        # TODO: check max n_images and len(indices) !!!
        # TODO: check valid CAM_ID if not None (int value, list... )
        
        # TODO: 
        # - Option to stack to ImageID ... and then sample that ... title will report TripletID and CAM ID
        
        # TODO: add stuff of mm @Jacobo did 
        # add option to report mm ... 
        # pix_size = self.triplet.loc[index].pix_size*1e3  # mm
        # xsize = out.shape[0]*pix_size # mm
        # ysize = out.shape[1]*pix_size # mm
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If CAM_ID is an integer (instead of list length 1 ... the CAM_ID dimension is dropped)
        da_subset = self.da.isel(TripletID = indices, CAM_ID = CAM_ID).transpose(...,'TripletID')
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            # Get list of images 
            l_imgs = [da_subset.isel(TripletID=i,).values for i in range(len(indices))]
            # Zoom a list of image 
            l_zoomed = [_get_zoomed_image(img) for img in l_imgs]
            # Ensure same shape across all zoomed images 
            l_shapes = [img.shape for img in l_zoomed]
            r_max, c_max = (max(n) for n in zip(*l_shapes)) # Get max number of row an columns 
            l_zoomed = [_center_image(img, nrow=r_max, ncol=c_max) for img in l_zoomed]
            # Reassign to da_subset 
            da_subset = da_subset.isel(x=slice(0, c_max), y=slice(0,r_max))
            new_arr = np.stack(l_zoomed, axis=-1) 
            da_subset.values = new_arr
        #--------------------------------------------------.
        # Plot triplet(s)
        row = "TripletID" if len(indices) > 1 else None
        p = da_subset.plot(x='x',y='y', 
                           row="TripletID", col_wrap=3, 
                           aspect=1, 
                           cmap='gray', add_colorbar=False, 
                           vmin=0, vmax=255,
                           **kwargs)
        # Nice layout 
        for i, ax in enumerate(p.axes.flat):
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_axis_off() 
        p.fig.subplots_adjust(wspace=0.01, hspace=0.1)
        #--------------------------------------------------. 
        return p  
     


 