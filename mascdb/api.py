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
import copy
# TODO 
# - RANDOM TO FALSE 

# TODO: add stuff of mm @Jacobo did 
# add option to report mm ... 
# pix_size = self.triplet.loc[index].pix_size*1e3  # mm
# xsize = out.shape[0]*pix_size # mm
# ysize = out.shape[1]*pix_size # mm

#-----------------------------------------------------------------------------.
### In future ###
## filter
# mascdb = mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter(triplet.Dmax > 2 and cam1.Xhi > 2)
# mascdb.sel(label)
#-----------------------------------------------------------------------------.

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
### Checks 

def _check_CAM_ID(CAM_ID): 
    "Return CAM_ID integer."
    if not isinstance(CAM_ID, (int, np.int64, list)):
        raise TypeError("'CAM_ID', if specified, must be an integer or list (of length 1).")
    if isinstance(CAM_ID, list):
        if len(CAM_ID) != 1:
           raise ValueError("Expecting a single value for 'CAM_ID'.")
        CAM_ID = int(CAM_ID[0])
    # Check value validity
    if CAM_ID not in [0,1,2]:
        raise ValueError("Valid values of 'CAM_ID' are [0,1,2].")
    # Return integer CAM_ID
    return CAM_ID

def _check_index(index, vmax):
    "Return index integer."
    if not isinstance(index, (int, np.int64, list)):
        raise TypeError("'index', if specified, must be an integer or list (of length 1).")
    if isinstance(index, list):
        if len(index) != 1:
           raise ValueError("Expecting a single value for 'index'.")
        index = int(index[0])
    # Check value validity
    if index < 0:
        raise ValueError("'index' must be a positive integer.")
    if index > vmax:
        raise ValueError("The largest 'index' can be {}".format(vmax))
    # Return integer index
    return index

def _check_indices(indices, vmax):
    "Return index integer."
    if not isinstance(indices, (int, np.int64, list)):
        raise TypeError("'indices', if specified, must be an integer or list (of length 1).")
    if isinstance(indices, (int, np.int64)):
        indices = [indices]
    # Check value validity
    if any([not isinstance(index, (int, np.int64)) for index in indices]):
        raise ValueError("All 'indices' values must be integers.")
    if any([index < 0 for index in indices]):
        raise ValueError("All 'indices' values must be positive integers.")
    if any([index > vmax for index in indices]):
        raise ValueError("The largest 'indices' value can be {}".format(vmax))
    # Return list of indices
    return indices

def _check_n_triplets(n_triplets, vmax):
    if not isinstance(n_triplets, int):
        raise TypeError("'n_triplets' must be an integer.")
    if n_triplets < 1:
        raise ValueError("'n_triplets' must be at least 1.")    
    if n_triplets > vmax: 
        raise ValueError("'n_triplets' must be maximum {}.".format(vmax))
    if n_triplets > 10: 
        raise ValueError("It's not recommended to plot more than 10 triplets of images.")       

def _check_n_images(n_images, vmax):
    if not isinstance(n_images, int):
        raise TypeError("'n_images' must be an integer.")
    if n_images < 1:
        raise ValueError("'n_images' must be at least 1.")    
    if n_images > vmax: 
        raise ValueError("'n_images' must be maximum {}.".format(vmax))
    if n_images > 25: 
        raise ValueError("It's not recommended to plot more than 25 images.")
        
def _check_random(random): 
    if not isinstance(random, bool):
        raise TypeError("'random' must be either True or False.")
        
def _check_zoom(zoom): 
    if not isinstance(zoom, bool):
        raise TypeError("'zoom' must be either True or False.")
        
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
        
        # Number of triplets 
        self.n_triplets = len(self.triplet)
        
    ##------------------------------------------------------------------------.
    #################
    ## Subsetting ###
    #################
    def isel(self, idx): 
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # TODO: Check valid idx 
        if isinstance(idx, pd.Series):
            if idx.dtype.name == 'bool':
                idx = np.where(idx)[0]
            else: 
                idx = np.array(idx)
        # if list of bool --> np.where 
        idx = list(idx)

        ##--------------------------------------------------------------------.
        # Subset all datasets 
        self.da = self.da.isel(TripletID=idx)
        if isinstance(idx[0], bool):
            self.cam0 = self.cam0[idx]  
            self.cam1 = self.cam1[idx]  
            self.cam2 = self.cam2[idx]  
            self.triplet = self.triplet[idx]  
        else: 
            self.cam0 = self.cam0.iloc[idx]
            self.cam1 = self.cam1.iloc[idx]
            self.cam2 = self.cam2.iloc[idx]
            self.triplet = self.triplet.iloc[idx]
        ##--------------------------------------------------------------------.
        # Update number of triplets 
        self.n_triplets = len(self.triplet)
        return self 
    
    def sample_n(self, n=10):
        idx = np.random.choice(self.n_triplets, n) 
        return self.isel(idx)
    
    def first_n(self, n=10):     
        idx = np.arange(n)
        return self.isel(idx)
    
    def head(self, n=10):     
        idx = np.arange(n)
        return self.isel(idx)
    
    def tail(self, n=10):     
        idx = np.arange(self.n_triplets-1,self.n_triplets-n-1, step=-1)
        return self.isel(idx)
    
    def last_n(self, n=10):     
        idx = np.arange(self.n_triplets-1,self.n_triplets-n-1, step=-1)
        return self.isel(idx)
    
   
    ##------------------------------------------------------------------------.
    ################
    ### Getters ####
    ################
    @property
    def env(self):
        columns = list(self.triplet.columns)
        env_variables = [column for column in columns if column.startswith("env_")]
        env_db = self.triplet[[*env_variables]]
        env_db.columns = [column.strip("env_") for column in env_variables]
        return env_db
    
    @property
    def bs(self):
        columns = list(self.triplet.columns)
        bs_variables = [column for column in columns if column.startswith("bs_")]
        bs_db = self.triplet[[*bs_variables]]
        bs_db.columns = [column.strip("bs_") for column in bs_variables]
        return bs_db
    
    @property
    def gan3d(self):
        # TODO: to modify 3dgan_ to gan3d_ in future
        columns = list(self.triplet.columns)
        gan3d_variables = [column for column in columns if column.startswith("3dgan_")]
        gan3d_db = self.triplet[[*gan3d_variables]]
        gan3d_db.columns = [column.strip("3dgan_") for column in gan3d_variables]
        return gan3d_db
    
    @property
    def full_db(self):
        triplet = self.triplet.copy()
        # Add CAM_ID to each cam db 
        l_cams = [self.cam0.copy(), self.cam1.copy(), self.cam2.copy()]
        for i, cam in enumerate(l_cams):
             cam['CAM_ID'] = i
        # Merge cam(s) db into 
        full_db = pd.concat(l_cams)
        # Add triplet variables to fulldb 
        # TODO: to remove/modify in future
        labels_vars = ['riming_deg_level', 'riming_id',  
                       'melting_id',
                       'melting_prob',  
                       'label_name',       
                       'label_id',
                       'label_id_prob']
        vars_not_add = ['pix_size','Xhi','n_roi', 'Dmax'] + labels_vars
        triplet = triplet.drop(columns=vars_not_add)
        full_db = full_db.merge(triplet, how="left")
        return full_db
    
    ##------------------------------------------------------------------------.
    ###########################
    ### Plotting routines #####
    ###########################

    def plot_triplets(self, indices=None, random = True, n_triplets = 1, zoom=True, **kwargs):
        #--------------------------------------------------.
        # Retrieve number of valid index
        n_idxs = len(self.triplet.index)
        if n_idxs == 0: 
            raise ValueError("No data to plot.") 
        #--------------------------------------------------.
        # Check args
        _check_random(random)
        _check_zoom(zoom)
        _check_n_triplets(n_triplets, vmax=n_idxs)
     
        #--------------------------------------------------.
        # Define index if is not provided
        if indices is None:
            if random: # (use df.sample)
               indices = list(np.random.choice(n_idxs, n_triplets))
            else: 
               indices = list(np.arange(0,n_triplets))
        #--------------------------------------------------.       
        # Check validity of indices
        indices = _check_indices(indices, vmax=n_idxs-1)
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
            new_arr = np.stack(l_stack, axis=-1)  
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
        _check_random(random)
        _check_zoom(zoom)
        #--------------------------------------------------.
        # Retrieve number of valid index
        n_idxs = len(self.triplet.index)
        if n_idxs == 0: 
            raise ValueError("No data to plot.")
        #--------------------------------------------------.
        # Define index if is not provided
        if index is None:
           if random:         
               index = list(np.random.choice(n_idxs, 1))[0]
           else:             
               index = 0   
        if CAM_ID is None: 
            if random:    
               CAM_ID = np.random.choice([0,1,2], 1)
            else:
               CAM_ID = 1
        #--------------------------------------------------.       
        # Check validty of CAM_ID and index 
        CAM_ID = _check_CAM_ID(CAM_ID)
        index = _check_index(index, vmax=n_idxs-1)
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If CAM_ID is an integer (instead of list of length 1), then the CAM_ID dimension is dropped)
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
        # Plot single image 
        # - TODO: 'aspect' cannot be specified without 'size'
        p = da_img.plot(x='x',y='y', 
                        cmap='gray', add_colorbar=False, 
                        vmin=0, vmax=255,
                        **kwargs)
        #--------------------------------------------------. 
        return p  

    def plot_flakes(self, CAM_ID=None, indices=None, random = True, 
                    n_images = 9, zoom=True, 
                    col_wrap = 3, **kwargs):
        # Retrieve number of valid index
        n_idxs = len(self.triplet.index) # TODO 
        # TODO: 
        # - Option to stack to ImageID ... and then sample that ... title will report TripletID and CAM ID

        #--------------------------------------------------.
        # Check args
        _check_random(random)
        _check_zoom(zoom)
        _check_n_images(n_images, vmax=n_idxs)
        #--------------------------------------------------
        # Define index if is not provided
        if indices is None:
            if random:  
               indices = list(np.random.choice(n_idxs, n_images))
            else: 
               indices = list(np.arange(0,n_images))
        if CAM_ID is None: 
            if random:    
                CAM_ID = np.random.choice([0,1,2], 1)[0]
            else:
                CAM_ID = 1
        #-------------------------------------------------- 
        # Check validity of indices and CAM_ID
        indices = _check_indices(indices, vmax=n_idxs-1)
        CAM_ID = _check_CAM_ID(CAM_ID)
        #-------------------------------------------------- 
        if len(indices) == 1: 
            print("It's recommended to use 'plot_flake()' to plot a single image.")
            return self.plot_flake(index=indices[0], CAM_ID=CAM_ID, random=random, zoom=zoom, *kwargs)
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
                           row=row, col_wrap=col_wrap, 
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
     
    ##------------------------------------------------------------------------.

 