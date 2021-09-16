#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:56:51 2021

@author: ghiggi
"""
import os
import shutil
import copy
import xarray as xr
import pandas as pd
import numpy as np
import mascdb.pd_sns_accessor 
from mascdb.utils_event import _define_event_id
from mascdb.utils_event import _get_timesteps_duration
 
from mascdb.aux import get_snowflake_class_name_dict
from mascdb.aux import get_riming_class_name_dict
from mascdb.aux import get_melting_class_name_dict
  
from mascdb.utils_img import xri_zoom
from mascdb.utils_img import xri_contrast_stretching 
from mascdb.utils_img import xri_hist_equalization 
from mascdb.utils_img import xri_local_hist_equalization 
  
#-----------------------------------------------------------------------------.
# TODO: add stuff of mm @Jacobo did 
# add option to report mm ... 
# pix_size = self._triplet.loc[index].pix_size*1e3  # mm
# xsize = out.shape[0]*pix_size # mm
# ysize = out.shape[1]*pix_size # mm

# zarr as a zarr store ... with ID ordered 
#-----------------------------------------------------------------------------.
### In future ###
## filter
# mascdb = mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter(triplet.Dmax > 2 and cam1.Xhi > 2)
# mascdb.sel(label)
#-----------------------------------------------------------------------------.

#-----------------------------------------------------------------------------.
##############
### Checks ###
##############
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

def _check_enhancement(enhancement): 
    if not isinstance(enhancement, (type(None), str)): 
        raise TypeError("'enhancement' must be a string (or None).")
    if isinstance(enhancement, str): 
        valid_enhancements=["histogram_equalization", "contrast_stretching", "local_equalization"]
        if enhancement not in valid_enhancements:
            raise ValueError("{!r} is not a valid enhancement. "
                             "Specify one of {}".format(enhancement, valid_enhancements))
            
def _check_isel_idx(idx, vmax):
    # Return a numpy array of positional idx 
    if not isinstance(idx, (int,list, pd.Series, np.ndarray)):
        raise ValueError("Expecting isel idx to be int, list/pd.Series/np.array of int or boolean.")
    # Reformat all types to unique format 
    if isinstance(idx, np.ndarray):
        if idx.dtype.name == 'bool':
            idx = np.where(idx)[0]
        elif idx.dtype.name == 'int64':
            idx = np.array(idx)
        else:
            raise ValueError("Expecting idx np.array to be of 'bool' or 'int64' type.")
    if isinstance(idx, pd.Series):
        if idx.dtype.name == 'bool':
            idx = np.where(idx)[0]
        else: 
            idx = np.array(idx)
    if isinstance(idx, int): 
        idx = np.array([idx])
    if isinstance(idx, list): 
        idx = np.array(idx)
        if idx.dtype.name == 'bool':
            idx = np.where(idx)[0]
        elif idx.dtype.name == 'int64':
            idx = np.array(idx)
        else:
            raise ValueError("Expecting values in the idx list to be of 'bool' or 'int' type.")
    #--------------------------------------------------------------------.
    # Check idx validity 
    if np.any(idx > vmax):
        raise ValueError("The maximum positional idx is {}".format(vmax))
    if np.any(idx < 0):
        raise ValueError("The positional idx must be positive integers.")
    #--------------------------------------------------------------------.
    # Return idx 
    return idx 

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
        
        self._da = xr.open_zarr(zarr_store_fpath)['data']
        self._da.name = "MASC Images"

        # Read data into dataframes
        self._cam0    = pd.read_parquet(cam0_fpath)
        self._cam1    = pd.read_parquet(cam1_fpath)
        self._cam2    = pd.read_parquet(cam2_fpath)
        self._triplet = pd.read_parquet(triplet_fpath)
        
        # Number of triplets 
        self._n_triplets = len(self._triplet)
        
        # Save dir_path info 
        self._dir_path = dir_path
    
    def __len__(self):
        return self._n_triplets
    
    ##------------------------------------------------------------------------.
    ####################
    ## Print method  ###
    ####################
    def __str__(self):
        print("MASCDB data structure:")
        print("-------------------------------------------------------------------------")
        print("- mascdb.da:")
        print(self._da)
        print("-------------------------------------------------------------------------")
        print("- mascdb.cam0, mascdb.cam1, mascdb.cam2:")
        print(self._cam0) 
        print("-------------------------------------------------------------------------")
        print("- mascdb.triplet")
        print(self._triplet)
        print("-------------------------------------------------------------------------")
        print("- mascdb.env")
        print(self.env)
        print("-------------------------------------------------------------------------")
        print("- mascdb.bs")
        print(self.bs)
        print("-------------------------------------------------------------------------")
        print("- mascdb.gan3d")
        print(self.gan3d)
        print("-------------------------------------------------------------------------")
        return "" 
    
    def __repr__(self):
        return self.__str__()
        
    ##------------------------------------------------------------------------.
    #################
    ## Save db    ###
    #################
    def save(self, dir_path, force=False):
        # - Check there are data to save
        if self._n_triplets == 0: 
            raise ValueError("Nothing to save. No data left in the MASCDB.")
        # - Check dir_path
        if dir_path == self._dir_path:
            if force: 
                print("- Overwriting existing 'source' MASCDB at {}".format(dir_path))
                shutil.rmtree(dir_path)
            else:
                raise ValueError("If you want to overwrite the existing MASCDB at {},"
                                 "please specify force=True".format(dir_path))
        if os.path.exists(dir_path): 
            if force: 
                print("- Replacing content of directory {}.".format(dir_path))
                shutil.rmtree(dir_path)  
            else: 
                raise ValueError("A directory already exists at {}."
                                 "Please specify force=True if you want to overwrite it.".format(dir_path))
        #---------------------------------------------------------------------.
        # - Create directory 
        os.makedirs(dir_path)
        
        # - Define fpath of databases 
        zarr_store_fpath = os.path.join(dir_path,"MASCdb.zarr")
        cam0_fpath = os.path.join(dir_path, "MASCdb_cam0.parquet")
        cam1_fpath = os.path.join(dir_path, "MASCdb_cam1.parquet")
        cam2_fpath = os.path.join(dir_path, "MASCdb_cam2.parquet")
        triplet_fpath = os.path.join(dir_path, "MASCdb_triplet.parquet")
        
        # - Write databases 
        ds = self._da.to_dataset(name="data") 
        ds.to_zarr(zarr_store_fpath)
        self._cam0.to_parquet(cam0_fpath)
        self._cam1.to_parquet(cam1_fpath)
        self._cam2.to_parquet(cam2_fpath)
        self._triplet.to_parquet(triplet_fpath)
        #---------------------------------------------------------------------.
        return None 
        
    ##------------------------------------------------------------------------.
    #################
    ## Subsetting ###
    #################
    def isel(self, idx): 
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # Check valid idx 
        idx = _check_isel_idx(idx, vmax=self._n_triplets-1)
        ##--------------------------------------------------------------------.
        # Subset all datasets 
        self._da = self._da.isel(TripletID=idx)
        if isinstance(idx[0], bool):
            self._cam0 = self._cam0[idx]  
            self._cam1 = self._cam1[idx]  
            self._cam2 = self._cam2[idx]  
            self._triplet = self._triplet[idx]  
        else: 
            self._cam0 = self._cam0.iloc[idx]
            self._cam1 = self._cam1.iloc[idx]
            self._cam2 = self._cam2.iloc[idx]
            self._triplet = self._triplet.iloc[idx]
        ##--------------------------------------------------------------------.
        # Update number of triplets 
        self._n_triplets = len(self._triplet)
        return self 
    
    def sample_n(self, n=10):
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))      
        idx = np.random.choice(self._n_triplets, n) 
        return self.isel(idx)
        
    def first_n(self, n=10):  
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))
        idx = np.arange(n)
        return self.isel(idx)
    
    def last_n(self, n=10): 
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))
        idx = np.arange(self._n_triplets-1,self._n_triplets-n-1, step=-1)
        return self.isel(idx)
    
    def head(self, n=10): 
        n = min(self._n_triplets, n)
        idx = np.arange(n)
        return self.isel(idx)
    
    def tail(self, n=10):
        n = min(self._n_triplets, n)
        idx = np.arange(self._n_triplets-1,self._n_triplets-n-1, step=-1)
        return self.isel(idx)
     
    ##------------------------------------------------------------------------.
    ############ 
    ### Sort ###
    ############
    def arrange(self, expression, decreasing=True):
        # Check expression type 
        if not isinstance(expression, str):
            raise TypeError("'expression' must be a string.")
        #------------------------------.
        # Retrieve db name and column 
        split_expression = expression.split(".")
        db_name = split_expression[0]
        db_column = split_expression[1]
        # Check valid format 
        if len(split_expression) != 2:
            raise ValueError("An unvalid 'expression' has been specified.\n"
                             "The expected format is <cam*/triplet/env/bs>.<column_name>.")
        # Check valid db 
        valid_db = ['cam0', 'cam1','cam2','triplet','bs','env','gan3d']
        if db_name not in valid_db:
            raise ValueError("The first component must be one of {}".format(valid_db))
        #------------------------------.
        # Get db 
        db = getattr(self, db_name)
        # Check valid column 
        valid_columns = list(db.columns)
        if db_column not in valid_columns:
            raise ValueError("{!r} is not a column of {!r}. Valid columns are {}".format(db_column, db_name, valid_columns))
        #------------------------------.
        # Retrieve sorting idx 
        idx = db[db_column].to_numpy().argsort()
        if decreasing: 
            idx = idx[::-1]
        #------------------------------.
        # Return sorted object
        return self.isel(idx)
    
    def select_max(self, expression, n=10):
        return self.arrange(expression, decreasing=True).isel(np.arange(min(n, self._n_triplets)))

    def select_min(self, expression, n=10):
        return self.arrange(expression, decreasing=False).isel(np.arange(min(n, self._n_triplets)))   
    
    ##------------------------------------------------------------------------.
    ################
    ### Filters ####
    ################
    def from_campaign(self, campaign):
        if not isinstance(campaign, (list, str)): 
            raise TypeError("'campaign' must be a string or a list of strings.")
        if isinstance(campaign, str): 
            campaign = [campaign]
        # Convert to numpy array with str type (not object...)
        campaign = np.array(campaign).astype(str)
        campaigns_arr = self.triplet['campaign'].values.astype(str)
        valid_campaigns = np.unique(campaigns_arr)
        unvalid_campaigns_arg = campaign[np.isin(campaign, valid_campaigns, invert=True)]
        if len(unvalid_campaigns_arg) > 0: 
            raise ValueError("{} is not a campaign of the current mascdb. "
                             "Valid campaign names are {}".format(unvalid_campaigns_arg.tolist(),
                                                                  valid_campaigns.tolist()))
        idx = np.isin(campaigns_arr, campaign)
        return self.isel(idx) 
    
    def exclude_campaign(self, campaign):
        if not isinstance(campaign, (list, str)): 
            raise TypeError("'campaign' must be a string or a list of strings.")
        if isinstance(campaign, str): 
            campaign = [campaign]
        # Convert to numpy array with str type (not object...)
        campaign = np.array(campaign).astype(str)
        campaigns_arr = self.triplet['campaign'].values.astype(str)
        valid_campaigns = np.unique(campaigns_arr)
        unvalid_campaigns_arg = campaign[np.isin(campaign, valid_campaigns, invert=True)]
        if len(unvalid_campaigns_arg) > 0: 
            raise ValueError("{} is already not a campaign of the current mascdb. "
                             "Current mascdb has campaign names {}".format(unvalid_campaigns_arg.tolist(),
                                                                           valid_campaigns.tolist()))
        idx = np.isin(campaigns_arr, campaign, invert=True)
        return self.isel(idx) 
     
    def select_snowflake_class(self, values, method='Praz2017'): 
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_snowflake_class_name_dict().values()) # id
            column = 'snowflake_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_snowflake_class_name_dict().keys())   # name
            column = 'snowflake_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve triplet column values 
        arr = self.triplet[column].values
        # Check values are valid 
        unvalid_values = values[np.isin(values, valid_names, invert=True)]
        if len(unvalid_values) > 0: 
            raise ValueError("{} is not a {} of the current mascdb. "
                             "Current mascdb has {} values {}".format(unvalid_values.tolist(),
                                                                      column, column,
                                                                      valid_names))
        #---------------------------------------------------------------------.
        # Subset the mascdb
        idx = np.isin(arr, values)
        return self.isel(idx) 

    def select_riming_class(self, values, method='Praz2017'): 
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_riming_class_name_dict().values()) # id
            column = 'riming_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_riming_class_name_dict().keys())   # name
            column = 'riming_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve triplet column values 
        arr = self.triplet[column].values
        # Check values are valid 
        unvalid_values = values[np.isin(values, valid_names, invert=True)]
        if len(unvalid_values) > 0: 
            raise ValueError("{} is not a {} of the current mascdb. "
                             "Current mascdb has {} values {}".format(unvalid_values.tolist(),
                                                                      column, column,
                                                                      valid_names))
        #---------------------------------------------------------------------.
        # Subset the mascdb
        idx = np.isin(arr, values)
        return self.isel(idx)     
 
    def select_melting_class(self, values, method='Praz2017'): 
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_melting_class_name_dict().values()) # id
            column = 'melting_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_melting_class_name_dict().keys())   # name
            column = 'melting_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve triplet column values 
        arr = self.triplet[column].values
        # Check values are valid 
        unvalid_values = values[np.isin(values, valid_names, invert=True)]
        if len(unvalid_values) > 0: 
            raise ValueError("{} is not a {} of the current mascdb. "
                             "Current mascdb has {} values {}".format(unvalid_values.tolist(),
                                                                      column, column,
                                                                      valid_names))
        #---------------------------------------------------------------------.
        # Subset the mascdb
        idx = np.isin(arr, values)
        return self.isel(idx)  
  
    ##------------------------------------------------------------------------.
    ################
    ### Getters ####
    ################  
    # The following properties are used to avoid accidental modification in place by the user
    @property
    def da(self):
        return self._da.copy()
    
    @property
    def cam0(self):
        return self._cam0.copy()
    
    @property
    def cam1(self):
        return self._cam1.copy()
    
    @property
    def cam2(self):
        return self._cam2.copy()
    
    @property
    def triplet(self):
        return self._triplet.copy()
    
    # The following properties are just utils
    @property
    def env(self):
        columns = list(self._triplet.columns)
        env_variables = [column for column in columns if column.startswith("env_")]
        env_db = self.triplet[[*env_variables]]
        env_db.columns = [column.strip("env_") for column in env_variables]
        return env_db
    
    @property
    def bs(self):
        columns = list(self._triplet.columns)
        bs_variables = [column for column in columns if column.startswith("bs_")]
        bs_db = self.triplet[[*bs_variables]]
        bs_db.columns = [column.strip("bs_") for column in bs_variables]
        return bs_db
    
    @property
    def gan3d(self):
        columns = list(self._triplet.columns)
        gan3d_variables = [column for column in columns if column.startswith("gan3d_")]
        gan3d_db = self.triplet[[*gan3d_variables]]
        gan3d_db.columns = [column.strip("gan3d_") for column in gan3d_variables]
        return gan3d_db
    
    @property
    def full_db(self):
        # Add CAM_ID to each cam db 
        l_cams = [self.cam0, self.cam1, self.cam2]
        for i, cam in enumerate(l_cams):
             cam['CAM_ID'] = i
        # Merge cam(s) db into 
        full_db = pd.concat(l_cams)
        # Add triplet variables to fulldb 
        # TODO: to remove/modify in future
        labels_vars = ['riming_deg_level', 'riming_id','riming_id_prob', 
                       'melting_id',
                       'melting_prob',  
                       'snowflake_class_name',       
                       'snowflake_class_id',
                       'snowflake_class_id_prob']
        vars_not_add = ['pix_size','quality_xhi_flake','n_roi', 'Dmax_flake'] + labels_vars
        triplet = self.triplet.drop(columns=vars_not_add)
        full_db = full_db.merge(triplet, how="left")
        return full_db
    
    def ds_images(self, CAM_ID = None, campaign=None, img_id='img_id'):
        #----------------------------------------------------------------------.
        # Subset by campaign 
        if campaign is not None:
            if isinstance(campaign, str): 
                campaign = [campaign]
            campaign = np.array(campaign).astype(str)
            db_campaigns = self.triplet['campaign'].values.astype(str)
            valid_campaigns = np.unique(db_campaigns)
            unvalid_campaigns_arg = campaign[np.isin(campaign, valid_campaigns, invert=True)]
            if len(unvalid_campaigns_arg) > 0: 
                raise ValueError("{} is not a campaign of the current mascdb. "
                                 "Valid campaign names are {}".format(unvalid_campaigns_arg.tolist(),
                                                                      valid_campaigns.tolist()))
            idx = np.isin(db_campaigns, campaign)
            da = self.isel(idx).da
        else: 
            da = self.da
        #----------------------------------------------------------------------.
        # Subset cam images
        if CAM_ID is not None:
            da = da.isel(CAM_ID = CAM_ID)
        #----------------------------------------------------------------------.
        ### Retrieve dimensions to eventually stack along a new third dimension 
        dims = list(da.dims)
        unstacked_dims = list(set(dims).difference(["x","y"]))
        # - If only x and y, add third dimension img_id
        if len(unstacked_dims) == 0: 
            stack_dict = {}
            da_stacked = da.expand_dims(img_id, axis=-1)
        # - If there is already a third dimension, transpose to the last 
        elif len(unstacked_dims) == 1:
            da = da.rename({unstacked_dims[0]: img_id})
            da_stacked = da.transpose(..., img_id) 
        # - If there is more than 3 dimensions, stack it all into a new third dimension
        elif len(unstacked_dims) > 1: 
            stack_dict = {img_id: unstacked_dims}
            # Stack all additional dimensions into a 3D array with all img_id in the last dimension 
            da_stacked = da.stack(stack_dict).transpose(..., img_id)
        else: 
            raise NotImplementedError()
        return da_stacked   
    
    ##------------------------------------------------------------------------.
    ##############################
    ### Datetime/Event utils #####
    ##############################
    def _add_event_n_triplets(self):
        if 'event_id' not in list(self._triplet.columns):
            raise ValueError("First define 'event_id' using mascdb.define_event_id().")
        event_triplets = self.triplet.groupby('event_id').size()
        event_triplets.name = "event_n_triplets"
        self._triplet = self._triplet.merge(event_triplets, on="event_id")
        
    def _add_event_duration(self):
        if 'event_id' not in list(self._triplet.columns):
            raise ValueError("First define 'event_id' using mascdb.define_event_id().")
        event_durations = self.triplet.groupby('event_id').apply(lambda x: _get_timesteps_duration(x.datetime, unit="m"))
        event_durations.name = "event_duration"
        self._triplet = self._triplet.merge(event_durations, on="event_id")
        
    def define_event_id(self, timedelta_thr): 
        # - Extract relevant columns from triplet db
        db = self.triplet[['campaign','datetime']]
        # - Retrieve campaign_ids 
        campaign_ids = np.unique(db['campaign'])
        # - Retrieve event_id column 
        db['event_id'] = -1
        max_event_id = 0
        for campaign_id in campaign_ids:
            # - Retrieve row index of specific campaign 
            idx_campaign = db['campaign'] == campaign_id
            # - Define event_ids for the campaign 
            campaign_event_ids = _define_event_id(timesteps = db.loc[idx_campaign,'datetime'], 
                                                 timedelta_thr=timedelta_thr)
            # - Add offset to ensure having an unique event_id across all campaigns 
            campaign_event_ids = campaign_event_ids + max_event_id  
            # - Add event_id to the campaign subset of the database 
            db.loc[idx_campaign, 'event_id'] = campaign_event_ids
            # - Update the current maximum event_id
            max_event_id = max(campaign_event_ids) + 1
         
        # - Add event_id column to all databases 
        self._cam0['event_id'] = db['event_id']
        self._cam1['event_id'] = db['event_id']
        self._cam2['event_id'] = db['event_id']
        self._triplet['event_id'] = db['event_id']
        
        # - Add also duration and n_triplets for each event 
        self._add_event_duration()
        self._add_event_n_triplets()
        
    ##------------------------------------------------------------------------.
    ###########################
    ### Plotting routines #####
    ###########################

    def plot_triplets(self, indices=None, random = False, n_triplets = 1,
                      enhancement="histogram_equalization",
                      zoom=True, **kwargs):
        #--------------------------------------------------.
        # Retrieve number of valid index
        n_idxs = len(self._triplet.index)
        if n_idxs == 0: 
            raise ValueError("No data to plot.") 
        #--------------------------------------------------.
        # Check args
        _check_random(random)
        _check_zoom(zoom)
        _check_enhancement(enhancement)
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
        da_subset = self._da.isel(TripletID = indices).transpose(...,'CAM_ID','TripletID')
        
        #--------------------------------------------------.
        # Apply enhancements 
        if enhancement is not None: 
            if enhancement == "histogram_equalization":
               da_subset =  xri_hist_equalization(da_subset, adaptive=True)
            elif enhancement == "contrast_stretching":
               da_subset =  xri_contrast_stretching(da_subset, pmin=2, pmax=98)
            elif enhancement == "local_equalization":
               da_subset = xri_local_hist_equalization(da_subset)
               
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            da_subset = xri_zoom(da_subset, squared=False)
            
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
            
    def plot_flake(self, CAM_ID=None, index=None, random = False,
                   enhancement="histogram_equalization",
                   zoom=True, ax=None, **kwargs):
        # Check args
        _check_random(random)
        _check_zoom(zoom)
        _check_enhancement(enhancement)
        
        #--------------------------------------------------.
        # Retrieve number of valid index
        n_idxs = len(self._triplet.index)
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
        da_img = self._da.isel(TripletID = index, CAM_ID = CAM_ID) 
        
        #--------------------------------------------------.
        # Apply enhancements 
        if enhancement is not None: 
            if enhancement == "histogram_equalization":
               da_img =  xri_hist_equalization(da_img, adaptive=True)
            elif enhancement == "contrast_stretching":
               da_img =  xri_contrast_stretching(da_img, pmin=2, pmax=98)
            elif enhancement == "local_equalization":
               da_img = xri_local_hist_equalization(da_img)
           
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            da_img = xri_zoom(da_img, squared=False)
            
        #--------------------------------------------------.
        # Plot single image 
        # - TODO: 'aspect' cannot be specified without 'size'
        p = da_img.plot(x='x',y='y', ax=ax,
                        cmap='gray', add_colorbar=False, 
                        vmin=0, vmax=255,
                        **kwargs)
        #--------------------------------------------------. 
        return p  

    def plot_flakes(self, CAM_ID=None, indices=None, random = False, 
                    n_images = 9, col_wrap = 3,
                    enhancement="histogram_equalization",
                    zoom=True, 
                    **kwargs):
        # Retrieve number of valid index
        n_idxs = len(self._triplet.index) # TODO 
        # TODO: 
        # - Option to stack to ImageID ... and then sample that ... title will report TripletID and CAM ID

        #--------------------------------------------------.
        # Check args
        _check_random(random)
        _check_zoom(zoom)
        _check_enhancement(enhancement)
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
        # If a single flake is specified, plot it with plot_flake 
        if len(indices) == 1: 
            print("It's recommended to use 'plot_flake()' to plot a single image.")
            return self.plot_flake(index=indices[0], CAM_ID=CAM_ID, random=random,
                                   enhancement=enhancement, zoom=zoom, *kwargs)
        
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If CAM_ID is an integer (instead of list length 1 ... the CAM_ID dimension is dropped)
        da_subset = self._da.isel(TripletID = indices, CAM_ID = CAM_ID).transpose(...,'TripletID')
        
        #--------------------------------------------------.
        # Apply enhancements 
        if enhancement is not None: 
            if enhancement == "histogram_equalization":
               da_subset =  xri_hist_equalization(da_subset, adaptive=True)
            elif enhancement == "contrast_stretching":
               da_subset =  xri_contrast_stretching(da_subset, pmin=2, pmax=98)
            elif enhancement == "local_equalization":
               da_subset = xri_local_hist_equalization(da_subset)
               
        #--------------------------------------------------.
        # Zoom all images to same extent 
        if zoom:
            da_subset = xri_zoom(da_subset, squared=False)
        
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

 