#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:56:51 2021

@author: ghiggi
"""
import os
import shutil
import copy
import dask
import xarray as xr
import pandas as pd
import numpy as np
import mascdb.pd_sns_accessor  # this is required to add pandas sns accessor
from mascdb.utils_event import _define_event_id
from mascdb.utils_event import _get_timesteps_duration
 
from mascdb.aux import get_snowflake_class_name_dict
from mascdb.aux import get_riming_class_name_dict
from mascdb.aux import get_melting_class_name_dict
from mascdb.aux import get_precip_class_name_dict
from mascdb.aux import get_vars_class
from mascdb.aux import var_units
from mascdb.aux import var_explanations
  
from mascdb.utils_img import _compute_2Dimage_descriptors
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

## add flake_id as cam_df and triplet_df index !!!

## check add_columns* 

#-----------------------------------------------------------------------------.
### In future ###
## filter
# mascdb = mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter("cam0.Dmax > 2 and cam1.Xhi > 2")
# mascdb.filter(triplet.Dmax > 2 and cam1.Xhi > 2)
# mascdb.sel(label)
#-----------------------------------------------------------------------------.

####-----------------------------------------------------------------------------.
###############
#### Checks ###
###############
def _check_cam_id(cam_id): 
    "Return cam_id integer."
    if not isinstance(cam_id, (int, np.int64, list)):
        raise TypeError("'cam_id', if specified, must be an integer or list (of length 1).")
    if isinstance(cam_id, list):
        if len(cam_id) != 1:
           raise ValueError("Expecting a single value for 'cam_id'.")
        cam_id = int(cam_id[0])
    # Check value validity
    if cam_id not in [0,1,2]:
        raise ValueError("Valid values of 'cam_id' are [0,1,2].")
    # Return integer cam_id
    return cam_id

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
    if not isinstance(idx, (int,list, slice,  pd.Series, np.ndarray)):
        raise ValueError("isel expect a slice object, an integer or a list/pd.Series/np.array of int or boolean.")
    # Reformat all types to unique format (numpy array of integers)
    if isinstance(idx, int): 
        idx = np.array([idx])
    if isinstance(idx, list): 
        idx = np.array(idx)
        if idx.dtype.name not in ['bool','int64']:
            raise ValueError("Expecting values in the idx list to be of 'bool' or 'int' type.")
    if isinstance(idx, slice): 
        idx = np.arange(idx.start, idx.stop, idx.step)
        
    if isinstance(idx, pd.Series):
        if idx.dtype.name in ['bool','boolean']:
            idx = np.where(idx.values)[0]
        else: 
            idx = np.array(idx.values)
            
    if isinstance(idx, np.ndarray):
        if idx.dtype.name in ['bool','boolean']:
            idx = np.where(idx)[0]
        if idx.dtype.name != 'int64':
            raise ValueError("Expecting idx np.array to be of 'bool' or 'int64' type.")

    #--------------------------------------------------------------------.
    # Check idx validity 
    if np.any(idx > vmax):
        raise ValueError("The maximum positional idx is {}".format(vmax))
    if np.any(idx < 0):
        raise ValueError("The positional idx must be positive integers.")
        
    #--------------------------------------------------------------------.
    # Return idx 
    return idx 

def _check_sel_ids(ids, valid_ids):
    # Return a numpy array of str 
    if not isinstance(ids, (str, list, pd.Series, np.ndarray)):
        raise ValueError("sel expect a string or a list/pd.Series/np.array of str")
    if valid_ids.dtype.name == "object":
        valid_ids = valid_ids.astype(str)
    # Reformat all types to unique format (numpy array of strings)
    if isinstance(ids, str): 
        ids = np.array([ids])
    if isinstance(ids, list): 
        ids = np.array(ids)
        if ids.dtype.name == 'object':
            ids = ids.astype(str)
        if not ids.dtype.name.startswith('str'):
            raise ValueError("Expecting values in the list to be strings.")
    if isinstance(ids, pd.Series):
        ids = np.array(ids.values)
    if isinstance(ids, np.ndarray):
        if ids.dtype.name == 'object':
            ids = ids.astype(str)
        if not ids.dtype.name.startswith('str'):
            raise ValueError("Expecting values in the np.array to be strings.")
    #--------------------------------------------------------------------.
    # Check idx validity 
    unvalid_ids = ids[np.isin(ids, valid_ids, invert=True)]
    if len(unvalid_ids) > 0:
        raise ValueError("The following ids are not valid: {}".format(unvalid_ids.tolist()))
    if len(ids) == 0: 
        raise ValueError("THIS SHOULD NOT OCCUR ...")
    #--------------------------------------------------------------------.    
    return ids 


def _check_timedelta(timedelta): 
    if not isinstance(timedelta, (pd.Timedelta, np.timedelta64)):
        raise TypeError("'timedelta' must be a pd.Timedelta or np.timedelta64 instance.")
    return timedelta 

def _check_df(df, name=None):
    if not isinstance(df, pd.DataFrame):
        if name is not None: 
            raise TypeError("{} is not a pd.DataFrame".format(name))
        else: 
           raise TypeError("Expecting a pd.DataFrame".format(name)) 
    return df 

def _count_occurence(x): 
    return [dict(zip(list(x.value_counts().keys()), list(x.value_counts())))]

def _convert_object_to_string(df):
    idx_object = df.dtypes.values == "object"
    columns = df.columns[np.where(idx_object)]
    for column in columns:
        df[column] = df[column].astype('string')
    return df 

####-----------------------------------------------------------------------------.
class MASC_DB:
    """
    Masc database class to read and manipulate the 4 databases of 
    descriptors (one for each cam) as well TODO

    """
    #####################
    #### Read MASCDB ###
    ##################### 
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
        self._da['flake_id'] = self._da['flake_id'].astype(str)
        self._da.name = "MASC Images"
  
        # Read data into dataframes
        self._cam0    = pd.read_parquet(cam0_fpath).set_index('flake_id', drop=False)
        self._cam1    = pd.read_parquet(cam1_fpath).set_index('flake_id', drop=False)
        self._cam2    = pd.read_parquet(cam2_fpath).set_index('flake_id', drop=False)
        self._triplet = pd.read_parquet(triplet_fpath).set_index('flake_id', drop=False)
        # - Ensure categorical/object columns are encoded as string 
        self._cam0 = _convert_object_to_string(self._cam0)
        self._cam1 = _convert_object_to_string(self._cam1)
        self._cam2 = _convert_object_to_string(self._cam2)
        self._triplet = _convert_object_to_string(self._triplet)
        
        # Number of triplets 
        self._n_triplets = len(self._triplet)
        
        # Save dir_path info 
        self._dir_path = dir_path
        
        # Add default events
        self._define_events(maximum_interval_without_images = np.timedelta64(4,'h'),
                            unit="ns")
        print(self._triplet)
    
    ####----------------------------------------------------------------------.
    #########################
    #### Builtins methods ###
    #########################
    def __len__(self):
        return self._n_triplets
    
  
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
        
    ####----------------------------------------------------------------------.
    #####################
    #### Write MASCDB ###
    ##################### 
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
        
    ####----------------------------------------------------------------------.
    ###################
    #### Subsetting ###
    ###################
    def isel(self, idx): 
        #---------------------------------------------------------------------.
        # Check valid (integer) idx 
        idx = _check_isel_idx(idx, vmax=self._n_triplets-1)
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        ### Subset all datasets 
        # - DataArray
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            self._da = self._da.isel(flake_id=idx)
        
        # - Dataframes
        # if isinstance(idx[0], bool):
        #     self._cam0 = self._cam0[idx]  
        #     self._cam1 = self._cam1[idx]  
        #     self._cam2 = self._cam2[idx]  
        #     self._triplet = self._triplet[idx]  
        # else: 
        self._cam0 = self._cam0.iloc[idx]
        self._cam1 = self._cam1.iloc[idx]
        self._cam2 = self._cam2.iloc[idx]
        self._triplet = self._triplet.iloc[idx]
        ##--------------------------------------------------------------------.
        # Update number of triplets 
        self._n_triplets = len(self._triplet)
        ##--------------------------------------------------------------------.
        return self 
        
    def sel(self, flake_ids): 
        #---------------------------------------------------------------------.
        # Check valid flake_ids 
        valid_flake_ids = self._da['flake_id'].values
        flake_ids = _check_sel_ids(flake_ids, valid_ids = valid_flake_ids)
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        ### Subset all datasets 
        # - DataArray
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            self._da = self._da.sel(flake_id=flake_ids)
        
        # - Dataframes
        self._cam0 = self._cam0.loc[flake_ids]
        self._cam1 = self._cam1.loc[flake_ids]
        self._cam2 = self._cam2.loc[flake_ids]
        self._triplet = self._triplet.loc[flake_ids]
        ##--------------------------------------------------------------------.
        # Update number of triplets 
        self._n_triplets = len(self._triplet)
        ##--------------------------------------------------------------------.
        return self 
    
    def sample_n(self, n=10):
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))      
        idx = np.random.choice(self._n_triplets, n) 
        return self.isel(idx)
        
    def first(self, n=1):  
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))
        idx = np.arange(n)
        return self.isel(idx)
    
    def last(self, n=1): 
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
     
    ####----------------------------------------------------------------------.
    ################ 
    #### Sorting ###
    ################
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
        idx = db[db_column].values.argsort()
        if decreasing: 
            idx = idx[::-1]
        #------------------------------.
        # Return sorted object
        return self.isel(idx)  
    
    def select_max(self, expression, n=10):
        return self.arrange(expression, decreasing=True).isel(np.arange(min(n, self._n_triplets)))

    def select_min(self, expression, n=10):
        return self.arrange(expression, decreasing=False).isel(np.arange(min(n, self._n_triplets)))   
    
    ####----------------------------------------------------------------------.
    ##########################
    #### Data explanation ####
    ##########################
    def get_var_units(self,varname):
        """
        Get units of a given variable

        Parameters
        ----------
        varname : str
            String containing the name (must be one of the columns of MASCDB dataframes)
            of a MASCDB variable
            
        Returns
        -------
        str
            Abbreviated units of the variable

        """
        if not isinstance(varname, str):
            raise TypeError("'varname' must be a string")
        units = var_units()
        if varname in units.keys():
            return units[varname]
        else:
            raise ValueError("{} units are not currently available. "
                             "Units are available for {}".format(varname, list(units.keys())))
    
    def get_var_explanation(self,varname):
        """
        Get verbose explanation of a given variable, including DOI of reference paper
        whenever relevant

        Parameters
        ----------
        varname : str
            String containing the name (must be one of the columns of MASCDB dataframes)
            of a MASCDB variable

        Returns
        -------
        str
            Verbose explanation of the variable
            
        """
        
        if not isinstance(varname, str):
            raise TypeError("'varname' must be a string")
        explanations = var_explanations()
        if varname in explanations.keys():
            return explanations[varname]
        else:
            raise ValueError("{} verbose explanation is not currently available. "
                             "Explanations are available for {}".format(varname, list(explanations.keys())))
    
    
    ####----------------------------------------------------------------------.
    #################
    #### Filters ####
    #################
    def select_campaign(self, campaign):
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
    
    def discard_campaign(self, campaign):
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
     
    def select_snowflake_class(self, values, method='Praz2017', invert = False): 
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
        idx = np.isin(arr, values, invert=invert)
        return self.isel(idx) 

    def select_riming_class(self, values, method='Praz2017', invert=False): 
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
        idx = np.isin(arr, values, invert=invert)
        return self.isel(idx)     
 
    def select_melting_class(self, values, method='Praz2017', invert=False): 
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
        idx = np.isin(arr, values, invert=invert)
        return self.isel(idx)  
    
    def select_precip_class(self, values, method='Schaer2020', invert=False): 
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_precip_class_name_dict().values()) # id
            column = 'bs_precip_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_precip_class_name_dict().keys())   # name
            column = 'bs_precip_class_name'
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
        idx = np.isin(arr, values, invert=invert)
        return self.isel(idx)  
       
    def discard_snowflake_class(self, values, method='Praz2017'):
        return self.select_snowflake_class(values=values, method=method, invert = True) 
    
    def discard_melting_class(self, values, method='Praz2017'):
        return self.select_melting_class(values=values, method=method, invert = True) 
    
    def discard_riming_class(self, values, method='Praz2017'):
        return self.select_riming_class(values=values, method=method, invert = True) 
      
    def discard_precip_class(self, values, method='Schaer2020'):
        return self.select_precip_class(values=values, method=method, invert = True)
    
    ####----------------------------------------------------------------------.
    #################
    #### Getters ####
    #################  
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
    def event(self):
        # columns = list(self._triplet.columns)
        # event_columns = [column for column in columns if column.startswith("event_")]
        event_columns = ['event_id', 'event_duration', 'event_n_triplets',
                         'campaign', 'datetime', 'latitude', 'longitude','altitude']
        event_db = self.triplet[event_columns].groupby('event_id').first().reset_index()
        event_db['start_time'] = self.triplet[["event_id","datetime"]].groupby('event_id').min()
        event_db['end_time'] = self.triplet[["event_id","datetime"]].groupby('event_id').max()
        event_db['month'] = event_db.datetime.dt.month
        event_db['year'] = event_db.datetime.dt.year
        _ = event_db.drop(columns="datetime", inplace=True)
        return event_db 
    
    ## BUG !!!
    # df_event = event_db
    # df_event[['start_time',c_id]].groupby(c_id).min()  
    # df_event1 = self.event
    
    @property
    def campaign(self): 
        #----------------------------------------------.
        # Retrieve data 
        df_event = self.event     
        df_triplet = self.triplet    
        c_id = 'campaign'
        #----------------------------------------------.
        # Compute location info  
        columns = ['latitude','longitude', 'altitude', 'campaign']
        info_location = df_event[columns].groupby('campaign').first().reset_index()
        
        #----------------------------------------------.
        # Compute number of triplets 
        # - n_triplet  (sum)
        n_triplets = df_triplet.groupby(c_id)['flake_id'].count()
        n_triplets.name = "n_triplets"
        
        #----------------------------------------------.
        ## Compute event summary
        n_events = df_event[['event_id',c_id]].groupby(c_id).count()
        n_events.columns = ["n_events"]
        start_time = df_event[['start_time',c_id]].groupby(c_id).min()  
        end_time = df_event[['end_time',c_id]].groupby(c_id).max()      
        event_duration_stats = df_event.groupby(c_id).agg({'event_duration': ['min','mean','max','sum']})
        event_duration_stats.columns = ["event_duration_min", "event_duration_mean", "event_duration_max", "total_event_duration"]
        
        #----------------------------------------------.
        ### Compute class occurence    
        snowflake_class_counts = df_triplet[['snowflake_class_name',c_id]].groupby(c_id)['snowflake_class_name'].apply(_count_occurence).apply(lambda x: x[0])     
        riming_class_counts = df_triplet[['riming_class_name',c_id]].groupby(c_id)['riming_class_name'].apply(_count_occurence).apply(lambda x: x[0])     
        melting_class_counts = df_triplet[['melting_class_id',c_id]].groupby(c_id)['melting_class_id'].apply(_count_occurence).apply(lambda x: x[0])     
        precipitation_class_counts = df_triplet[['bs_precip_class_name',c_id]].groupby(c_id)['bs_precip_class_name'].apply(_count_occurence).apply(lambda x: x[0])     
        snowflake_class_counts.name = "snowflake_class"
        riming_class_counts.name = "riming_class"
        melting_class_counts.name = "melting_class"
        precipitation_class_counts.name = "precipitation_class"
        
        ## TODO replace 
        # snowflake_class_counts = df_triplet[['snowflake_class_name',c_id]].groupby(c_id)['snowflake_class_name'].apply(_count_occurence).apply(lambda x: x[0]) 
        # riming_class_counts = df_triplet[['riming_class_name',c_id]].groupby(c_id)['riming_class_name'].apply(_count_occurence).apply(lambda x: x[0]) 
        # melting_class_counts = df_triplet[['melting_class_name',c_id]].groupby(c_id)['melting_class_name'].apply(_count_occurence).apply(lambda x: x[0]) 
        # precipitation_class_counts = df_triplet[['bs_precipitation_class',c_id]].groupby(c_id)['bs_precipitation_class'].apply(_count_occurence).apply(lambda x: x[0]) 
        # snowflake_class_counts.name = "snowflake_class"
        # riming_class_counts.name = "riming_class"
        # melting_class_counts.name = "melting_class"
        # precipitation_class_counts.name = "precipitation_class"
        
        #----------------------------------------------.
        ## Compute other time infos 
        # months_list = df_triplet.groupby(c_id)['datetime'].apply(lambda x: list(np.unique(x.dt.month_name())))
        # months_list.name = "months"
        # years_list = df_triplet.groupby(c_id)['datetime'].apply(lambda x: list(np.unique(x.dt.year)))
        # years_list.name = "years"
        # years_months_list = df_triplet.groupby(c_id)['datetime'].apply(lambda x: list(np.unique(x.dt.strftime('%Y-%m'))))
        # years_months_list.name = "years_months"
        
        #----------------------------------------------.
        # Define summary dataframe
        summary = pd.merge(start_time, end_time, right_index=True, left_index=True)
        summary = summary.join(info_location)
        summary = summary.join(n_triplets)
        summary = summary.join(n_events)
        summary = summary.join(event_duration_stats)
        summary = summary.join(snowflake_class_counts)
        summary = summary.join(riming_class_counts)
        summary = summary.join(melting_class_counts)
        summary = summary.join(precipitation_class_counts)
        #----------------------------------------------.
        return summary 

    @property
    def full_db(self):
        # TODO: check same order as ds_images ... maybe add cam_id  and campaign args 
        # Add cam_id to each cam db 
        l_cams = [self.cam0, self.cam1, self.cam2]
        for i, cam in enumerate(l_cams):
             cam['cam_id'] = i
        # Merge cam(s) db into 
        full_db = pd.concat(l_cams)
        # Add triplet variables to fulldb 
        labels_vars = get_vars_class()
        vars_not_add = ['flake_quality_xhi','flake_n_roi', 'flake_Dmax','flake_id'] + labels_vars
        triplet = self.triplet.drop(columns=vars_not_add)
        full_db = full_db.merge(triplet, how="left")
        return full_db
    
    def ds_images(self, cam_id = None, campaign=None, img_id='img_id'):
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
        if cam_id is not None:
            da = da.isel(cam_id = cam_id)
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
    
    ####----------------------------------------------------------------------.
    ###################### 
    #### Event utils #####
    ###################### 
    def _add_event_n_triplets(self):
        if 'event_id' not in list(self._triplet.columns):
            raise ValueError("First define 'event_id' using mascdb.define_event_id().")
        event_triplets = self._triplet.groupby('event_id').size()
        event_triplets.name = "event_n_triplets"
        # - Remove existing column to remerge 
        if "event_n_triplets" in list(self._triplet.columns):
            _ = self._triplet.drop(columns="event_n_triplets", inplace=True)
        # - Add column 
        self._triplet = self._triplet.merge(event_triplets, on="event_id").set_index('flake_id', drop=False)
        return None
        
    def _add_event_duration(self, unit="ns"):
        if 'event_id' not in list(self._triplet.columns):
            raise ValueError("First define 'event_id' using mascdb.define_event_id().")
        event_durations = self._triplet.groupby('event_id').apply(lambda x: _get_timesteps_duration(x.datetime, unit=unit))
        event_durations.name = "event_duration"
        # - Remove existing column to remerge
        if "event_duration" in list(self._triplet.columns):
            _ = self._triplet.drop(columns="event_duration", inplace=True)
        # - Add column 
        self._triplet = self._triplet.merge(event_durations, on="event_id").set_index('flake_id', drop=False)
        return None
    
    def _define_events(self, 
                      maximum_interval_without_images = np.timedelta64(4,'h'),
                      unit="ns"):   
        # This function modify in place       
        #----------------------------------------------------------.
        # - Extract relevant columns from triplet db
        db = self.triplet[['campaign','datetime']]
        # - Retrieve campaign_ids 
        campaign_ids = np.unique(db['campaign'])
        #----------------------------------------------------------.
        # - Retrieve event_id column 
        db['event_id'] = -1
        max_event_id = 0
        for campaign_id in campaign_ids:
            # - Retrieve row index of specific campaign 
            idx_campaign = db['campaign'] == campaign_id
            # - Define event_ids for the campaign 
            campaign_event_ids = _define_event_id(timesteps = db.loc[idx_campaign,'datetime'], 
                                                  maximum_interval_without_timesteps = maximum_interval_without_images)
            # - Add offset to ensure having an unique event_id across all campaigns 
            campaign_event_ids = campaign_event_ids + max_event_id  
            # - Add event_id to the campaign subset of the database 
            db.loc[idx_campaign, 'event_id'] = campaign_event_ids
            # - Update the current maximum event_id
            max_event_id = max(campaign_event_ids) + 1
        #----------------------------------------------------------.
        # - Add event_id column to all databases 
        self._cam0['event_id'] = db['event_id']
        self._cam1['event_id'] = db['event_id']
        self._cam2['event_id'] = db['event_id']
        self._triplet['event_id'] = db['event_id']
        #----------------------------------------------------------.
        # - Add duration and n_triplets for each event to triplet db 
        self._add_event_duration(unit=unit)
        self._add_event_n_triplets()
    
    ##--------------------------------------------------
    ## Filtering events utils 
    def select_events_with_more_triplets_than(self, n):
        # Note: also used in define_events
        if not isinstance(n, int): 
            raise TypeError('Expects an integer')
        if n < 1: 
            raise ValueError("'n' must be equal or larger than 1.")
        # If n=1, nothing to subset 
        if n == 1:
            return self
        # Else subset 
        df_event = self.event 
        idx_event_ids = df_event['event_n_triplets'] >= n
        subset_event_ids = df_event.loc[idx_event_ids, 'event_id'].values  
        idx_db_subset = np.isin(self._triplet['event_id'].values, subset_event_ids)
        return self.isel(idx_db_subset)
    
    def select_events_with_less_triplets_than(self, n):
        if not isinstance(n, int): 
            raise TypeError('Expects an integer')
        if n < 1: 
            raise ValueError("'n' must be equal or larger than 1.")
        df_event = self.event 
        idx_event_ids = df_event['event_n_triplets'] <= n
        subset_event_ids = df_event.loc[idx_event_ids, 'event_id'].values  
        idx_subset = np.isin(self._triplet['event_id'].values, subset_event_ids)
        return self.isel(idx_subset)
    
    def discard_events_with_more_triplets_than(self, n):
        return self.select_events_with_less_triplets_than(n) 
    
    def discard_events_with_less_triplets_than(self, n):
        return self.select_events_with_more_triplets_than(n) 
        
    def select_events_with_min_duration(self, timedelta):
        timedelta = _check_timedelta(timedelta)
        df_event = self.event
        idx_event_ids = df_event['event_duration'] >= timedelta
        subset_event_ids = df_event.loc[idx_event_ids, 'event_id'].values  
        idx_subset = np.isin(self._triplet['event_id'].values, subset_event_ids)
        return self.isel(idx_subset)
    
        # TODO: BUG !!!
        # df_event[['start_time',c_id]].groupby(c_id).min()  
        # df_event1 = self.event
        # df_event1[['start_time',c_id]].groupby(c_id).min() 
        # a = self.isel(idx_subset)
        # df_event2 = a.event
        # df_event2[['start_time',c_id]].groupby(c_id).min() 
        
    def select_events_with_max_duration(self, timedelta):
        timedelta = _check_timedelta(timedelta)
        df_event = self.event
        idx_event_ids = df_event['event_duration'] <= timedelta
        subset_event_ids = df_event.loc[idx_event_ids, 'event_id'].values  
        idx_subset = np.isin(self._triplet['event_id'].values, subset_event_ids)
        return self.isel(idx_subset)
    
    def select_longest_events(self, n=1):
        longest_event_ids = self.arrange('triplet.event_duration', decreasing=True)._triplet['event_id'].iloc[0:n]
        idx_longest_events = np.isin(self._triplet['event_id'].values, longest_event_ids)
        return self.isel(idx_longest_events)
    
    def select_shortest_events(self, n=1):
        shortest_event_ids = self.arrange('triplet.event_duration', decreasing=False)._triplet['event_id'].iloc[0:n]
        idx_shortest_events = np.isin(self._triplet['event_id'].values, shortest_event_ids)
        return self.isel(idx_shortest_events)
        
    def discard_events_with_max_duration(self, timedelta):
        return self.select_events_with_min_duration(timedelta) 
    
    def discard_events_with_min_duration(self, timedelta):
        return self.select_events_with_max_duration(timedelta)
    
    ##--------------------------------------------------
    ## Redefine events utils  
    def redefine_events(self, 
                        maximum_interval_without_images = np.timedelta64(4,'h'),
                        minimum_duration = None,
                        minimum_n_triplets = None,
                        unit="ns"): 
        # Copy new instance 
        self = copy.deepcopy(self)
        # Define event_id 
        self._define_events(maximum_interval_without_images=maximum_interval_without_images)
        #----------------------------------------------------------.
        # Select only events with min_n_triplets and minimum_event_duration
        if minimum_n_triplets is not None: 
            self = self.select_events_with_more_triplets_than(minimum_n_triplets)
        if minimum_duration is not None: 
            self = self.select_events_with_min_duration(minimum_duration)
        #----------------------------------------------------------.
        # Return the object 
        return self 
        
    ####------------------------------------------------------------------------.
    #################################
    #### Image plotting routines ####
    #################################

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
        da_subset = self._da.isel(flake_id = indices).transpose(...,'cam_id','flake_id')
        
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
        row = "flake_id" if len(indices) > 1 else None
        p = da_subset.plot(x='x',y='y',
                           col="cam_id",
                           row=row,
                           aspect=1,
                           yincrease=False,
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
            
    def plot_flake(self, cam_id=None, index=None, random = False,
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
        if cam_id is None: 
            if random:    
               cam_id = np.random.choice([0,1,2], 1)
            else:
               cam_id = 1
        #--------------------------------------------------.       
        # Check validty of cam_id and index 
        cam_id = _check_cam_id(cam_id)
        index = _check_index(index, vmax=n_idxs-1)
        
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If cam_id is an integer (instead of list of length 1), then the cam_id dimension is dropped)
        da_img = self._da.isel(flake_id = index, cam_id = cam_id) 
        
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
                        yincrease=False,
                        cmap='gray', add_colorbar=False, 
                        vmin=0, vmax=255,
                        **kwargs)
        #--------------------------------------------------. 
        return p  

    def plot_flakes(self, cam_id=None, indices=None, random = False, 
                    n_images = 9, col_wrap = 3,
                    enhancement="histogram_equalization",
                    zoom=True, 
                    **kwargs):
        # Retrieve number of valid index
        n_idxs = len(self._triplet.index) # TODO 
        # TODO: 
        # - Option to stack to ImageID ... and then sample that ... title will report flake_id and CAM ID

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
        if cam_id is None: 
            if random:    
                cam_id = np.random.choice([0,1,2], 1)[0]
            else:
                cam_id = 1
        #-------------------------------------------------- 
        # Check validity of indices and cam_id
        indices = _check_indices(indices, vmax=n_idxs-1)
        cam_id = _check_cam_id(cam_id)
        #-------------------------------------------------- 
        # If a single flake is specified, plot it with plot_flake 
        if len(indices) == 1: 
            print("It's recommended to use 'plot_flake()' to plot a single image.")
            return self.plot_flake(index=indices[0], cam_id=cam_id, random=random,
                                   enhancement=enhancement, zoom=zoom, *kwargs)
        
        #--------------------------------------------------.
        # Subset triplet(s) images 
        # - If cam_id is an integer (instead of list length 1 ... the cam_id dimension is dropped)
        da_subset = self._da.isel(flake_id = indices, cam_id = cam_id).transpose(...,'flake_id')
        
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
        row = "flake_id" if len(indices) > 1 else None
        p = da_subset.plot(x='x',y='y', 
                           row=row, col_wrap=col_wrap, 
                           aspect=1, 
                           yincrease=False,
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
     
    ####-----------------------------------------------------------------------.
    ####################### 
    #### MASCDB Updates ###
    ####################### 
    def compute_2Dimage_descriptors(self, fun, labels, fun_kwargs = None, force=False,
                                    dask = "parallelized"):
        #---------------------------------------------------------------------.
        # Check if specified labels are already columns of mascdb.cam* 
        existing_cam_columns = list(self._cam0.columns)
        overwrited_columns = list(np.array(existing_cam_columns)[np.isin(existing_cam_columns, labels)])
        if not force: 
            if len(overwrited_columns) > 0:
                raise ValueError("Columns {} would be overwrited. "
                                 "Specify force=True if you want to overwrite existing columns.".format(overwrited_columns))
                
        #---------------------------------------------------------------------.
        # Compute descriptors 
        da_descriptors = _compute_2Dimage_descriptors(da = self.da, 
                                                      fun = fun,
                                                      labels = labels, 
                                                      fun_kwargs = fun_kwargs,
                                                      dask = dask)
        #---------------------------------------------------------------------.
        # Retrieve cam dataframes 
        cam0 = da_descriptors.isel(cam_id = 0).to_dataset('descriptor').to_pandas().drop(columns='cam_id')
        cam1 = da_descriptors.isel(cam_id = 1).to_dataset('descriptor').to_pandas().drop(columns='cam_id')
        cam2 = da_descriptors.isel(cam_id = 2).to_dataset('descriptor').to_pandas().drop(columns='cam_id')
        
        #---------------------------------------------------------------------.
        # Attach to new mascdb instance
        new_mascdb = self.add_cam_columns(cam0=cam0, cam1=cam1, cam2=cam2, force=force, complete=True)
        # Return new mascdb instance
        return new_mascdb 


    def add_cam_columns(self, cam0, cam1, cam2, force=False, complete=True): 
        """
        Method allowing to safely add columns to cam dataframes of MASCDB.
        
        Parameters
        ----------
        cam0 : pd.DataFrame
            pd.DataFrame with index 'flake_id' .
        cam1 : pd.DataFrame
            pd.DataFrame with index 'flake_id' .
        cam2 : pd.DataFrame
            pd.DataFrame with index 'flake_id' .
        force : bool, optional
            Wheter to overwrite existing column of mascdb. The default is False.
        complete : vool, optional
            Wheter to merge only when the cam dataframes have same 'flake_id' of 
            the current mascdb. The default is True.
    
        Returns
        -------
        MASCDB class instance
        
        """
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # Check all cam* are pd.DataFrame 
        cam0 = _check_df(cam0, name='cam0')
        cam1 = _check_df(cam1, name='cam1')
        cam2 = _check_df(cam2, name='cam2')
        
        # Check length is the same across all dataframes 
        n_cam0 = len(cam0)
        n_cam1 = len(cam1)
        n_cam2 = len(cam2) 
        if not n_cam0 == n_cam1: 
            raise ValueError("cam0 has {} rows, while cam1 has {} rows.".format(n_cam0, n_cam1))
        if not n_cam0 == n_cam2: 
            raise ValueError("cam0 has {} rows, while cam2 has {} rows.".format(n_cam0, n_cam2))
          
        # Check columns are the same across all dataframes
        cam0_columns = np.sort(list(cam0.columns))
        cam1_columns = np.sort(list(cam1.columns))
        cam2_columns = np.sort(list(cam2.columns))
        if not np.array_equal(cam0_columns, cam1_columns): 
            raise ValueError("cam0 and cam1 does not have the same column names.")
        if not np.array_equal(cam0_columns, cam2_columns): 
            raise ValueError("cam0 and cam2 does not have the same column names.")
        
        #---------------------------------------------------------------------.
        ### - Check cam flake_id match each others
        cam0_flake_ids = np.sort(cam0.index.values)
        cam1_flake_ids = np.sort(cam0.index.values)
        cam2_flake_ids = np.sort(cam0.index.values)
        if not np.array_equal(cam0_flake_ids, cam1_flake_ids): 
            raise ValueError("cam0 and cam1 does not have the same 'flake_id' index.")
        if not np.array_equal(cam0_flake_ids, cam2_flake_ids): 
            raise ValueError("cam0 and cam2 does not have the same 'flake_id' index.")
        
        #---------------------------------------------------------------------.
        # Check if column names already exist in mascdb.cam* 
        existing_cam_columns = list(self._cam0.columns)
        overwrited_columns = list(np.array(existing_cam_columns)[np.isin(existing_cam_columns, cam0_columns)])
        if not force: 
            if len(overwrited_columns) > 0:
                raise ValueError("Columns {} would be overwrited. "
                                 "Specify force=True if you want to overwrite existing columns.".format(overwrited_columns))
        
        #---------------------------------------------------------------------.
        # Drop columns that must be overwrited 
        if len(overwrited_columns) > 0:
            _ = self._cam0.drop(columns=overwrited_columns, inplace=True)
            _ = self._cam1.drop(columns=overwrited_columns, inplace=True)
            _ = self._cam2.drop(columns=overwrited_columns, inplace=True)
            
        #---------------------------------------------------------------------.
        # Ensure columns order is the same across all dataframes 
        cam0 = cam0[cam0_columns]
        cam1 = cam1[cam1_columns]
        cam1 = cam2[cam2_columns]
           
        #---------------------------------------------------------------------.
        ### - Check flake_id match between mascdb and provided cam dataframes 
        new_flake_ids = cam0_flake_ids
        existing_flake_ids = self._cam0['flake_id'].values
        
        missing_flake_ids = existing_flake_ids[np.isin(existing_flake_ids, new_flake_ids, invert=True)]
        matching_flake_ids = new_flake_ids[np.isin(new_flake_ids, existing_flake_ids)]
        non_matching_flake_ids = new_flake_ids[np.isin(new_flake_ids, existing_flake_ids, invert=True)]
       
        #----------------------------------------------------------------------. 
        # Check at least 1 flake_id match 
        if len(matching_flake_ids) == 0: 
            raise ValueError("No matching 'flake_id' between current mascdb and provided cam dataframes.")
        
        #---------------------------------------------------------------------. 
        # Check flake_id index and number of rows of new cam correspond to existing one 
        if complete: 
            # Check that there are the same flake_id 
            if len(missing_flake_ids) > 0: 
                msg = ("There are {} flake_id missing in the provided cam dataframes. \n " 
                       "If you want to still merge the new columns, specify complete=False. \n "  
                      "New columns with non-matching rows will be filled by NaN.".format(len(missing_flake_ids)))
                raise ValueError(msg)
            # Check number of rows 
            if self._n_triplets != n_cam0: 
                msg = ("The provided cam dataframes have {} rows, while "  
                      "the current mascdb has {} rows. \n "  
                      "If you want to still merge the new columns, specify complete=False. \n "  
                      "New columns with non-matching rows will be filled by NaN.".format(n_cam0, self._n_triplets))
                raise ValueError(msg)
                
        #---------------------------------------------------------------------. 
        # Print a message if some flake_id does not have a match         
        if len(non_matching_flake_ids) > 0:
            msg = ("There are {} flake_id in the provided cam dataframes which "
                   "will not be merged to the mascdb because of non-matching flake_id.".format(len(non_matching_flake_ids)))
            print(msg)
            
        #---------------------------------------------------------------------.
        # Join data     
        self._cam0 = self._cam0.merge(cam0, how="left", on='flake_id').set_index('flake_id', drop=False)
        self._cam1 = self._cam1.merge(cam1, how="left", on='flake_id').set_index('flake_id', drop=False)
        self._cam2 = self._cam2.merge(cam2, how="left", on='flake_id').set_index('flake_id', drop=False)
        
        #---------------------------------------------------------------------.
        # Return the new mascdb 
        return self 
    
    def add_triplet_columns(self, df, force=False, complete=True): 
        """
        Method allowing to safely add columns to cam dataframes of MASCDB.
        
        Parameters
        ----------
        df : pd.DataFrame
            pd.DataFrame with index 'flake_id' .
        force : bool, optional
            Wheter to overwrite existing column of mascdb. The default is False.
        complete : vool, optional
            Wheter to merge only when the provided dataframe has the same 'flake_id' of 
            the current mascdb triplet dataframe. The default is True.
    
        Returns
        -------
        MASCDB class instance
        
        """
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # Check df is pd.DataFrame 
        df = _check_df(df, name='df')
            
        # Check length is the same across all dataframes 
        n_df = len(df)
        
        # Check columns are the same across all dataframes
        df_columns = list(df.columns)
       
        #---------------------------------------------------------------------.
        ### - Check df flake_id match each others
        df_flake_ids = np.sort(df.index.values)
                
        #---------------------------------------------------------------------.
        # Check if column names already exist in mascdb.triplet 
        existing_triplet_columns = list(self._triplet.columns)
        overwrited_columns = list(np.array(existing_triplet_columns)[np.isin(existing_triplet_columns, df_columns)])
        if not force: 
            if len(overwrited_columns) > 0:
                raise ValueError("Columns {} would be overwrited. "
                                 "Specify force=True if you want to overwrite existing columns.".format(overwrited_columns))
        
        #---------------------------------------------------------------------.
        # Drop columns that must be overwrited 
        if len(overwrited_columns) > 0:
            _ = self._triplet.drop(columns=overwrited_columns, inplace=True)
                            
        #---------------------------------------------------------------------.
        ### - Check flake_id match between mascdb and provided dataframes 
        new_flake_ids = df_flake_ids
        existing_flake_ids = self._triplet['flake_id'].values
        
        missing_flake_ids = existing_flake_ids[np.isin(existing_flake_ids, new_flake_ids, invert=True)]
        matching_flake_ids = new_flake_ids[np.isin(new_flake_ids, existing_flake_ids)]
        non_matching_flake_ids = new_flake_ids[np.isin(new_flake_ids, existing_flake_ids, invert=True)]
       
        #---------------------------------------------------------------------. 
        # Check at least 1 flake_id match 
        if len(matching_flake_ids) == 0: 
            raise ValueError("No matching 'flake_id' between current mascdb and provided dataframe.")
        
        #---------------------------------------------------------------------. 
        # Check flake_id index and number of rows of new cam correspond to existing one 
        if complete: 
            # Check that there are the same flake_id 
            if len(missing_flake_ids) > 0: 
                msg = ("There are {} flake_id missing in the provided dataframe. \n " 
                       "If you want to still merge the new columns, specify complete=False. \n "  
                      "New columns with non-matching rows will be filled by NaN.".format(len(missing_flake_ids)))
                raise ValueError(msg)
            # Check number of rows 
            if self._n_triplets != n_df: 
                msg = ("The provided dataframe have {} rows, while "  
                      "the current mascdb has {} rows. \n "  
                      "If you want to still merge the new columns, specify complete=False. \n "  
                      "New columns with non-matching rows will be filled by NaN.".format(n_df, self._n_triplets))
                raise ValueError(msg)
                
        #---------------------------------------------------------------------. 
        # Print a message if some flake_id does not have a match         
        if len(non_matching_flake_ids) > 0:
            msg = ("There are {} flake_id in the provided df dataframe which "
                   "will not be merged to the mascdb because of non-matching flake_id.".format(len(non_matching_flake_ids)))
            print(msg)
        #---------------------------------------------------------------------.
        # Join data     
        self._triplet = self._triplet.merge(df, how="left", on='flake_id').set_index('flake_id', drop=False)
    
        #---------------------------------------------------------------------.
        # Return the new mascdb 
        return self 
      
    
    
    
    