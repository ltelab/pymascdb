#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:56:51 2021

TODO: 
    
- Add description of instance variables in the header of 
 MASCdb (for Sphinx)
    
- Add to the docs that mascdb.events and mascdb.campaigns
  should not be used to derive idxs for mascdb.isel() JGR
  
- Maybe modify something in the event code so that when 
  mascdb.redefine_events() it's called with min_duration or
  min_n_triplets arguments ... if some event_ids are removed,
  we redefine event_id to range from 0 to the number of events.

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
  
####--------------------------------------------------------------------------.
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
    if not isinstance(df, (pd.DataFrame, pd.Series)):
        if name is not None: 
            raise TypeError("{} is not a pd.DataFrame or pd.Series".format(name))
        else: 
           raise TypeError("Expecting a pd.DataFrame or pd.Series") 
    if isinstance(df, pd.Series): 
        df = df.to_frame()
    return df 

def _check_columns(columns):
    if not isinstance(columns, (str, list, np.ndarray)):
        raise TypeError("'columns' must be a string or list/np.array of strings.")
    if isinstance(columns, str): 
        columns = [columns]
    if isinstance(columns, list): 
        columns = np.array(columns)
    if isinstance(columns, np.ndarray):
        if columns.dtype.name == 'object':
            columns = columns.astype(str)
        if not columns.dtype.name.startswith('str'):
            raise ValueError("Expecting columns in the np.array to be strings.")
    return columns

def _check_df_source(df_source):
    if not isinstance(df_source, str): 
        raise TypeError("'df_source' must be a string. Either 'triplet', 'cam0','cam1', 'cam2'.")
    valid_source = ['triplet', 'cam0', 'cam1', 'cam2']
    if not df_source in valid_source:
        raise ValueError("Valid 'df_source' values are {}".format(valid_source))

def _get_df_values(self, column, df_source='triplet'):
    if df_source == 'triplet':
        return self._triplet[column].to_numpy()
    elif df_source == 'cam0':
        return self._cam0[column].to_numpy()
    elif df_source == 'cam1':
        return self._cam1[column].to_numpy()
    elif df_source == 'cam2':
        return self._cam2[column].to_numpy()
    else:
        raise ValueError("Unvalid 'source'.")
        
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
    Read MASCDB database from a specific directory.
    
    Parameters
    ----------
    dir_path : str
        Filepath to a directory storing a MASCDB.
        5 files are expected in the directory:
          - MASCdb_cam<0/1/2>.parquet   
          - MASCdb_triplet.parquet 
          - MASCdb.zarr
    
    Returns
    -------
    MASCDB class instance.
    
    """
    
    #####################
    #### Read MASCDB ###
    ##################### 
    def __init__(self, dir_path):
        """
        Initialize MASC_DB object.
        It reads 4 parquet databases as well as the zarr database of MASC greyscale images.

        Returns
        -------
        MASCDB class instance.

        """
        zarr_store_fpath = os.path.join(dir_path,"MASCdb.zarr")
        cam0_fpath = os.path.join(dir_path, "MASCdb_cam0.parquet")
        cam1_fpath = os.path.join(dir_path, "MASCdb_cam1.parquet")
        cam2_fpath = os.path.join(dir_path, "MASCdb_cam2.parquet")
        triplet_fpath = os.path.join(dir_path, "MASCdb_triplet.parquet")
        
        # - Check if the Zarr DirectoryStore has not been unzipped
        if not os.path.exists(zarr_store_fpath):
            zarr_zipstore_fpath = zarr_store_fpath + ".zip"
            if os.path.exists(zarr_zipstore_fpath):
                raise ValueError("You need to unzip {}".format(zarr_zipstore_fpath))

        # - Read image dataset 
        self._da = xr.open_zarr(zarr_store_fpath)['data']
        self._da['flake_id'] = self._da['flake_id'].astype(str)
        self._da.name = "MASC Images"
  
        # - Read dataframes
        self._cam0    = pd.read_parquet(cam0_fpath).set_index('flake_id', drop=False)
        self._cam1    = pd.read_parquet(cam1_fpath).set_index('flake_id', drop=False)
        self._cam2    = pd.read_parquet(cam2_fpath).set_index('flake_id', drop=False)
        self._triplet = pd.read_parquet(triplet_fpath).set_index('flake_id', drop=False)
       
        # - Ensure categorical/object columns are encoded as string 
        self._cam0 = _convert_object_to_string(self._cam0)
        self._cam1 = _convert_object_to_string(self._cam1)
        self._cam2 = _convert_object_to_string(self._cam2)
        self._triplet = _convert_object_to_string(self._triplet)
        
        #------------   
        # - Temporary solution for timedelta issue in parquets/arrow
        if "event_duration" in list(self._triplet.columns):
            self._triplet["event_duration"] = self._triplet["event_duration"].astype("timedelta64[ns]") # astype("m8[ns]")
             
        #------------    
        # - Define number of triplets 
        self._n_triplets = len(self._triplet)
        
        # - Save source MASCDB directory 
        self._dir_path = dir_path
        
        # - Add default events
        self._define_events(max_interval_without_images = np.timedelta64(4,'h'), unit="ns")
    
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
        """
        Save MASCDB object to disk into 4 parquet files and one Zarr store.

        Parameters
        ----------
        dir_path : str
            Directory path where to save the current MASCDB database.
        force :  Bool, default False
           If dir_path is the same as the source path of MASCDB object, 
           force=True should allows to overwrite the original source database.
           
        """
        # TODO: if overwriting, put all DataArray in memory first... otherwise deleting on disk remove lazy loaded data
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
        
        # - Ensure "correct" chunks of DataArray 
        da = self.da
        new_chunks = [max(chunk) for chunk in da.chunks]
        da = da.chunk(new_chunks)
        
        # - Write databases
        #------------
        # Temporary solution because timedelta cannot be saved 
        #   currently to parquet: https://issues.apache.org/jira/browse/ARROW-6780
        # - event_duration timedelta is converted to int 
        # - It assume no other timedelta columns are present in dataframes
        triplet = self.triplet 
        triplet["event_duration"] = triplet["event_duration"].astype("timedelta64[ns]").view(int)
        #------------        
        ds = da.to_dataset(name="data") 
        ds.to_zarr(zarr_store_fpath)
        self._cam0.to_parquet(cam0_fpath, engine="auto")
        self._cam1.to_parquet(cam1_fpath, engine="auto")
        self._cam2.to_parquet(cam2_fpath, engine="auto")
        triplet.to_parquet(triplet_fpath, engine="auto")
        #---------------------------------------------------------------------.
        return None 
        
    ####----------------------------------------------------------------------.
    ###################
    #### Subsetting ###
    ###################
    def isel(self, idx): 
        """
        Positional-index subsetting of MASCDB DataArray and MASCDB DataFrames.

        Parameters
        ----------
        idx : (np.ndarray, list, int)
            List or np.ndarray of integer/boolean values used as positional indices for subsetting.

        Returns
        -------
        MASCDB class instance subsetted (or index-based reordered).

        """
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
        """
        Subset MASCDB based on specified flake_ids.

        Parameters
        ----------
        flake_ids : (np.ndarray, list, str)
            List or np.ndarray of string specifing flake_id values to subset.

        Returns
        -------
        MASCDB class instance subsetted.

        """
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
        """
        Sample randomly 'n' flakes in the current MASC_DB object.

        Parameters
        ----------
        n : int,float; optional
             Number of samples to extract The default is 10.


        Returns
        -------
        MASC_DB object with n sampled flakes.

        """
        
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))      
        idx = np.random.choice(self._n_triplets, n) 
        return self.isel(idx)
        
    def first(self, n=1):  
        """
        Extract first 'n' flakes in the database.

        Parameters
        ----------
        n : int,float; optional
             Number of samples to extract The default is 1.


        Returns
        -------
        MASC_DB object containing only the n first flakes of the current database

        """
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))
        idx = np.arange(n)
        return self.isel(idx)
    
    def last(self, n=1): 
        """
        Extract last 'n' flakes in the database.

        Parameters
        ----------
        n : int,float; optional
             Number of samples to extract The default is 1


        Returns
        -------
        MASC_DB object containing only the n last flakes of the current database.

        """
        if n > len(self): 
            raise ValueError("The MASCDB instance has currently only {} triplets.".format(len(self)))
        idx = np.arange(self._n_triplets-1,self._n_triplets-n-1, step=-1)
        return self.isel(idx)
    
    def head(self, n=10): 
        """
        Extract first 'n' flakes in the database or less if the database contains
        less rows than n.

        Parameters
        ----------
        n : int,float; optional
             Number of samples to extract The default is 10.


        Returns
        -------
        MASC_DB object containing only the n first flakes of the current database.

        """
        
        n = min(self._n_triplets, n)
        idx = np.arange(n)
        return self.isel(idx)
    
    def tail(self, n=10):
        """
        Extract last 'n' flakes in the database or less if the database contains.
        less rows than n.

        Parameters
        ----------
        n : int,float; optional
             Number of samples to extract The default is 10.


        Returns
        -------
        MASC_DB object containing only the n last flakes of the current database.

        """
        n = min(self._n_triplets, n)
        idx = np.arange(self._n_triplets-1,self._n_triplets-n-1, step=-1)
        return self.isel(idx)
     
    ####----------------------------------------------------------------------.
    ################ 
    #### Sorting ###
    ################
    def arrange(self, expression, decreasing=True):
        """
        Reorder the MASCDB based on the DataFrame column values specified with expression.

        Parameters
        ----------
        expression : str
            Expression specifying the DataFrame and column used to sort the MASCDB.
            The expression must have the following pattern '<df_name>.<column_name>' .
            Valid df_names are : ['cam0', 'cam1','cam2','triplet','bs','env','gan3d','flake','labels'] .
        decreasing : bool, optional
            Whether to sort MASCDB by increasing or decreasing values of the DataFrame colum.
            The default is True.

        Returns
        -------
        MASCDB class instance sorted.

        """
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
        valid_db = ['cam0', 'cam1','cam2','triplet','bs','env','gan3d','flake','labels']
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
    
    ####----------------------------------------------------------------------.
    ##########################
    #### Data explanation ####
    ##########################
    def get_var_units(self,varname):
        """
        Get units of a given variable.

        Parameters
        ----------
        varname : str
            String specifying a single column of the cam or triplet dataframe.
            
        Returns
        -------
        str
            Abbreviated units of the variable.

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
        Get verbose explanation of a given variable.
        
        It includes DOI of reference paper whenever relevant.

        Parameters
        ----------
        varname : str
            String specifying a single column of the cam or triplet dataframe.

        Returns
        -------
        str
            Verbose explanation of the variable.
            
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
        """
        Select MASCDB data of specific campaigns.
        
        Parameters
        ----------
        campaign : (str, list)
            String or list of string specifying MASCDB campaigns to select.
        
        Returns
        -------
        MASCDB class instance with data of specific campaigns.
        
        """
        if not isinstance(campaign, (list, str)): 
            raise TypeError("'campaign' must be a string or a list of strings.")
        if isinstance(campaign, str): 
            campaign = [campaign]
        # Convert to numpy array with str type (not object...)
        campaign = np.array(campaign).astype(str)
        campaigns_arr = self._triplet['campaign'].to_numpy().astype(str)
        valid_campaigns = np.unique(campaigns_arr)
        unvalid_campaigns_arg = campaign[np.isin(campaign, valid_campaigns, invert=True)]
        if len(unvalid_campaigns_arg) > 0: 
            raise ValueError("{} is not a campaign of the current mascdb. "
                             "Valid campaign names are {}".format(unvalid_campaigns_arg.tolist(),
                                                                  valid_campaigns.tolist()))
        idx = np.isin(campaigns_arr, campaign)
        return self.isel(idx) 
    
    def discard_campaign(self, campaign):
        """
        Discard MASCDB data from specific campaigns.

        Parameters
        ----------
        campaign : (str, list)
            String or list of string specifying MASCDB campaigns to discard.

        Returns
        -------
        MASCDB class instance with data of specific campaigns.
    
        """
        if not isinstance(campaign, (list, str)): 
            raise TypeError("'campaign' must be a string or a list of strings.")
        if isinstance(campaign, str): 
            campaign = [campaign]
        # Convert to numpy array with str type (not object...)
        campaign = np.array(campaign).astype(str)
        campaigns_arr = self._triplet['campaign'].to_numpy().astype(str)
        valid_campaigns = np.unique(campaigns_arr)
        unvalid_campaigns_arg = campaign[np.isin(campaign, valid_campaigns, invert=True)]
        if len(unvalid_campaigns_arg) > 0: 
            raise ValueError("{} is already not a campaign of the current mascdb. "
                             "Current mascdb has campaign names {}".format(unvalid_campaigns_arg.tolist(),
                                                                           valid_campaigns.tolist()))
        idx = np.isin(campaigns_arr, campaign, invert=True)
        return self.isel(idx) 
     
    def select_snowflake_class(self, values, method='Praz2017', invert = False, df_source="triplet"): 
        """
        Select MASCDB data with specific snowflake classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the snowflake classes to select.
            If integers, it assumes snowflake_class_id.
            If strings, it assumes snowflake_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_snowflake_class_name_dict(method)'.
        method : str, optional
            Method used to determine snowflake_class. The default is 'Praz2017'.
        invert : bool, optional
            If True, instead of selecting it discard the specified snowflake_class.
            The default is False.
        df_source: str, optional 
            The dataframe from which retrieve the class. 
            Either 'cam0', 'cam1', 'cam2' or 'triplet'.
            The  default is 'triplet'.

        Returns
        -------
        MASCDB class instance with specific snowflake classes.
    
        """
        #---------------------------------------------------------------------.
        ## Check default args 
        if not isinstance(invert, bool):
            raise TypeError("'invert' must be either True or False'.")
        _check_df_source(df_source)
        #---------------------------------------------------------------------.
        ## Check values 
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_snowflake_class_name_dict(method=method).values()) # id
            column = 'snowflake_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_snowflake_class_name_dict(method=method).keys())   # name
            column = 'snowflake_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve column values (by default from triplet df) 
        arr = _get_df_values(self, df_source=df_source, column=column) 
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

    def select_riming_class(self, values, method='Praz2017', invert=False, df_source="triplet"): 
        """
        Select MASCDB data with specific riming classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the riming classes to select.
            If integers, it assumes riming_class_id.
            If strings, it assumes riming_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_riming_class_name_dict(method)'.
        method : str, optional
            Method used to determine riming_class. The default is 'Praz2017'.
        invert : bool, optional
            If True, instead of selecting it discard the specified riming_class.
            The default is False.
        df_source: str, optional 
            The dataframe from which retrieve the class. 
            Either 'cam0', 'cam1', 'cam2' or 'triplet'.
            The  default is 'triplet'.
            
        Returns
        -------
        MASCDB class instance with specific riming classes.
    
        """
        #---------------------------------------------------------------------.
        ## Check default args 
        if not isinstance(invert, bool):
            raise TypeError("'invert' must be either True or False'.")
        _check_df_source(df_source)
        #---------------------------------------------------------------------.
        ## Check values
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_riming_class_name_dict(method=method).values()) # id
            column = 'riming_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_riming_class_name_dict(method=method).keys())   # name
            column = 'riming_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve column values (by default from triplet df) 
        arr = _get_df_values(self, df_source=df_source, column=column) 
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
 
    def select_melting_class(self, values, method='Praz2017', invert=False, df_source="triplet"): 
        """
        Select MASCDB data with specific melting classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the melting classes to select.
            If integers, it assumes melting_class_id.
            If strings, it assumes melting_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_melting_class_name_dict(method)'.
        method : str, optional
            Method used to determine melting_class. The default is 'Praz2017'.
        invert : bool, optional
            If True, instead of selecting it discard the specified melting_class_id.
            The default is False.
        df_source: str, optional 
            The dataframe from which retrieve the class. 
            Either 'cam0', 'cam1', 'cam2' or 'triplet'.
            The  default is 'triplet'.
            
        Returns
        -------
        MASCDB class instance with specific melting classes.
    
        """
        #---------------------------------------------------------------------.
        ## Check default args 
        if not isinstance(invert, bool):
            raise TypeError("'invert' must be either True or False'.")
        _check_df_source(df_source)
        #---------------------------------------------------------------------.
        ## Check values
        if not isinstance(values,(int, str, list, np.ndarray)):
            raise TypeError("'values' must be either (list of) integers (for class ids) or str (for class names).")
        # Convert to numpy array object 
        if isinstance(values, (int,str)):
            values = np.array([values])
        else: 
            values = np.array(values)
        # If values are integers --> Assume it provide the class id
        if isinstance(values[0].item(), int):
            valid_names = list(get_melting_class_name_dict(method=method).values()) # id
            column = 'melting_class_id'    
        # If values are str --> Assume it provide the class name
        elif isinstance(values[0].item(), str):
            valid_names = list(get_melting_class_name_dict(method=method).keys())   # name
            column = 'melting_class_name'
        else:
            raise TypeError("'values' must be either integers (for class ids) or str (for class names).")  
        #---------------------------------------------------------------------.
        # Retrieve column values (by default from triplet df) 
        arr = _get_df_values(self, df_source=df_source, column=column) 
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
        """
        Select MASCDB data with specific precipitation types.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the precipitation classes to select.
            If integers, it assumes bs_precip_class_id.
            If strings, it assumes bs_precip_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_precip_class_name_dict(method)'.
        method : str, optional
            Method used to determine bs_precip_class. The default is 'Schaer2020'.
        invert : bool, optional
            If True, instead of selecting it discard the specified bs_precip_class.
            The default is False.

        Returns
        -------
        MASCDB class instance with specific precipitation classes.
    
        """
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
        arr = self._triplet[column].to_numpy()
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
       
    def discard_snowflake_class(self, values, method='Praz2017', df_source="triplet"):
        """
        Discard MASCDB data with specific snowflake classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the snowflake classes to discard.
            If integers, it assumes snowflake_class_id.
            If strings, it assumes snowflake_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_snowflake_class_name_dict(method)'.
        method : str, optional
            Method used to determine snowflake_class. The default is 'Praz2017'.
        df_source: str, optional 
            The dataframe from which retrieve the class. 
            Either 'cam0', 'cam1', 'cam2' or 'triplet'.
            The  default is 'triplet'.

        Returns
        -------
        MASCDB class instance with specific snowflake classes.
    
        """
        return self.select_snowflake_class(values=values, method=method, invert = True, df_source=df_source)  
    
    def discard_melting_class(self, values, method='Praz2017', df_source="triplet"):
        """
        Discard MASCDB data with specific melting classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the melting classes to discard.
            If integers, it assumes melting_class_id.
            If strings, it assumes melting_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_melting_class_name_dict(method)'.
        method : str, optional
            Method used to determine melting_class. The default is 'Praz2017'.
        df_source: str, optional 
            The dataframe from which retrieve the class. 
            Either 'cam0', 'cam1', 'cam2' or 'triplet'.
            The  default is 'triplet'.

        Returns
        -------
        MASCDB class instance with specific melting classes.
    
        """
        return self.select_melting_class(values=values, method=method, invert = True, df_source=df_source)  
    
    def discard_riming_class(self, values, method='Praz2017', df_source="triplet"):
        """
        Discard MASCDB data with specific riming classes.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the riming classes to discard.
            If integers, it assumes riming_class_id.
            If strings, it assumes riming_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_riming_class_name_dict(method)'.
        method : str, optional
            Method used to determine riming_class. The default is 'Praz2017'.
        df_source: str, optional 
                The dataframe from which retrieve the class. 
                Either 'cam0', 'cam1', 'cam2' or 'triplet'.
                The  default is 'triplet'.
                
        Returns
        -------
        MASCDB class instance with specific riming classes.
    
        """
        return self.select_riming_class(values=values, method=method, invert = True, df_source=df_source) 
      
    def discard_precip_class(self, values, method='Schaer2020'):
        """
        Discard MASCDB data with specific precipitation types.

        Parameters
        ----------
        values : (str, int, list)
            Values specifying the precipitation classes to discard.
            If integers, it assumes bs_precip_class_id.
            If strings, it assumes bs_precip_class_name.
            Valid values can be retrieved by calling 'mascdb.aux.get_precip_class_name_dict(method)'.
        method : str, optional
            Method used to determine bs_precip_class. The default is 'Schaer2020'.

        Returns
        -------
        MASCDB class instance with specific precipitation classes.
    
        """
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
        env_db = self._triplet[[*env_variables]].copy()
        env_db.columns = [column.strip("env_") for column in env_variables]
        return env_db
    
    @property
    def bs(self):
        columns = list(self._triplet.columns)
        bs_variables = [column for column in columns if column.startswith("bs_")]
        bs_db = self._triplet[[*bs_variables]].copy()
        bs_db.columns = [column.strip("bs_") for column in bs_variables]
        return bs_db
    
    @property
    def gan3d(self):
        columns = list(self._triplet.columns)
        gan3d_variables = [column for column in columns if column.startswith("gan3d_")]
        gan3d_db = self._triplet[[*gan3d_variables]].copy()
        gan3d_db.columns = [column.strip("gan3d_") for column in gan3d_variables]
        return gan3d_db
    
    @property
    def flake(self):
        columns = list(self._triplet.columns)
        flake_variables = [column for column in columns if column.startswith("flake_")]
        flake_db = self._triplet[[*flake_variables]].copy()
        flake_db.columns = [column.strip("flake_") for column in flake_variables]
        return flake_db
    
    @property
    def labels(self):
        labels_variables = get_vars_class()
        labels_db = self._triplet[[*labels_variables]].copy()
        return labels_db

    @property
    def event(self):
        # columns = list(self._triplet.columns)
        # event_columns = [column for column in columns if column.startswith("event_")]
        event_columns = ['event_id', 'event_duration', 'event_n_triplets',
                         'campaign', 'datetime', 'latitude', 'longitude','altitude']
        event_db = self._triplet[event_columns].groupby('event_id').first().reset_index()
        # Compute month and year 
        event_db['month'] = event_db['datetime'].dt.month
        event_db['year'] = event_db['datetime'].dt.year
        # Compute start_time and end_time 
        start_time = self._triplet[["event_id","datetime"]].groupby('event_id').min()
        start_time.columns = ['start_time']
        end_time = self._triplet[["event_id","datetime"]].groupby('event_id').max()
        end_time.columns = ['end_time']   
        _ = event_db.drop(columns="datetime", inplace=True)
        event_db = event_db.merge(start_time, left_on="event_id", right_index=True)
        event_db = event_db.merge(end_time, left_on="event_id", right_index=True)
        return event_db 
        
    @property
    def campaign(self): 
        #----------------------------------------------.
        # Retrieve data 
        df_event = self.event     
        df_triplet = self._triplet    
        c_id = 'campaign'
        #----------------------------------------------.
        # Compute location info  
        columns = ['latitude','longitude', 'altitude', 'campaign']
        info_location = df_event[columns].groupby('campaign').first()
        
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
        melting_class_counts = df_triplet[['melting_class_name',c_id]].groupby(c_id)['melting_class_name'].apply(_count_occurence).apply(lambda x: x[0]) 
        precipitation_class_counts = df_triplet[['bs_precip_class_name',c_id]].groupby(c_id)['bs_precip_class_name'].apply(_count_occurence).apply(lambda x: x[0]) 
        snowflake_class_counts.name = "snowflake_class"
        riming_class_counts.name = "riming_class"
        melting_class_counts.name = "melting_class"
        precipitation_class_counts.name = "precipitation_class"
        
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
        triplet = self._triplet.drop(columns=vars_not_add)
        full_db = full_db.merge(triplet, how="left")
        return full_db
    
    def ds_images(self, cam_id = None, campaign=None, img_id='img_id'):
        #----------------------------------------------------------------------.
        # Subset by campaign 
        if campaign is not None:
            if isinstance(campaign, str): 
                campaign = [campaign]
            campaign = np.array(campaign).astype(str)
            db_campaigns = self._triplet['campaign'].to_numpy().astype(str)
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
                      max_interval_without_images = np.timedelta64(4,'h'),
                      unit="ns"):   
        # This function modify in place       
        #----------------------------------------------------------.
        # - Extract relevant columns from triplet db
        db = self._triplet[['campaign','datetime']].copy()    
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
                                                  maximum_interval_without_timesteps = max_interval_without_images)
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
    def select_events_with_n_triplets(self, min=0, max=np.inf):
        """
        Select events with number of triplets between min and max. 

        Parameters
        ----------
        min : int, optional
            Minimum number of triplets. The default is 0.
        max : int, optional
            Maximum number of triplets. The default is np.inf.

        Returns
        -------
        MASCDB class instance
        
        """
        ## Check min and max values validity
        if not isinstance(min, int): 
            raise TypeError("'min' must be an integer.")
        if not isinstance(max, (int, float)):
            raise TypeError("'max' must be an integer (or np.inf).")
        if min < 0: 
            raise ValueError("'min' must be an integer larger or equal to 0.")
        if max < 1: 
            raise ValueError("'max' must be an integer larger or equal to 1.")
        if isinstance(max, float):
            if max != np.inf:
                raise ValueError("'max' must be an integer (or np.inf).")
        if min > max:
            raise ValueError("'min' must be smaller than 'max'.")
        if max < min:
            raise ValueError("'max' must be larger than 'min'.")
        #---------------------------------------------------------------------.
        # Retrieve subset index 
        df_event = self.event 
        idx_event_ids = (df_event['event_n_triplets'] >= min) & (df_event['event_n_triplets'] <= max) 
        event_ids_subset = df_event.loc[idx_event_ids, 'event_id'].to_numpy()  
        idx_bool_subset = np.isin(self._triplet['event_id'].to_numpy(), event_ids_subset)
        #---------------------------------------------------------------------.
        # Subset the data and return 
        return self.isel(idx_bool_subset)
    
    ## ------------------------------------------------------------------------.
    def select_events_with_duration(self, min=np.timedelta64(0,'ns'), max=np.timedelta64(365,'D')):
        """
        Select events with duration between min and max. 

        Parameters
        ----------
        min : (np.timedelta64, pd.Timedelta), optional
            Minimum duration. The default is 0 ns.
        max : (np.timedelta64, pd.Timedelta), optional
            Maximum duration. The default is 1 year.

        Returns
        -------
        MASCDB class instance
        
        """
        # Check min and max values validity
        if not isinstance(min, (np.timedelta64, pd.Timedelta)): 
            raise TypeError("'min' must be a np.timedelta64 or pd.Timedelta object.")
        if not isinstance(max, (np.timedelta64, pd.Timedelta)):
            raise TypeError("'max' must be a np.timedelta64 or pd.Timedelta object.")
        if isinstance(min, pd.Timedelta): 
            min = min.to_numpy()
        if isinstance(max, pd.Timedelta): 
            max = max.to_numpy()
        if min < np.timedelta64(0,'ns'): 
            raise ValueError("'min' must be a positive timedelta object")
        if max < np.timedelta64(0,'ns'): 
            raise ValueError("'max' must be a positive timedelta object (larger than 0).")
        if min > max:
            raise ValueError("'min' must be smaller than 'max'.")
        if max < min:
            raise ValueError("'max' must be larger than 'min'.")
        #---------------------------------------------------------------------.
        # Retrieve subset index 
        df_event = self.event 
        idx_event_ids = (df_event['event_duration'] >= min) & (df_event['event_duration'] <= max) 
        subset_event_ids = df_event.loc[idx_event_ids, 'event_id'].to_numpy()  
        idx_subset = np.isin(self._triplet['event_id'].to_numpy(), subset_event_ids)
        # Subset the data and return 
        return self.isel(idx_subset)
            
    def select_events_longest(self, n=1):
        """
        Select MASCDB data corresponding to the 'n' events with longest duration.

        Parameters
        ----------
        n : int, optional
            The number of events to retrieve. The default is 1.

        Returns
        -------
        MASCDB class instance

        """
        longest_event_ids = self.arrange('triplet.event_duration', decreasing=True)._triplet['event_id'].iloc[0:n]
        idx_longest_events = np.isin(self._triplet['event_id'].to_numpy(), longest_event_ids)
        return self.isel(idx_longest_events)
    
    def select_events_shortest(self, n=1):
        """
        Select MASCDB data corresponding to the 'n' events with shortest duration.

        Parameters
        ----------
        n : int, optional
            The number of events to retrieve. The default is 1.

        Returns
        -------
        MASCDB class instance

        """
        shortest_event_ids = self.arrange('triplet.event_duration', decreasing=False)._triplet['event_id'].iloc[0:n]
        idx_shortest_events = np.isin(self._triplet['event_id'].to_numpy(), shortest_event_ids)
        return self.isel(idx_shortest_events) 
    
    ##--------------------------------------------------
    ## Redefine events utils  
    def redefine_events(self, 
                        max_interval_without_images = np.timedelta64(4,'h'),
                        min_duration = None, max_duration = None,
                        min_n_triplets = None, max_n_triplets = None,
                        unit="ns"): 
        """
        Enable selection and custom definition of an 'event'.
        
        If <min/max>_<duration/n_triplets> are specified, the MASCDB will likely be subsetted.

        Parameters
        ----------
        max_interval_without_images : (np.timedelta64, pd.Timedelta), optional
            Maximum interval of time without images to consider 
            consecutive images to belong the same event. 
            The default is np.timedelta64(4,'h').
        min_duration : (np.timedelta64, pd.Timedelta), optional
            Minimum duration of an event to retained. The default is np.timedelta64(0,'ns').
        max_duration : (np.timedelta64, pd.Timedelta), optional
            Maximum duration of an event to retained. The default is np.timedelta64(365,'D').
        min_n_triplets : int, optional
            Minimum number of triplets within an event to retain the event. The default is 0.
        max_n_triplets : int, optional
            Maximum number of triplets within an event to retain the event.. The default is Inf.
        unit : str, optional
            Unit of timedelta to consider for events definition. 
            The default is "ns".

        Returns
        -------
        MASCDB class instance with the custom event definition.

        """
        # Copy new instance 
        self = copy.deepcopy(self)
        # Define event_id 
        self._define_events(max_interval_without_images=max_interval_without_images,unit=unit)
        #----------------------------------------------------------.
        EVENT_FILTERING = False
        # Select only events with specific min/max n_triplets and duration
        if (min_n_triplets is not None) or (max_n_triplets is not None): 
            EVENT_FILTERING = True
            if min_n_triplets is None: 
                min_n_triplets = 0
            if max_n_triplets is None: 
                max_n_triplets = np.inf
            self = self.select_events_with_n_triplets(min=min_n_triplets, max=max_n_triplets)
        if (min_duration is not None) or (max_duration is not None): 
            EVENT_FILTERING = True
            if min_duration is None: 
                min_duration = np.timedelta64(0,'ns')
            if max_duration is None: 
                max_duration = np.timedelta64(365,'D')
            self = self.select_events_with_duration(min=min_duration, max=max_duration)
        #----------------------------------------------------------.
        # Ensure event_id incremental order (0,1,..,n_events) if filtering out events 
        if EVENT_FILTERING: 
            event_ids = self._triplet['event_id'].values
            event_ids_new = np.unique(event_ids, return_inverse=True)[1]
            self._cam0['event_id'] = event_ids_new
            self._cam1['event_id'] = event_ids_new
            self._cam2['event_id'] = event_ids_new
            self._triplet['event_id'] = event_ids_new
        # Return the object 
        return self 
        
    ####------------------------------------------------------------------------.
    #################################
    #### Image plotting routines ####
    #################################
    def plot_triplets(self, indices=None, random = False, n_triplets = 1,
                      enhancement="histogram_equalization",
                      zoom=True, squared=True, 
                      wspace=0.01, hspace=0.01,
                      **kwargs):
        """
        Plotting routine to display specific triplets of MASC snowflake images.
        
        By default:
        - images are enhanced with histogram_equalization and zoomed.
        - 'n_triplets' and 'random' are effective only if 'indices' are not specified.
        - If indices are unspecified, the chosen triplets correspond to the first 'n_triplets' of MASCDB.
        
        Parameters
        ----------
        indices : (int, list), optional
            Integer list of rows to display. The default is None.
        random : bool, optional
           Specify if the displayed MASCDB triplets must be choosen randomly.
           It's effective only if 'indices' are not specified.
           The default is False.
        n_triplets : int, optional
           Specify the number of MASCDB triplets to be displayed.   
           It's effective only if 'indices' are not specified.
           The default is 1.
        enhancement : str, optional
            Type of enhancement to use to improve the image quality. 
            Valid enhancements are : [None, "histogram_equalization", "contrast_stretching", "local_equalization"]
            The default is "histogram_equalization".
        zoom : bool, optional
            Specify if zooming close to the snowflake bounding box. 
            The image shape is defined by selecting the smallest possible shapes 
            across all the snowflakes to be plotted
            The default is True.
        squared : bool, optional
            Specify if the zoomed images must have equal height,width. 
            The default is True.
        hspace : float 
            Define the space across images in the vertical dimension.
            The default is 0.01.
        wspace : float 
            Define the space across images in the horizontal dimension.
            The default is 0.01.
        **kwargs : dict
            Optional arguments to be passed to DataArray.plot.

        Returns
        -------
        xarray.plot.facetgrid.FacetGrid object for additional customization

        """
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
        p.fig.subplots_adjust(wspace=wspace, hspace=hspace)   
        #--------------------------------------------------. 
        return p       
            
    def plot_flake(self, cam_id=None, index=None, random = False,
                   enhancement="histogram_equalization",
                   zoom=True, squared=True, ax=None, **kwargs):
        """
        Plotting routine to display a specific MASC snowflake image.
        
        By default:
        - The image is enhanced with histogram_equalization and zoomed.
        - 'random' is effective only if 'index' is not specified.
        - If index is unspecified, it plot an image of the first MASCDB triplet.
        
        Parameters
        ----------
        cam_id : int, optional
            The camera from which display the snowflake image.
            If not specified, the camera is randomly chosen. 
            Valid cam_id values are 0, 1 and 2.
            The default is None.
        index : int, optional
            Row index of the MASCDB triplet image to display. 
            The default is None.
        random : bool, optional
           Specify if the displayed MASCDB image must be choosen randomly.
           It's effective only if 'index' is not specified.
           The default is False.
        enhancement : str, optional
            Type of enhancement to use to improve the image quality. 
            Valid enhancements are : [None, "histogram_equalization", "contrast_stretching", "local_equalization"]
            The default is "histogram_equalization".
        zoom : bool, optional
            Specify if zooming close to the snowflake bounding box. 
            The image shape is defined by selecting the smallest possible shape 
            to include the entire snowflake.
            The default is True.
        squared : bool, optional
            Specify if the zoomed images must have equal height,width. 
            The default is True.
        ax: matplotlib axis, optional
            Optional matplotlib axis on which to plot the image.
            The default is None.
        **kwargs : dict
            Optional arguments to be passed to DataArray.plot.
            
        Returns
        -------
        p : TYPE
            DESCRIPTION.

        """
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
               cam_id = list(np.random.choice([0,1,2], 1))
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
            da_img = xri_zoom(da_img, squared=squared)
            
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
                    zoom=True, squared=True, 
                    hspace=0.1, wspace=0.1,
                    **kwargs):
        """
        Plotting routine to display MASC snowflake images.
        
        By default:
        - images are enhanced with histogram_equalization and zoomed.
        - 'n_images' and 'random' are effective only if 'indices' are not specified.
        - If indices are unspecified:
        * If cam_id is unspecified: it displays the first 'n_images' from a randomly selected camera of MASCDB.
        * If cam_id specify 1 camera: it displays the first 'n_images' of the specified camera of MASCDB.
        * If cam_id specifies more than 1 camera: it displays the first 'n_images' of each of the specified camera of MASCDB.
             
        Parameters
        ----------
        cam_id : (int, list), optional 
            The camera(s) from which display the snowflake images.
            If not specified, a single camera is randomly chosen. 
            If specified, it can be any subset of the 3 camera.
            Valid cam_id values are 0, 1 and 2.
            The default is None.
         indices : (int, list), optional
            Integer list of rows to display. The default is None.
        random : bool, optional
           Specify if the displayed MASCDB images must be choosen randomly.
           It's effective only if 'indices' are not specified.
           The default is False.
        n_images : int, optional
           Specify the number of MASCDB images to be displayed for each camera.   
           It's effective only if 'indices' are not specified.
           The default is 1.
        enhancement : str, optional
            Type of enhancement to use to improve the image quality. 
            Valid enhancements are : [None, "histogram_equalization", "contrast_stretching", "local_equalization"]
            The default is "histogram_equalization".
        zoom : bool, optional
            Specify if zooming close to the snowflake bounding box. 
            The image shape is defined by selecting the smallest possible shapes 
            across all the snowflakes to be plotted.
            The default is True.
        squared : bool, optional
            Specify if the zoomed images must have equal height,width. 
            The default is True.
        hspace : float 
            Define the space across images in the vertical dimension.
            The default is 0.1.
        wspace : float 
            Define the space across images in the horizontal dimension.
            The default is 0.1.
        **kwargs : dict
            Optional arguments to be passed to DataArray.plot.

        Returns
        -------
        xarray.plot.facetgrid.FacetGrid object for additional customization

        """
        #--------------------------------------------------
        # Check args
        n_idxs = len(self)
        _check_random(random)
        _check_zoom(zoom)
        _check_enhancement(enhancement)
        #--------------------------------------------------
        # Check cam_id 
        if cam_id is None: 
           cam_id = np.random.choice([0,1,2], 1).tolist()
           
        if isinstance(cam_id, int):
            cam_id = [cam_id]
        #--------------------------------------------------    
        # Define indices if is not provided
        _check_n_images(n_images, vmax=n_idxs)
        if indices is None:
            if random:  
               indices = list(np.random.choice(n_idxs, n_images))
            else: 
               indices = list(np.arange(0,n_images))
        
        #--------------------------------------------------        
        # Check indices and recompute n_images 
        indices = _check_indices(indices, vmax=n_idxs-1)
        if isinstance(indices, int): 
            indices = [indices]
        n_images = len(indices)*len(cam_id)
        
        #-------------------------------------------------- 
        # If a single flake is specified, plot it with plot_flake 
        if len(indices) == 1 and len(cam_id) == 1: 
            print("It's recommended to use 'plot_flake()' to plot a single image.")
            return self.plot_flake(index=indices[0], cam_id=cam_id, random=random,
                                   enhancement=enhancement, zoom=zoom, *kwargs)
       
        #--------------------------------------------------.
        # Retrieve DataArray and subset cam_id   
        da = self.da
        da = da.isel(cam_id=cam_id, flake_id = indices)
        #--------------------------------------------------.
        # If more than 1 between cam_id and flake_id dimensions are present --> Stack 
        dims = list(da.dims)
        unstacked_dims = list(set(dims).difference(["x","y"]))
        # - If only x and y, do nothing
        if len(unstacked_dims) == 0: 
            stack_dict = {}
            da_stacked = da
            n_idxs = 1
        # - If there is already a third dimension, transpose to the last 
        elif len(unstacked_dims) == 1: 
              img_id = unstacked_dims[0]
              da_stacked = da
              n_idxs = len(da_stacked[img_id])
        #     stack_dict = {}
        #     da_stacked = da.stack(stack_dict).transpose(..., img_id) 
        #     n_idxs = len(da_stacked[img_id])
        # - If there is more than 3 dimensions, stack it all into a new third dimension
        elif len(unstacked_dims) > 1: 
            img_id = "img_id"
            stack_dict = {img_id: ("cam_id", "flake_id")}
            # stack_dict = {img_id: unstacked_dims}
            # Stack all additional dimensions into a 3D array with all img_id in the last dimension 
            da_stacked = da.stack(stack_dict).transpose(..., img_id)
            n_idxs = len(da_stacked[img_id])
        else: 
            raise NotImplementedError()

        #--------------------------------------------------.
        da_subset = da_stacked 
        
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
            da_subset = xri_zoom(da_subset, squared=squared)
        
        #--------------------------------------------------.
        # Retrieve title from stacked dimension
        xr_indexes = da_subset.img_id.xindexes[img_id]
        FLAG_MULTI_INDEX = False
        if isinstance(xr_indexes,  xr.core.indexes.PandasMultiIndex):
            FLAG_MULTI_INDEX = True
            pd_indexes = xr_indexes.to_pandas_index()
            names = list(pd_indexes.names)
            titles = []
            for i in range(n_images):
                tmp_str_list = ", ".join([str(pd_indexes.get_level_values(name)[i]) for name in names])
                # tmp_str_list = ", ".join([name + ": " + str(pd_indexes.get_level_values(name)[i]) for name in names])
                titles.append(tmp_str_list)
                
        #--------------------------------------------------.
        # Plot flakes(s)
        row = img_id #  if len(indices) > 1 else None
        p = da_subset.plot(x='x',y='y', 
                           row=row, col_wrap=col_wrap, 
                           aspect=1, 
                           yincrease=False,
                           cmap='gray', add_colorbar=False, 
                           vmin=0, vmax=255,
                           #**kwargs,
                          )
        # Nice layout 
        for i, ax in enumerate(p.axes.flat):
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_axis_off() 
            if FLAG_MULTI_INDEX and i < len(titles):
                ax.set_title(titles[i])
        p.fig.subplots_adjust(wspace=wspace,hspace=hspace)
        #--------------------------------------------------. 
        return p     
 
    ####-----------------------------------------------------------------------.
    ####################### 
    #### MASCDB Updates ###
    ####################### 
    def compute_2Dimage_descriptors(self, fun, labels, fun_kwargs = None, force=False,
                                    dask = "parallelized"):
        """
        Compute user-specific image descriptors with parallelized computations
        and add it to the cam dataframes. 
        
        It requires the specification of a function ('fun') expecting the image 2D array 
        and returning the descriptor(s) value(s).
        It also require the specification of the expected descriptors names ('labels'). 

        Parameters
        ----------
        fun : callable
            A function computing the descriptor(s) of a 2D image.
            The function must expects a grayscale 2D array and return the descriptor(s) value(s).
        labels : (str, list)
            String or list of string specifying the descriptor names computed by 'fun'.
            These labels will become the colums added to cam dataframe.
        fun_kwargs : dict, optional
            Optional arguments to be passed to 'fun'. The default is None.
        force : bool, optional
            force=True enable to overwrite existing descriptors present in the cam dataframes.
            The default is False.
        dask : str, optional
            Option to be passed to xr.apply_u_func. 
            The default is "parallelized".

        Returns
        -------
        MASCDB class instance with new descriptors in cam dataframes.

        """
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
        cam0_flake_ids = np.sort(cam0.index.values.astype(str))
        cam1_flake_ids = np.sort(cam0.index.values.astype(str))
        cam2_flake_ids = np.sort(cam0.index.values.astype(str))
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
        existing_flake_ids = self._cam0.index.values.astype(str)
        
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
        self._cam0 = self._cam0.merge(cam0, left_index=True, right_index=True, how="left")
        self._cam1 = self._cam1.merge(cam1, left_index=True, right_index=True, how="left")
        self._cam2 = self._cam2.merge(cam2, left_index=True, right_index=True, how="left")
        
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
        complete : bool, optional
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
        df_flake_ids = np.sort(df.index.values.astype(str))
                
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
        existing_flake_ids = self._triplet.index.values.astype(str)
        
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
        self._triplet = self._triplet.merge(df, left_index=True, right_index=True, how="left")
        
        #---------------------------------------------------------------------.
        # Return the new mascdb 
        return self 
      
    def drop_cam_columns(self, columns): 
        """
        Method allowing to safely remove columns from all cam dataframes of MASCDB.
        
        Parameters
        ----------
        columns : list 
            List with column names of MASCDB cam dataframes to be removed
       
        Returns
        -------
        MASCDB class instance
        
        """
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # Check columns 
        columns = _check_columns(columns)  
        # Check columns are valid columns 
        current_columns = np.array(list(self._cam0.columns))
        unvalid_columns = np.array(columns)[np.isin(columns, current_columns, invert=True)]
        if len(unvalid_columns) > 0: 
            raise ValueError("{} are not columns of cam dataframes.".format(unvalid_columns.tolist()))
        #---------------------------------------------------------------------.
        # - Remove columns 
        columns = columns.tolist()
        _ = self._cam0.drop(columns=columns, inplace=True)
        _ = self._cam1.drop(columns=columns, inplace=True)
        _ = self._cam2.drop(columns=columns, inplace=True)
        #---------------------------------------------------------------------.
        return self   
        
    def drop_triplet_columns(self, columns): 
        """
        Method allowing to safely remove columns from the MASCDB triplet dataframe.
        
        Parameters
        ----------
        columns : list 
            List with column names of cam dataframes of MASCDB to be removed
       
        Returns
        -------
        MASCDB class instance
        
        """
        #---------------------------------------------------------------------.
        # Copy new instance 
        self = copy.deepcopy(self)
        #---------------------------------------------------------------------.
        # Check columns 
        columns = _check_columns(columns)  
        # Check columns are valid columns 
        current_columns = np.array(list(self._triplet.columns))
        unvalid_columns = np.array(columns)[np.isin(columns, current_columns, invert=True)]
        if len(unvalid_columns) > 0: 
            raise ValueError("{} are not columns of cam dataframes.".format(unvalid_columns.tolist()))
        #---------------------------------------------------------------------.
        # - Remove columns 
        columns = columns.tolist()
        _ = self._triplet.drop(columns=columns, inplace=True)
        #---------------------------------------------------------------------.
        return self  
    