"""
Helper functions for basic manipulation of the MASC database
of particle descriptors (parquet) and images (zarr archive)
"""

from typing import TypedDict
import pandas as pd
import numpy as np
import zarr

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from copy import deepcopy

import seaborn as sns

def wet_bulb_t(t,rh):

    return t * np.arctan(0.152*(rh+8.3136)**(1/2)) + np.arctan(t + rh) - np.arctan(rh - 1.6763) + 0.00391838 *(rh)**(3/2) * np.arctan(0.0231 * rh) - 4.686



class masc_database:
    """
    Masc database class to read and manipulate the 4 databases of 
    descriptors (one for each cam) as well TODO

    """
    def __init__(self,parquet_path,parquet_id='MASCdb',zarr_id='MASCdb'):
        """
        TO DO
        
        """

        # Filenames 
        fn_cam0 = parquet_path+parquet_id+'_cam0.parquet'
        fn_cam1 = parquet_path+parquet_id+'_cam1.parquet'
        fn_cam2 = parquet_path+parquet_id+'_cam2.parquet'
        fn_triplet = parquet_path+parquet_id+'_triplet.parquet'

        self.db_files   = [fn_cam0,fn_cam1,fn_cam2,fn_triplet]        
        self.zarr_files = parquet_path+zarr_id+'.zarr'

        # Read data into dataframes
        self.cam0    = pd.read_parquet(fn_cam0)
        self.cam1    = pd.read_parquet(fn_cam1)
        self.cam2    = pd.read_parquet(fn_cam2)
        self.triplet = pd.read_parquet(fn_triplet)

    def keep_precip_type(self,type='blowing_snow',mix_level=0.5):
        """
        Filter all dataframes according to precipitation type
        based on the paper of TODO

        type: "blowing_snow" (default) --> keep pure blowing snow
                "precip"                 --> keep pure precip
                "blowing_snow_and_mixed" --> keep blowing_snow and mixed cases (mixing larger than mix_level)
                "precip_and_mixed"       --> keep pure precip and mixed cases  (mixing lower than mix_level)
                "mixed"                  --> keep only mixed (any mixing ratio)
                
                TODO
        """
        
        masc_db = deepcopy(self)

        if type == 'blowing_snow':
            cond = masc_db.triplet.bs_nor_angle > 0.881
        elif type == 'precip':
            cond = masc_db.triplet.bs_nor_angle < 0.193
        elif type == 'mixed':
            cond = (masc_db.triplet.bs_mix_ind >= 0) & (masc_db.triplet.bs_mix_ind <= 1)
        elif type == 'blowing_snow_and_mixed':
            cond = (masc_db.triplet.bs_nor_angle > 0.881) | (masc_db.triplet.bs_mix_ind >= mix_level)
        elif type == 'precip_and_mixed':
            cond = (masc_db.triplet.bs_nor_angle < 0.193) | (masc_db.triplet.bs_mix_ind <= mix_level)
        else: 
            print("Doing nothing, unrecongnized filter type: "+type)

        # Apply condition
        masc_db.cam0 = masc_db.cam0[cond]
        masc_db.cam1 = masc_db.cam1[cond]
        masc_db.cam2 = masc_db.cam2[cond]
        masc_db.triplet = masc_db.triplet[cond]

        return(masc_db)

    def plot_flake(self, index=None, random = True, info = False,zoom=True):
        """
        
            index:   the index of the snowflake to plot
            random: if true, a random snowflake is plotted from the class 
            zoom: TODO
            info: TODO
        
        """
        if random: # (use df.sample)
            index = self.triplet.sample().index[0]

        if self.zarr_files is None:
            print("Zarr archive not set, cannot plot flakes")
            return None

        zz = zarr.open(self.zarr_files,mode='r')

        # PLot the snowflake
        dat = zz[index,:,:,:]

        if zoom:
            true_points = np.argwhere(dat)
            # take the smallest points and use them as the top left of your crop
            top_left = true_points.min(axis=0)
            # take the largest points and use them as the bottom right of your crop
            bottom_right = true_points.max(axis=0)
            out = dat[top_left[0]:bottom_right[0]+1,  # plus 1 because slice isn't
                top_left[1]:bottom_right[1]+1,:]  # inclusive
        else: 
            out = dat

        pix_size = self.triplet.loc[index].pix_size*1e3  # mm
        xsize = out.shape[0]*pix_size # mm
        ysize = out.shape[1]*pix_size # mm

        fig = plt.figure( tight_layout=True)
        gs = gridspec.GridSpec(ncols=3,nrows=1,hspace=0.01,wspace=0.05)

        for j in range(3):
            plt.subplot(gs[0,j])
            plt.imshow(out[:,:,j],cmap='gray',extent=[0,xsize,0,ysize],aspect=1)        
            if j==0:
                plt.gca().tick_params(left=True, labelleft=True)
                plt.ylabel("[mm]")
            else:
                plt.gca().tick_params(left=False, labelleft=False)

        gs.tight_layout(fig)

        return fig
















