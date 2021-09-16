#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 11:26:36 2021

@author: ghiggi
"""
##########################################
### MASCDB Data Manipulation tutorial ####
##########################################
#-----------------------------------------------------------------------------.
import os
os.chdir("/home/ghiggi/Projects/pymascdb")
# os.chdir("/home/grazioli/CODES/python/pymascdb")
import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt
import mascdb.api
from mascdb.api import MASC_DB

dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
#dir_path = "/data/MASC_DB"
 
##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Display structure
print(mascdb)
len(mascdb)

##----------------------------------------------------------------------------.
### General 
# mascdb.<mascdb_method> 
# mascdb.da.<xarray.DataArray methods>
# mascdb.<cam*,triplet,env,bs,gan3d,full_db>.<pandas.DataFrame methods> # Do not modify original objects !
# mascdb.<cam*,triplet,env,bs,gan3d,full_db>.sns.<seaborn plot methods>

##----------------------------------------------------------------------------.
## Properties 
## Property class mascdb
mascdb.env   
mascdb.bs     
mascdb.gan3d  # Error ... strippa via ... r_g --> r

mascdb.full_db  #  slow !

mascdb.ds_images()
mascdb.ds_images(CAM_ID=[0,1])
mascdb.ds_images(campaign=['Valais-2016', 'PLATO-2019'])
mascdb.ds_images(campaign=['Valais-2016', 'PLATO-2020'])
##----------------------------------------------------------------------------.
# Property plots 
mascdb.env.sns.scatterplot(x="T", y="DD", hue="P")

mascdb.full_db.sns.scatterplot(x="Dmax", y="perim", hue="CAM_ID")

##----------------------------------------------------------------------------.
# It's not possible to assign values to mascdb dataframe directly 
cam0 = mascdb.cam0
cam0.loc[0] = 1
print(cam0.loc[0])        # Here the assignement works

mascdb.cam0.loc[0] = 1
print(mascdb.cam0.loc[0]) # Here nothing is modified

##----------------------------------------------------------------------------.
# Filtering 
idx = (mascdb.cam0['Dmax'] > 0.02) & (mascdb.env['RH'] > 50)
mascdb1 = mascdb.isel(idx)  
mascdb1.cam0.sns.scatterplot(x="Dmax", y="solidity")

mascdb.isel(idx).cam0.sns.scatterplot(x="Dmax", y="solidity") 
mascdb.isel((mascdb.cam0['Dmax'] > 0.02)).cam0.sns.scatterplot(x="Dmax", y="solidity") 


mascdb.select_max('cam0.Dmax', n=5)
mascdb.select_min('cam0.Dmax', n=5)

mascdb.from_campaign(campaign='PLATO-2020')
mascdb.from_campaign(campaign='PLATO-2019')
mascdb.from_campaign(campaign=['Valais-2016','PLATO-2019'])
mascdb.from_campaign(campaign=['Valais-2016','PLATO-2020'])
mascdb.exclude_campaign(campaign='PLATO-2019')   
mascdb.exclude_campaign(campaign='PLATO-2020')   

mascdb.select_flakes('aggregate') # select_class ?
mascdb.select_rimed('medium')
mascdb.select_melting(0)
   
mascdb.select_rimed('medium').select_max(cam0.Dmax, n=5).plot_triplets()

"""
 - Example to filter on blowing snow / precip
 - Example to filter quality_xhi
 
"""
##----------------------------------------------------------------------------.
# Filtering and sorting 
idx = mascdb.cam0['Dmax'] > 0.02
mascdb_largeD = mascdb.isel(idx) 
mascdb_largeD = mascdb_largeD.arrange('cam0.Dmax', decreasing=True)  
mascdb_largeD.plot_triplets(n_triplets = 3, zoom=False)
mascdb_largeD.plot_triplets(n_triplets = 3, zoom=True)    

##----------------------------------------------------------------------------.
### Define events (timedelta_thr define the allowed time interval without images)
timedelta_thr = pd.Timedelta(2,'h')
timedelta_thr = np.timedelta64(2, 'h')
mascdb.define_event_id(timedelta_thr=timedelta_thr)

# See some stats on event durations and n_triplets per events 
mascdb.triplet.sns.scatterplot(x="event_duration", y="event_n_triplets")

# Retrieve longest event 
mascdb.arrange('triplet.event_duration', decreasing=True).triplet
longest_event_id = mascdb.arrange('triplet.event_duration', decreasing=True).triplet['event_id'].iloc[0]
idx_longest_event = mascdb.cam0['event_id'] == longest_event_id

# Retrieve event leading to largest Dmax 
mascdb.arrange('cam0.Dmax', decreasing=True) # TODO CHECK: WHY THIS WARNING 
largest_Dmax_event_id = mascdb.arrange('cam0.Dmax', decreasing=True).cam0['event_id'].iloc[0]
idx_largest_Dmax_event = mascdb.cam0['event_id'] == largest_Dmax_event_id

# - Subset event and plot some stuffs
mascdb_event = mascdb.isel(idx_largest_Dmax_event).arrange('triplet.datetime', decreasing=False)
mascdb_event = mascdb.isel(idx_longest_event).arrange('triplet.datetime', decreasing=False)
print(mascdb_event)

mascdb_event.cam0.sns.scatterplot(x="datetime", y="Dmax")
mascdb_event.full_db.sns.scatterplot(x="datetime", y="env_T")
mascdb_event.full_db.sns.scatterplot(x="datetime", y="Dmax", hue="env_T")

mascdb_event.cam0.sns.kdeplot(x="datetime", y="Dmax", shade=True)
mascdb_event.cam0.sns.kdeplot(x="datetime", y="Dmax", fill=True, thresh=0, levels=100, cmap="mako")

# mascdb_event.full_db.sns.kdeplot(x="datetime", y="Dmax", shade=True)
# mascdb_event.full_db.sns.kdeplot(x="datetime", y="3dgan_vol_ch", fill=True, thresh=0, levels=100, cmap="mako")

##----------------------------------------------------------------------------.
### Selecting specific dataframe variables 
mascdb.cam0[['Dmax','perim','pix_size']].sns.pairsplot()

##----------------------------------------------------------------------------.
### Sample triplets
mascdb.sample_n(10).cam0
mascdb.first_n().cam0
mascdb.last_n().cam0
mascdb.tail().cam0
mascdb.head().cam0