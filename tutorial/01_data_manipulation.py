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
#os.chdir("/home/ghiggi/Projects/pymascdb")
os.chdir("/home/grazioli/CODES/python/pymascdb")

import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt
import mascdb.api
from mascdb.api import MASC_DB

#dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
dir_path = "/data/MASC_DB"
 
##----------------------------------------------------------------------------.
###########################
#### MASCDB Introduction ##
###########################
### MASCDB General Structure  
# mascdb.<mascdb_method> 
# mascdb.da.<xarray.DataArray methods>
# mascdb.<cam*,triplet,env,bs,gan3d,event,campaign, full_db>.<pandas.DataFrame methods>  
# mascdb.<cam*,triplet,env,bs,gan3d,event,campaign, full_db>.sns.<seaborn plot methods>

### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Display structure
print(mascdb)
len(mascdb)

# Display attributes and methods 
dir(mascdb)

##----------------------------------------------------------------------------.
#########################
#### MASCDB Attributes ##
#########################
# They are implemented as 'properties' in order to return a copy of the original data !
        
## Core dataset of MASCDB 
# - DataArray with MASC images 
mascdb.da 
# - Dataframes of image descriptors and classes for each camera 
mascdb.cam0
mascdb.cam1
mascdb.cam2 
# - Dataframe with characteristics of each image acquisition
mascdb.triplet 

## Additional dataframes with application-specific informations (contained in mascdb.triplet) 
mascdb.flake   # flake avg properties across cam images 
mascdb.labels  # class labels 
mascdb.env     # atmospheric environment 
mascdb.bs      # blowing snow estimation using Schaer et al., 2019 method
mascdb.gan3d   # mass & volume estimation using GAN algorithm (Leinonen et al., 2021)

# - Get units and explanation of variables/columns
mascdb.get_var_units('flake_Dmax')
mascdb.get_var_explanation('flake_Dmax')

mascdb.get_var_units('snowflake_mood')       # this returns an error
mascdb.get_var_explanation('flake_religion') # also this

# Event summary information 
event_summary = mascdb.event     
print(event_summary)
print(event_summary.columns)

# Campaign summary information 
campaign_summary = mascdb.campaign  
print(campaign_summary)
print(campaign_summary.columns)
campaign_summary[["start_time", "end_time"]]
campaign_summary[["n_triplets", "n_events", 'total_event_duration']]
campaign_summary[['event_duration_min','event_duration_mean', 'event_duration_max', ]]
campaign_summary[['snowflake_class']]
campaign_summary[['riming_class']]
campaign_summary[['melting_class']]
campaign_summary[['precipitation_class']]
 
## Other database that can extracted (STILL IN DEVELOPMENT)
mascdb.full_db  #  slow !

mascdb.ds_images()
mascdb.ds_images(cam_id=[0,1])
mascdb.ds_images(campaign=['Valais-2016', 'PLATO-2019'])
mascdb.ds_images(campaign=['Valais-2016', 'PLATO-2020'])

##----------------------------------------------------------------------------.
###################
#### Assignement ##
###################
# Warning: currently it is not possible to add additional MASC images to a MASCDB instance

# It is possible to extract a copy of a dataframe and add columns or modify values directly
cam0 = mascdb.cam0
cam0.loc[0] = 1
print(cam0.loc[0])              # Here the assignement works

cam0['new_column'] = 2
print(cam0["new_column"])       # Here a column is added 

# It's not possible to add columns or modify values directly to a MASCDB dataframe   
mascdb.cam0.loc[0] = 1
print(mascdb.cam0.loc[0])       # Here nothing is assigned

mascdb.cam0['new_column'] = 1
print(mascdb.cam0['new_column']) # Here nothing is added

### Add columns to a MASCDB dataframe the following functions must be used  
# ---> mascdb.add_triplet_columns()
# ---> mascdb.add_cam_columns()

# - Calculate wet bulb temperature and add it as a column to triplet dataframe
from mascdb.utils_env import wet_bulb_t # Util to generate wet bulb temp

df = mascdb.triplet
df['env_Twb'] = wet_bulb_t(df["env_T"].to_numpy(), df["env_RH"].to_numpy())

new_mascdb = mascdb.add_triplet_columns(df['env_Twb'], force=False, complete=True)
new_mascdb.triplet

# - Add columns only for a partial subset of flake_id 
df_subset = df.iloc[0:100]
new_mascdb = mascdb.add_triplet_columns(df_subset['env_Twb'], force=False, complete=True)
new_mascdb = mascdb.add_triplet_columns(df_subset['env_Twb'], force=False, complete=False)
new_mascdb.triplet

##----------------------------------------------------------------------------.
########################
### Dropping columns ###
########################
# - Remove columns from triplet dataframe 
new_mascdb.drop_triplet_columns("env_Twb").triplet

# - Remove columns from triplet dataframe 
mascdb.cam0
mascdb.drop_cam_columns("event_id").cam0

##----------------------------------------------------------------------------.
##################
#### Subsetting ##
##################

# Sample randomly n rows
mascdb.sample_n(10).cam0

# Select first(s) 
mascdb.first().cam0
mascdb.head().cam0

# Select last(s)
mascdb.last().cam0
mascdb.tail().cam0

# Subset specific rows with boolean indices
# - Boolean indices must have same length of the database 
idx = mascdb.cam0['Dmax'] > 0.02
print(idx)
mascdb_largeD = mascdb.isel(idx) 

# Subset specific rows with positional indices
mascdb.isel(1) 
mascdb.isel([1,2,3,4]) 

# Subset with slice
# - Remember that the last element of slice is not included ! 
mascdb.isel(slice(100,125))   

# Unvalid ways to subsettig mascdb
mascdb.isel(-1)  
mascdb.isel(1.2) 
mascdb.isel(False)
mascdb.isel(True)
mascdb.isel(None) 
mascdb.isel(np.nan) 
mascdb.isel(np.inf) 
mascdb.isel([np.nan]) 
mascdb.isel([np.inf]) 

# Subsetting with boolean indices with length < len(mascdb)
# - Proceed with caution using such subsetting procedure !  
mascdb.isel([False])  # Empty DBs
mascdb.isel([True])   # First row (index=0)
mascdb.isel([True, False])  # First row
mascdb.isel([False, True])  # Second row 

idx = mascdb.cam0['Dmax'] > 0.02
mascdb.isel(idx[0:2]) 
mascdb.isel(idx[0:100000]) 
 
# Subsetting dataframe columns
mascdb.cam0['Dmax']
mascdb.cam0[['Dmax','perim','pix_size']] 

# Subsetting dataframe rows by position index 
mascdb.cam0.iloc[slice(0,10)]
mascdb.cam0.iloc[[1,2,5]]
mascdb.cam0.iloc[1]   # return a pd.Series
mascdb.cam0.iloc[[1]] # return a pd.Dataframe

# Subsetting dataframe rows by flake_id 
mascdb.cam0.loc['2015.02.10_11.55.15_flake_9']   
mascdb.cam0.loc[['2015.02.10_11.55.15_flake_9','2015.02.10_11.55.16_flake_10']]   

mascdb.cam0.loc[0] # This clearly does not work 
mascdb.cam0.loc['0'] # This clearly does not work 

# Subsettings mascdb triplets by flake_id 
mascdb.sel('2015.02.10_11.55.15_flake_9')   
mascdb.sel(['2015.02.10_11.55.15_flake_9','2015.02.10_11.55.16_flake_10'])
mascdb.sel('0')     # This clearly does not work 
mascdb.sel('')      # This clearly does not work 
mascdb.sel([1,2])   # This clearly does not work 
mascdb.sel(['0','2015.02.10_11.55.15_flake_9'])  # This clearly does not work 

##----------------------------------------------------------------------------.
################
#### Sorting  ##
################
# Sort by snowflake diameter 
mascdb_largeD = mascdb.arrange('cam0.Dmax', decreasing=True)     
mascdb_largeD.cam0['Dmax']
mascdb_largeD.plot_triplets(n_triplets = 3)

mascdb_smallD = mascdb.arrange('cam0.Dmax', decreasing=False) 
mascdb_smallD.cam0['Dmax']
mascdb_smallD.plot_triplets(n_triplets = 3, zoom=True)
mascdb_smallD.plot_triplets(n_triplets = 3, zoom=False)

##----------------------------------------------------------------------------.
#################
#### Filtering ##
################# 
############################ 
#### - Descriptors-based  ##
############################ 
# expression = 'triplet.env_T

mascdb.select_max('cam0.Dmax', n=5)    
mascdb.select_min('cam0.Dmax', n=5)

mascdb.select_min('triplet.env_T', n=5)   
mascdb.select_min('env.T', n=5)

mascdb.select_max('triplet.bs_normalized_angle', n=5)
mascdb.select_min('triplet.bs_mixing_ind', n=5)
mascdb.select_min('bs.mixing_ind', n=5)


mascdb.select_min('gan3d.volume', n=5)  
mascdb.select_min('gan3d.mass', n=5)
 
# Select based on complex conditions 
idx = (mascdb.cam0['Dmax'] > 0.02) & (mascdb.env['RH'] > 50)
mascdb1 = mascdb.isel(idx)  
mascdb1.cam0.sns.scatterplot(x="Dmax", y="solidity")

mascdb.isel(idx).cam0.sns.scatterplot(x="Dmax", y="solidity") 

mascdb.isel((mascdb.cam0['Dmax'] > 0.02)).cam0.sns.scatterplot(x="Dmax", y="solidity") 

#####################
#### - Time-based ###
#####################
dir(mascdb.cam0.datetime.dt)

# Subset specific month  
mascdb.isel(mascdb.cam0.datetime.dt.month == 2).cam0["datetime"] 
mascdb.isel(mascdb.cam0.datetime.dt.month_name() == "February").cam0["datetime"] 

# Subset multiple months
mascdb.isel(mascdb.cam0.datetime.dt.month.isin([2,4])).cam0["datetime"] 
mascdb.isel(mascdb.cam0.datetime.dt.month_name().isin(["February","April"])).cam0["datetime"] 

# Select specific years 
mascdb.isel(mascdb.cam0.datetime.dt.year == 2016).cam0["datetime"] 

# Select morning aquisition 
mascdb.isel(mascdb.cam0.datetime.dt.hour.isin(np.arange(6,12))).cam0["datetime"] 

# Select night-time aquisition 
mascdb.isel(mascdb.cam0.datetime.dt.hour.isin(np.append(np.arange(0,6),np.arange(20,24)))).cam0["datetime"] 

######################## 
#### - Campaign-based ##
######################## 
mascdb.select_campaign('PLATO-2020')
mascdb.select_campaign('PLATO-2019')
mascdb.select_campaign(['Valais-2016','PLATO-2019'])
mascdb.select_campaign(['Valais-2016','PLATO-2020'])
mascdb.discard_campaign('PLATO-2019')   
mascdb.discard_campaign('PLATO-2020')   

#####################  
#### - Class-based ## 
##################### 
mascdb.select_melting_class('dry')  
mascdb.select_melting_class('melted')  

mascdb.select_riming_class('medium')
mascdb.select_riming_class('rimed').triplet['riming_class_name']  
mascdb.select_riming_class(['rimed', 'densely_rimed']).triplet['riming_class_name'] 
mascdb.select_riming_class([2,3]).triplet['riming_class_id']
mascdb.select_riming_class(2).triplet['riming_class_id']
mascdb.select_riming_class(5)
 
mascdb.select_snowflake_class("rimed")
mascdb.select_snowflake_class(10) 
mascdb.select_snowflake_class(1).triplet['snowflake_class_id']

mascdb.select_snowflake_class('aggregate').triplet['snowflake_class_name']
mascdb.select_snowflake_class(['columnar_crystal', 'planar_crystal']).triplet['snowflake_class_name']
   
mascdb.select_riming_class('medium').select_max(cam0.Dmax, n=5).plot_triplets()

mascdb.select_precip_class('rain')
mascdb.select_precip_class('precip')
mascdb.discard_precip_class('blowing_snow')
mascdb.discard_precip_class(['blowing_snow','undefined'])
mascdb.discard_precip_class(['blowing_snow','rain'])

##----------------------------------------------------------------------------.
####################### 
#### - Quality-based ## 
####################### 
# Select images with  average high quality and plot it 
idx = (mascdb.triplet['flake_quality_xhi'] > 9)
mascdb_high_quality_img = mascdb.isel(idx).arrange('triplet.flake_quality_xhi', decreasing=True)
mascdb_high_quality_img.plot_triplets(n_triplets = 3, zoom=True)

# Select images from specific camera with high quality and plot it 
idx = (mascdb.cam0['quality_xhi'] > 9)
mascdb_cam0_high_quality = mascdb.isel(idx).arrange('cam0.quality_xhi', decreasing=True)
mascdb_cam0_high_quality.plot_flakes(cam_id=0, n_images=9, col_wrap=3, zoom=True)
mascdb_cam0_high_quality.plot_triplets(n_triplets = 3, zoom=True)

##----------------------------------------------------------------------------.
#######################
#### - Event-based   ##
#######################
# - By default, events are defined as periods of time with at least 1 
#   image acquisition within a time interval of 4 hours 

# - The event information columns are: 
#    - event_duration  
#    - event_n_triplets
#    - event_id 
# - event_id is present in triplet and cam* dataframes 
# - event_n_triplets and event_duration are present only in the triplet dataframe
 
# The default arguments of mascdb.define_events are: 
# - maximum_interval_without_images = 4h
# - minimum_duration = None   # include everything
# - minimum_n_triplets = 0    # include everything
# --> If minimum_event_duration or minimum_n_triplets > 0, filtering is performed !!!
# --> 'event_id' change each time define_events is called !

##-----------------------------------------
# Watch event information 
mascdb.cam0['event_id'] 
mascdb.triplet[['event_id','event_duration','event_n_triplets']]
mascdb.event

##-----------------------------------------
# (Re)define default events 
# mascdb = mascdb.redefine_events() # as default event definition 
# mascdb.cam0['event_id'] 
# mascdb.triplet[['event_id','event_duration','event_n_triplets']]
# mascdb.event

##-----------------------------------------
# Filter by number of n_triplet images within an event   
# - Select events with more than n triplets
new_mascdb = mascdb.select_events_with_n_triplets(min=50)
print(new_mascdb.triplet[['event_id','event_n_triplets']])
print(new_mascdb.event)
print(np.unique(new_mascdb.event['event_n_triplets']))

# - Select events with less than n triplets
new_mascdb = mascdb.select_events_with_n_triplets(max=1)
print(new_mascdb.triplet[['event_id','event_n_triplets']])
print(new_mascdb.event)
print(np.unique(new_mascdb.event['event_n_triplets']))
print(np.unique(new_mascdb.event['event_duration'])) # When only 1 triplet, event_duration is 0 

# - Select events with more than n triplets and less than n triplets 
new_mascdb = mascdb.select_events_with_n_triplets(min=1, max=10)
print(new_mascdb.triplet[['event_id','event_n_triplets']])
print(new_mascdb.event)
print(np.unique(new_mascdb.event['event_n_triplets']))

##-----------------------------------------
# Filter by event duration
# - Less than 10 minutes
new_mascdb = mascdb.select_events_with_duration(max=np.timedelta64(10,'m'))
print(new_mascdb.triplet[['event_id','event_duration']])
print(np.unique(new_mascdb.event['event_n_triplets'])) # The shortest the duration, the less the number of triplets
 
# - Less than 10 seconds
# --> When there is only 1 triplet, event_duration is 0 
new_mascdb = mascdb.select_events_with_duration(max=np.timedelta64(10,'s'))
print(new_mascdb.triplet[['event_id','event_duration']])
print(np.unique(new_mascdb.event['event_n_triplets'])) # The shortest the duration, the less the number of triplets
 
# - More than 2 hours duration 
new_mascdb = mascdb.select_events_with_duration(min=np.timedelta64(2,'h'))
print(new_mascdb.triplet[['event_id','event_duration']])
print(np.unique(new_mascdb.event['event_n_triplets']))
 
# - More than 1 day duration 
new_mascdb = mascdb.select_events_with_duration(min=np.timedelta64(1,'D'))
print(new_mascdb.triplet[['event_id','event_duration']])
print(np.unique(new_mascdb.event['event_n_triplets']))

##-----------------------------------------
# Select longest or shortest events 
mascdb.select_events_longest()
mascdb.select_events_longest(1)

mascdb.select_events_longest(2)

mascdb.select_events_shortest(10)

##-----------------------------------------
### Redefine events with custom thresholds  
# - This can filter out images from events that do not match min/max n_triplets and duration
max_interval_without_images = pd.Timedelta(2,'h')
max_interval_without_images = np.timedelta64(2, 'h')
min_duration = np.timedelta64(2, 'm')
min_n_triplets = 10 

new_mascdb = mascdb.redefine_events(max_interval_without_images = max_interval_without_images,
                                    min_duration = min_duration,
                                    min_n_triplets = min_n_triplets) 
print(new_mascdb.cam0['event_id'])
print(new_mascdb.campaign[["start_time", "end_time"]])
print(new_mascdb.triplet[['event_id','event_duration','event_n_triplets']])
print(new_mascdb.event)

##----------------------------------------------------------------------------.
################### 
#### Image plots ##
################### 

##----------------------------------------------------------------------------.
####################### 
#### Dataframe plots ##
####################### 
# Property plots 
mascdb.env.sns.scatterplot(x="T", y="DD", hue="P")

# mascdb.full_db.sns.scatterplot(x="Dmax", y="perim", hue="CAM_ID")

##----------------------------------------------------------------------------.
#############################
#### Data analysis example ##
#############################
### Blowing snow analysis 
## Plot precip_type for event  with duration less 1 minute (should be blowing snow no?)
mascdb_short_events = mascdb.select_events_with_duration(max=np.timedelta64(60,'s'))
print(mascdb_short_events.event)
mascdb_short_events.triplet.sns.boxplot(x="bs_precip_class_name", y="event_n_triplets") # This make sense

### Investigate blowing snow occurence 
mascdb_few_triplets = mascdb.select_events_with_n_triplets(max=10)
print(mascdb_few_triplets.event)

mascdb_few_triplets.triplet.sns.boxplot(x="bs_precip_class_name", y="event_n_triplets") # This is strange 
 
# See some stats on event durations and n_triplets per events 
mascdb.triplet.sns.scatterplot(x="event_duration", y="event_n_triplets")

### Dmax analysis  
# Retrieve event leading to largest Dmax 
mascdb.arrange('cam0.Dmax', decreasing=True)  
largest_Dmax_event_id = mascdb.arrange('cam0.Dmax', decreasing=True).cam0['event_id'].iloc[0]
idx_largest_Dmax_event = mascdb.cam0['event_id'] == largest_Dmax_event_id
mascdb_event = mascdb.isel(idx_largest_Dmax_event).arrange('triplet.datetime', decreasing=False)
 
mascdb_event.cam0.sns.scatterplot(x="datetime", y="Dmax")
mascdb_event.full_db.sns.scatterplot(x="datetime", y="env_T")
mascdb_event.full_db.sns.scatterplot(x="datetime", y="Dmax", hue="env_T")

mascdb_event.cam0.sns.kdeplot(x="datetime", y="Dmax", shade=True)
mascdb_event.cam0.sns.kdeplot(x="datetime", y="Dmax", fill=True, thresh=0, levels=100, cmap="mako")

# mascdb_event.full_db.sns.kdeplot(x="datetime", y="Dmax", shade=True)
# mascdb_event.full_db.sns.kdeplot(x="datetime", y="3dgan_vol_ch", fill=True, thresh=0, levels=100, cmap="mako")



##----------------------------------------------------------------------------.
