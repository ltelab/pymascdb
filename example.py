#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 22:50:08 2021

@author: ghiggi
"""
#-----------------------------------------------------------------------------.
import os
os.chdir("/home/grazioli/CODES/python/pymascdb")
#os.chdir("/home/grazioli/CODES/python/pymascdb")
import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt
import mascdb.api
from mascdb.api import MASC_DB

dir_path = "/data/MASC_DB"
#dir_path = "/data/MASC_DB/"

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

##----------------------------------------------------------------------------.
### Image plots 
# - By default the 'enhancement' is (adaptive) histogram_equalization
mascdb.plot_flake(CAM_ID=0, random = True, zoom=True)
mascdb.plot_flake(CAM_ID=0, random = False, zoom=True)
mascdb.plot_flake(CAM_ID=0, index=0, random = True, zoom=True)
mascdb.plot_flake(CAM_ID=[0], index=[0], random = True, zoom=True)

mascdb.plot_flakes(CAM_ID=0, indices=[0], random = True, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=0, random = True, zoom=True) # ERROR
mascdb.plot_flakes(CAM_ID=0, indices=None, random = True, n_images = 1, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=None, random = True, n_images = 2, zoom=True)
mascdb.plot_flakes(CAM_ID=[0], indices=None, random = True, n_images = 9, zoom=True) 
mascdb.plot_flakes(CAM_ID=0, indices=[1,2], random = True, n_images = 2, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=[1,2], random = True, n_images = 9, zoom=True)  

mascdb.plot_triplets(indices=[0,10], random = True, n_triplets = 1, zoom=True)
mascdb.plot_triplets(indices=None, random = True, n_triplets = 1, zoom=True)
mascdb.plot_triplets(indices=None, random = True, n_triplets = 3, zoom=True)

##----------------------------------------------------------------------------.
## Added accessor for seaborn plot routines ;) 
# - hue: categorical group over which to compute stats 
# - https://seaborn.pydata.org/examples/index.html: Gallery from which to draw examples :) 

mascdb.cam0.sns.boxplot(x="Dmax")   
mascdb.cam0.sns.boxenplot(x="Dmax")   
mascdb.cam0.sns.violinplot(x="Dmax") 

mascdb.cam0.sns.violinplot(x="Dmax", y="label_name") 
mascdb.cam0.sns.stripplot(x="Dmax", y="label_name") 

mascdb.cam0.sns.violinplot(x="perim", y="label_name", hue="riming_id") 

mascdb.cam0.sns.scatterplot(x="Dmax", y="perim")

mascdb.cam0.sns.barplot(x="label_id", y="riming_id")
mascdb.cam0.sns.catplot(x="label_id", y="riming_id") # problems with '' or None ? 

mascdb.cam0.sns.histplot(x="label_id", y="riming_id")
mascdb.cam0.sns.histplot(x="label_id", hue="riming_id")

# For complex stuff ... first filter to a small number 
mascdb.cam0.sns.jointplot(x="Dmax", y="perim", kind="kde", hue="label_name") 
 
#-----------------------------------------------------------------------------.
### EDA tools 
# All this plotting method are implemented ;) Let's then do example in an EDA.ipnyb
# --> For some plotting method ... we might create a check to not provide too much data or it stucks
violinplot
boxplot
boxenplot
stripplot
pointplot
swarmplot # heavy computation ? 
 
lmplot
pairplot # select first few variables 

scatterplot
displot # kind
catplot
barplot

histplot
jointplot # heavy computation ? 
kdeplot   # heavy computation ? 

relplot(data, x='datetime')  # https://seaborn.pydata.org/examples/timeseries_facets.html
lineplot(data, x='datetime') # https://seaborn.pydata.org/examples/wide_data_lineplot.html
relplot(data, x='datetime')  # https://seaborn.pydata.org/examples/faceted_lineplot.html

## TODO: 
corrplot # https://seaborn.pydata.org/examples/many_pairwise_correlations.html
kde_ridgeplot  # https://seaborn.pydata.org/examples/kde_ridgeplot.html
kde_marginals  # https://seaborn.pydata.org/examples/smooth_bivariate_kde.html

# https://seaborn.pydata.org/examples/horizontal_boxplot.html
# https://seaborn.pydata.org/examples/pairgrid_dotplot.html
# https://seaborn.pydata.org/examples/palette_generation.html           

#-----------------------------------------------------------------------------.
