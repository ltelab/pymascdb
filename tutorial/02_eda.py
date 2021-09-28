#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 11:28:36 2021

@author: ghiggi
"""
##################################################
### MASCDB Exploratory Data Analysis Tutorial ####
##################################################
#-----------------------------------------------------------------------------.
import os
#os.chdir("/home/ghiggi/Projects/pymascdb")
os.chdir("/home/grazioli/CODES/python/pymascdb")
import numpy as np
import pandas as pd 
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt
import mascdb.api
from mascdb.api import MASC_DB

#dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
dir_path = "/data/MASC_DB"
 
##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Display structure
print(mascdb)
len(mascdb)

##----------------------------------------------------------------------------.
### Image plots 
# - By default the 'enhancement' is (adaptive) histogram_equalization
# - By default random is set to False (plot ordered)
# - By defaulr zoom is enabled (zoom=True) and image is squared (squared=True)
# - n_images/n_triplets and random are used only if indices is not provided 
# 
# - Plot flakes triplet(s)
mascdb.plot_triplets(random = True, zoom=True) # by default plot a single triplet 
mascdb.plot_triplets(indices=[0,10], n_triplets = 1, zoom=True) # --> n_triplet has no effect
mascdb.plot_triplets(indices=[0,10], n_triplets = 1, zoom=True, enhancement=None)
mascdb.plot_triplets(random = True, n_triplets = 1, zoom=True) 
mascdb.plot_triplets(random = True, n_triplets = 3, zoom=True)

# - Plot single flake 
mascdb.plot_flake(random = True, zoom=True)  # random cam, random indices   
mascdb.plot_flake(cam_id=0, random = True,  zoom=True)
mascdb.plot_flake(cam_id=0, random = False, zoom=True)
mascdb.plot_flake(cam_id=0,   index=0,   random = True, zoom=True)
mascdb.plot_flake(cam_id=[0], index=[0], random = True, zoom=True)

mascdb.plot_flakes(cam_id=0,   indices=0, random = True, zoom=True)
mascdb.plot_flakes(cam_id=[0], indices=0, random = True, zoom=True)  
mascdb.plot_flakes(indices=[0], random = True, zoom=True)  

# - Plot multiple flakes from a single camera 
mascdb.plot_flakes(random = True, zoom=True)  # By default: single random camera and n_images = 9
mascdb.plot_flakes(random = True, n_images = 9, zoom=True)  # single random camera and n_images = 9

mascdb.plot_flakes(cam_id=0, indices=[1,2], n_images = 1, zoom=True) # --> 'indices' has high priority over 'n_images' 
mascdb.plot_flakes(cam_id=0, indices=[1,2], n_images = 9, zoom=True) # --> 'indices' has high priority over 'n_images' 
mascdb.plot_flakes(cam_id=0,   random = True, n_images = 2, zoom=True)
mascdb.plot_flakes(cam_id=[0], random = True, n_images = 9, zoom=True) 
 
# - Plot multiple flakes for each specific camera
mascdb.plot_flakes(cam_id=[0,1], random = True, n_images = 2, col_wrap=2, zoom=True)
mascdb.plot_flakes(cam_id=[0,1], indices=[0], random = True, zoom=True)   
mascdb.plot_flakes(cam_id=[0,1], indices=[0,1], col_wrap=2, zoom=True)  
mascdb.plot_flakes(cam_id=[1,0], indices=[0,1], col_wrap=2, zoom=True)  
mascdb.plot_flakes(cam_id=[0,1], indices=[0,1,2], col_wrap=3, zoom=True)  
mascdb.plot_flakes(cam_id=[0,1,2], indices=[0,1,2], col_wrap=3, zoom=True) 

##----------------------------------------------------------------------------.
## Added accessor for seaborn plot routines ;) 
# - hue: categorical group over which to compute stats 
# - https://seaborn.pydata.org/examples/index.html: Gallery from which to draw examples :) 

mascdb.cam0.sns.boxplot(x="Dmax")   
mascdb.cam0.sns.boxenplot(x="Dmax")   
mascdb.cam0.sns.violinplot(x="Dmax") 

# Physically-sound example 1 (temperature and snowflake properties)
# - plot temperature and Dmax
mascdb.triplet.sns.jointplot(x="env_T",y="flake_Dmax") # All data --> blowing snow contaminates
# - let's try the same plot but keeping only "precip"
mascdb1=mascdb.select_precip_class('precip')  
mascdb1.triplet.sns.jointplot(x="env_T",y="flake_Dmax") # Better 

# Physically-sound example 2 (wind and Dmax)
# - We know that high wind conditions hampers the measurement of large snowflakes
mascdb.triplet.sns.jointplot(x="env_FF",y="flake_Dmax") # Indeed

# Physically-sound example 3 (size and complexity)
# - Let's take cam0 as example
mascdb1 = mascdb1.select_snowflake_class('aggregate') # let's get aggregates 
mascdb1.cam0.sns.histplot(x='Dmax',y='complexity')




mascdb.cam0.sns.violinplot(x="Dmax", y="snowflake_class_name") 
mascdb.cam0.sns.stripplot(x="Dmax", y="snowflake_class_name") 

mascdb.cam0.sns.violinplot(x="perim", y="snowflake_class_name", hue="riming_class_id") 

mascdb.cam0.sns.scatterplot(x="Dmax", y="perim")

mascdb.cam0.sns.barplot(x="snowflake_class_name", y="riming_class_id")
mascdb.cam0.sns.catplot(x="snowflake_class_name", y="riming_class_id") # problems with '' or None ? 

mascdb.cam0.sns.histplot(x="snowflake_class_name", y="riming_class_id")
mascdb.cam0.sns.histplot(x="snowflake_class_name", hue="riming_class_id")

# For complex stuff ... first filter to a small number 
mascdb1 = mascdb.isel(slice(0,4000))
mascdb1.cam0.sns.jointplot(x="Dmax", y="perim", kind="kde", hue="snowflake_class_name") 
 

cam_descriptors = ['n_roi', 'area','perim','Dmax','area_porous','compactness',
                   'bbox_width','bbox_len','solidity','nb_holes','complexity']
 
mascdb1.cam0.sns.pairplot(vars = cam_descriptors[0:5])

mascdb1.cam0.sns.corrplot(vars = cam_descriptors[0:5],
                         vmin = -1, vmax=1, center=0,
                         cbar_kws={"shrink": .5},  
                         linewidths=.5)

mascdb1.cam0.sns.kdeplot(x="Dmax",y="perim", 
                        cmap="rocket")

sns.set_theme(style="white")
mascdb1.cam0.sns.kde_marginals(x="Dmax",y="perim", 
                              # xlim=(), ylim=(), 
                              space=0, thresh=0, 
                              levels=100, cmap="rocket",
                              hist_color = "#03051A", hist_alpha=1,hist_bins=25)

pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
mascdb1.cam0.sns.kde_ridgeplot(x = "compactness",
                              group = "snowflake_class_name",
                              linewidth = 2, 
                              pal = pal, bw_adjust=.5, height=.7,
                              aspect=15, hspace=-.0) # increase hspace to superpose


#-----------------------------------------------------------------------------.
### EDA tools 
# All this plotting method are implemented ;) Let's then do example in an EDA.ipnyb
# --> For some plotting method ... we might create a check to not provide too much data or it stucks
 
pointplot
swarmplot # heavy computation ? 
 
lmplot

scatterplot
displot # kind
 
relplot(data, x='datetime')  # https://seaborn.pydata.org/examples/timeseries_facets.html
lineplot(data, x='datetime') # https://seaborn.pydata.org/examples/wide_data_lineplot.html
relplot(data, x='datetime')  # https://seaborn.pydata.org/examples/faceted_lineplot.html

# Nice examples
# https://seaborn.pydata.org/examples/horizontal_boxplot.html
# https://seaborn.pydata.org/examples/pairgrid_dotplot.html
# https://seaborn.pydata.org/examples/palette_generation.html           

#-----------------------------------------------------------------------------.