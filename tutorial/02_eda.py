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
os.chdir("/home/ghiggi/Projects/pymascdb")
# os.chdir("/home/grazioli/CODES/python/pymascdb")
import numpy as np
import pandas as pd 
import xarray as xr
import seaborn as sns
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
### Image plots 
# - By default the 'enhancement' is (adaptive) histogram_equalization
mascdb.plot_flake(CAM_ID=0, random = True, zoom=True)
mascdb.plot_flake(CAM_ID=0, random = False, zoom=True)
mascdb.plot_flake(CAM_ID=0, index=0, random = True, zoom=True)
mascdb.plot_flake(CAM_ID=[0], index=[0], random = True, zoom=True)

mascdb.plot_flakes(CAM_ID=0, indices=[0], random = True, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=0, random = True, zoom=True)  
mascdb.plot_flakes(CAM_ID=0, indices=None, random = True, n_images = 1, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=None, random = True, n_images = 2, zoom=True)
mascdb.plot_flakes(CAM_ID=[0], indices=None, random = True, n_images = 9, zoom=True) 
mascdb.plot_flakes(CAM_ID=0, indices=[1,2], random = True, n_images = 2, zoom=True)
mascdb.plot_flakes(CAM_ID=0, indices=[1,2], random = True, n_images = 9, zoom=True)  

mascdb.plot_triplets(indices=[0,10], random = True, n_triplets = 1, zoom=True)
mascdb.plot_triplets(indices=[0,10], random = True, n_triplets = 1, zoom=True, enhancement=None)
mascdb.plot_triplets(indices=None, random = True, n_triplets = 1, zoom=True) #
mascdb.plot_triplets(indices=None, random = True, n_triplets = 3, zoom=True)

##----------------------------------------------------------------------------.
## Added accessor for seaborn plot routines ;) 
# - hue: categorical group over which to compute stats 
# - https://seaborn.pydata.org/examples/index.html: Gallery from which to draw examples :) 

mascdb.cam0.sns.boxplot(x="Dmax")   
mascdb.cam0.sns.boxenplot(x="Dmax")   
mascdb.cam0.sns.violinplot(x="Dmax") 

mascdb.cam0.sns.violinplot(x="Dmax", y="snowflake_class_name") 
mascdb.cam0.sns.stripplot(x="Dmax", y="snowflake_class_name") 

mascdb.cam0.sns.violinplot(x="perim", y="snowflake_class_name", hue="riming_class_id") 

mascdb.cam0.sns.scatterplot(x="Dmax", y="perim")

mascdb.cam0.sns.barplot(x="snowflake_class_name", y="riming_class_id")
mascdb.cam0.sns.catplot(x="snowflake_class_name", y="riming_class_id") # problems with '' or None ? 

mascdb.cam0.sns.histplot(x="snowflake_class_name", y="riming_class_id")
mascdb.cam0.sns.histplot(x="snowflake_class_name", hue="riming_class_id")

# For complex stuff ... first filter to a small number 
mascdb = mascdb.isel(slice(0,4000))
mascdb.cam0.sns.jointplot(x="Dmax", y="perim", kind="kde", hue="snowflake_class_name") 
 

cam_descriptors = ['n_roi', 'area','perim','Dmax','area_porous','compactness',
                   'bbox_width','bbox_len','solidity','nb_holes','complexity']
 
mascdb.cam0.sns.pairplot(vars = cam_descriptors[0:5])

mascdb.cam0.sns.corrplot(vars = cam_descriptors[0:5],
                         vmin = -1, vmax=1, center=0,
                         cbar_kws={"shrink": .5},  
                         linewidths=.5)

mascdb.cam0.sns.kdeplot(x="Dmax",y="perim", 
                        cmap="rocket")

sns.set_theme(style="white")
mascdb.cam0.sns.kde_marginals(x="Dmax",y="perim", 
                              # xlim=(), ylim=(), 
                              space=0, thresh=0, 
                              levels=100, cmap="rocket",
                              hist_color = "#03051A", hist_alpha=1,hist_bins=25)

pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
mascdb.cam0.sns.kde_ridgeplot(x = "compactness",
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