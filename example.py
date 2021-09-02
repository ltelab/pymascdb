#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 22:50:08 2021

@author: ghiggi
"""
#-----------------------------------------------------------------------------.
import os
os.chdir("/home/ghiggi/Projects/pymascdb")
#os.chdir("/home/grazioli/CODES/python/pymascdb")
import mascdb.api
from mascdb.api import MASC_DB

dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
#dir_path = "/data/MASC_DB/"

##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

##----------------------------------------------------------------------------.
## Properties 
## Property class mascdb

# mascdb.full_db  #  slow !
mascdb.env   
mascdb.bs     
mascdb.gan3d  # Error ... strippa via ... r_g --> r

##----------------------------------------------------------------------------.
# Property plots 
mascdb.env.sns.scatterplot(x="T", y="DD", hue="P")
mascdb.full_db.sns.scatterplot(x="Dmax", y="perim", hue="CAM_ID")

##----------------------------------------------------------------------------.
# Filtering 
idx = (mascdb.cam0['Dmax'] > 0.02) & (mascdb.env['RH'] > 50)
idx = mascdb.cam0['Dmax'] > 0.02

mascdb1 = mascdb.isel(idx)  
mascdb1.cam0.sns.scatterplot(x="Dmax", y="solidity")

mascdb.isel(idx).cam0.sns.scatterplot(x="Dmax", y="solidity") 

mascdb.isel((mascdb.cam0['Dmax'] > 0.02)).cam0.sns.scatterplot(x="Dmax", y="solidity") 

##----------------------------------------------------------------------------.
### Selecting specific dataframe variables 
mascdb.cam0[['Dmax','perim','pix_size']].sns.pairsplot()

##----------------------------------------------------------------------------.
## Sample 
mascdb.sample_n(10).cam0
mascdb.first_n().cam0
mascdb.last_n().cam0
mascdb.tail().cam0
mascdb.head().cam0

##----------------------------------------------------------------------------.
## Image plots 

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
