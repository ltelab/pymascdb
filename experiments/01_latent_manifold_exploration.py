#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 21:15:36 2021

@author: ghiggi
"""
###########################################
### MASCDB Latent Manifold Exploration ####
###########################################
#-----------------------------------------------------------------------------.
import os

os.chdir("/home/ghiggi/Projects/pymascdb")
# os.chdir("/home/grazioli/CODES/python/pymascdb")
from pathlib import Path
import umap
import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib as mpl 
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

import mascdb.api
from mascdb.api import MASC_DB

dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
figs_path = "/home/ghiggi/Projects/pymascdb/figs/LatentManifold"

#dir_path = "/data/MASC_DB"
 
##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

##----------------------------------------------------------------------------.
## TODO 
# - Dimension reduction for specific snowflake class ... to see patterns 
# - Dimension reduction for specific rimed class .... to see patterns 
 
##----------------------------------------------------------------------------.

from mascdb.aux import get_vars_cam_descriptors
from mascdb.aux import get_vars_class_ids
from mascdb.aux import get_vars_class_names

from mascdb.aux import get_snowflake_class_name_colors_dict
from mascdb.aux import get_riming_class_name_colors_dict
from mascdb.aux import get_snowflake_class_id_colors_dict
from mascdb.aux import get_riming_class_id_colors_dict
from mascdb.aux import get_campaign_colors_dict

# Define descriptors to use for dimension reduction
cam_descriptors = get_vars_cam_descriptors()
cam_descriptors = ['n_roi', 'area','perim','Dmax','area_porous','compactness',
                   'bbox_width','bbox_len','solidity','nb_holes','complexity']
 
#------------------------------------------------------------------------------
### - Define colors for snow particles type
riming_class_id_colors_dict = get_riming_class_id_colors_dict()
snowflake_class_id_colors_dict = get_snowflake_class_id_colors_dict()
 
#----------------------------------------------------------------------------- 
### - Define numeric descriptors you want to visualize 
viz_descriptors = ['perim','Dmax',"roundness","compactness",'nb_holes',] 

### - Define categorical class you want to visualize 
viz_classes = ["snowflake_class_name","campaign"] # [*get_vars_class_names(),"campaign"]

### - Define color dictionary for each categorical class 
colors_dict = {'snowflake_class_name': get_snowflake_class_name_colors_dict(),
               'campaign': get_campaign_colors_dict(),
              }

####--------------------------------------------------------------------------. 
###############################
#### Data preprocessing data ##
###############################
# - Subset mascdb for faster dimension reduction example 
subset_campaigns = ['APRES3-2016', 'Davos-2015', 'Jura-2019', 'ICEGENESIS-2021']
mascdb_subset = mascdb.select_campaign(subset_campaigns)
# mascdb_subset = mascdb 

### Explore the manifold using cam0 image descriptors 
# - Retrieve descriptors to use 
X = mascdb_subset.cam0[cam_descriptors]
# - Coerce columns to be of same datatype  
X = X.astype(float)
# - Retrieve row where there are missing values in the descriptors
idx_nan = np.isnan(X).any(axis=1)
# - Remove rows with NaN values 
X = X.loc[~idx_nan]

### Retrieve variables you want to visualize 
Y_viz_class = mascdb_subset.triplet[viz_classes].loc[~idx_nan]
X_viz_numeric = mascdb_subset.cam0[viz_descriptors].loc[~idx_nan]

### Standardize data matrix used for dimension reduction 
scaler = MinMaxScaler()
scaler.fit(X)
X_std = scaler.transform(X)

####--------------------------------------------------------------------------.
#################################
#### Dimensionality reduction ###
#################################
### - PCA
pca = PCA(n_components=2)
pca.fit(X_std)
pca_embedding = pca.transform(X_std)

### - UMAP
reducer = umap.UMAP(n_neighbors=30, min_dist=0.1, verbose=True)
umap_embedding = reducer.fit_transform(X_std)

####--------------------------------------------------------------------------.
#############################################
#### Dimension reduction plots with PCA  ####
#############################################
from mascdb.utils_figs import cm2inch
from mascdb.utils_figs import get_c_cmap_from_color_dict
from mascdb.utils_figs import get_legend_handles_from_colors_dict
from mascdb.utils_figs import minmax
 
algorithms = ["PCA", "UMAP"]
  
for algorithm in algorithms:
    print(" - Generating manifold plots for:", algorithm)
    #---------------------------------------------------------------.  
    # Retrieve latent codes 
    if algorithm == "PCA":
        x_latent = pca_embedding[:, 0]
        y_latent = pca_embedding[:, 1]
        xlabel = "PC1"
        ylabel = "PC2"
        title = 'PCA projection'
    elif algorithm == "UMAP": 
        x_latent = umap_embedding[:, 0]
        y_latent = umap_embedding[:, 1]
        xlabel = "L1"
        ylabel = "L2"
        title = 'UMAP projection'
    else: 
        raise NotImplementedError
    #---------------------------------------------------------------.
    # Define directory and filepath
    folderpath =  os.path.join(figs_path, algorithm)
    if not os.path.exists(folderpath):
        os.makedirs(folderpath)
    #---------------------------------------------------------------.  
    # Define xlim and ylim 
    xlim = minmax(x_latent)
    ylim = minmax(y_latent)
    #---------------------------------------------------------------.  
    #### - Plot numeric variables 
    for column in viz_descriptors:
        #---------------------------------------------------------------.  
        # Define image filepath
        tmp_filename = algorithm + "_" + str(column) + ".png"
        tmp_filepath = os.path.join(folderpath, tmp_filename)
        #---------------------------------------------------------------.  
        # - Define title
        cbar_title = column.title().replace("_", " ")
        #---------------------------------------------------------------.
        # - Define figure layout 
        plt.figure(figsize=cm2inch(12,12), dpi=400)
        plt.style.use('dark_background')
        #---------------------------------------------------------------.
        # - Plot scatterplot 
        plt.scatter(x_latent,
                    y_latent,  
                    c = X[column].values, 
                    cmap = 'Spectral',  
                    marker = '.', 
                    s = 0.8, 
                    edgecolors = 'none')
        #---------------------------------------------------------------.
        # - Set axis labels, limits and title 
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title, fontsize=12)
        # plt.gca().set_aspect('equal', 'datalim')
        #---------------------------------------------------------------.
        # - Plot colorbar 
        cbar = plt.colorbar() # aspect= ... for width
        cbar.set_label(cbar_title, rotation=270)
        cbar.ax.get_yaxis().labelpad = 15   
        #---------------------------------------------------------------.
        # - Save figure 
        plt.savefig(tmp_filepath)
        plt.close()    # to not display in IPython console
        
        #---------------------------------------------------------------.
    #-------------------------------------------------------------------------.  
    #### - Plot categorical variables 
    for column in viz_classes:
        # Define image filepath
        tmp_filename = algorithm + "_" + str(column) + ".png"
        tmp_filepath = os.path.join(folderpath, tmp_filename)
        #---------------------------------------------------------------.  
        # - Define title
        legend_title = column.title().replace("_", " ")
        #---------------------------------------------------------------.
        # - Retrieve marker colors 
        color_dict = colors_dict[column]
        c, cmap = get_c_cmap_from_color_dict(color_dict, labels=Y_viz_class[column])  
        
        #---------------------------------------------------------------.
        # - Define figure layout 
        fig, ax = plt.subplots(1, 1, figsize= cm2inch(14,10), dpi=400)
        plt.style.use('dark_background')
        #---------------------------------------------------------------.
        # - Plot scatterplot 
        ax.scatter(x_latent,
                    y_latent,  
                    c = c, cmap = cmap,  
                    marker = '.', 
                    s = 0.8, 
                    edgecolors = 'none')
        #---------------------------------------------------------------.
        # - Set axis labels, limits and title 
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize = 12)
        # ax.set_aspect('equal', 'datalim')
        #---------------------------------------------------------------.
        # - Create space for legend 
        box = ax.get_position() # Shrink axis on the bottom to create space for the legend
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height *0.9]) # left, bottom, width, height
        # - Retrieve legend handles and make pretty labels 
        handles = get_legend_handles_from_colors_dict(color_dict)
        labels = [handle.get_label() for handle in handles]
        labels = [label.title().replace("_", " ") for label in labels]
        _ = [handle.set_label(label) for label, handle in zip(labels, handles)]      
        # - Plot legend  
        ax.legend(handles = handles,  
                   title = legend_title, 
                   loc='upper center',
                   ncol = 3,
                   bbox_to_anchor=(0.5, -0.12))
        
        #---------------------------------------------------------------.
        # - Save figure  
        fig.tight_layout()
        fig.savefig(tmp_filepath)
        plt.close()  # to not display in IPython console
        
        #---------------------------------------------------------------------.     