#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 21:15:36 2021

@author: ghiggi
"""
##----------------------------------------------------------------------------.
## TODO 
# It could be interesting to:
# - Color the manifold for riming degree level
# - Color the manifold by Campaign to see if the feature manifold (and occurence) is different between campaigns
# - Derive a manifold only for rimed particles and then use the riming_deg_level to see patterns

##----------------------------------------------------------------------------.
import os
from pathlib import Path
import matplotlib as mpl 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap

### - Set filepaths  
dir_path = "/data/MASC_DB/"
figs_path = "/home/grazioli/tmp/Figs/"
feature_cam0_fpath = os.path.join(dir_path, "CH_cam0.parquet")
feature_cam1_fpath = os.path.join(dir_path, "CH_cam1.parquet")
feature_cam2_fpath = os.path.join(dir_path, "CH_cam2.parquet")

### - Set columns to remove (not useful for manifold analysis) 
columns_discarded = ['label_id_prob',
                     "label_prob_1",
                     "label_prob_2",
                     "label_prob_3",
                     "label_prob_4",
                     "label_prob_5",
                     "label_prob_6",
                     "riming_prob_1",
                     "riming_prob_2",
                     "riming_prob_3",
                     "riming_prob_4",
                     "riming_prob_5",
                     "cam",
                     "datetime",
                     "flake_id",
                     "pix_size",
                     'melting_id',
                     'melting_prob',
                     "Campaign",
                     "label_id",
                     "label_name",
                     "riming_deg_level",
                     'fallspeed', # per evitare di eliminare troppe obs (xke nan)
                     ]

#------------------------------------------------------------------------------
### - Define colors for snow particles type
particles_colors_dict = {'agg':'forestgreen', \
                         'capcol': 'darkblue',  
                         'col': 'red',  
                         'dend': 'orange',  
                         'graup': 'yellow', 
                         'smal': 'gray'
                         }
#----------------------------------------------------------------------------- 
### - Define feature/columns we want to visualize 
variables_viz = ['perim','area','Dmax'] 

#------------------------------------------------------------------------------. 
#########################
### - Create Database ###
#########################
cam0 = pd.read_parquet(feature_cam0_fpath)
cam1 = pd.read_parquet(feature_cam1_fpath)
cam2 = pd.read_parquet(feature_cam2_fpath)

df = cam0   # you might want to also add cam1 and cam2

#-----------------------------------------------------------------------------. 
###################### 
### Preprocessing ####
###################### 
# - Retrieve columns used for dimension reduction
columns_used_bool = np.isin(df.columns.tolist(), columns_discarded, invert=True)
columns_used = df.columns[columns_used_bool].tolist()

### - Retrieve features we would like to predict (or investigate)
campaign_name = df['Campaign'].values
label_id = df['label_id'].values
label_name = df['label_name'].values
riming_id = df['riming_id'].values
riming_deg_level = df['riming_deg_level'].values

### - Drop columns not used for manifold retrieval 
df_manifold = df.drop(columns=columns_discarded)
df_manifold.dtypes # different dtypes !

### - Coerce columns to same datatype (float64)
df_manifold = df_manifold.astype(float)
data_matrix = df_manifold.values

### - Remove missing values 
idx_NAN = np.isnan(data_matrix).any(axis=1)

df = df[~idx_NAN]
data_matrix = data_matrix[~idx_NAN,:]

campaign_name = campaign_name[~idx_NAN]
label_id = label_id[~idx_NAN]
label_name = label_name[~idx_NAN]
riming_id = riming_id[~idx_NAN]
riming_deg_level = riming_deg_level[~idx_NAN]

### - Standardize data 
scaler = StandardScaler()
scaler.fit(data_matrix)
data_matrix_std = scaler.transform(data_matrix)

#-----------------------------------------------------------------------------.
#################################
### Dimensionality reduction ####
#################################
### - PCA
pca = PCA(n_components=2)
pca.fit(data_matrix_std)
pca_embedding = pca.transform(data_matrix_std)

### - UMAP
reducer = umap.UMAP(n_neighbors=30, min_dist=0.1,)
umap_embedding = reducer.fit_transform(data_matrix_std)

#-----------------------------------------------------------------------------.
######################################
### Define functions for analysis ####
######################################
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
    
def get_c_cmap_from_color_dict(color_dict, labels): 
    """
    # Retrieve c and cmap argument for plt.scatter provided a custom color dictionary 
    # assign_dict_colors = lambda x : campaign_colors_dict[x]
    # c_names = list(map(assign_dict_colors, experiments))
    """
    c_names = [color_dict[x] for x in labels]
    # Retrieve c integer values 
    c, c_unique_name = pd.factorize(c_names, sort=False)
    # Create cmap
    cmap = mpl.colors.ListedColormap(c_unique_name)
    # Return object 
    return[c, cmap]
    
def get_legend_handles_from_colors_dict(colors_dict, marker='o'):
    """
    Retrieve 'handles' for the matplotlib.pyplot.legend
    # marker : "s" = filled square, 'o' = filled circle
    # marker : "PATCH" for filled large rectangles 
    """
    import matplotlib as mpl
    if (marker == 'PATCH'):
        # PATCH ('filled large rectangle')
        handles = []
        for key in colors_dict:
            data_key = mpl.patches.Patch(facecolor=colors_dict[key], edgecolor=colors_dict[key], label=key)
            handles.append(data_key)    
    else:
        # Classical Markers
        handles = []
        for key in colors_dict:
            data_key = mpl.lines.Line2D([0], [0], linewidth=0, \
                                        marker=marker, label=key, \
                                        markerfacecolor=colors_dict[key], \
                                        markeredgecolor=colors_dict[key], \
                                        markersize=3)
            handles.append(data_key)    
    return(handles)  


def minmax(x):
   return([np.min(x),np.max(x)])       

#------------------------------------------------------------------------------.
#####################################
## Dimension reduction with PCA  ####
#####################################

### - Plot continuous variables 
for variable in variables_viz:
    tmp_folderpath =  Path(figs_path,"PCA")
    tmp_folderpath.mkdir(exist_ok=True) # create dir if does not exist
    tmp_filename = "PCA_" + str(variable) + ".png"
    tmp_filepath = Path(tmp_folderpath, tmp_filename)
    
    # Plot scatter 
    plt.figure(figsize=(12,12), dpi=800)
    plt.style.use('dark_background')
    plt.scatter(pca_embedding[:, 0], pca_embedding[:, 1],  
                c=df[[variable]].to_numpy().squeeze(), cmap='Spectral',  
                marker='.', s=0.8, edgecolors='none')
    plt.xlim(minmax(pca_embedding[:, 0]))
    plt.ylim(minmax(pca_embedding[:, 1]))
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    # plt.gca().set_aspect('equal', 'datalim')
    plt.title('PCA projection', fontsize=12)
    cbar = plt.colorbar() # aspect= ... for width
    cbar.set_label(variable, rotation=270)
    cbar.ax.get_yaxis().labelpad = 15         
    plt.savefig(tmp_filepath)
    plt.close()    # to not display in IPython console
    

### - Plot categorical variables 
cat_variable = "label_name"
color_dict = particles_colors_dict
c, cmap = get_c_cmap_from_color_dict(color_dict, labels=df[cat_variable])  

plt.figure(figsize=cm2inch(12,12), dpi=800)
plt.style.use('dark_background')
plt.scatter(pca_embedding[:, 0], pca_embedding[:, 1],  
            c=c, cmap=cmap,  
            marker='.', s=0.8, edgecolors='none')
plt.xlim(minmax(pca_embedding[:, 0]))
plt.ylim(minmax(pca_embedding[:, 1]))
plt.xlabel("PC1")
plt.ylabel("PC2")
# plt.gca().set_aspect('equal', 'datalim')
box = plt.gca().get_position() # Shrink axis on the left to create space for the legend
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(handles=get_legend_handles_from_colors_dict(color_dict),  
           title="Snow particles",  
           bbox_to_anchor=(1,0.5), loc="center left")
plt.title('PCA projection', fontsize=12)       
plt.savefig(Path(tmp_folderpath, "PCA_Particles.png"))
plt.close()

#------------------------------------------------------------------------------.
###########################################
## Dimension reduction with UMAP in 2D ####
###########################################
### - Plot continuous variables 
for variable in variables_viz:
    tmp_folderpath =  Path(figs_path,"UMAP")
    tmp_folderpath.mkdir(exist_ok=True) # create dir if does not exist
    tmp_filename = "UMAP_" + str(variable) + ".png"
    tmp_filepath = Path(tmp_folderpath, tmp_filename)
    
    # Plot scatter 
    plt.figure(figsize=(12,12), dpi=800)
    plt.style.use('dark_background')
    plt.scatter(umap_embedding[:, 0], umap_embedding[:, 1],  
                c=df[[variable]].to_numpy().squeeze(), cmap='Spectral',  
                marker='.', s=0.8, edgecolors='none')
    plt.xlim(minmax(umap_embedding[:, 0]))
    plt.ylim(minmax(umap_embedding[:, 1]))
    plt.xlabel("L1")
    plt.ylabel("L2")
    # plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection', fontsize=12)
    cbar = plt.colorbar() # aspect= ... for width
    cbar.set_label(variable, rotation=270)
    cbar.ax.get_yaxis().labelpad = 15         
    plt.savefig(tmp_filepath)
    plt.close()    # to not display in IPython console

### - Plot categorical variables 
cat_variable = "label_name"
color_dict = particles_colors_dict
c, cmap = get_c_cmap_from_color_dict(color_dict, labels=df[cat_variable])  

plt.figure(figsize=cm2inch(12,12), dpi=800)
plt.style.use('dark_background')
plt.scatter(pca_embedding[:, 0], pca_embedding[:, 1],  
            c=c, cmap=cmap,  
            marker='.', s=0.8, edgecolors='none')
plt.xlim(minmax(umap_embedding[:, 0]))
plt.ylim(minmax(umap_embedding[:, 1]))
plt.xlabel("L1")
plt.ylabel("L2")
# plt.gca().set_aspect('equal', 'datalim')
box = plt.gca().get_position() # Shrink axis on the left to create space for the legend
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(handles=get_legend_handles_from_colors_dict(color_dict),  
           title="Snow particles",  
           bbox_to_anchor=(1,0.5), loc="center left")
plt.title('UMAP projection', fontsize=12)       
plt.savefig(Path(tmp_folderpath, "UMAP_Particles.png"))
plt.close()

 