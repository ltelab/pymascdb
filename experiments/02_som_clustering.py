#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:03:56 2019

@author: ghiggi
"""
import os
os.chdir("/home/ghiggi/Projects/pymascdb")
#os.chdir("/home/grazioli/CODES/python/pymascdb")
import collections
import somoclu
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

from mascdb.api import MASC_DB
from mascdb.utils_img import xri_zoom
from mascdb.utils_img import _get_zoomed_image
from mascdb.utils_figs import cm2inch
from mascdb.utils_figs import get_colors_from_cmap
from mascdb.aux import get_vars_cam_descriptors
from mascdb.aux import get_vars_cam_descriptors
from mascdb.aux import get_vars_class_ids
from mascdb.aux import get_vars_class_names
from mascdb.aux import get_snowflake_class_name_colors_dict
from mascdb.aux import get_riming_class_name_colors_dict
from mascdb.aux import get_snowflake_class_id_colors_dict
from mascdb.aux import get_riming_class_id_colors_dict
from mascdb.aux import get_campaign_colors_dict

##----------------------------------------------------------------------------.
### Specify MASCDB directory 
dir_path = "/media/ghiggi/New Volume/Data/MASCDB"
dir_path = "/ltenas3/MASC_DB/"

figs_path = "/home/ghiggi/Projects/pymascdb/figs/SOM_Clustering"
if not os.path.exists(figs_path):
    os.makedirs(figs_path)
    
##----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Define snoflake descriptors to use for clustering 
# - They could be latent codes from an InfoGAN ...
cam_descriptors = get_vars_cam_descriptors()
cam_descriptors = ['n_roi', 'area','perim','Dmax','area_porous','compactness',
                   'bbox_width','bbox_len','solidity','nb_holes','complexity']

# Retrieve snowflake descriptors  
X = mascdb.cam0[cam_descriptors]  
 
### Standardize data matrix used for clustering
scaler = MinMaxScaler()
scaler.fit(X)
X_std = scaler.transform(X)

#-----------------------------------------------------------------------------.
#################### 
## SOM training ####
####################
# Define SOM grid 
n_rows, n_columns = (20,20)

# Sample some initial codebooks 
# - To speed up training and convergence
n_codebooks = n_rows*n_columns
n_features = X_std.shape[1]
sample_idxs = np.random.choice(len(X_std), size=n_codebooks)
initial_codebook = X_std[sample_idxs,:] 
initial_codebook = initial_codebook.reshape(n_rows, n_columns, n_features)
initial_codebook.shape

# Define SOM structure 
# - maptype : "toroid","planar"
# - gridtype . "hexagonal", "rectangular"
som = somoclu.Somoclu(n_columns=n_columns, n_rows=n_rows,  
                      gridtype='rectangular', maptype='planar',
                      initialcodebook=initial_codebook) 

# Train SOM
som.train(data = X_std,
          epochs = 50, # 50,  
          radius0 = 0, radiusN = 1,  
          scale0 = 0.5, scaleN = 0.001)

#-----------------------------------------------------------------------------.
## Retrieve Umatrix 
# Umatrix = som.umatrix 
# Umatrix.shape

# # Retrieve the activation map  
# # - Euclidean distance between input and codebooks  
# ActivationMatrix = som.activation_map
# # ActivationMatrix.shape  # (nobs x n_nodes)

# ## Retrieve position of input obs in the map  (BMUs = Best Matching Units) 
# node_assignement = som.bmus 
# node_assignement.shape

# ## Get distance between new data and codebooks 
# # newdata = X_std
# # ActivationMatrix_new = som.get_surface_state(data=newdata) # ! This is slow 

# ## Get BMUs indexes of the activation map 
# # new_data_assignement = som.get_bmus(ActivationMatrix_new)

# ##-----------------------------------------------------------------------------.
# ### Provide new data and update the som (--> online learning)
# # som.update_data(data=newdata)
# # som.train()

# ##-----------------------------------------------------------------------------.
# ## Display the node centroid/codebook values for each of the dimension
# som.view_component_planes([0]) # cam_descriptors[0]

# ## Display activation map (of one obs) (aka distance to codebooks)
# som.view_activation_map(data_index=1) 

# ## Display U-matrix
# colors = ["red"] * 50
# colors.extend(["green"] * 50)
# colors.extend(["blue"] * 50)
# som.view_umatrix(bestmatches=True, bestmatchcolors=colors)

## Display the similarity matrix 
# # som.view_similarity_matrix()  # Require few data --> len(data)²/2 array 

##----------------------------------------------------------------------------.
############################################################
## Plot MASC grid with neighbour distance color on axis ####
############################################################
# Define custom distance
def dist(x,y):   
    return np.sqrt(np.sum((x-y)**2))

# Retrieve distance between SOM neighbours using SOM codebooks
# - Average distance across all features ... maybe extend to specific columns of codebooks ... 

def retrieve_nearest_neighbor_distance(codebooks, row, col, dist):
    dim_codebooks = codebooks.shape  # (nrow, ncol, dim_variables)
    # ----------------------------------------------------------------------.
    # Above (top)
    tmp_row = row - 1
    tmp_col = col
    if tmp_row < 0 or row == 0:
        above_dist = float('NaN')  # np.nan
    else: 
        above_dist = dist(codebooks[row,col,:], codebooks[tmp_row, tmp_col,:])
    # ----------------------------------------------------------------------.
    # Below (bottom)
    tmp_row = row + 1
    tmp_col = col
    if tmp_row > (dim_codebooks[0] - 1) or row == (dim_codebooks[0]-1):
        below_dist = float('NaN')
    else: 
        below_dist = dist(codebooks[row,col,:], codebooks[tmp_row, tmp_col,:])
    # ----------------------------------------------------------------------.
    # Right
    tmp_col = col + 1 
    tmp_row = row
    if tmp_col > (dim_codebooks[1] - 1) or col == (dim_codebooks[1]-1):
        right_dist = float('NaN')
    else: 
        right_dist = dist(codebooks[row,col,:], codebooks[tmp_row, tmp_col,:])
    # ----------------------------------------------------------------------.
    # Left 
    tmp_col = col - 1
    tmp_row = row
    if tmp_col < 0 or col == 0:
        left_dist = float('NaN')
    else: 
        left_dist = dist(codebooks[row,col,:], codebooks[tmp_row, tmp_col,:])
    # ----------------------------------------------------------------------.    
    d = dict([("top", above_dist), ("bottom",below_dist),("right", right_dist),("left",left_dist)])    
    return(d)
 
#------------------------------------------------------------------------------.
# Retrieve SOM node to which each image is associated  
node_assignement = som.bmus 
node_assignement.shape

# Retrieve image example for each node 
nodes_selected, img_idx_sample = np.unique(node_assignement, return_index=True, axis=0)

# Retrieve feature example for each node (nodes centroids)
codebooks = som.codebook
codebooks.shape

# Retrieve images 
images = mascdb.da.isel(flake_id=img_idx_sample, cam_id=0)
images.shape
images = images.compute()

# Zoom images (to the same zoom level)
imgs_zoomed = xri_zoom(images, squared=True)
imgs_zoomed = imgs_zoomed.values

# Define colors limits for distance 
vmin = 0
vmax = 0.08

# l_values = []
# for (k, (i,j)) in enumerate(nodes_selected):
#     l_values.append(list(retrieve_nearest_neighbor_distance(codebooks,row=i, col=j, dist=dist).values()))
# plt.hist(np.array(l_values).flatten())

# Retrieve distance between SOM neighbours and assign colors
dist_color_dict = collections.defaultdict(lambda : collections.defaultdict(dict))
for (k, (i,j)) in enumerate(nodes_selected):
    dist_color_dict[i][j]['dist'] = retrieve_nearest_neighbor_distance(codebooks,row=i, col=j, dist=dist)
    dist_color_dict[i][j]['color'] = get_colors_from_cmap(dist_color_dict[i][j]['dist'], 
                                                          cmap_name='Spectral', 
                                                          vmin=vmin, vmax=vmax,
                                                          nan_color = 'black' )      
    
## Define grid layout 
gs = gridspec.GridSpec(n_rows,n_columns)
gs.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 

## Plot SOM maps with example snowflakes 
fig = plt.figure(figsize=cm2inch(15, 15), dpi=400)

for (k, (i,j)) in enumerate(nodes_selected):
    # Plot the snowflakes 
    ax = fig.add_subplot(gs[i,j])
    ax.imshow(imgs_zoomed[k,:,:], cmap="gray", vmin=0, vmax=255)
    
    # Set the color of the axis 
    for position in ['bottom','left','top','right']:
        ax.spines[position].set_color(dist_color_dict[i][j]['color'][position])
    
    # Disable axis 
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False) 
  
# Save the figure     
fig.savefig(os.path.join(figs_path, 'MASC_SOM_Cluster.png'))

#------------------------------------------------------------------------------. 
################################
## Clustering the SOM nodes ####
################################
# from sklearn.cluster import DBSCAN
# algorithm = DBSCAN()
# som.cluster(algorithm=algorithm)

# plt.imshow(som.clusters)
# som.view_umatrix(bestmatches=True)

