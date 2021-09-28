#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 09:52:51 2021

@author: ghiggi
"""
##----------------------------------------------------------------------------.
############ 
### JGR ####
############
### CAM DB 
# - Appropriate rounding of descriptors  --> DONE. To be checked for errors

## PLOTS SNOWFLAKES
# -- Add flake_ids to plot_flake, plot_flakes, plot_triplets? 
# -- Add legend |--| 1 mm over the image 
# -- Add value of a descriptor in image corner 

# add option to report mm ... 
# pix_size = self._triplet.loc[index].pix_size*1e3  # mm
# xsize = out.shape[0]*pix_size # mm
# ysize = out.shape[1]*pix_size # mm

### TO DECIDE: 
# get_vars_class --> get_vars_labels ? # include prob and riming_deg_level?
# mascdb.labels  # this return columns of get_vars_class()

# Default event definition: currently time interval of 4h without images ... to reduce ??
               
##----------------------------------------------------------------------------.
############ 
### GG  ####
############
# Working example of computing descriptors on full dataset 

##----------------------------------------------------------------------------.
##################
### TUTORIALS ####
##################
# - To be finalized !

#################
### Optional ####
#################
##----------------------------------
## GETTERS
# Ensure the same order 
mascdb.get_ds_images() # Done 
mascdb.get_full_db  (cam_id , campaign)

ds_image, df_descriptors, df_class_ids = mascdb.get_triplet_descriptors_set(sample=1)

##----------------------------------
## PLOTS SNOWFLAKES
# -- Minimum size of zoomed image ? To be set? 
# -- Enhancement on zoomed image? Now on full 1024x1024
#    If applied on zoomed image, define minimum size if to apply ... 

##----------------------------------------------------------------------------.







 

 
 
