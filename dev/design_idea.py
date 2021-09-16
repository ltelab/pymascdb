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
#### Consistency 
# - Zarr dimensions
flake_id     # LOWERCASE  (flake_id?)
cam_id       # LOWERCASE  (0,1,2)

# triplet_id currently do not have dimension values !

# flake_id <--> Triplet_ID ? 

### AUX 
# - Built in function to get verbose explanation of variables

#### Triplets DB 
# - Rounding to 3 decimals of triplets  
# quality_xhi_flake
# melting_prob
# n_roi       # senza decimali 
# Dmax_flake  
# fallspeed   
# bs_nor_angle
# bs_mix_ind

#-  Rename in triplet ... flake_*  for best guess values 
# flake_quality_xhi 
# flake_Dmax

# - Remove pix_size from triplet 

### CAM DB 
# - Appropriate rounding of descriptors  

### Add to cam
# - columns with 0-1 id for image used for manual classification? 
# - hl_snowflake
# - hl_riming
# - hl_melting

##----------------------------------------------------------------------------. 
############ 
### GG #####
############ 
##----------------------------------
### EVENTS 
# Number of events per campaign (with minimum duration threshold ... )
# Total duration of events (with minimum duration threshold ... )(it sums event_id durations)

##----------------------------------
## SETTERS 
# - Add columns to df 
mascdb.add_cam_column 
mascdb.add_cam_columns
mascdb.add_triplet_column 
mascdb.add_triplet_columns

# - Compute descriptors based on fun applied to each image 
mascdb.compute_image_descriptor(fun, fun_kwargs, force=False)
mascdb.compute_image_descriptors(fun, fun_kwargs, force=False)

##----------------------------------------------------------------------------.
##################
### TUTORIALS ####
##################
# - To be finalized 

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
## PLOTS
# -- Minimum size of zoomed image ? To be set? 
# -- Enhancement on zoomed image? Now on full 1024x1024
#    If applied on zoomed image, define minimum size if to apply ... 
# -- Add legend |--| 1 mm over the image 
# -- Add value of a descriptor in image corner 
# -- Add sns wrapper 


##----------------------------------------------------------------------------.
####################
### TODO CHECKS ####
####################
### 1
# Retrieve the event leading to largest Dmax 
timedelta_thr = np.timedelta64(2, 'h')
mascdb.define_event_id(timedelta_thr=timedelta_thr)
mascdb.arrange('cam0.Dmax', decreasing=True)        # TODO CHECK WHY THIS RAISE WARNING

event_id = mascdb.arrange('cam0.Dmax', decreasing=True).cam0['event_id']
idx_event = mascdb.cam0['event_id'] == event_id[0]
mascdb_event = mascdb.isel(idx_event).arrange('cam0.datetime', decreasing=False)
print(mascdb_event)

##----------------------------------------------------------------------------.




 

 
 