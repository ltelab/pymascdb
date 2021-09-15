#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 09:52:51 2021

@author: ghiggi
"""
##----------------------------------------------------------------------------.
#### Consistency 

# Zarr dimensions
triplet_id   # LOWERCASE  (flake_id?)
cam_id       # LOWERCASE  (0,1,2)

# triplet_id currently do not have dimension values !

# flake_id <--> triplet_id ? 




##----------------------------------------------------------------------------.
np.unique(db['bs_precip_type']) # ! Do not match with https://github.com/jacgraz/pymascdb/blob/master/reader/database_reader.py#L57 
                                # What ['', means? --> JGR answer: indeed they do not match (the function in reader is older)
                                # '' means that for any reason the blowing snow estimation is not available. Not a good choice. Better None?
## Triplet_ID 
'flake_id'

### Predicted variables (average of the three cam?)
'riming_deg_level', # what is that? 
'riming_id',        # create riming_name ? or riming_degree?
'melting_id',
'melting_prob',  
'label_name',      # --> snowflake_label ?
'label_id',
'label_id_prob',

# Other descriptors (single estimate based on 3 images?)
'fallspeed' # not in cam DB 

 
# This also present in cam DB  ... same? average? else? 
'pix_size',
'Xhi',
'n_roi', 
'Dmax'    

# TODO: get_feature_list (without predicted classes, datetime, flake_id) 
# TODO: get_classes_list (riming_*, label_*, melting_*) 
 
# TODO:
# - columns with 0-1 id for image used for manual classification? 
# get_riming_name_dict(name: id) TODO: add riming_name in dataset as column!

### TODO: check images are not all 0 !!!!

##----------------------------------------------------------------------------.
#######################
### TODO IMPLEMENT ####
#######################
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

##----------------------------------
## GETTERS
"""Possibility to put together data from different cameras into a unique set"""
mascdb.ds_images() # Done 

ds_image, df_descriptors, df_class_ids = mascdb.get_triplet_descriptors_set(sample=1)


##----------------------------------
## PLOTS
# -- Minimum size of zoomed image ? To be set? 
# -- Enhancement on zoomed image? Now on full 1024x1024
#    If applied on zoomed image, define minimum size if to apply ... 
# -- Add relevant quantitative info to the plot (you mean mm? )
# -- Add sns wrapper 
##----------------------------------
## AUX 
# - Built in function to get verbose explanation of variables

##----------------------------------
# CAM statistics   
# - Histograms and statistics of variables --> Summary function? min mean max std?

##----------------------------------
### Filtering 
mascdb.exclude_campaign(...) # without_campaign
mascdb.from_campaign(...)

mascdb.select_flakes('aggregate') # select_class ?
mascdb.select_rimed('medium')

mascdb.select_max(cam0.Dmax, n=10) # arrange in the background and subset first
mascdb.select_min(cam0.Dmax, n=10)
"""
 - Functions to filter on blowing snow / precip
 - Filter on quality (Xhi for example)
 - Functions to filter on any parameter or combination of parameters --> isel() is enough? 
"""

##----------------------------------
### TUTORIALS 
# - To be finalized 

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




 

 
 