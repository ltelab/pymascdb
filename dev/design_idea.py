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

### AUX 

### CAM DB 
# - Appropriate rounding of descriptors  --> DONE. To be checked for errors

### Parquet dataframe savings 
# ... strings column are currently saved as 'object' 
# ... reading then requires always conversion 

##----------------------------------------------------------------------------. 
############ 
### GG #####
############ 

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
# -- Add legend |--| 1 mm over the image 
# -- Add value of a descriptor in image corner 

## PLOTS DATAFRAME
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




 

 
 