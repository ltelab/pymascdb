#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 09:52:51 2021

@author: ghiggi
"""

import masc_api as masc




mascdb = masc.open(fpath, lazy=True)
 
mascdb.plot(camera=1, histogram_equalization=False)
mascdb.plot_triplet()


mascdb.plt_feature static? 
mascdb.eda? 

# Filter
# - 
mascdb.filter(...)
mascdb.filter_rimed( )
mascdb.crystal_set

# Get image triplet and features 
image_ds, feature_db = mascdb.image_feature_set(sample=1)