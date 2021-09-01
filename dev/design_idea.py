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

# Filter (jgr ideas)
"""
 - Functions to filter on blowing snow / precip
 - Select field campaign(s) or exclude field campaign(s)
 - Filter on quality (Xhi for example)
 - Functions to filter on any parameter or combination of parameters
"""


# Get image triplet and features 
image_ds, feature_db = mascdb.image_feature_set(sample=1)

# (jgr ideas)
"""
- Possibility to put together data from different cameras into a unique set
- Histograms and statistics of variables
"""


# Auxiliary functions / support
"""
 - Built in function to get units
 - Built in function to get verbose explanation of variables

"""

# Plotting
"""
- Display snowflakes (with / without zoom)
 -- Display triplet or individual cam
 -- Add relevant quantitative info to the plot
 -- Image enhancement?
- Display in combination with filter (i.e. the translation of queries like "show me a nice quality aggregate at least 5 mm long" )
"""
