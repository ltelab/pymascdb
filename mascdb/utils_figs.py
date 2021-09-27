#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 22:13:58 2021

@author: ghiggi
"""
#-----------------------------------------------------------------------------.
import numpy as np
import pandas as pd 
import matplotlib as mpl 
import matplotlib.pyplot as plt

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def minmax(x):
   return([np.min(x),np.max(x)])    

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

def get_colors_from_cmap(x, cmap_name='Spectral', vmin=None, vmax=None, nan_color=None): 
   # # - mpl.cm.<Spectral> : colormaps are defined by default between 0 and 255  mpl.cm.Spectral(257)
   flag_dict = False
   # Preprocess x if dictionary
   if (isinstance(x, dict)):
       flag_dict = True
       keys = list(x.keys())
       x = np.asarray(list(x.values()))
   # Get index with NaN 
   idx_nan = np.isnan(x)  
   # Retrieve colormap 
   cmap = mpl.cm.get_cmap(cmap_name)
   # Rescale x and assign colormap values 
   rgb_val = cmap(((x - vmin)/vmax))
   # norm_fun = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
   # rgb_val = cmap(norm_fun(x)) 
   #----------------------------------------------------------------.
   ######################
   ## Convert to hex ####
   ######################
   # If only a single value
   if isinstance(rgb_val, tuple):
      rgb_hex = mpl.colors.rgb2hex(rgb_val)
   # If multiple values (matrix with rows of rgb values)
   else:
      rgb_hex = []
      for i in range(rgb_val.shape[0]):
          rgb_hex.append(mpl.colors.rgb2hex(rgb_val[i,]))
   #----------------------------------------------------------------.
   # Set back NaN
   rgb_hex = np.array(rgb_hex)
   if np.sum(idx_nan) > 0:
     rgb_hex[idx_nan] = 'NaN'
   #----------------------------------------------------------------.
   if nan_color is not None: 
       rgb_hex[idx_nan] = nan_color    
   #################################################
   # If x is dictionary --> Recreate dictionary ####
   #################################################
   if flag_dict:
      rgb_hex = dict(zip(keys, rgb_hex))
   #---------------------------------------------------------------------.
   return(rgb_hex)     

#### Test get_colors_from_cmap  
#get_colors_from_cmap(2, cmap_name="Spectral", vmin=0, vmax = 8)
#get_colors_from_cmap(np.array([2,np.nan,5]), cmap_name="Green", vmin=0, vmax = 8)
#get_colors_from_cmap(np.array([0,2,3,5]), cmap_name="Wistia", vmin=0, vmax = 5)