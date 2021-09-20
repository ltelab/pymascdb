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