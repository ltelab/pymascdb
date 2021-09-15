#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 11:53:06 2021

@author: ghiggi
"""
import pandas as pd 
import numpy as np 

def _define_event_id(timesteps, timedelta_thr):
    # Check type validity 
    if not isinstance(timesteps, (list, pd.Series,np.ndarray)):
        raise TypeError("'timesteps' must be a list, pd.Series or np.array with datetime values.")
    if isinstance(timesteps, list): 
        timesteps = np.array(timesteps)
        if not np.issubdtype(timesteps.dtype, np.datetime64):
            raise TypeError("'timesteps' must have datetime values") 
    if isinstance(timesteps, pd.Series):
        timesteps = timesteps.to_numpy()
        if not np.issubdtype(timesteps.dtype, np.datetime64):
            raise TypeError("'timesteps' must have np.datetime64 dtype")
    if isinstance(timesteps, np.ndarray):
        if not np.issubdtype(timesteps.dtype, np.datetime64):
            raise TypeError("'timesteps' must have np.datetime64 dtype")
            
    if not isinstance(timedelta_thr, (np.timedelta64, pd.Timedelta)):
        raise TypeError("'timedelta_thr' must be a np.timedelta64 or pd.Timedelta object.")
    #-------------------------------------------------------------------------.
    # Check there are timesteps 
    if len(timesteps) == 0: 
        raise ValueError("No timesteps provided.")
    #-------------------------------------------------------------------------.
    # Retrieve event id 
    if len(timesteps) == 1:
        event_ids = np.array([0])
    else:
        cont_groups = np.diff(timesteps)  < timedelta_thr
        cont_groups = np.insert(cont_groups, 0, True)
        event_id = 0
        l_event_id = [] 
        is_previous_True = True
        for is_True in cont_groups:
            if not is_True:
                event_id += 1 
                l_event_id.append(event_id)
                is_previous_True = is_True
            elif is_True and not is_previous_True:
                event_id += 1 
                l_event_id.append(event_id) 
                is_previous_True = is_True
            else: # is_True and is_previous_True:
                l_event_id.append(event_id)
                is_previous_True = is_True
        event_ids = np.array(l_event_id)
    return event_ids      

def _get_timesteps_duration(timesteps, unit="s"):
    duration = np.max(timesteps) - np.min(timesteps)
    return duration.to_numpy().astype('timedelta64['+ unit + "]")
