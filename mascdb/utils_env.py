#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 11:23:36 2021

@author: grazioli
"""

import numpy as np

def wet_bulb_t(t,rh):
    """
    Returns Wet bulb temperature estimated from T and RH

    Parameters
    ----------
    t : float/int scalar, list or numpy.ndarray
        Temperature in degree Celsius
    rh : float/int scalar, list or numpy.ndarray
        Relative humidity in percentage

    Returns
    -------
    TYPE
        Wet bulb temperature in °C, same data type as as t/rh

    """
    
    if not isinstance(t,(int, float, list, np.ndarray)):
        raise TypeError("'t' must be either (list of) integers or float.")
        
    if not isinstance(rh,(int, float, list, np.ndarray)):
        raise TypeError("'rh' must be either (list of) integers or float.")

    if isinstance(t,(list,np.ndarray)):
        if len(t) != len(rh):
            raise TypeError("'t' and 'rh' must have the same length")

    t_wetbulb = (t * np.arctan(0.152*(rh+8.3136)**(1/2)) + np.arctan(t + rh) - 
            np.arctan(rh - 1.6763) + 0.00391838 *(rh)**(3/2) * np.arctan(0.0231 * rh) - 4.686
            )

    return t_wetbulb

