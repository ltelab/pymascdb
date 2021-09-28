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
    
    tt = np.asarray(t)
    rr = np.asarray(rh)

 
    t_wetbulb = (tt * np.arctan(0.152*(rr+8.3136)**(1/2)) + np.arctan(tt + rr) - 
            np.arctan(rr - 1.6763) + 0.00391838 *(rr)**(3/2) * np.arctan(0.0231 * rr) - 4.686
            )

    return t_wetbulb

