#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 11:23:36 2021

@author: grazioli
"""

import numpy as np

def wet_bulb_t(t,rh):
    """
    

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    rh : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    


    return (t * np.arctan(0.152*(rh+8.3136)**(1/2)) + np.arctan(t + rh) - 
            np.arctan(rh - 1.6763) + 0.00391838 *(rh)**(3/2) * np.arctan(0.0231 * rh) - 4.686
            )

