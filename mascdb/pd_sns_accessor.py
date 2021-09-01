#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:55:56 2021

@author: ghiggi
"""
import pandas as pd
import seaborn as sns
    
def _get_sns_fun_names(): 
    sns_fun_names = ["boxplot",
                     "violinplot",
                     "boxplot",
                     "boxenplot",
                     "swarmplot",
                     "stripplot",
                     "pointplot",
                      
                     "lmplot",
                     "pairplot",

                     "scatterplot",
                     "relplot",
                     "lineplot",
                    
                     "displot",
                     "catplot",
                     "barplot",

                     "histplot",
                     "jointplot",
                     "kdeplot"]
    return sns_fun_names

def _sns_fun_factory(sns_fun_name):
    fun = getattr(sns, sns_fun_name)
    def method_fun(self, **kwargs):
        # Remove 'data' argument if specified   
        _ = kwargs.pop('data', None)
        # Plot 
        return fun(data = self._obj, **kwargs)
    return method_fun

@pd.api.extensions.register_dataframe_accessor("sns")
class SeabornAccessor:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        # Add all methods dynamically 
        sns_fun_names = _get_sns_fun_names()
        for sns_fun_name in sns_fun_names:
             setattr(type(self), sns_fun_name, _sns_fun_factory(sns_fun_name))

