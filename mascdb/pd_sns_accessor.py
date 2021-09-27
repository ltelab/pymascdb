#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 21:55:56 2021

@author: ghiggi
"""
import numpy as np 
import matplotlib.pyplot as plt 
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
             
    def corrplot(self, vars=None, vmin=-.3, vmax=.3, center=0,
                 cbar_kws={"shrink": .5},  linewidths=.5):
        # https://seaborn.pydata.org/examples/many_pairwise_correlations.html
        if vars is not None:
            df = self._obj.loc[:, vars]
        else:
            df = self._obj
        # Compute the correlation matrix
        corr = df.corr()

        # Generate a mask for the upper triangle
        mask = np.triu(np.ones_like(corr, dtype=bool))

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11, 9))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(230, 20, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(corr, mask=mask,
                    cmap=cmap, vmin=vmin, vmax=vmax, center=center,
                    cbar_kws = cbar_kws,
                    square=True, 
                    linewidths=linewidths)
        return f 
    
    def kde_marginals(self, x, y, xlim=None, ylim=None, space=0, thresh=0, 
                      levels=100, cmap="rocket",
                      hist_color = "#03051A", hist_alpha=1,hist_bins=25):
        # https://seaborn.pydata.org/examples/smooth_bivariate_kde.html
        df = self._obj
        g = sns.JointGrid(data=df, x=x, y=y, space=space)
        g.plot_joint(sns.kdeplot,
                     fill=True, 
                     clip=(xlim, ylim),
                     thresh=thresh, levels=levels, cmap=cmap)
        g.plot_marginals(sns.histplot, bins=hist_bins,
                         color=hist_color, alpha=hist_alpha)
        return g
    
    def kde_ridgeplot(self, x, group, pal=None, bw_adjust=.5, height=.5,
                      aspect=15, hspace=-.25, linewidth=2):
        # https://seaborn.pydata.org/examples/kde_ridgeplot.html
        # sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
 
        # Reorder dataframe by group 
        df = self._obj
        df = df[[x, group]]
        # m = df[group].map(ord)
        # df["x"] += m

        # Initialize the FacetGrid object
        if pal is None: 
            pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
        
        g = sns.FacetGrid(df, row=group, hue=group, 
                          aspect=aspect, height=height, palette=pal)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, x,
              bw_adjust=bw_adjust,
              linewidth=linewidth,
              clip_on=False,
              fill=True, alpha=1,
              )
           
        g.map(sns.kdeplot, x,
              bw_adjust=bw_adjust, lw=2,
              clip_on=False, color="w", 
              )
        
        # Passing color=None to refline() uses the hue mapping
        g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

        # Define and use a simple function to label the plot in axes coordinates
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, .2, label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes)

        g.map(label, x)

        # Set the subplots to overlap
        g.figure.subplots_adjust(hspace=hspace)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[], ylabel="")
        g.despine(bottom=True, left=True)
        
        # Return the object 
        return g 