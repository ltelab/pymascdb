#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SCRIPT providing step-by step examples to perform investigations of the orientation
habit of falling snowflakes. 

With this script, part of the results presented in Grazioli et al, GRL
can be reproduced 


@author: grazioli
"""

# Adapt paths according to local users
import os
os.chdir("/home/grazioli/CODES/python/pymascdb") 

import sys
sys.path.insert(1,'/home/grazioli/CODES/python/py-masc-gan-3Deval/')

dir_path  = '/home/grazioli/tmp/MASC_DB/'
out_path  = '/home/grazioli/Documents/Publications/Grazioli_GRL_2024/Raw_images/'


from statsmodels.graphics.gofplots import qqplot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib import rc
from matplotlib import gridspec
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


import seaborn as sns
from scipy.stats import norm
import mascdb.api
from mascdb.api import MASC_DB
from scipy import stats


#-----------------------------------------------------------------------------
# Some utility functions before the main script

def plot_kde_ori(x,
                  hue=None,
                  title='',
                  xlabel='Orientation [°]'):
    
    # Create a figure and axis
    plt.figure(figsize=(10, 6))
    
    # Customize the grid
    sns.set_style('whitegrid')

    # Plot KDE for each column in the DataFrame with different line styles
    sns.kdeplot(x=x,
                 hue=hue,
                 fill=0,
                 common_norm=False)            

    # Set axis limits and labels
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel('Hist. Density [-]', fontsize=14)

    # Add a legend
    plt.legend(fontsize=14)

    # Add a title
    plt.title(title, fontsize=16)

    # Show the plot
    plt.show()

def plot_hist_ori(x,
              hue=None,
              kde=False,
              filename='Histogram_Orientation.png',
              xlabel='Orientation [°]',
              binwidth=5,
              hue_order=None,
              plot_gaussian=False): 

    # Set style and font sizes
    sns.set(rc={'axes.labelsize': 16, 
                'axes.titlesize': 16,
                'xtick.labelsize': 16,
                'ytick.labelsize': 16})

    sns.set_style('whitegrid')

    # Create a figure with specified size and layout
    plt.figure(figsize=(8, 6), tight_layout=True)

    # Create the histogram with customized aesthetics
    ax = sns.histplot(x=x, 
                      common_norm=False,
                      stat='density',
                      kde=kde,
                      hue=hue,
                      binwidth=binwidth,
                      hue_order=hue_order,
                      line_kws={'lw': 2.5, 'linestyle': '-'},
                      fill=0,
                      element='step',
                      lw=1)

    # Set titles and labels
    plt.xlabel('Orientation [°]')
    plt.ylabel('Density')
    plb.setp(ax.get_legend().get_texts(), fontsize='14') # for legend text

    # Customize grid lines
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Customize tick parameters
    ax.tick_params(axis='both', which='major', length=6, width=1)
    ax.tick_params(axis='both', which='minor', length=4, width=0.5)

    # Customize spine visibility
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim())
    
    if plot_gaussian:
        x0, x1 = ax.get_xlim()  # extract the endpoints for the x-axis
        x_pdf = np.linspace(x0, x1, 100)

        for i in range(5,40,5):
            y_pdf = norm.pdf(x_pdf,loc=np.mean(x),scale=i)
            if i == 5:
                plt.plot(x_pdf, y_pdf, 'gray', lw=1,alpha=0.6,
                         label='Gaussian, \u03C3 [5,40]')   
            else:
                plt.plot(x_pdf, y_pdf, 'gray', lw=1,alpha=0.6) 
    
        plt.legend()

    # Show the plot
    plt.savefig(out_path+filename,dpi=450)

    plt.show()
    
def get_wind_std(x,wind):
    
    
    wind_th  = np.asarray(range(100))*0.5+0.5
    x_std    = np.asarray(range(100))+np.nan
    x_number = np.asarray(range(100))+np.nan
    
    i=0
    for th in wind_th:
        idx = np.where(wind < th)
        x_number[i] = len(idx[0])
        x_std[i]    = np.std(x[idx])
        i=i+1
        
    return x_std, wind_th, x_number


def get_ori_a_r(mascdb,method = 'mean',absolute=False):
    
    
    if method == 'min':
        vector = [np.abs(mascdb.cam0.ell_fit_ori),
                  np.abs(mascdb.cam1.ell_fit_ori),
                  np.abs(mascdb.cam2.ell_fit_ori)]

        vector_sign =[mascdb.cam0.ell_fit_ori,
                  mascdb.cam1.ell_fit_ori,
                  mascdb.cam2.ell_fit_ori]
    
        vector_a_r = [mascdb.cam0.ell_fit_a_r,
                  mascdb.cam1.ell_fit_a_r,
                  mascdb.cam2.ell_fit_a_r]
    
        # Take minimum absolute orientation but preserve the sign
        argmin = np.argmin(vector,axis=0)
        x   = [vector_sign[argmin[z]][z] for z in range(len(argmin))  ]
        y   = [vector_a_r[argmin[z]][z] for z in range(len(argmin))  ] 
    elif method == 'ellipsoid':
        x = mascdb.triplet.ellipsoid_fit_pitch
        y = mascdb.triplet.ellipsoid_fit_c/mascdb.triplet.ellipsoid_fit_a
    elif method == 'cam0':
        x = mascdb.cam0.ell_fit_ori
        y = mascdb.cam0.ell_fit_a_r
    elif method == 'cam1':
        x = mascdb.cam1.ell_fit_ori
        y = mascdb.cam1.ell_fit_a_r
    elif method == 'cam2':
        x = mascdb.cam2.ell_fit_ori
        y = mascdb.cam2.ell_fit_a_r
    elif method == 'mean':
        vector = [mascdb.cam0.ell_fit_ori,
                  mascdb.cam1.ell_fit_ori,
                  mascdb.cam2.ell_fit_ori]
    
        vector_a_r = [mascdb.cam0.ell_fit_a_r,
                  mascdb.cam1.ell_fit_a_r,
                  mascdb.cam2.ell_fit_a_r]
        
        # Take minimum absolute orientation but preserve the sign
        x = np.mean(vector,axis=0)
        y = np.mean(vector_a_r,axis=0)
    elif method == 'median':
        vector = [mascdb.cam0.ell_fit_ori,
                  mascdb.cam1.ell_fit_ori,
                  mascdb.cam2.ell_fit_ori]
    
        vector_a_r = [mascdb.cam0.ell_fit_a_r,
                  mascdb.cam1.ell_fit_a_r,
                  mascdb.cam2.ell_fit_a_r]
        
        # Take minimum absolute orientation but preserve the sign
        x = np.median(vector,axis=0)
        y = np.median(vector_a_r,axis=0)
        
    if absolute: 
        x = np.abs(x)
        
    return x, y

#----------------------------------------------------------------------------.
### Create MASC_DB instance 
mascdb = MASC_DB(dir_path=dir_path)

# Keep only data at negative temperatures
idx = mascdb.triplet['env_T'] < 0
print(idx)
mascdb = mascdb.isel(idx) 

# Remove blowing snow
mascdb = mascdb.discard_precip_class(['blowing_snow'])

# Remove small particles
idx = mascdb.cam0['snowflake_class_name'] != 'small_particle'
mascdb = mascdb.isel(idx)
idx = mascdb.cam1['snowflake_class_name'] != 'small_particle'
mascdb = mascdb.isel(idx)
idx = mascdb.cam2['snowflake_class_name'] != 'small_particle'
mascdb = mascdb.isel(idx)

mascdb = mascdb.discard_snowflake_class('small_particle')


#-----------------------------------------------------------------------------
# 1. Plot histogram of all data

# Get orientation
x, y = get_ori_a_r(mascdb)

# Set style and font sizes
sns.set(rc={'axes.labelsize': 20, 
            'axes.titlesize': 22,
            'xtick.labelsize': 22,
            'ytick.labelsize': 22})

sns.set_style('whitegrid')

# Create a figure with specified size and layout
fig, axs = plt.subplots(2, 1, figsize=(8, 10), 
                        sharex=True,
                        gridspec_kw={'height_ratios': [1, 1]})

# Create the histogram in the first subplot
ax = sns.histplot(x=x, color='grey',
                  stat='density',
                  kde=False,
                  line_kws={'lw': 4, 'linestyle': '--'},
                  alpha=0.5,
                  binwidth=5,
                  fill=0,
                  element='bars',
                  lw=3,
                  label='Observed',
                  ax=axs[0])

# Add a Gaussian to the histogram plot
x0, x1 = ax.get_xlim()  
x_pdf = np.linspace(x0, x1, 100)
y_pdf = norm.pdf(x_pdf, loc=np.mean(x), scale=np.std(x))
ax.plot(x_pdf, y_pdf, 'grey', lw=4, label='Gaussian',linestyle='dotted') 

# Customize legend
legend = ax.legend(fontsize=16)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')

# Add kde plot in the second subplot
x, y = get_ori_a_r(mascdb)
hue = np.asarray(["Unsheltered"  for i in range(len(x))])

idx = np.where(mascdb.triplet['campaign'] == 'Davos-2015')
hue[idx]= "DFIR"

idx = np.where(mascdb.triplet['campaign'] == 'ICEPOP-2017')
hue[idx]= "DFIR"


axs[1] = sns.histplot(x=x, hue=hue,
                      lw=2,
                      binwidth=5,
                      stat='density',
                      fill=0,
                      kde=False,
                      element='step',
                      palette={'DFIR':'blue','Unsheltered':'red'},
                      common_norm=False,
                      legend=True)


# Add Gaussians
x_in = x[np.where(hue == 'DFIR')]
y_pdf = norm.pdf(x_pdf, loc=np.mean(x_in), scale=np.std(x_in))
axs[1].plot(x_pdf, y_pdf, 'blue', lw=4,
            linestyle='dotted') 

x_in = x[np.where(hue == 'Unsheltered')]
y_pdf = norm.pdf(x_pdf, loc=np.mean(x_in), scale=np.std(x_in))
axs[1].plot(x_pdf, y_pdf, 'red', lw=4,
            linestyle='dotted') 

plt.setp(axs[1].get_legend().get_texts(), fontsize='16') # for legend text
plt.xlabel('Orientation [°]')


# Show the plot
plt.tight_layout()
plt.savefig(out_path+'00_AllData_orientation_hist.png', dpi=450)
plt.savefig(out_path+'00_AllData_orientation_hist.pdf')
plt.show()


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 2. Plot an investigation of the effect of horizonatl winds on distribution

idx = mascdb.triplet['env_FF'] > 0
mascdb = mascdb.isel(idx) 

x, y = get_ori_a_r(mascdb)


# Binned statistics for unshletered campaigns
mascdb_in = mascdb.discard_campaign(['Davos-2015','ICEPOP-2018'])  
x, y = get_ori_a_r(mascdb_in)

bin_sd_un, bin_edges, binnumber = stats.binned_statistic(
                        np.asarray(mascdb_in.triplet.env_FF),
                        np.asarray(x),
                        statistic='std',
                        range=(0,20),
                        bins=30)

bin_count_un, bin_edges, binnumber = stats.binned_statistic(
                        np.asarray(mascdb_in.triplet.env_FF),
                        np.asarray(x),
                        statistic='count',
                        range=(0,20),
                        bins=30)

bins = [np.mean(bin_edges[i:i+2]) for i in range(len(bin_edges)-1)]


# Binned statistics for sheltered campaigns
mascdb_in = mascdb.select_campaign(['Davos-2015','ICEPOP-2018']) 
x, y = get_ori_a_r(mascdb_in)

bin_sd_sh, bin_edges, binnumber = stats.binned_statistic(
                        np.asarray(mascdb_in.triplet.env_FF),
                        np.asarray(x),
                        statistic='std',
                        range=(0,20),
                        bins=30)

bin_count_sh, bin_edges, binnumber = stats.binned_statistic(
                        np.asarray(mascdb_in.triplet.env_FF),
                        np.asarray(x),
                        statistic='count',
                        range=(0,20),
                        bins=30)


sns.set_style('whitegrid')

# Create a figure with specified size and layout
plt.figure(figsize=(8, 8))

# set height ratios for subplots
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.5]) 
ax0 = plt.subplot(gs[0])

ax0.plot(bins, bin_sd_sh, color='b', marker='o', markersize=8, label='DFIR', lw=2)
ax0.plot(bins, bin_sd_un, color='r', marker='v', markersize=8, label='Unsheltered', lw=2)
ax0.set_ylabel('$\sigma$($\Theta$) [$^\circ$]')
ax0.legend(fontsize=16)
ax0.set_ylim(25,45)


ax1 = plt.subplot(gs[1], sharex=ax0)
ax1.plot(bins, bin_count_sh, color='b', marker='o', markersize=8, label='DFIR', lw=2)
ax1.plot(bins, bin_count_un, color='r', marker='v', markersize=8, label='Unsheltered', lw=2)
ax1.set_yscale("log")
ax1.set_ylabel('Points per bin [-]')
ax1.set_xlabel('Wind speed on-site [ms$^{-1}$]')
ax0.set_xticks(np.arange(0, 21, 2))  # Adjust as needed


yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

# Hide x-axis for top subplot
plt.setp(ax0.get_xticklabels(), visible=False)

# Adjust layout to remove space between subplots
plt.subplots_adjust(hspace=0.1)

plt.savefig(out_path+'02_sigma_wind.png',dpi=450)
plt.savefig(out_path+'02_sigma_wind.pdf')
plt.show()

# ----------------------------------------------------------------------------
# 3. Plot only the sheltered for low wind, stratified by hydro type
idx = mascdb.triplet['env_FF'] < 4
mascdb_in = mascdb.isel(idx) 
mascdb_in = mascdb_in.discard_snowflake_class('columnar_planar_combination')


mascdb_in=mascdb_in.select_campaign(['Davos-2015','ICEPOP-2018'])

x, y = get_ori_a_r(mascdb_in)

# Get the hue on axis ratio
# 0.65 is roughly Q50%
hue = np.asarray(['a_r > 0.65'  for i in range(len(x))])
y1  = np.asarray(y)

idx = np.where(y1 < 0.65)
hue[idx]='a_r < 0.65'


rc('text', usetex=False)

plt.figure(figsize=(10, 8))
ax = sns.violinplot(x=x,
                    y=mascdb_in.triplet.snowflake_class_name,
                    split=True,
                    orient='y',
                    hue=hue,
                    inner='quart',
                    inner_kws=dict(color=".8",linewidth=2))

# Add a title and adjust ylabel
ax.set_title('DFIR, wind < 4 m/s', fontsize=25)
ax.set_xlabel('Orientation [°]', fontsize=25)
ax.set_xlim(-91,91)
ax.set_ylabel('')

# Adjust tick parameters for better visibility
ax.tick_params(axis='both', which='major', labelsize=25, width=2, length=6)

# Add grid lines
plt.grid(True, linestyle='--', linewidth=1)

# Customize legend if needed
plt.legend(fontsize='20', loc='upper right')

# Show the plot
plt.tight_layout()  # Adjust layout to prevent clipping of labels


plt.savefig(out_path+'02_hydro_axis_ratio_DFIR.png',dpi=450)
plt.savefig(out_path+'02_hydro_axis_ratio_DFIR.pdf')
plt.show()

#----------------------------------------------------------------------------
# 4. Plot only the sheltered for high wind, stratified by hydro type
idx = mascdb.triplet['env_FF'] > 8
mascdb_in = mascdb.isel(idx) 
mascdb_in = mascdb_in.discard_snowflake_class('columnar_planar_combination')


mascdb_in=mascdb_in.select_campaign(['Davos-2015','ICEPOP-2018'])

x, y = get_ori_a_r(mascdb_in)

# Get the hue on axis ratio
# 0.7 is roughly Q50%
hue = np.asarray(['a_r > 0.7'  for i in range(len(x))])
y1  = np.asarray(y)

idx = np.where(y1 < 0.7)
hue[idx]='a_r < 0.7'


rc('text', usetex=False)

plt.figure(figsize=(10, 8))
ax = sns.violinplot(x=x,
                    y=mascdb_in.triplet.snowflake_class_name,
                    split=True,
                    orient='y',
                    hue=hue,
                    inner='quart',
                    inner_kws=dict(color=".8",linewidth=2),
                    order=['graupel','aggregate','planar_crystal','columnar_crystal'])

# Add a title and adjust ylabel
ax.set_title('DFIR, wind > 8 m/s', fontsize=25)
ax.set_xlabel('Orientation [°]', fontsize=25)
ax.set_xlim(-91,91)
ax.set_ylabel('')

# Adjust tick parameters for better visibility
ax.tick_params(axis='both', which='major', labelsize=25, width=2, length=6)

# Add grid lines
plt.grid(True, linestyle='--', linewidth=0.5)

# Customize legend if needed
plt.legend(fontsize='20', loc='upper right')

# Show the plot
plt.tight_layout()  # Adjust layout to prevent clipping of labels


plt.savefig(out_path+'02_hydro_axis_ratio_DFIR_highwind.png',dpi=450)
plt.savefig(out_path+'02_hydro_axis_ratio_DFIR_highwind.pdf')
plt.show()