"""
Various function to read weather data from different sources, to be eventually added to the database
"""


from numpy.core.arrayprint import dtype_is_implied
import pandas as pd
import csv
import fnmatch
import os
import netCDF4 
import numpy as np


def files_in_dir_recursive(top, pattern="*", include_dir=True):
    for (root, dirs, files) in os.walk(top):
        match_files = (fn for fn in files if 
            fnmatch.fnmatchcase(fn, pattern))
        if include_dir:
            match_files = (os.path.join(root,fn) for fn in match_files)
        for fn in match_files:
            yield fn

def concatenate_wprof(data_dir,out_file):

    """
    Function to concatenate environmental data contained in Wprof files
    into a single dataframe and store it as a pickle

    Input:
        data_dir: Wprof base path data
        out_file: output storage filename

    """

    time= np.asarray([])    # Unix time
    T   = np.asarray([])    # Temperature
    RH  = np.asarray([])    # RH [%]
    FF  = np.asarray([])    # Wind speed [m/s]
    P   = np.asarray([])    # Pressure [hPa]
    DD  = np.asarray([])    # Wind dir [°]
    #LWP = np.asarray([])    # Liquid W path [g/m3]

    for fn_full in files_in_dir_recursive(
        data_dir, pattern="*LV0.nc"):

        print(fn_full)
        ncobj = netCDF4.Dataset(fn_full)
        ncvars = ncobj.variables

        # Put variables together
        time=np.concatenate([time, np.asarray(ncvars['Time'])]).astype(int)
        T   =np.concatenate([T,  np.asarray(ncvars['Environment-temp'])-273.15])
        RH  =np.concatenate([RH, np.asarray(ncvars['Rel-humidity'])])
        FF  =np.concatenate([FF, np.asarray(ncvars['Wind-speed'])/3.6])
        P   =np.concatenate([P,  np.asarray(ncvars['Barometric-pressure'])]) 
        DD  =np.concatenate([DD, np.asarray(ncvars['Wind-direction'])])
        #LWP =np.concatenate([LWP, np.asarray(ncvars['Liquid-water-path'])])

        # Close NC file (GAN)   
        ncobj.close()
    #  Make output dataframe
    df = pd.DataFrame({'time':pd.to_datetime(time,unit='s'),
                        'T':T,
                        'RH':RH,
                        'P':P,
                        'FF':FF,
                        'DD':DD}) 

    # Save pickle
    df = df.sort_values(by='time')        # sometimes order is not respected
    df.drop_duplicates(inplace=True)      # drop eventual duplicates
    df = df.set_index(["time"])
    df.to_pickle(out_file)


class idaweb:
    """
    Reads some modified IDAWEB data 

    Modifications are done by hand in the header in order to have the variables named as
    time, T, P, RH, WS, WD

    """
    def __init__(self,fname,res='10min'):

        if res == '10min':
            format="%Y%m%d%H%M"
        elif res == '1h':
            format="%Y%m%d%H"

        self.df = pd.read_csv(fname,sep=";",parse_dates=["time"],
            date_parser=lambda x: pd.to_datetime(x, format=format),na_values='-')
        self.df = self.df.drop('stn',axis='columns')

        self.df['time']=self.df['time'].astype('str')
        self.df['time']=self.df['time'].astype('datetime64[ns]')
        self.df = self.df.set_index(["time"])

    def resample(self,sampling='60S',method='nearest'): # downsampling
        resampled = self.df.resample(sampling,label='right',closed='right')
        interpolated = resampled.interpolate(method=method)
        self.df = interpolated

class blowingsnow:
    def __init__(self,fname, dtype='csv'):
        if dtype == 'csv':
            self.df = pd.read_csv(fname)
        else:
            print("Datatype expected CSV")

        self.df['date_vec_unique']=pd.to_datetime(self.df['date_vec_unique'], dayfirst=True)
        self.df['date_vec_unique']=self.df['date_vec_unique'].astype('datetime64[ns]') 

        self.flake_uid = self.df.date_vec_unique.apply(lambda x: x.strftime('%Y%m%d%H%M%S'))+'_'+self.df.ID.apply(str)





#class Wprof:
#concatenate_wprof('/home/grazioli/tmp/Y2017/','/data/MASC_DB/rawinput/Valais-2016/Weather/Valais-2016_wprof_weather.pickle')


file='/data/MASC_DB/rawinput/Jura-2023/Weather/Jura_2023_IDAWEB.txt'
ss=idaweb(file,res='10min')

print("Hi")

#
#i = np.argmin(np.abs(df.index.to_pydatetime() - date))
#i =minute_out.index.searchsorted(date)
