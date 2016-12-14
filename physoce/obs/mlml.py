# -*- coding: utf-8 -*-

'''
Tools for working with data from MLML public data portal:
http://pubdata.mlml.calstate.edu
'''

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
import re
import os
import sys
import numpy as np
from glob import glob
from datetime import datetime, timedelta
from physoce import util
try:
    import pandas as pd
except ImportError:
    pass
try:
    import xarray as xr
except ImportError:
    pass

def make_netcdf(station_dir,netcdf_file,station,download=False,overwrite=False):
    """
Create a netcdf file containing MLML historical seawater or weather data. The file will be created from csv and readme files already on disk, or they can be downloaded.

INPUT:
station_dir - string specifying the location of csv files (e.g. '/home/username/data/')
netcdf_file - string specifying the location and name of netcdf file to be created (e.g. '/home/username/data/mlml_seawater.nc')
station     - either 'seawater' or 'weather' (default: 'seawater')
download    - boolean specifying whether to download new files
              (default: False)
overwrite   - boolean specifying whether to overwrite the existing files, only used if downloading new data (default: False)
    """
    
    # download new data, if specified    
    if download == True:    
        download_station_data(station_dir,station,overwrite)
    
    # read data in csv files to xarray dataset
    d = read_csv_data(station_dir,format='dataset')
    
    # specify location of readme file and add metadata to dataset
    readme_file = station_dir + '1_README.TXT'
    _add_metadata_xarray(d,station,readme_file)
    d.attrs['history'] = d.attrs['history'] + 'netcdf file created using physoce.obs.mlml.make_netcdf(station_dir'+station_dir+',netcdf_file='+netcdf_file+',station='+station+'download='+str(download)+',overwrite='+str(overwrite)+'): ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', '
    
    # create netcdf file
    d.to_netcdf(netcdf_file,mode='w')
        
def download_station_data(station_dir,station='seawater',overwrite=True):
    '''
Download all historical csv files for the MLML seawater intake or weather station. A latest version of the readme file is also downloaded. It is highly recommended to use different directories for seawater and weather, since the readme files have the same name. By default, new files are downloaded and existing files are overwritten.

INPUT:
station_dir - string specifying the local directory where you want to put 
              the data files
station     - either 'seawater' or 'weather' (default: 'seawater')
overwrite   - boolean specifying whether to overwrite the existing files 
              (default: False)
    
    '''
    # remote directories
    base_url = 'http://pubdata.mlml.calstate.edu/mlml_last/'    
    station_url = base_url + '/' + station + '/'

    # local directory
    station_dir = station_dir + '/'    
    
    # check whether a directory exists for this station
    if os.path.isdir(station_dir) == False:
        os.makedirs(station_dir)
    
    # find names of csv files that exist on the web and create a list
    # the csv filenames are in format yyyy-mm.csv
    urlpath =urlopen(station_url)
    html_string = urlpath.read().decode()
    urlpath.close()    
    file_pattern = '[0-9][0-9][0-9][0-9]-[0-9][0-9].csv'
    csv_list = re.findall(file_pattern,html_string)
    
    # get updated readme file
    urlr = urlopen(station_url + '1_README.TXT')
    fr = open(station_dir + '1_README.TXT','w')
    fr.write(str(urlr.read()))
    fr.close()
    urlr.close()
    
    # loop through each remote csv file and download if:
    # the file does not exist, the overwrite option is True or it is the last
    # file in the list (there may be new data in that file)
    for csv_name in csv_list:
        print('downloading ' + station + ': ' + csv_name)
        remote_file = station_url + csv_name
        local_file = station_dir + csv_name
        write_conditions = [os.path.exists(local_file) == False,
                            overwrite == True,
                            csv_name == csv_list[-1]]
        if any(write_conditions):
            urlfile = urlopen(remote_file)
            f = open(local_file,'w')
            filebytes = urlfile.read()
            f.write(filebytes.decode('utf8'))
            f.close()
            urlfile.close()
            
def read_csv_data(data_dir,format='dict'):
    '''
Read historical text data (.csv files) from the MLML seawater intake or weather station. The data must be stored locally, and can be downloaded automatically with the download_station_data() function.

Inputs:
data_dir - Specifies the directory where the data files are located. All files with the format yyyy-mm.csv in this directory will be read.

Options:
    format: output format
        format = 'dict' (default): dictionary
        format = 'dataframe': pandas DataFrame
        format = 'dataset': xarray DataSet

Output: dictionary, pandas DataFrame or xarray DataSet with keys/variable names taken from column headers
    '''
    
    file_list = glob(data_dir+'*.csv')
    
    # get list of variable names from header of first file
    f = open(file_list[0],'r')
    header = f.readline()
    f.close()
    header = header.strip('\r\n')
    varnames = header.split(',')
    
    #initialize dictionary with key and empty list for each variable
    d = dict()
    for ii,var in enumerate(varnames[0:]):
        d[varnames[ii]] = []
    
    # specify which columns contain numeric data
    floatcols = range(2,len(varnames))
    allcols = range(0,len(varnames))
    strcols = list(set(allcols)-set(floatcols))
    
    for file_name in file_list:
        print('reading ' + file_name)        
        
        # get numeric data, with missing values as NaN
        datamasked = np.genfromtxt(file_name,
                             skip_header=1,
                             delimiter=',',
                             missing_values='-99999',
                             usemask=True)
        data = datamasked.filled(np.nan)
        
        # get string data
        datastr = np.genfromtxt(file_name,
                         skip_header=1,
                         delimiter=',',
                         usecols=tuple(strcols),
                         dtype=str)
        
        # append data variables    
        if data.size != 0:    
            for col in floatcols:
                vname = varnames[col]
                d[vname] = np.append(d[vname],data[:,col])
            for si,col in enumerate(strcols):
                vname = varnames[col]
                d[vname] = np.append(d[vname],datastr[:,si])
    
    # create date variables
    # put in a numpy array for easy indexing
    # new variable for datetime
    dtime = np.array(util.list2date(d['utc_time'],'%Y-%m-%dT%H:%M:%SZ'))    

    # remove duplicate times
    ii = np.where(np.diff(dtime) > timedelta(0.))[0]
    dtime = dtime[ii]
    for var in varnames:
        d[var] = d[var][ii]
        
    # Try loading in pandas or xarray format if specified, default to dictionary format
    if format == 'dataset':
        if 'xarray' not in sys.modules:
            format = 'dataframe'
            print("Warning: xarray not installed, loading MLML data in pandas dataframe format instead")  
    if format == 'dataframe':
        if 'pandas' not in sys.modules:
            format = 'dict'
            print("Warning: pandas not installed, loading MLML data in dictionary format instead")
    
    if format is 'dataframe':
        # turn dictionary into pandas dataframe
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'
    elif format is 'dataset':
        # turn dictionary in xarray dataset, using dataframe as intermediate format
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'        
        d = xr.Dataset(d)
        d.attrs['history'] = 'dataset created using physoce.obs.mlml.read_csv_data: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
    else:
        # default format: dictionary containing numpy arrays
        d['dtime'] = []    
        d['dtime'] = dtime
    
    return d
    
def _add_metadata_xarray(d,station,readme_file):
    """
Add metadata to xarray dataset. Currently this adds lat and lon coordinates and puts the contents of the readme in an attribute. For the weather data, the anemometer height is also added as a coordinate.
    """    
    
    if station is 'seawater':
        d.coords['lon'] = -121.7915
        d.coords['lat'] = 36.8025
    elif station is 'weather':
        d.coords['lon'] = -121.78842
        d.coords['lat'] = 36.80040
        d.coords['z'] = 3.3
        d.coords['z'].attrs['name'] = 'anemometer height'
        d.coords['z'].attrs['units'] = 'meters'
        
    with open(readme_file) as f:
        contents = f.read()
        d.attrs['readme'] = contents
        
    d.attrs['history'] = d.attrs['history'] + 'attributes added to dataset using physoce.obs.mlml._add_metadata_xarray: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
            