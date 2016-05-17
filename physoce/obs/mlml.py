# -*- coding: utf-8 -*-

'''
Tools for working with data from MLML public data portal:
http://pubdata.mlml.calstate.edu
'''

from urllib2 import urlopen
import re
import os
import numpy as np
from glob import glob
from physoce import util
try:
    import pandas as pd
    import xarray as xr
except ImportError:
    pass

def download_station_data(station_dir,station='seawater',overwrite=True):
    '''
Download all historical csv files for the MLML seawater intake or weather station. A latest version of the readme file is also downloaded. It is highly recommended to use different directories for seawater and weather, since the readme files have the same name. By default, new files are downloaded and existing files are overwritten.

INPUT:
station_dir - string specifying the local directory where you want to put 
              the data files
station     - either 'seawater' or 'weather' (default: 'seawater')
overwrite   - boolean specifying whether to overwrite the existing files 
              (default: 'False')
    
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
    fr.write(urlr.read())
    fr.close()
    urlr.close()
    
    # loop through each remote csv file and download if:
    # the file does not exist, the overwrite option is True or it is the last
    # file in the list (there may be new data in that file)
    for csv_name in csv_list:
        remote_file = station_url + csv_name
        local_file = station_dir + csv_name
        write_conditions = [os.path.exists(local_file) == False,
                            overwrite == True,
                            csv_name == csv_list[-1]]
        if any(write_conditions):
            urlfile = urlopen(remote_file)
            f = open(local_file,'w')
            f.write(urlfile.read())
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
                         dtype='S')
        
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

    # Try loading in pandas or xarray format if specified, default to dictionary format
    if format == 'dataset':
        try:
            reload(xr)
        except NameError:
            format = 'dataframe'
            print "Warning: xarray not installed, loading MLML data in pandas dataframe format instead"  
    if format == 'dataframe':
        try:
            reload(pd)
        except NameError:
            format = 'dict'
            print "Warning: pandas not installed, loading MLML data in dictionary format instead"
    
    if format is 'dataframe':
        # turn dictionary into pandas dataframe
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'
    elif format is 'dataset':
        # turn dictionary in xarray dataset, using dataframe as intermediate format
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'        
        d = xr.Dataset(d)
    else:
        # default format: dictionary containing numpy arrays
        d['dtime'] = []    
        d['dtime'] = dtime
    
    return d
    
def add_metadata_xarray(ds,station,readme_file):
    
    if station is 'seawater':  
        ds.coords['lon'] = -121.7915
        ds.coords['lat'] = 36.8025
    
        ds['temp'].attrs['name'] = 'temperature'
        ds['temp'].attrs['units'] = 'degrees Celcius'
        ds['temp_flg'].attrs['name'] = 'temperature flag'       
        ds['temp_flg'].attrs['sensor range test'] = '-5<temp<35'
        ds['temp_flg'].attrs['user range test'] = '4<temp<22'
        ds['temp_flg'].attrs['spike test'] = 'mean(temp_1.42hours)+/-2.4*stdev'
        
        ds['otemp'].attrs['name'] = 'optode temperature'
        ds['otemp'].attrs['units'] = 'degrees Celcius'
        ds['otemp_flg'].attrs['name'] = 'optode temperature flag'       
        ds['otemp_flg'].attrs['sensor range test'] = '-5<otemp<35'
        ds['otemp_flg'].attrs['user range test'] = '4<otemp<22'
        ds['otemp_flg'].attrs['spike test'] = 'mean(otemp_1.42hours)+/-2.4*stdev'
        
        ds['cond'].attrs['name'] = 'conductivity'
        ds['cond'].attrs['units'] = 'Siemens per meter'
        ds['cond_flg'].attrs['name'] = 'conductivity flag'        
        ds['cond_flg'].attrs['sensor range test'] = '0<cond<7'
        ds['cond_flg'].attrs['user range test'] = '0.75<cond<1.1'
        ds['cond_flg'].attrs['spike test'] = 'mean(cond_1.42hours)+/-2.4*stdev'            
                
        ds['sal'].attrs['name'] = 'practical salinity'
        ds['sal'].attrs['units'] = ''
        ds['sal'].attrs['scale'] = 'PSU-78'
        ds['sal_flg'].attrs['name'] = 'salinity flag'   
        ds['sal_flg'].attrs['sensor range test'] = '20<sal<50'
        ds['sal_flg'].attrs['user range test'] = '32<sal<35'
        ds['sal_flg'].attrs['spike test'] = 'mean(sal_1.42hours)+/-2.4*stdev'
        
        ds['fluor'].attrs['name'] = 'fluorescence'
        ds['fluor'].attrs['units'] = 'micrograms per liter'    
        ds['fluor_flg'].attrs['name'] = 'fluorescence flag' 
        ds['fluor_flg'].attrs['sensor range test'] = '0.03<fluor<75'
        ds['fluor_flg'].attrs['user range test'] = '0<fluor<20'
        ds['fluor_flg'].attrs['spike test'] = 'mean(fluor_1.42hours)+/-2.4*stdev'            
        
        ds['ba'].attrs['name'] = 'beam attenuation'
        ds['ba'].attrs['units'] = 'meters^-1'    
        ds['ba_flg'].attrs['name'] = 'beam attenuation flag' 
        ds['ba_flg'].attrs['sensor range test'] = '0<ba<100'
        ds['ba_flg'].attrs['user range test'] = '0<ba<50'
        ds['ba_flg'].attrs['spike test'] = 'mean(ba_1.42hours)+/-2.4*stdev'            

        ds['trans'].attrs['name'] = 'transmission'
        ds['trans'].attrs['units'] = 'percent (%)'    
        ds['trans_flg'].attrs['name'] = 'transmission flag' 
        ds['trans_flg'].attrs['sensor range test'] = '0<trans<100'
        ds['trans_flg'].attrs['user range test'] = '0<trans<100'
        ds['trans_flg'].attrs['spike test'] = 'mean(ba_1.42hours)+/-2.4*stdev'            
        
        ds['do2'].attrs['name'] = 'dissolved oxygen'
        ds['do2'].attrs['units'] = 'micromoles per liter'    
        ds['do2_flg'].attrs['name'] = 'dissolved oxygen flag' 
        ds['do2_flg'].attrs['sensor range test'] = '0<do2<100'
        ds['do2_flg'].attrs['user range test'] = 'do2<100'
        ds['do2_flg'].attrs['spike test'] = 'mean(do2_1.42hours)+/-2.4*stdev'                    

        ds['osat'].attrs['name'] = 'oxygen saturation'
        ds['osat'].attrs['units'] = 'percent (%)'    
        ds['osat_flg'].attrs['name'] = 'oxygen saturation flag' 
        ds['osat_flg'].attrs['sensor range test'] = '0<osat<100'
        ds['osat_flg'].attrs['user range test'] = '0<trans<100'
        ds['osat_flg'].attrs['spike test'] = 'mean(ba_1.42hours)+/-2.4*stdev'            

    with open(readme_file) as f:
        contents = f.read()
        for n, line in enumerate(contents):
            n
            