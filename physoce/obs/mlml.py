# -*- coding: utf-8 -*-

'''
Tools for working with data from MLML public data portal:
http://pubdata.mlml.calstate.edu
'''

from urllib2 import urlopen
import re
import os

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