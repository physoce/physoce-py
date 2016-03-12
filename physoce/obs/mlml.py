# -*- coding: utf-8 -*-

from urllib2 import urlopen
import re
import os

def download_station_data(data_root_dir,station='seawater',overwrite=False):
    # remote directories
    base_url = 'http://pubdata.mlml.calstate.edu/mlml_last/'    
    station_url = base_url + '/' + station + '/'
    
    # local directory
    station_dir = data_root_dir + '/' + station + '/'
    
    # check whether a directory exists for this station
    if os.path.isdir(station_dir) == False:
        os.makedirs(station_dir)
    
    # find names of csv files that exist on the web and create a list
    # the csv filenames are in format yyyy-mm.csv
    urlpath =urlopen(station_url)
    html_string = urlpath.read().decode()
    pattern = re.compile('[0-9][0-9][0-9][0-9]-[0-9][0-9].csv')
    csv_list = pattern.findall(html_string)
    urlpath.close()
    
    # get updated readme file
    urlr = urlopen(station_url + '1_README.TXT')
    fr = open(station_dir + '1_README.TXT','w')
    fr.write(urlr.read())
    fr.close()
    urlr.close()
    
    # loop through each remote csv file and download if necessary
    for csv_name in csv_list:
        remote_file = station_url + csv_name
        local_file = station_dir + csv_name
        if (os.path.exists(local_file) == False) or (overwrite == True):
            urlfile = urlopen(remote_file)
            f = open(local_file,'w')
            f.write(urlfile.read())
            f.close()
            urlfile.close()