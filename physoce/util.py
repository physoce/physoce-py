# -*- coding: utf-8 -*-
"""
General purpose functions.
"""
from datetime import datetime
import numpy as np

def list2date(datestr_list,fmt='%a %b %d %H:%M:%S %Y'):
    '''Convert a list of date strings to datetime format.
    
    INPUT:
    datestr_list: a list of strings that represent dates
    fmt: format of the date string, as would be input to strftime() or strptime()
    
    see https://docs.python.org/library/datetime.html#strftime-and-strptime-behavior
    
    OUTPUT:
    list of datetimes
    '''
    datetime_list = [datetime.strptime(datestr, fmt) for datestr in datestr_list]
    return datetime_list
    
def compass2polar(theta):
    '''
    Convert an angle from compass direction (clockwise from true North) to direction in polar coordinates (counter-clockwise from x-axis, pointing East).
    
    INPUT: 
    - theta - compass direction (degrees), numpy array
        
    OUTPUT:
    - direction in polar coordinates (degrees), numpy array of same size
    '''
    
    theta = theta*np.pi/180 # convert to radians
    x = -np.sin(-theta)
    y = np.cos(-theta)
    theta_pol = np.arctan2(y,x)
    theta_pol = theta_pol*180/np.pi # convert back to degrees
    return theta_pol
    
def polar2compass(theta):
    '''
    Convert an angle from direction in polar coordinates (counter-clockwise from x-axis, pointing East) to compass direction (clockwise from true North).
    
    INPUT: 
    - theta - compass direction (degrees), numpy array
        
    OUTPUT:
    - direction in polar coordinates (degrees), numpy array of same size
    '''
    
    theta = theta*np.pi/180 # convert to radians
    x = np.cos(theta)
    y = np.sin(theta)
    theta_comp = -np.arctan2(x,-y)
    theta_comp = theta_comp*180/np.pi # convert back to degrees
    return theta_comp
