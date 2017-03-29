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
    
def compasstransform(theta):
    '''
    Converts angles between compass direction (clockwise from true North) to direction in polar coordinates (counter-clockwise from x-axis, pointing East). 
    
    Note that regardless of which way the conversion is being done (compass -> polar, or polar -> compass), the output of this function will be the same for a given input.
    
    INPUT: 
    - theta: direction (degrees), numpy array
        
    OUTPUT:
    - converted direction (degrees), numpy array of same size
        (0 to 360 degrees)
    '''
    theta = np.array(theta)
    theta = theta*np.pi/180. # convert to radians
    x = -np.sin(-theta)
    y = np.cos(-theta)
    theta_out = np.arctan2(y,x)
    theta_out = theta_out*180/np.pi # convert back to degrees
    neg = theta_out < 0
    theta_out[neg] = theta_out[neg]+360
    return theta_out
