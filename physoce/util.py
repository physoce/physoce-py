# -*- coding: utf-8 -*-
"""
General purpose functions.
"""
from datetime import datetime
import numpy as np
from matplotlib import path

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

def matlab2datetime64(datenum,unit='s'):
    '''
    Convert Matlab serial date number to NumPy datetime64 format.

    INPUTS:
    datenum - Matlab serial date number, can be array
    unit - time unit of datetime64 output (default 's')

    OUTPUT:
    array of datetime64 objects
    '''
    origin = np.datetime64('0000-01-01 00:00:00', unit) - np.timedelta64(1, 'D')
    daylength = int(np.timedelta64(1,'D')/np.timedelta64(1, unit))
    dt64 = datenum * np.timedelta64(daylength, unit) + origin
    return dt64

def haversine(lat, lon):
    """
    Calculate the great circle distances between a set of points
    on the earth (specified in decimal degrees).

    INPUTS:
    lat,lon -  1D arrays of coordinate pairs (length N)

    OUTPUT:
    great circle distances between consecutive points in km (length N-1)
    """

    lat1 = lat[:-1]
    lon1 = lon[:-1]
    lat2 = lat[1:]
    lon2 = lon[1:]

    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371 * c
    return km

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

def inside(xpt,ypt,xpoly,ypoly):
    """
    Check whether a set of (x,y) points are within a polygon.

    INPUTS:
    xpt,ypt - define a set of points (arrays of length N)
    xpoly,ypoly - define the vertices of an arbitrary polygon (arrays of length M)

    OUTPUT:
    returns an array of Booleans, length N, True where (xpt,ypt) is inside the polygon
    """

    mask_poly = np.array([0,0])
    for xx,yy in zip(xpoly,ypoly):
        mask_poly = np.vstack((mask_poly,np.array([xx,yy])))
    mask_poly = mask_poly[1:,:]
    mask_path = path.Path(mask_poly)
    points = np.transpose(np.vstack((xpt,ypt)))
    mask = mask_path.contains_points(points)
    return mask
