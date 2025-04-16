'''Functions for working with meteorological data'''

import numpy as np
from airsea import windstress as ws

def uv_from_spddir(wspd, wdir):
    '''Convert wind speed and direction to eastward and northward components.
    
    Inputs
    wspd - wind speed
    wdir - wind direction (wind is blowing FROM this direction, 
                           clockwise from true north)
    
    Output:
    u - eastward velocity component (TOWARDS the east)
    v - northward velocity component (TOWARDS the north)
    '''
    
    theta = np.array(wdir) # direction CW from true north
    theta = theta*np.pi/180. # convert to radians
    x = -np.sin(-theta)
    y = np.cos(-theta)
    theta_cart = np.arctan2(y, x) # direction CCW from east (Cartesian)
    u = -wspd*np.cos(theta_cart) # eastward component
    v = -wspd*np.sin(theta_cart) # northward component
    
    return u,v

def tauxy_from_uv(u, v, z=3., drag='largepond', Ta=10.):
    '''Calculate wind stress components from wind velocity components. Requires
    the python-airsea package: https://github.com/pyoceans/python-airsea
    
    Inputs:
    u (float, array-like): x-component of wind velocity
    v (float, array-like): y-component of wind velocity
    z (float): anemometer height
    drag (str): neutral drag formulation
           'largepond' <-- default
           'smith'
           'vera'
    Ta (float, array_like): optional for drag='smith', air temperature [deg C]
    
    Output:
    taux - eastward velocity component (TOWARDS the east)
    tauy - northward velocity component (TOWARDS the north)
    '''

    spd = np.sqrt(u**2 + v**2)
    taumag = ws.stress(spd, z)
    theta = np.arctan2(v, u)
    taux = taumag*np.cos(theta)
    tauy = taumag*np.sin(theta)
    
    return taux, tauy