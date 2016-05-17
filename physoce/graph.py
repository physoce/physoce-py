# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from matplotlib import dates
import seawater as sw

def plts(date,y):
    """
Plot a multi-year time series y as a function of time of year (rather than absolute time). The time series from each separate year is plotted in a different color, with months on the x-axis. Useful for visualizing seasonal patterns.
    
Inputs:
date - a list of datetime objects
y - a numpy array
    """
    
    sdate = []
    sy = []
    ii = 0
    for nn,t in enumerate(date):
        sdate.append(datetime(1980,t.month,t.day,t.hour,t.minute,t.second,t.microsecond))
        sy.append(y[nn])       
        
        # make sure this is not the last index in the date list
        # and check whether year is about to change
        if nn<len(date)-2: 
            if date[nn].year!=date[nn+1].year:
                plt.plot(sdate,sy)
                sdate=[]
                sy=[]      
        ii=ii+1
    # Set major x ticks on months
    ax = plt.gca()
    ax.xaxis.set_major_locator(dates.MonthLocator())
    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    
def TS_contours(SP_range,T_range,sigma_levels,**kwargs):
    ''' 
Plot contours of density anomaly (sigma) on a T-S plot. Uses EOS-80 equation of state. If the T_range input is in-situ T, then the sigma-t values are contoured. If the temperature input is potential temperature (theta), then the sigma-theta values are contoured (see Stewart 2005 for definitions).
    
INPUTS:
SP_range: Practical salinity, minimum and maximum values
T_range: In-situ or potential temperature [C], minimum and maximum values
sigma_levels: density anomaly values to contour
**kwargs: these will be passed to the contour function, see http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contour

RETURNS:
cs: the matplotlib.contour.QuadContourSet object returned by contour (can 
    be used as input to pyplot.clabel function, for example).

REQUIRES:
Seawater toolbox: https://pypi.python.org/pypi/seawater

Reference:
Stewart, R. H. (2005) Introduction to Physical Oceanography. http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_05.htm
    '''
    
    smin = SP_range[0]
    smax = SP_range[1]
    tmin = T_range[0]
    tmax = T_range[1]
    
    sgrid = np.linspace(smin,smax,101)
    tgrid = np.linspace(tmin,tmax,101)
    
    sigma = np.nan*np.zeros((len(tgrid),len(sgrid)))   
    
    # Loop to fill in grid with densities
    for i,s in enumerate(sgrid):
        for j,t in enumerate(tgrid):
            sigma[j,i]=sw.dens(s,t,0)-1000 # sigma-t
    
    # contour and return contour object
    cs = plt.contour(sgrid,tgrid,sigma,sigma_levels,**kwargs)
    return cs
    
if __name__ == '__main__':
    
    # demonstrate T_S_contours function
    plt.figure()
    cs = TS_contours([30,34],[10,20],np.arange(20,26,0.6),colors='k')
    plt.clabel(cs,fmt='%1.1f',fontsize=10)
    plt.show()