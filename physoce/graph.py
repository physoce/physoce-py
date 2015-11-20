# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
from datetime import datetime
from matplotlib import dates

def plts(date,y):
    """
    Plot a multi-year time series y as a function of time of year (rather than absolute time).
    The time series from each separate year is plotted in a different color, with months on the 
    x-axis. Useful for visualizing seasonal patterns.
    
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