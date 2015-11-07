import numpy as np
from scipy import stats

def fillgapwithnan(x,date):
    """
    (newx,newdate) = fillgapwithnan(x,date)
    ---------------------------------------
    Fill in missing data with NaN values. This is intended for a regular time 
    series that has gaps where data or even missing values are reported.
    
    Inputs:
    x - a numpy array of data
    date - datetime values that correspond to x
    
    note: should work if x is a 2D array but this has not been tested yet
    Tom Connolly (tconnolly@mlml.calstate.edu)
    """
          
    xnd = np.ndim(x) # remember original number of dims in x
    x = np.array(x,ndmin=2) # make array 2D for general use
    
    # find most common timedelta
    alldeltat = np.diff(date)
    deltat = stats.mode(alldeltat)[0][0] # stats.mode returns (value, number of occurences) in arrays
    
    # build new date array
    t = date[0]
    newdate = [date[0]]
    newx = np.array(x[:,0],ndmin=2)
    cnt = 1 # counter for the indices in the input array x (along the time dimension)
    while(t<=date[-1]):
        # if there is an actual date near t, use this date and corresponding data
        if abs(date[cnt]-t) < deltat/2:
            newdate.append(date[cnt])
            newx = np.hstack((newx,np.array(x[:,cnt],ndmin=2)))
            cnt = cnt+1
        # if there is no date near t, create a date and fill data with NaN
        else:
            newdate.append(t)
            newx = np.hstack((newx,np.nan*np.ones(np.shape(np.array(x[:,cnt],ndmin=2)))))
        t = t+deltat
        
    # put data back into original form
    if xnd == 1:
        newx = newx[0] # reduce back to 1D array if necessary
           
    return (newx,newdate)