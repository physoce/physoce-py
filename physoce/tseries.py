import numpy as np
from scipy import stats

def fillgapwithnan(x,date):
    """
    (newx,newdate) = fillgapwithnan(x,date)
    ---------------------------------------
    Fill in missing data with NaN values. This is intended for a regular time 
    series that has gaps where no data are reported. 
    
    Although this function is inteneded for regularly-sampled time series (with gaps), 
    it does allow for some slight irregularity (e.g. 13:00,14:00,15:00,15:59,16:59...)
    
    Inputs:
    x - a numpy array of data
    date - datetime values that correspond to x

    Returns:
    (newx,newdate)
    newx - new array of data
    newdate = new datetime values
    
    Tom Connolly (tconnolly@mlml.calstate.edu)
    """
          
    xnd = np.ndim(x) # remember original number of dims in x  
    x = np.array(x,ndmin=2) # make array 2D for general use
    
    # ensure that x is oriented correctly and has one dimension with same length as date
    flipdim = False    
    if np.shape(x)[0] != np.shape(date)[0]:
        x = x.transpose()
        flipdim = True
        if np.shape(x)[0] != np.shape(date)[0]: 
            raise Exception('one dimension of x must have same length as date')
    
    # find most common timedelta
    alldeltat = np.diff(date)
    deltat = stats.mode(alldeltat)[0][0] # stats.mode returns (value, number of occurences) in arrays
    
    gapi = np.where(alldeltat > deltat*3/2)[0]    
    
    # build new arrays
    newdate = date[0:gapi[0]+1]
    newx = np.array(x[0:gapi[0]+1,:],ndmin=2)
    cnt = 0 # counter for looping through the gaps
    for ii in gapi:
        tdiff = date[ii+1]-date[ii] # size of gap
        nstep = int(round((tdiff.total_seconds()/deltat.total_seconds())))-1 # number of new values needed to fill gap
        for step in np.arange(nstep):
            t = newdate[-1]+deltat
            newdate.append(t)
        gapnans = np.nan*np.ones((nstep,np.shape(x)[1]))
        newx = np.vstack((newx,gapnans))
        if ii!=gapi[-1]:
            i1 = ii+1
            i2 = gapi[cnt+1]
            newdate = newdate+date[i1:i2+1]
            newx = np.vstack((newx,x[i1:i2+1,:]))
        else:
            newdate = newdate + date[ii+1:]
            newx = np.vstack((newx,x[ii+1:,:]))
        cnt=cnt+1

    if flipdim:
        newx = newx.transpose()        
        
    if xnd == 1:
        newx = newx.flatten() # reduce back to 1D array if necessary

    return (newx,newdate)