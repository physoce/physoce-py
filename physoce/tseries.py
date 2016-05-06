import numpy as np
from scipy import stats, signal

def hanning(x,N):
    """ 
Filter a time series x with a Hanning window of length N
    
Inputs:
x - a numpy array to be filtered
N - width of window
    
Output: numpy array of filtered time series
    """
    
    wts = signal.hann(N) # filter weights
    xf = _filt(x,wts)
    return xf
    
def lancz(x,cutoff=40,dt=1):
    """ 
Filter a time series x with cosine-Lanczos filter

The default cutoff (half power period) of 40 hours corresponds to a frequency of 0.6 cpd. A cutoff of 34.29h corresponds to 0.7 cpd. The 40 hour cutoff is more effective at reducing diurnal-band variability but shifts periods of variability in low passed time series to >2 days.
    
Inputs:
x - a numpy array to be filtered
cutoff - half-power period (hours), default = 40
dt - sample interval (hours), default = 1
    
Output: numpy array of filtered time series, same size as input with ends NaN values at start and end.

Reference: Emery and Thomson, 2004, Data Analysis Methods in Physical Oceanography. 2nd Ed., pp. 539-540. Section 5.10.7.4 - The Hanning window.
    """

    cph = 1./dt   # samples per hour
    nwts = int(120*cph) # number of weights
    
    # create the filter weights
    wts = signal.firwin(nwts, 
                        1./cutoff, 
                        window='hanning', 
                        nyq=cph/2.)  
                        
    xf = _filt(x,wts)
    return xf
    
def _filt(x,wts):
    """
Private function to filter a time series and pad the ends of the filtered time series with NaN values. For N weights, N/2 values are padded at each end of the time series. The filter weights are normalized so that the sum of weights = 1.
   
Inputs: 

x - the time series    
wts - the filter weights

Output: the filtered time series

    """
    
    wtsn = wts/sum(wts) # normalize weights so sum = 1
    xf = signal.convolve(x,wtsn,mode='same')

    # pad ends of time series
    nwts = len(wts) # number of filter weights
    npad = np.ceil(0.5*nwts) 
    xf[:npad] = np.nan
    xf[-npad:] = np.nan
    return xf

def fillgapwithnan(x,date):
    """
Fill in missing data with NaN values. This is intended for a regular time series that has gaps where no data are reported. 
    
Although this function is intended for regularly-sampled time series (with gaps), it does allow for some slight irregularity (e.g. 13:00,14:00,15:00,15:59,16:59...)
    
Inputs:
x - a numpy array of data (1 or 2 dimensions)
date - datetime values that correspond to x

Returns:
(newx,newdate)
newx - new array of data
newdate = new datetime values
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
    # stats.mode returns (value, number of occurences) in arrays
    deltat = stats.mode(alldeltat)[0][0] 
    
    gapi = np.where(alldeltat > deltat*3/2)[0]    
    
    # build new arrays
    newdate = date[0:gapi[0]+1]
    newx = np.array(x[0:gapi[0]+1,:],ndmin=2)
    cnt = 0 # counter for looping through the gaps
    for ii in gapi:
        tdiff = date[ii+1]-date[ii] # size of gap
        # number of new values needed to fill gap
        nstep = int(round((tdiff.total_seconds()/deltat.total_seconds())))-1 
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