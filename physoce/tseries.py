import numpy as np
from scipy import stats, signal

def hanning(x,N):
    """ 
Filter a time series x with a Hanning window of length N. If x is 2D, the time series will be filtered along columns.
    
Inputs:
x - a numpy array to be filtered
N - width of window
    
Output: numpy array of filtered time series
    """
    
    wts = signal.hann(N) # filter weights
    xf = _filt(x,wts)
    return xf
    
def lancz(x,dt=1,T=40):
    """ 
Filter a time series x with cosine-Lanczos filter. If x is 2D, the time series will be filtered along columns.

The default half amplitude period of 40 hours corresponds to a frequency of 0.6 cpd. A half amplitude period of 34.29h corresponds to 0.7 cpd. The 40 hour half amplitude period is more effective at reducing diurnal-band variability but shifts periods of variability in low passed time series to >2 days.
    
Inputs:
x - a numpy array to be filtered. 
dt - sample interval (hours), default = 1
T - half-amplitude period (hours), default = 40
    
Output: numpy array of filtered time series, same size as input with ends NaN values at start and end.

Reference: Emery and Thomson, 2004, Data Analysis Methods in Physical Oceanography. 2nd Ed., pp. 539-540. Section 5.10.7.4 - The Hanning window.
    """

    cph = 1./dt   # samples per hour
    nwts = int(round(120*cph)) # number of weights
    
    # create filter weights
    wts = signal.firwin(nwts, 
                        1./T, 
                        window='hanning', 
                        nyq=cph/2.)  
                        
    xf = _filt(x,wts)
    return xf
    
def pl66(x,dt=1,T=33):
    """
Filter a time series x with the PL66 filter. If x is 2D, the time series will be filtered along columns.
    
Inputs:
x - a numpy array to be filtered. 
dt - sample interval (hours), default = 1
T - half-amplitude period (hours), default = 33
    
Output: numpy array of filtered time series, same size as input with ends NaN values at start and end.    
    
Reference: Rosenfeld (1983), WHOI Technical Report 85-35
Matlab code: http://woodshole.er.usgs.gov/operations/sea-mat/bobstuff-html/pl66tn.html
    """
    
    Tn=float(T)/dt # normalized cutoff period
    fqn=1./Tn # normalized cutoff frequency
    nw = int(round(2.*T/dt)) # number of weights on one side
    
    # create filter weights
    j = np.arange(1,nw)
    tn = np.pi*j
    den=fqn*fqn*tn**3
    wts = (2*np.sin(2*fqn*tn)-np.sin(fqn*tn)-np.sin(3*fqn*tn))/den 
    
    # make symmetric
    wts = np.hstack((wts[::-1],2*fqn,wts))
    
    xf = _filt(x,wts)
    return xf
    
def _filt(x,wts):
    """
Private function to filter a time series and pad the ends of the filtered time series with NaN values. For N weights, N/2 values are padded at each end of the time series. The filter weights are normalized so that the sum of weights = 1.
   
Inputs: 

x - the time series (may be 2d, will be filtered along columns)
wts - the filter weights

Output: the filtered time series
    """
    
    # convert to 2D array if necessary (general case)
    ndims = np.ndim(x)
    if ndims == 1:
        x = np.expand_dims(x,axis=1)
    
    # normalize weights and convolve
    wtsn = wts*sum(wts)**-1 # normalize weights so sum = 1
    xf = signal.convolve(x,wtsn[:,np.newaxis],mode='same')  
    
    # note: np.convolve may be faster 
    # http://scipy.github.io/old-wiki/pages/Cookbook/ApplyFIRFilter
    
    # pad ends of time series
    nwts = len(wts) # number of filter weights
    npad = int(np.ceil(0.5*nwts))
    xf[:npad,:] = np.nan
    xf[-npad:,:] = np.nan
    
    # return array with same number of dimensions as input
    if ndims == 1:
        xf = xf.flatten()    
    return xf

def princax(u,v=None):
    '''
Principal axes of a vector time series.

Usage:
theta,major,minor = princax(u,v) # if u and v are real-valued vector components
    or
theta,major,minor = princax(w)   # if w is a complex vector

Input:
u,v - 1-D arrays of vector components (e.g. u = eastward velocity, v = northward velocity)
    or
w - 1-D array of complex vectors (u + 1j*v)

Output:
theta - angle of major axis (math notation, e.g. east = 0, north = 90)
major - standard deviation along major axis
minor - standard deviation along minor axis

Reference: Emery and Thomson, 2001, Data Analysis Methods in Physical Oceanography, 2nd ed., pp. 325-328.
Matlab function: http://woodshole.er.usgs.gov/operations/sea-mat/RPSstuff-html/princax.html
    '''  
    
    # if one input only, decompose complex vector
    if v is None:
        w = np.copy(u)
        u = np.real(w)
        v = np.imag(w)
        
    # Make sure inputs are numpy arrays
    u = np.array(u)
    v = np.array(v)  
    
    # only use finite values for covariance matrix
    ii = np.isfinite(u+v)
    uf = u[ii]
    vf = v[ii]    
    
    # compute covariance matrix
    C = np.cov(uf,vf)
    
    # calculate principal axis angle (ET, Equation 4.3.23b)
    theta = 0.5*np.arctan2(2.*C[0,1],(C[0,0] - C[1,1])) * 180/np.pi
    
    # calculate variance along major and minor axes (Equation 4.3.24)
    term1 = C[0,0] + C[1,1]
    term2 = ((C[0,0] - C[1,1])**2 + 4*(C[0,1]**2))**0.5
    major = np.sqrt(0.5*(term1 + term2))
    minor = np.sqrt(0.5*(term1 - term2))
    
    return theta,major,minor
    
def rot(u,v,theta):
    """
Rotate a vector counter-clockwise OR rotate the coordinate system clockwise. 

Usage:
ur,vr = rot(u,v,theta)

Input: 
u,v - vector components (e.g. u = eastward velocity, v = northward velocity)
theta - rotation angle (degrees)

Output:
ur,vr - rotated vector components

Example: 
rot(1,0,90) returns (0,1)
    """

    # Make sure inputs are numpy arrays
    u = np.array(u)
    v = np.array(v)    
    
    w = u + 1j*v            # complex vector
    ang = theta*np.pi/180   # convert angle to radians
    wr = w*np.exp(1j*ang)  # complex vector rotation
    ur = np.real(wr)        # return u and v components
    vr = np.imag(wr)
    return ur,vr

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
    
if __name__ == '__main__':
    
    # Test princax function
    u = np.array([1,2,4,5,np.nan])
    v = np.array([1,2,3,5,6])
    theta,major,minor = princax(u,v)
    theta,major,minor = princax(u+1j*v)
    mat_theta = 43.0138 # From Matlab output
    mat_major = 2.4763
    mat_minor = 0.3432
    test = np.isclose(np.array([theta,major,minor]),
                      np.array([mat_theta,mat_major,mat_minor]),
                          atol = 1e-4)             
    if test.all():
        print('princax test: passed')
    else:
        raise ValueError('princax test: failed')
        
    # Test rot function
    x = [1,0,-1]
    y = [0,1,0]
    xr,yr = rot(x,y,90)
    test1 = np.isclose(xr,[0,-1,0])
    test2 = np.isclose(yr,[1,0,-1])
    if test1.all() & test2.all():
        print('rot test: passed')
    else:
        raise ValueError('rot test: failed')