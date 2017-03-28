import numpy as np
import numpy.ma as ma
from scipy import stats as st

def nancorr(x,y):
    """ 
    r = nancorr(x,y)
    Calculate correlation matrix, treating NaN values as missing data
    """	
    x_msk = ma.masked_invalid(x)
    y_msk = ma.masked_invalid(y)
    r = ma.corrcoef(x_msk,y_msk)
    return r

def maxcorr(x,y,**options):
    """
    (rmax,lag,ind) = maxcorr(x,y,**'maxlag'=int(len(x)/4)):
    Calculate the maximum lagged correlation between two 1D arrays
    Inputs:
    x,y are 1D arrays
    Options
    'maxlag' the maximum number of lagged correlations to calculate (default: 1/4 of array length)
    Output:
    r is the correlation coefficient with the maximum absolute value
    lag is the lag of the maximum correlation (positive: y lags x)
    """
    
    nrows = len(x)
    maxlag = int(np.floor(nrows/4))
    if ('maxlag' in options):
        maxlag = options['maxlag']   
    
    # use masked arrays (mask NaNs)  
    x = ma.masked_invalid(x)
    y = ma.masked_invalid(y)
    
    lags = np.arange(-maxlag,maxlag+1)
    rs = np.zeros(np.shape(lags))
    for ni, lag in enumerate(lags):
        lag = lags[ni]
        if lag < 0:
            rs[ni] = ma.corrcoef(x[-lag:],y[:lag])[0,1]
        elif lag > 0:
            rs[ni] = ma.corrcoef(x[:-lag],y[lag:])[0,1]
        else:
            rs[ni] = ma.corrcoef(x,y)[0,1]
            
    ind = ma.argmax(np.abs(rs))
    rmax = rs[ind] 
    lag = lags[ind]
        
    return (rmax,lag,ind)
	
def rsig(r,nu):
    """ p-value for correlation coefficient r, and degrees of freedom nu

    INPUTS:
    r - correlation coefficient
    nu - degrees of freedom (N-2)
    
    OUTPUT:
    p - significance level/p-value
 
    significance levels of 0.05 and 0.01 correspond with Appendix E in
    Emery and Thomson (2004) Data Analysis Methods in Physical 
    Oceanography
    """

    # t value
    t = abs(r)*np.sqrt(nu)/np.sqrt(1-r**2)

    # significance level, using the "survival function" (1-cdf)
    p = 2*(st.t.sf(t,nu))

    return p
    
def rcrit(nu,sig=0.05):
    """
    Critical r (correlation coefficient), given significance level
    and degrees of freedom.
    
    INPUTS:
    nu - degrees of freedom (N-2)        
    sig - significance level (default 0.05)
    
    OUTPUT:
    rcrit - critical r value
    
    Values for 0.05 and 0.01 correspond with Appendix E in
    Emery and Thomson (2004) Data Analysis Methods in Physical 
    Oceanography
    """
    
    # critical t value (this is equivalent to Matlab tinv function)
    t = st.t.ppf(1 - sig/2,nu)
    
    # critical r value
    rc = t/np.sqrt(t**2+nu)
    
    return rc