import numpy as np
import numpy.ma as ma

def nancorr(x,y):
# r = nancorr(x,y)
#
# calculate correlation matrix, treating NaN values as missing data
	
	x_msk = ma.masked_invalid(x)
	y_msk = ma.masked_invalid(y)
	r = ma.corrcoef(x_msk,y_msk)
	return r
	
def maxcorr(x,y,**options):
# (rmax,lag,ind) = maxcorr(x,y,**'maxlag'=int(len(x)/8)):
#
# Calculate the maximum lagged correlation between two 1D arrays
# INPUTS
# x,y are 1D arrays
# OPTIONS
# 'maxlag' the maximum number of lagged correlations to calculate (default: 1/4 of array length)
# OUTPUT
# r is the maximum correlation coefficient
# lag is the lag of the maxiumum correlation (positive: y lags x)

    nrows = len(x)
    maxlag = int(np.floor(nrows/8))
    masknan = True
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
            
    ind = ma.argmax(rs)
    rmax = rs[ind] 
    lag = lags[ind]
        
    return (rmax,lag,ind)
	