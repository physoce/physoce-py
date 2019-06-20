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
    nwts = int(np.round(120*cph)) # number of weights

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
    nw = int(np.round(2.*T/dt)) # number of weights on one side

    # create filter weights
    j = np.arange(1,nw)
    tn = np.pi*j
    den=fqn*fqn*tn**3
    wts = (2*np.sin(2*fqn*tn)-np.sin(fqn*tn)-np.sin(3*fqn*tn))/den

    # make symmetric
    wts = np.hstack((wts[::-1],2*fqn,wts))

    xf = _filt(x,wts)
    return xf

def pl64(x,dt=1,T=33):
    """
Filter a time series x with the PL64 filter. If x is 2D, the time series will be filtered along columns.

Inputs:
x - a numpy array to be filtered.
dt - sample interval (hours), default = 1
T - half-amplitude period (hours), default = 33

Output: numpy array of filtered time series, same size as input with ends NaN values at start and end.

Reference: CODE-2: Moored Array and Large-Scale Data Report, WHOI 85-35
    """

    Tn=float(T)/dt # normalized cutoff period
    fqn=1./Tn # normalized cutoff frequency
    nw = int(np.round(64/dt)) # number of weights on one side

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

    # normalize weights
    wtsn = wts*sum(wts)**-1 # normalize weights so sum = 1

    # Convolve using 'direct' method. In older versions of scipy, this has to
    # be specified because the default 'auto' method could decide to use the
    # 'fft' method, which does not work for time series with NaNs. In newer
    # versions, there is no method option.
    try:
        xf = signal.convolve(x,wtsn[:,np.newaxis],mode='same',method='direct')
    except:
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
            newdate = np.append(newdate,t)
        gapnans = np.nan*np.ones((nstep,np.shape(x)[1]))
        newx = np.vstack((newx,gapnans))
        if ii!=gapi[-1]:
            i1 = ii+1
            i2 = gapi[cnt+1]
            newdate = np.append(newdate,date[i1:i2+1])
            newx = np.vstack((newx,x[i1:i2+1,:]))
        else:
            newdate = np.append(newdate,date[ii+1:])
            newx = np.vstack((newx,x[ii+1:,:]))
        cnt=cnt+1

    if flipdim:
        newx = newx.transpose()

    if xnd == 1:
        newx = newx.flatten() # reduce back to 1D array if necessary

    return (newx,newdate)

def depthavg(x,z,h,ssh=None,surface='mixed',bottom='zero'):
    '''
Compute depth average of each row in 2D array x, with corresponding depths z.

Designed to accomodate upward looking ADCP data, with moving sea surface and
blank bins with no data near surface. If no sea surface height is specified,
it is assumed to be at z=0 for all times.

x: variable to be depth-averaged, 2D array with shape N rows, M columns
z: measurement depths (z=0 is surface, z=-h is bottom), array of length M
h: bottom depth (positive value)
ssh: sea surface height (optional, set to zero if None or where value is undefined)
surface: boundary condition for surface
        'mixed' (default) or 'extrap'
bottom:  boundary condition for bottom
        'zero' (default),'mixed' or 'extrap'
    '''

    # ensure that inputs are arrays of floats
    x = np.array(x).astype('float')
    z = np.array(z).astype('float')

    if np.ndim(x) == 1:
        x = x[np.newaxis]

    ni,nj = np.shape(x)

    # If SSH not specified, create an array of zeros
    if ssh is None:
        ssh = np.zeros(ni)

    ssh = np.array(ssh)
    if np.ndim(ssh) == 0:
        ssh = ssh[np.newaxis]

    # ssh in 2D column array
    ssh2 = np.array([ssh]).T

    # depths in 2D array
    sorti = np.argsort(z)
    zs = z[sorti]
    zs2 = np.tile(zs,[ni,1])
    xs2 = x[:,sorti] # sort data in same manner as depths

    # water depth in 2D column array
    h2 = np.tile(h,[ni,1])

    # new 2D x and z arrays to work with, with bottom and surface included
    zmat = np.hstack([-h2,zs2,ssh2])
    nans2 = np.nan*np.ones([ni,1])
    xmat = np.hstack([nans2,xs2,nans2])

    # only do calculations for rows where finite data exist
    fini = np.isfinite(xmat)
    ii, = np.where(np.sum(fini,axis=1) > 0)

    # bottom calculation
    if bottom == 'zero':
        xmat[ii,0] = 0.
    elif bottom == 'mixed':
        xmat[:,0] = xmat[:,1]
    elif bottom == 'extrap':
        xmat[:,0] = (xmat[:,2]-xmat[:,1])*(zmat[:,0]-zmat[:,1]) \
                    /(zmat[:,2]-zmat[:,1]) \
                    + xmat[:,1]
    else:
        raise ValueError('depthavg: bottom condition not understood (should be \'mixed\', \'extrap\' or \'zero\')')

    # find where depths are higher than sea surface or where there is no data,
    # mask with NaN Values
    xmatz = np.copy(xmat)
    xmatz[ii,-1] = 0.
    msk = (zmat > ssh2) | np.isnan(xmatz)
    zmatnan = np.copy(zmat)
    if np.any(msk):
        zmatnan[msk] = np.nan

    # sort each row of arrays by depth
    sj = np.argsort(zmatnan)
    si = np.arange(np.shape(zmat)[0])[:,np.newaxis]
    zmats = zmatnan[si,sj]
    xmats = xmat[si,sj]

    # column index of surface in each row where data exists
    jj = (np.sum(np.isfinite(zmats),axis=1)-1)[ii]

    # calculate surface value
    if surface == 'mixed':
        xmats[ii,jj] = xmats[ii,jj-1]
    elif surface == 'extrap':
        xmats[ii,jj] = (xmats[ii,jj-1]-xmats[ii,jj-2])*(zmats[ii,jj]-zmats[ii,jj-2]) \
                     /(zmats[ii,jj-1]-zmats[ii,jj-2]) \
                     + xmats[ii,jj-2]
    else:
        raise ValueError('depthavg: surface condition not understood (should be \'mixed\' or \'extrap\')')

    # integrate vertically using trapezoidal rule
    xm = 0.5*(xmats[:,:-1]+xmats[:,1:])
    dz = np.diff(zmats,axis=1)
    xint = np.nansum(xm*dz,axis=1)

    # divide by instantaneous water depth to compute depth average
    xda = np.nan*ssh
    xda[ii] = xint[ii]/(ssh[ii] + h)

    return xda

def surface_transport(u,z,ssh=None,surface='mixed'):
    '''
Compute surface transport Us from velocity component u, with corresponding depths z.

Integration is performed from the sea surface to the first zero crossing in the u profile.

If u has units of m/s, Us will have units of m^2/s.

Designed to accomodate upward looking ADCP data, with moving sea surface and
blank bins with no data near surface. If no sea surface height is specified,
it is assumed to be at z=0 for all times.

INPUT:
x: variable to be depth-averaged, 2D array with shape N rows, M columns
z: measurement depths (z=0 is surface, negative below surface), array of length M
ssh: sea surface height (optional, set to zero if None or where value is undefined)
surface: boundary condition for surface
        'mixed' (default) or 'extrap'

OUTPUT:
Us - transport integrated from surface to first zero crossing of velocity profiled
zs - depth of first zero crossing
    '''

    # ensure that inputs are arrays of floats
    x = np.array(u).astype('float')
    z = np.array(z).astype('float')

    if np.ndim(x) == 1:
        x = x[np.newaxis]

    ni,nj = np.shape(x)

    # If SSH not specified, create an array of zeros
    if ssh is None:
        ssh = np.zeros(ni)

    ssh = np.array(ssh)
    if np.ndim(ssh) == 0:
        ssh = ssh[np.newaxis]

    # ssh in 2D column array
    ssh2 = np.array([ssh]).T

    # depths in 2D array
    sorti = np.argsort(z)
    zsort = z[sorti]
    zs2 = np.tile(zsort,[ni,1])
    xs2 = x[:,sorti] # sort data in same manner as depths

    # make sure that there is at least one zero crossing
    # just below the deepest level
    eps = np.finfo(float).eps
    zend = zs2[:,0] - np.sqrt(eps)
    xend = -np.sign(xs2[:,0])*np.sqrt(eps)

    zend = zend[:,np.newaxis]
    xend = xend[:,np.newaxis]

    # new 2D x and z arrays to work with, with bottom and surface included
    zmat = np.hstack([zend,zs2,ssh2])
    nans2 = np.nan*np.ones([ni,1])
    xmat = np.hstack([xend,xs2,nans2])

    # only do calculations for rows where finite data exist
    fini = np.isfinite(xmat)
    ii, = np.where(np.sum(fini,axis=1) > 0)

    # find where depths are higher than sea surface or where there is no data,
    # mask with NaN Values
    xmatz = np.copy(xmat)
    xmatz[ii,-1] = 0.
    msk = (zmat > ssh2) | np.isnan(xmatz)
    zmatnan = np.copy(zmat)
    if np.any(msk):
        zmatnan[msk] = np.nan

    # sort each row of arrays by depth
    sj = np.argsort(zmatnan)
    si = np.arange(np.shape(zmat)[0])[:,np.newaxis]
    zmats = zmatnan[si,sj]
    xmats = xmat[si,sj]

    # column index of surface in each row where data exists
    jj = (np.sum(np.isfinite(zmats),axis=1)-1)[ii]

    # calculate surface value
    if surface == 'mixed':
        xmats[ii,jj] = xmats[ii,jj-1]
    elif surface == 'extrap':
        xmats[ii,jj] = (xmats[ii,jj-1]-xmats[ii,jj-2])*(zmats[ii,jj]-zmats[ii,jj-2]) \
                     /(zmats[ii,jj-1]-zmats[ii,jj-2]) \
                     + xmats[ii,jj-2]
    else:
        raise ValueError('surface_transport: surface condition not understood (should be \'mixed\' or \'extrap\')')

    # Flip 2D array so that surface is first column
    xmatr = np.fliplr(xmats)
    zmatr = np.fliplr(zmats)

    # Find index just above first zero crossing
    zj = np.argmax(np.sign(xmatr[:,:-1]*xmatr[:,1:]) == -1,axis=1)

    # Compute depth of first zero crossing using linear interpolation
    i = range(ni)
    zs = zmatr[i,zj] - xmatr[i,zj]*(zmatr[i,zj]-zmatr[i,zj+1])/(xmatr[i,zj]-xmatr[i,zj+1])

    # create rectangular array, where columns are integers starting from zero
    jmat = np.arange(np.shape(xmatr)[1])*np.ones(np.shape(xmatr)[0])[:,np.newaxis]

    # use to create mask of depths between surface and zero crossing
    mask = np.less_equal(jmat,zj[:,np.newaxis])

    # Replace data below first zero crossing with zeros
    xmatr = mask*xmatr

    # Use depth of first zero crossing for zero velocity closest to surface
    zmatr[i,zj+1] = zs

    # integrate vertically using trapezoidal rule
    xm = 0.5*(xmatr[:,:-1]+xmatr[:,1:])
    dz = -np.diff(zmatr,axis=1) # negative because depths are decreasing
    Us = np.nansum(xm*dz,axis=1)

    badi, = np.where(np.isnan(np.nanmax(xm,axis=1)))
    Us[badi] = np.nan
    zs[badi] = np.nan

    return Us,zs

def surface_flux(u,z,tr,ssh=None,surface='mixed'):
    '''
Compute integrated surface-layer flux of a tracer tr, due to velocity component u with corresponding depths z.

Integration is performed from the sea surface to the first zero crossing in the u profile.

If u has units of m/s, and tr has units of C, result will have units of (C m^2)/s

Designed to accomodate upward looking ADCP data, with moving sea surface and
blank bins with no data near surface. If no sea surface height is specified,
it is assumed to be at z=0 for all times.

INPUT:
u: velocity, 2D array with shape N rows, M columns
z: measurement depths (z=0 is surface, negative below surface), array of length M
tr: tracer values (N rows, M columns corresponding to depths z)
ssh: sea surface height (optional, set to zero if None or where value is undefined)
surface: boundary condition for surface
        'mixed' (default) or 'extrap'

OUTPUT:
trflux - tracer flux integrated from surface to first zero crossing of velocity profile
tr0 - tracer value at first zero crossing of velocity profile
    '''

    # ensure that inputs are arrays of floats
    x = np.array(u).astype('float')
    z = np.array(z).astype('float')
    tr = np.array(tr).astype('float')

    if np.ndim(x) == 1:
        x = x[np.newaxis]
        tr = tr[np.newaxis]

    ni,nj = np.shape(x)

    # If SSH not specified, create an array of zeros
    if ssh is None:
        ssh = np.zeros(ni)

    ssh = np.array(ssh)
    if np.ndim(ssh) == 0:
        ssh = ssh[np.newaxis]

    # ssh in 2D column array
    ssh2 = np.array([ssh]).T

    # depths in 2D array
    sorti = np.argsort(z)
    zsort = z[sorti] # sort depths
    zs2 = np.tile(zsort,[ni,1])
    xs2 = x[:,sorti] # sort data in same manner as depths
    trs2 = tr[:,sorti] # sort data in same manner as depths

    # make sure that there is at least one zero crossing
    # just below the deepest level
    eps = np.finfo(float).eps
    zend = zs2[:,0] - np.sqrt(eps)
    xend = -np.sign(xs2[:,0])*np.sqrt(eps)
    trend = trs2[:,0]

    zend = zend[:,np.newaxis]
    xend = xend[:,np.newaxis]
    trend = trend[:,np.newaxis]

    # new 2D x and z arrays to work with, with bottom and surface included
    zmat = np.hstack([zend,zs2,ssh2])
    nans2 = np.nan*np.ones([ni,1])
    xmat = np.hstack([xend,xs2,nans2])
    trmat = np.hstack([trend,trs2,nans2])

    # only do calculations for rows where finite data exist
    fini = np.isfinite(xmat)
    ii, = np.where(np.sum(fini,axis=1) > 0)

    # find where depths are higher than sea surface or where there is no data,
    # mask with NaN Values
    xmatz = np.copy(xmat)
    xmatz[ii,-1] = 0.
    msk = (zmat > ssh2) | np.isnan(xmatz)
    zmatnan = np.copy(zmat)
    if np.any(msk):
        zmatnan[msk] = np.nan

    # sort each row of arrays by depth
    sj = np.argsort(zmatnan)
    si = np.arange(np.shape(zmat)[0])[:,np.newaxis]
    zmats = zmatnan[si,sj]
    xmats = xmat[si,sj]
    trmats = trmat[si,sj]

    # column index of surface in each row where data exists
    jj = (np.sum(np.isfinite(zmats),axis=1)-1)[ii]

    # calculate surface value
    if surface == 'mixed':
        xmats[ii,jj] = xmats[ii,jj-1]
        trmats[ii,jj] = trmats[ii,jj-1]
    elif surface == 'extrap':
        xmats[ii,jj] = (xmats[ii,jj-1]-xmats[ii,jj-2])*(zmats[ii,jj]-zmats[ii,jj-2]) \
                     /(zmats[ii,jj-1]-zmats[ii,jj-2]) \
                     + xmats[ii,jj-2]
        trmats[ii,jj] = (trmats[ii,jj-1]-trmats[ii,jj-2])*(zmats[ii,jj]-zmats[ii,jj-2]) \
                     /(zmats[ii,jj-1]-zmats[ii,jj-2]) \
                     + trmats[ii,jj-2]
    else:
        raise ValueError('surface_transport: surface condition not understood (should be \'mixed\' or \'extrap\')')

    # Flip 2D array so that surface is first column
    xmatr = np.fliplr(xmats)
    zmatr = np.fliplr(zmats)
    trmatr = np.fliplr(trmats)

    # Find index just above first zero crossing
    zj = np.argmax(np.sign(xmatr[:,:-1]*xmatr[:,1:]) == -1,axis=1)

    # Compute depth of first zero crossing using linear interpolation
    i = range(ni)
    zs = zmatr[i,zj] - xmatr[i,zj]*(zmatr[i,zj]-zmatr[i,zj+1])/(xmatr[i,zj]-xmatr[i,zj+1])

    # Compute tracer value at first zero crossing using linear interpolation
    trs = trmatr[i,zj] - xmatr[i,zj]*(trmatr[i,zj]-trmatr[i,zj+1])/(xmatr[i,zj]-xmatr[i,zj+1])

    # create rectangular array, where columns are integers starting from zero
    jmat = np.arange(np.shape(xmatr)[1])*np.ones(np.shape(xmatr)[0])[:,np.newaxis]

    # use to create mask of depths between surface and zero crossing
    mask = np.less_equal(jmat,zj[:,np.newaxis])

    # Replace data below first zero crossing with zeros
    xmatr = mask*xmatr
    trmatr = mask*trmatr

    # Use depth of first zero crossing for zero velocity closest to surface
    zmatr[i,zj+1] = zs
    trmatr[i,zj+1] = trs

    # integrate vertically using trapezoidal rule
    xm = 0.5*(xmatr[:,:-1]+xmatr[:,1:])
    trm = 0.5*(trmatr[:,:-1]+trmatr[:,1:])
    dz = -np.diff(zmatr,axis=1) # negative because depths are decreasing
    trflux = np.nansum(xm*trm*dz,axis=1)

    badi, = np.where(np.isnan(np.nanmax(xm,axis=1)))
    trflux[badi] = np.nan
    zs[badi] = np.nan

    return trflux, trs

def lombscargle(t,x,ofac=4,hifac=1,t0=None,return_onesided=True,return_zero=False, window='boxcar',scaling='classical'):
    '''
    Compute the discrete Fourier transform and periodogram for unevenly-spaced
    data using the Lomb-Scargle periodogram. Follows methods outlined in Scargle
    (1989).

    INPUTS

    t - array of numerical time values (length N)
    x - array of data values, may be complex (length N)

    RETURNS

    f - array of frequencies
    ftx - array of complex coefficients
          discrete Fourier transform of x
    px - periodogram of ftx, proportional to |ftx|**2
         with the default "classical" scaling used in Scargle (1989),
         px = (1/N)*|ftx|**2

    OPTIONAL PARAMETERS

    ofac - oversampling parameter
           ratio of number of frequencies used to number of samples in x
           (default 4)
    hifac - high frequency parameter
            ratio of highest frequency to pseudo-Nyquist frequency
            (default 1)
    t0 - time origin
         reference point for phase calculation
         (default - None, first value in t array is used)
    return_onesided - boolean for returning a one-sided spectrum
        If True, return a one-sided spectrum for real data.
        If False return a two-sided spectrum for real data.
        Note that for complex data, a two-sided spectrum is always returned.
        (default - True)
    return_zero - boolean for evaluating zero frequency
        If True, include zero frequency. If False, do not include zero frequency.
        Uses expressions for the limit as frequency approaches zero, following
        Scargle (1989).
        (default - False)
    window - String specifying desired window to use. See `scipy.signal.get_window` for
        a list of windows and required parameters.
    scaling - Selects between computing the classical periodogram used by Scargle
        ('classical') or the power spectral density ('density'). The classical
        periodogram  has units of x**2. The power  spectral density has units of x**2/f.
        The scaling determines how the periodogram px is calculated from the discrete
        Fourier transform ftx:
        'classical': px = (1/N)*|ftx|**2, where N is the number of samples
        'density': px = (deltat/N)*|ftx|**2, where deltat is the average time step

    REFERENCE

    Scargle, J.D. (1989) Studies in astronomical time series analysis III: Fourier transforms,
    autocorrelation functions, and cross-correlation functions of unevenly spaced data. The
    Astrophysical Journal, 343, 874-887
    '''

    i = 1j # square root of -1
    N = len(x) # number of samples

    wts = signal.get_window(window,N)
    wts = N*wts/np.sum(wts) # make sum of weights equal to N

    x = x*wts  # apply window

    intm = np.mean(np.diff(t))

    flo = ((intm)**-1)/(len(x)*ofac)  # lowest freq
    fhi = hifac*(2*intm)**-1          # highest freq

    f = np.arange(flo,fhi+flo,flo)

    if return_zero == True:
        f = np.append(0,f)

    # if complex, evaluate two-sided spectrum regardless of user choice
    if np.any(np.iscomplex(x)):
        return_onesided = False

    # two-sided spectrum
    if return_onesided == False:
        if return_zero == True:
            f = np.append(-f[1:][::-1],f)
        else:
            f = np.append(-f[::-1],f)

    # time origin (reference point for phase calculation)
    if t0 is None:
        t0 = t[0]

    # initialize DFT as array of complex numbers
    ftx = np.nan*np.ones(len(f)) + i*np.nan*np.ones(len(f))

    for k,fk in enumerate(f):
        wrun = 2*np.pi*fk # angular frequency

        if fk == 0:
            # use well-defined limit as frequency approaches zero
            tau = np.sum(t)/N
            ftx[k] = np.sum(x)/np.sqrt(N)

        else:
            Fo = ((N/2)**0.5)*np.exp(-i*wrun*t0)

            tau = np.arctan2(np.sum(np.sin(2*wrun*t)),np.sum(np.cos(2*wrun*t)))/(2*wrun)
            tprime = t - tau

            A = np.sum(np.cos(wrun*tprime)**2)**-0.5
            B = np.sum(np.sin(wrun*tprime)**2)**-0.5

            # Note apparent typo in Scargle (1989), which has a plus
            # sign (+) instead of a minus sign below. This only makes a
            # difference in the periodogram if the input values in x are complex.

            ftx[k] = Fo*np.sum(A*x*np.cos(wrun*tprime) - i*B*x*np.sin(wrun*tprime))

    if scaling == 'classical':
        px = np.abs(ftx)**2/N
    elif scaling == 'density':
        px = np.abs(ftx)**2/N*intm
    else:
        raise ValueError('Scaling argument not understood. Acceptable options are classical or density')

    return f, ftx, px

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

    # test case #1 - depth avg
    u1= np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
               np.nan,         np.nan,         np.nan, -0.0018506 ,  0.00057345,
       -0.00027954,  0.00304925,  0.0056888 ,  0.01057738, -0.00096978,
        0.00614675,  0.00302453, -0.00028928,  0.00077288, -0.00768713,
       -0.01823976,  0.00612571, -0.00397687, -0.00580832, -0.00833382,
        0.0017868 , -0.00530538, -0.01031236])
    u2 = np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
               np.nan,         np.nan,         np.nan,  0.00150274,  0.00745662,
        0.00235997, -0.00074239, -0.00298656,  0.00302234,  0.00318832,
       -0.00169822, -0.00439177, -0.00226204, -0.00400032, -0.01001337,
       -0.00913997, -0.00681736, -0.01331132, -0.00100251, -0.01532928,
       -0.01763108, -0.01194093, -0.01909814])
    u = np.vstack([u1,u2])
    z = np.array([ 1.49,  1.24,  0.99,  0.74,  0.49,  0.24, -0.01, -0.26, -0.51,
       -0.76, -1.01, -1.26, -1.51, -1.76, -2.01, -2.26, -2.51, -2.76,
       -3.01, -3.26, -3.51, -3.76, -4.01, -4.26, -4.51, -4.76, -5.01, -5.26])
    ubar = depthavg(u,z,7)

    ub = u[:,::-1]
    zb = z[::-1]
    ubarb = depthavg(ub,zb,7)

    test = (ubar == ubarb)
    if test.all():
        print('depthavg test: passed')
    else:
        raise ValueError('depthavg test: failed')

    # surface transport test case #1
    u1= np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
               np.nan,         np.nan,         np.nan, -0.0018506 ,  0.00057345,
       -0.00027954,  0.00304925,  0.0056888 ,  0.01057738, -0.00096978,
        0.00614675,  0.00302453, -0.00028928,  0.00077288, -0.00768713,
       -0.01823976,  0.00612571, -0.00397687, -0.00580832, -0.00833382,
        0.0017868 , -0.00530538, -0.01031236])
    u2 = np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
               np.nan,         np.nan,         np.nan,  0.00150274,  0.00745662,
        0.00235997, -0.00074239, -0.00298656,  0.00302234,  0.00318832,
       -0.00169822, -0.00439177, -0.00226204, -0.00400032, -0.01001337,
       -0.00913997, -0.00681736, -0.01331132, -0.00100251, -0.01532928,
       -0.01763108, -0.01194093, -0.01909814])
    u = np.vstack([u1,u2])
    z = np.array([ 1.49,  1.24,  0.99,  0.74,  0.49,  0.24, -0.01, -0.26, -0.51,
       -0.76, -1.01, -1.26, -1.51, -1.76, -2.01, -2.26, -2.51, -2.76,
       -3.01, -3.26, -3.51, -3.76, -4.01, -4.26, -4.51, -4.76, -5.01, -5.26])
    Us,zs = surface_transport(u,z)

    test = np.array([np.isfinite(Us).all(),np.isfinite(zs).all()])
    if test.all():
        print('surface transport/flux test #1: passed')
    else:
        raise ValueError('surface transport/flux test #1: failed')

    # surface transport/flux test case #2
    u = np.array([[-1, -1, -1,  1,  1,  1],
                  [ 1,  2, -1,  4, -1, -2]])

    z = np.array([-6, -5, -4, -3, -2, -1])

    tr = np.array([[1, 1, -1,  -1,  -1,  -1],
                  [ 1,  1, 1,  1, 1, 1]])

    Us2,zs2 = surface_transport(u,z)
    trflux2,trs2 = surface_flux(u,z,tr)

    test = np.array([np.isfinite(Us2).all(),np.isfinite(zs2).all(),
                     np.isfinite(trflux2).all(),np.isfinite(trs2).all(),
                     zs2[0]>-4,zs2[0]<-3,zs2[1]>-3,zs2[1]<-2,
                     trflux2[0]==-Us2[0],trflux2[1]==Us2[1]])
    if test.all():
        print('surface transport/flux test #2: passed')
    else:
        raise ValueError('surface transport/flux test #2: failed')

    # surface transport/flux test case #3 (no zero crossing)
    z3 = np.array([-2. , -2.5, -3. , -3.5, -4. , -4.5, -5. , -5.5, -6. , -6.5, -7. ,
       -7.5, -8. , -8.5, -9. , -9.5])
    u3 = np.array(
        [ 4.21860615e-04,  2.31305385e-03,  1.31356857e-03,
         1.94045833e-03,  2.89933335e-04,  2.86449031e-03,
         1.97513672e-03,  6.18211638e-03,  3.36873352e-03,
         5.39932390e-03,  5.50479300e-03,  4.40279194e-03,
         5.36942569e-03,  6.60713067e-03,  5.74801833e-03,
         3.50485563e-03])
    tr3 = np.array(
       [-0.08513565, -0.08006332, -0.0802531 , -0.07925762, -0.07745643,
        -0.07925762, -0.07736153, -0.07698196, -0.07726664, -0.07608137,
        -0.07608137, -0.07418528, -0.0730949 , -0.07499098, -0.07608137,
        -0.07328469])

    Us3,zs3 = surface_transport(u3,z3)
    trflux3,trs3 = surface_flux(u3,z3,tr3)

    test = np.array([np.isfinite(Us3),np.isfinite(zs3),np.isfinite(trflux3),np.isfinite(trs3)])
    if test.all():
        print('surface transport/flux test #3: passed')
    else:
        raise ValueError('surface transport/flux test #3: failed')

    # test lomb-scargle function
    # use example from Trauth - MATLAB Recipes for Earth Sciences (3rd Ed)
    np.random.seed(0)
    t = np.arange(0,1000,3)
    t = t + np.random.randn(len(t))
    t = np.sort(t)

    x1 = (0.5*np.sin(2*np.pi*t/100) +
          1.0*np.sin(2*np.pi*t/40) +
          0.5*np.sin(2*np.pi*t/20))

    x2 = (0.5*np.sin(2*np.pi*t/100) +
          0.5*np.sin(2*np.pi*t/20))

    x = x1
    x[149:] = x2[149:]
    x = x + 0.5*np.random.randn(len(x))

    f,ftx,px_dft = lombscargle(t,x)
    px_scipy = signal.lombscargle(t, x, f*2*np.pi)

    test = np.isclose(px_dft,px_scipy)
    if test.all():
        print('Lomb-Scargle test: passed')
    else:
        raise ValueError('Lomb-Scargle test: failed')
