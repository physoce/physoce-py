'''
Tools for working with xarray Datasets and DataArrays
'''

import numpy as np
import xarray as xr
import physoce.tseries as ts

def filt(da,dim='time',dt=None,T=33,filter_name='pl64'):
    '''
Low pass filter a DataArray along a temporal dimension

Inputs:
da  - the DataArray to be filtered
dim - name of the temporal dimension (default 'time')
dt  - time interval of DataArray, in units of hours (default, selected automatically), for data obtained at hourly intervals, dt=1
      note: this must be consistent throughout the dataset. Specifying a value may be useful if the interval is nearly consistent, but not exactly
filter_name - name of the low pass filter (functions in physoce.tseries), should be 'pl64', 'pl66' or 'lancz' (default 'pl64')

Returns:
- the DataArray, filtered along the temporal dimension
    '''

    # automatically find time step if not specified by user
    if dt == None:
        dt_array = da[dim].diff(dim=dim)/np.timedelta64(1,'h')
        dt_unique = np.unique(dt_array)

        if len(dt_unique) == 1:
            dt = float(dt_unique)
        else:
            raise ValueError('time step must be consistent, or specified manually as input')

    # get filter weights
    if filter_name == 'pl64':
        wts = ts.pl64(dt=dt,T=T,return_weights=True)
    elif filter_name == 'pl66':
        wts = ts.pl66(dt=dt,T=T,return_weights=True)
    elif filter_name == 'lancz':
        wts = ts.lancz(dt=dt,T=T,return_weights=True)
    else:
        raise NameError('filter_name not understood, must be pl64, pl66 or lancz')

    # filter along specified dimension
    # follows example at: http://xarray.pydata.org/en/stable/user-guide/computation.html#rolling-window-operations
    weight = xr.DataArray(wts,dims=['window'])
    daf = da.rolling({dim:len(weight)},center='True').construct(time='window').dot(weight)

    return daf
