'''
Functions related to oceanographic instrumentation.
'''

import numpy as np
import pandas as pd
import xarray as xr

def rditext_to_dataset(data_path,time_span=None):
    '''
    Load RDI ADCP text data file into xarray dataset.
    
    INPUTS:
    data_path: path to series text file created by RDI software (WinADCP)
    time_span (optional): list of two datetime64 values defining time span of data to load
    
    RETURNS
    xarray dataset
    
    EXAMPLE:
    import numpy as np
    data_path = 'RDI_ADCP_data.txt'
    t1 = np.datetime64('2017-11-10T10:00Z')
    t2 = np.datetime64('2018-10-17T10:00Z')  
    ds = rditext_to_dataset(data_path,[t1,t2])
    '''
    
    ### Note: tested on two Workhorse 4-beam datasets (bottom mounted and vessel mounted) ###
    
    nhead = 16
    headerlines = list()

    with open(data_path) as f: 
        for i in range(nhead):
            line = f.readline()
            line = line.rstrip('\n')
            line = line.replace('\"','')
            headerlines.append(line)

    # create list of variable names from header line
    var_names = headerlines[12].split('\t')
    for i,var in enumerate(var_names):
        var_names[i] = var_names[i].strip()
    var_units = headerlines[13].split('\t')

    df = pd.read_csv(data_path,
                     skip_blank_lines=False,
                     header=None,
                     names=var_names,
                     skiprows=nhead,
                     delimiter='\t',
                     skipinitialspace=True)

    # create datetime variable
    datestr = '20'+df['YR'].map(str)+'-'+df['MO'].map(str)+'-'+df['DA'].map(str)
    timestr = df['HH'].map(str)+':'+df['MM'].map(str)+':'+df['SS'].map(str)+'.'+df['HH.1'].map(str)
    df['datetime'] = pd.to_datetime(datestr+' '+timestr,utc=True)

    hrow, = np.where(['Pings/Ens' in s for s in headerlines])
    PingsPerEns = int(headerlines[hrow.squeeze()].split('\t')[1])

    hrow, = np.where(['Time/Ping' in s for s in headerlines])
    TimePerPing = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['First Ensemble Date' in s for s in headerlines])
    FirstEnsDate = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['First Ensemble Time' in s for s in headerlines])
    FirstEnsTime = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['Ensemble Interval' in s for s in headerlines])
    EnsInterval = float(headerlines[hrow.squeeze()].split(' = ')[1])

    hrow, = np.where(['1st Bin Range' in s for s in headerlines])
    BlankDist = float(headerlines[hrow.squeeze()].split(' = ')[1])

    hrow, = np.where(['Bin Size' in s for s in headerlines])
    BinSize = float(headerlines[hrow.squeeze()].split('\t')[1])

    binnumbers = np.unique(headerlines[14].split('\t'))[1:]
    nbins = np.max([int(binnum) for binnum in binnumbers])

    BinDist = np.arange(BlankDist,BlankDist+nbins*BinSize,BinSize)

    ntime = len(df['datetime'])

    ds = xr.Dataset(coords={'time': df['datetime'].values})
    ds['Ens'] = ('time', df['Ens'])
    
    # parse variables after date columns
    var_names_sub = var_names[8:]
    var_units_sub = var_units[8:]
    
    for varname in pd.unique(var_names_sub):
        
        var_match = np.array(var_names_sub)==varname 
        mi, = np.where(var_match)
        units = var_units_sub[int(mi[0])]
        count = np.sum(var_match)
        
        # fix cases where columns are separated by spaces instead of tabs
        var_list = []
        if df[varname].dtype == 'O':
            df_obj = df[varname].str.split(expand = True)
            var_names_obj = var.split()
            for i,varobj in enumerate(var_names_obj):
                df[varobj] = df_obj.iloc[i].astype(float)
                var_list.append(varobj)
        else:
            var_list.append(varname)
        
        for var in var_list:
            if (count == 1) & (np.sum(np.isfinite(df[var])) > 0):
                vardata = df[var].values
                if units == 'mm/s':
                    vardata = vardata/1000
                    units = 'm/s'
                ds[var] = ('time', vardata)
                ds[var].attrs['units'] = units
            elif (count == nbins) & (np.sum(np.isfinite(df[var])) > 0):
                vardata = np.nan*np.ones([ntime,nbins])
                for n in range(nbins):
                    if n == 0:
                        col = var
                    else:
                        col = var+'.'+str(n)
                    vardata[:,n] = df[col].values
                if units == 'mm/s':
                    vardata = vardata/1000
                    units = 'm/s'
                ds[var] = (['time','bin'],vardata)
                ds[var].attrs['units'] = units
            
    # distance between bins and instrument
    ds['BinDist'] = ('bin',BinDist)
    ds['BinDist'].attrs['units'] = 'm'
    ds['BinDist'].attrs['description'] = 'distance between bins and instrument'
    
    if time_span is not None:
        t1 = time_span[0]
        t2 = time_span[1]
        ii, = np.where((ds['time'] >= t1) & (ds['time'] <= t2))
        ds = ds.isel(time=ii)
    
    # add metadata
    ds.attrs['PingsPerEns'] = PingsPerEns
    ds.attrs['TimePerPing'] = TimePerPing    
    ds.attrs['First Ensemble Date'] = FirstEnsDate 
    ds.attrs['First Ensemble Time'] = FirstEnsTime
    ds.attrs['Ensemble Interval'] = EnsInterval
    ds.attrs['1st Bin Range'] = BlankDist
    ds.attrs['Bin Size'] = BinSize
    ds.attrs['RDI binary file'] = headerlines[1]
    ds.attrs['Instrument'] = headerlines[2]
            
    return ds