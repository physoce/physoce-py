# -*- coding: utf-8 -*-
import glob
from datetime import datetime
import numpy as np

# TO DO: 
# Finish load_met function
# Inlude a lookup table for site position (lat/lon)
# Convert wind speed and direction to velocity vector
# Anemometer height

def string_convert(string):
    """
    Convert strings to numbers, convert to nan if not a number
    """
    try:
        num = float(string)
    except:
        num = np.nan
    return num

def load_met(data_dir,siteid,flaglist=[],frequency=15):
    """
    Read NERR meteorological data files from multiple years into one dictionary.
    Numerical data are in NumPy arrays with missing data as NaN values.
    All other data are in lists of strings.    
    Variable names and units are kept the same as in the original files. See NERR
    documentation for more information.
    
    Additional variables added:    
    datetime - Python datetime object
    
    flaglist specifies which flags should be used to identify bad data to be returned as NaN
    example: load_met(data_dir,siteid,flaglist=['<-5>','<-4>','<-3>','<-2>','<1>'])
    See NERR data readme.rtf for more information
    
    Frequency 
    """
    
    # WORK IN PROGRESS #    
    
    file_list = glob.glob(data_dir + siteid + 'met*.csv')
    cnt = 0 # initialize variable that counts total number of lines
    d = dict() # create new dictionary (to be returned)
    
    # make sure flaglist is a list 
    # (if it is a string, all data may be interpreted as bad)
    if isinstance(flaglist, list) != True:  
        flaglist = [flaglist]
    
    # get list of variable names from header of first file
    f = open(file_list[0],'r')
    header = f.readline()
    varnames = header.split(',')
    for ii,var in enumerate(varnames[0:-1]):
        varnames[ii] = var.strip('"') # strip quotes from each variable name 
        d[varnames[ii]] = [] #initialize dictionary with key and empty list
    # last column is junk
    varnames = varnames[0:-1]
    f.close()
    
    d['datetime'] = []    
    
    datecol = "DateTimeStamp"
    

def load_wq(data_dir,siteid,flaglist=[]):
    """
    Read NERR water quality data files from multiple years into one dictionary.
    Numerical data are in NumPy arrays with missing data as NaN values.
    All other data are in lists of strings.    
    Variable names and units are kept the same as in the original files. See NERR
    documentation for more information.
    
    Additional variables added:    
    datetime - Python datetime object
    
    flaglist specifies which flags should be used to identify bad data to be returned as NaN
    example: load_wq(data_dir,siteid,flaglist=['<-5>','<-4>','<-3>','<-2>','<1>'])
    See NERR data readme.rtf for more information
    """
    
    file_list = glob.glob(data_dir + siteid + 'wq*.csv')
    cnt = 0 # initialize variable that counts total number of lines
    d = dict() # create new dictionary (to be returned)
    
    # make sure flaglist is a list 
    # (if it is a string, all data may be interpreted as bad)
    if isinstance(flaglist, list) != True:  
        flaglist = [flaglist]
    
    # get list of variable names from header of first file
    f = open(file_list[0],'r')
    header = f.readline()
    varnames = header.split(',')
    for ii,var in enumerate(varnames[0:-1]):
        varnames[ii] = var.strip('"') # strip quotes from each variable name 
        d[varnames[ii]] = [] #initialize dictionary with key and empty list
    # last column is junk
    varnames = varnames[0:-1]
    f.close()
    
    d['datetime'] = []    
    
    datecol = "DateTimeStamp"
    # specifiy columns with numeric data that should be converted to numpy arrays
    floatcols=["Temp","SpCond","Sal","DO_Pct","DO_mgl","Depth","cDepth","Level","cLevel","pH","Turb","ChlFluor"]
    
    for data_file in file_list:       
        # check to see if data file exists
        try:
            f = open(data_file,'r')
        except:
            print("Cannot find " + data_file)   
        
        # since first line is header, read lines after header
        for line in f.readlines()[1:]:
            cols =  line.split(',')         
            # if line is empty, skip line
            # (last lines in file sometimes empty)
            if len(cols) > 1:        
                for ii,var in enumerate(varnames):
                    if var == datecol:
                        d['datetime'].append(datetime.strptime(cols[ii], '%m/%d/%Y %H:%M'))
                    
                    if var in floatcols:
                        value = string_convert(cols[ii])
                        if value == -99:
                            value = np.nan
                        d[var].append(value)
                    else:
                        d[var].append(cols[ii])
                    
                cnt = cnt+1      
                
        f.close()
    
    # convert numeric data to NumPy arrays
    for fvar in floatcols:
        d[fvar] = np.asarray(d[fvar])
        
    # apply flags specified in "flaglist" input variable
    for var in varnames:
        if 'F_' in var and var != 'F_Record':
            datavar = var.strip('F_')
            for flag in flaglist:
                flagi = [i for i, s in enumerate(d[var]) if flag in s]
                d[datavar][flagi] = np.nan
   
    return d
        

def load_nut(data_dir,siteid,flaglist=[]):
    """
    Read NERR nutrient data files from multiple years into one dictionary.
    Numerical data are in NumPy arrays with missing data as NaN values.
    All other data are in lists of strings.    
    Variable names and units are kept the same as in the original files. See NERR
    documentation for more information.
    
    Additional variables added:    
    datetime - Python datetime object
    PO4F_uM,NH4F_uM,NO2F_uM,NO3F_uM,NO23F_uM - nutrient concentrations in uM (micromolar)
    
    flaglist specifies which flags should be used to identify bad data to be returned as NaN
    example: load_nut(data_dir,siteid,flaglist=['<-5>','<-4>','<-3>','<-2>','<1>'])
    See NERR data readme.rtf for more information
    """
    
    file_list = glob.glob(data_dir + siteid + 'nut*.csv')
    
    cnt = 0 # initialize variable that counts total number of lines
    d = dict() # create new dictionary (to be returned)

    # make sure flaglist is a list 
    # (if it is a string, all data may be interpreted as bad)
    if isinstance(flaglist, list) != True:  
        flaglist = [flaglist]
    
    # Initialize dictionary with keys and empty values
    d['StationCode'] = []
    d['isSWMP'] =  []
    d['DateTimeStamp'] = []
    d['datetime'] = []
    d['Historical'] = []
    d["ProvisionalPlus"] = []
    d["CollMethd"] = []
    d["REP"] = []
    d["F_Record"] = []
    d["PO4F"] = []
    d["F_PO4F"] = []
    d["NH4F"] = []
    d["F_NH4F"] = []
    d["NO2F"] = []
    d["F_NO2F"] = []
    d["NO3F"] = []
    d["F_NO3F"] = []
    d["NO23F"] = []
    d["F_NO23F"] = []
    d["CHLA_N"] = []
    d["F_CHLA_N"] = []
    
    for data_file in file_list:       
        # check to see if data file exists
        try:
            f = open(data_file,'r')
        except:
            print("Cannot find " + data_file)   
        
        # since first line is header, read lines after header
        for line in f.readlines()[1:]:
            cols =  line.split(',')            
            # if line is empty, skip line
            # (last lines in file sometimes empty)
            if len(cols) > 1:            
                # columns to variables
                # strip quotes and whitespace from strings
                d['StationCode'].append(cols[0].strip('"').strip())
                d['isSWMP'].append(cols[1].strip('"'))
                d['DateTimeStamp'].append(cols[2].strip('"'))
                d['datetime'].append(datetime.strptime(cols[2].strip('"'), '%m/%d/%Y %H:%M'))
                d['Historical'].append(cols[3])
                d["ProvisionalPlus"].append(cols[4])
                d["CollMethd"].append(cols[5])
                d["REP"].append(cols[6])
                d["F_Record"].append(cols[7])
                d["PO4F"].append(string_convert(cols[8]))
                d["F_PO4F"].append(cols[9])
                d["NH4F"].append(string_convert(cols[10]))
                d["F_NH4F"].append(cols[11])
                d["NO2F"].append(string_convert(cols[12]))
                d["F_NO2F"].append(cols[13])
                d["NO3F"].append(string_convert(cols[14]))
                d["F_NO3F"].append(cols[15])
                d["NO23F"].append(string_convert(cols[16]))
                d["F_NO23F"].append(cols[17])
                d["CHLA_N"].append(string_convert(cols[18]))
                d["F_CHLA_N"].append(cols[19])
                cnt = cnt + 1
            
        f.close()
            
    # convert numerical data to numpy arrays
    d["PO4F"] = np.asarray(d["PO4F"])
    d["NH4F"] = np.asarray(d["NH4F"])
    d["NO2F"] = np.asarray(d["NO2F"])
    d["NO3F"] = np.asarray(d["NO3F"])
    d["NO23F"] = np.asarray(d["NO23F"])
    d["CHLA_N"] = np.asarray(d["CHLA_N"])
    
    # apply flags specified in "flaglist" input variable
    varnames = ['F_PO4F','F_NH4F','F_NO2F','F_NO3F','F_NO23F','F_CHLA_N']
    for var in varnames:
        datavar = var[2:]
        for flag in flaglist:
            flagi = [i for i, s in enumerate(d[var]) if flag in s]
            d[datavar][flagi] = np.nan    
    
    # convert from mg/L to to uM
    # From the NERR nutrient metadata:
    # "All parameter values at Moss Landing Marine Laboratories (MLML) are 
    # calculated and reported in µM.  For purposes of consistency in the NERR 
    # System, Elkhorn Slough NERR calculates the concentrations as mg/ l-1 
    # based on atomic weights of 14.01, 30.97 for N and P respectively.  
    # Therefore, Elkhorn Slough NERR staff multiplies the concentrations 
    # reported by MLML by 0.01401, 0.03097 to yield concentrations in mg/L 
    # as N and P respectively."

    d["PO4F_uM"] = d["PO4F"]/0.03097
    d["NH4F_uM"] = d["NH4F"]/0.01401
    d["NO2F_uM"] = d["NO2F"]/0.01401
    d["NO3F_uM"] = d["NO3F"]/0.01401
    d["NO23F_uM"] = d["NO23F"]/0.01401
    
    return d