import numpy as np
import csv
from datetime import datetime
from matplotlib import dates
from physoce import tseries
import pandas as pd

def read_txt_data(data_file,format='dict',fillgap=False):     
    """
    read_txt_data(data_file,format='dict'):   
    Read in MBARI LOBO data file obtained in text format from http://www.mbari.org/lobo/getLOBOdata.htm
    
    Inputs:
        data_file: path to text data
    
    Options:
        format: output format
            format = 'dict' (default): dictionary
            format = 'dataframe': pandas DataFrame
    
    Output: dictionary (default) or pandas DataFrame with variable names taken from column headers
    """
    
    try:
        f = open(data_file)
    except:
        print("Cannot find " + data_file)        
    
    # see where the text starts in order to count header lines
    # (blank lines may be present at beginning of file)
    n = 0
    for line in f:
    	if line.find("LOBOVIZ") == 0:
    		break
    	n = n + 1
    	
    csvraw = list(csv.reader(f))
    hdrrow = n+2
    ncols = np.shape(csvraw[hdrrow+1])[0]
    
    # load numerical data into array 
    # first column with date is the only string, so this is skipped
    # also skip days since 1900 since this is redundant
    firstcol = 2
    data = np.genfromtxt(data_file, 
                         delimiter=',', 
                         skip_header=hdrrow+1, 
                         usecols=(np.arange(firstcol,ncols)))
    f.close()
    
    # list of variable names in the header
    # length should be the numer of columns in 'data' array
    varnames = csvraw[hdrrow-1]
    varnames = varnames[firstcol:]
    
    # get date strings from first column and put into datetime format
    date = []
    for row in csvraw[hdrrow:]:
        date.append(datetime.strptime(row[0], "%m/%d/%Y %H:%M")) 
    
    # convert huge values (missing data) to nan
    np.seterr(invalid='ignore')
    data[data > 1e+300] = np.nan
    
    if fillgap == True:
        (data,date) = tseries.fillgapwithnan(data,date)
    
    # date number  
    dnum = dates.date2num(date)
        
    ### Create a dictionary to be returned by this function 
    LoboDict = {}
    
    ### Loop through variables, adding keys and data to dictionary
    ii = 0
    for varstr in varnames:
        ui = varstr.find('[')
        varkey = varstr.lower()[0:ui]
        LoboDict[varkey] = data[:,ii]
        ii = ii+1
        
    ### add date number variable    
    LoboDict["dnum"] = dnum        
    
    # Check for variable-specific missing values
    if 'waterdepth' in LoboDict:
        missing = np.isclose(LoboDict['waterdepth'],-9.999)
        LoboDict['waterdepth'][missing] = np.nan
    
    if format=='dataframe':
        LoboDF = pd.DataFrame(LoboDict,index=date)
        return LoboDF
    else:
        LoboDict["date"] = date        
        return LoboDict