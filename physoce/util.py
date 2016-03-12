# -*- coding: utf-8 -*-
"""
General purpose functions.
"""
from datetime import datetime

def list2datetime(datestr_list,fmt='%a %b %d %H:%M:%S %Y'):
    '''Convert a list of date strings to datetime format.
    
    INPUT:
    datestr_list: a list of strings that represent dates
    fmt: format of the date string, as would be input to strftime() or strptime()
    
    see https://docs.python.org/library/datetime.html#strftime-and-strptime-behavior
    
    OUTPUT:
    list of datetimes
    '''
    datetime_list = [datetime.strptime(datestr, fmt) for datestr in datestr_list]
    return datetime_list