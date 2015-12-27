import scipy.io as spio
import numpy as np

def loadmat(filename):
    '''
    This function is an alternative to directly calling scipy.io.loadmat
    and cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check_keys() to cure all entries
    which are still mat-objects.
    
    Useful for loading comlicated structures with multiple nested levels.
    
    Source: Stack Overflow 
    http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    Author: mergen
    http://stackoverflow.com/users/887597/mergen
    License: Creative Commons Attribution-ShareAlike 3.0
    http://creativecommons.org/licenses/by-sa/3.0/
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)
    
def _check_keys(dict):
    '''
    Checks if entries in dictionary are mat-objects. If yes,
    _todict is called to change them to nested dictionaries
    
    Source: Stack Overflow 
    http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    Author: mergen
    http://stackoverflow.com/users/887597/mergen
    License: Creative Commons Attribution-ShareAlike 3.0
    http://creativecommons.org/licenses/by-sa/3.0/    
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries.
    
    Source: Stack Overflow 
    http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    Author: mergen
    http://stackoverflow.com/users/887597/mergen
    License: Creative Commons Attribution-ShareAlike 3.0
    http://creativecommons.org/licenses/by-sa/3.0/    
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def print_mat_nested(d, indent=0, nkeys=0):
    """
    Print nested structures in a dictionary, d, created with loadmat()

    Source: PyHOGS
    http://pyhogs.github.io/reading-mat-files.html   
    Author: JP Rinehimer
    
    Inspired by: 
    Stack Overflow 
    http://stackoverflow.com/questions/3229419/pretty-printing-nested-dictionaries-in-python
    Author: sth
    http://stackoverflow.com/users/56338/sth
    
    Modified by T. Connolly for Python 3 compatibility, as suggested by "Drunken Master" 
    (http://stackoverflow.com/users/4592067/drunken-master) in response to sth's SO post above
    
    License: Creative Commons Attribution-ShareAlike 3.0
    http://creativecommons.org/licenses/by-sa/3.0/
    """
    
    # Subset dictionary to limit keys to print.  Only works on first level
    if nkeys>0:
        d = {k: d[k] for k in d.keys()[:nkeys]}  # Dictionary comprehension: limit to first nkeys keys.

    if isinstance(d, dict):
        for key, value in d.items():         # iteritems loops through key, value pairs
          print('\t' * indent + 'Key: ' + str(key))
          print_mat_nested(value, indent+1)

    if isinstance(d,np.ndarray) and d.dtype.names is not None:  # Note: and short-circuits by default
        for n in d.dtype.names:    # This means it's a struct, it's bit of a kludge test.
            print('\t' * indent + 'Field: ' + str(n))
            print_mat_nested(d[n], indent+1)