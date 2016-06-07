import numpy as np
try:
    import shapefile
except ImportError:
    print('warning: shapefile package not installed')
try:
    import pyproj
except ImportError:
    print('warning: pyproj package not installed')

def shapefile_latlon(es_shapefile,thresh=500):
    """
    Return latitude and longitude points from a GIS shapefile downloaded from
    http://www.elkhornslough.org/gis/index.htm
    
    Inputs: 
    es_shapefile: file path/prefix (for exmaple, if the coastline data files,
                              (cz.dbf, cz.shp, cz.shx) are in a directory 
                              called CZ, this would be 'CZ/cz')
    thresh: defines a threshold for gaps between points (in m), gaps more than
            this distance apart separated by NaN values to make lines look
            better when they are plotted
            
    Output:
    A dictionary with keys 'lon' and 'lat'
    
    Required packages:
    pyproj - https://pypi.python.org/pypi/pyproj
    pyshp - https://pypi.python.org/pypi/pyshp
    """    
    
    """
    Tom Connolly, MLML    
    """
    
    sf = shapefile.Reader(es_shapefile)
    lons = np.array(np.nan)
    lats = np.array(np.nan)
    for shape in sf.shapes():
        points = np.asarray(shape.points)
        x = points[:,0]
        y = points[:,1]
        
        dist = (np.diff(x)**2+np.diff(y)**2)**0.5
        ii = np.where(dist>thresh)
    
        p = pyproj.Proj(proj="utm",zone=10,datum='WGS84')
        lon, lat = p(x,y,inverse=True)
    
        # if there are distances above threshold, loop through and insert nan values
        if np.shape(ii)[1]>0:
            for idx in ii[0]:
                lon = np.hstack((lon[0:idx],np.nan,lon[idx+1:]))
                lat = np.hstack((lat[0:idx],np.nan,lat[idx+1:]))
                
        lons = np.hstack((lons,np.nan,lon))[2:]
        lats = np.hstack((lats,np.nan,lat))[2:]
        
        
    # return dictionary with lon/lat
    lld = dict()
    lld['lon'] = lons[0:-1]
    lld['lat'] = lats[0:-1]
    return lld