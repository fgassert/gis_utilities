"""
flowaccumulator.py

A simple basin to basin flow accumulation script

Requires numpy

Based on work by ISciences L.L.C.
http://isciences.com/

Copyright 2013 World Resources Institute

Licensed under the Creative Commons Attribution 4.0 International Public License (CC-BY 4.0)
http://creativecommons.org/licenses/by/4.0/
"""

import numpy as np

# maximum levels to accumulate before breaking
MAXLEVELS = 1000
VERBOSE = False

def dprint(s):
    if VERBOSE:
        print s

def accumulate(ids, d_ids, f0, f, *args):
    """
    Accumulates values over basins based on downstream relationships
    
    Parameters:
        ids:numpy.array     basin ids
        d_ids:numpy.array   id of basin which is immediately downstream of the given basin
        f0:function         function to compute values for basins without upstream basins
                            f0 must take at least one parameter:
                                i:int, the index of given basin within <ids> (this index in independent of the basin id,
                                   to retrieve the basin id, pass the <ids> array as an addition argument. ids[i] -> basin id
                            f0 must return a single numeric value
        f:function          function to compute values for basins with upstream basins
                            f must take at least three parameters:
                                i:int, the index of the given basin within <ids>
                                idx:numpy.array, boolean vector of indicating the basins immediately upstream of the given basin
                                values:numpy.array, computed values for all basins 
                                    values[idx] -> values of basins immediately upstream of the given basin
                            for simple flow accumulation, f should return ( sum(values[idx]) + f0(i, *args) )
                            f must return a single numeric value
        *args:arguments     [optional] additional arguments to pass through to f0 and f.
                            Both f0 and f must be able to accept the same additional arguments.
        
    Returns:
        numpy.array         accumulated basin values
        
    """
    return accumulate_vector(ids, d_ids, f0, f, 1, *args)
    
def accumulate_vector(ids, d_ids, f0, f, _len=1, *args):
    """
    Same as accumulate, but values is an m*n array where
        m = len(ids)
        n = _len
    Allows accumulation of vectors of values rather than single values
    f0 and f must return a 1d array-like of length _len
    """
    x = len(ids)
    dprint ("ids: %s" % x)
    
    computed = np.zeros(x,dtype=bool)
    if _len > 1:
        values = np.zeros((x,_len))
    else:
        values = np.zeros(x)
    
    # build upstream index array
    up_idx = np.empty((x,x),dtype=bool);
    for i in range(x):
        up_idx[i,:] = d_ids==ids[i]
    
    # no basins are upstream of the given basin
    no_upstream = ~np.any(up_idx, 1)
    
    # compute values for the basins with no upstream basins
    for i in np.arange(x)[no_upstream]:
        values[i] = f0(i, *args)
    computed[no_upstream] = 1
    
    level = 0
    while ~np.all(computed) and level<MAXLEVELS:
        
        dprint ("computed: %s/%s" % (np.sum(computed),x))
        
        # none of the uncomputed basins are upstream of the given basin
        up_computed = ~np.any(up_idx[:,~computed], 1)
        # and the given basin hasn't been computed yet
        to_compute = up_computed & ~computed
        
        # just compute the ones we need to
        for i in np.arange(x)[to_compute]:
            values[i] = f(i, up_idx[i,:], values, *args)
        computed[to_compute] = 1
        
        level +=1
    
    dprint ("longest path: %s" % level)
    
    return values


#def upstream_ids(ids, d_ids):
#    """
#    returns a list of list
#    where the inner lists contain the ids of all upstream basins
#    """
#    x = len(ids)
#    dprint ("ids: %s" % x)
#    
#    computed = np.zeros(x,dtype=bool)
#    up_ids = [[] for i in range(x)]
#    
#    up_idx = np.empty((x,x),dtype=bool);
#    for i in range(x):
#        up_idx[i,:] = d_ids==ids[i]
#    
#    no_upstream = ~np.any(up_idx, 1)
#    computed[no_upstream] = 1
#    
#    level = 0
#    while ~np.all(computed) and level<MAXLEVELS:
#        
#        dprint ("computed: %s/%s" % (np.sum(computed),x))
#        up_computed = ~np.any(up_idx[:,~computed], 1)
#        to_compute = up_computed & ~computed
#        
#        for i in np.arange(x)[to_compute]:
#            up_ids[i]=list(ids[up_idx[i,:]]).extend(up_ids[up_idx[i,:]])
#        computed[to_compute] = 1
#        
#        level +=1
#    
#    dprint ("longest path: %s" % level)
#    
#    return up_ids


def test():
    """Test script"""
    import matplotlib.mlab as mlab
    import time
    import gen_merge

    BASINCSV = r"C:\Users\francis.gassert\Documents\ArcGIS\GISSync\global_maps\basins_15006.csv"
    BASINID = "basinid"
    DWNBASIN = "dwnbasinid"
    OUTCSV = r"C:\Users\francis.gassert\Documents\ArcGIS\GISSync\global_maps\bt_test.csv"
    runoffcsv = r"C:\Users\francis.gassert\Documents\ArcGIS\GISSync\global_maps\global-GLDAS-2.0_Noah-3.3_M.020-20121211-filled-20130821-RO.csv"
    
    basin_arr = mlab.csv2rec(BASINCSV)
    ids = basin_arr[BASINID]
    d_ids = basin_arr[DWNBASIN]
    r_arr = mlab.csv2rec(runoffcsv)
    r = r_arr["2010"]
    assert np.all(r_arr[BASINID]==ids)
    
    def f0( i, r ):
        return r[i]
    def f( i, idx, values, *args ):
        return np.sum(values[idx]) + f0(i, *args)
    
    
    time.clock()
    #id_dict = dict(zip(ids, upstream_ids(ids, d_ids)))
    #r2 = gen_merge.arrange_vector_by_ids(r, ids, np.arange(len(ids)+1))
    #out1 = np.array([np.sum(r2[id_dict[i]])+r2[i] for i in ids])
    #t1 = time.clock()

    out2 = accumulate(ids, d_ids, f0, f, r)
    t2 = time.clock()
    
    btcsv = r"C:\Users\francis.gassert\Documents\ArcGIS\GISSync\global_maps\global-GLDAS-2.0_Noah-3.3_M.020-20121211-filled-20130821-Bt.csv"
    bt_arr = mlab.csv2rec(btcsv)
    bt = bt_arr["2010"]

    #print ("time1: %s" % t1)
    print ("time2: %s" % t2)

    #print ("error1: %s " % (np.sum(out1-bt)/np.sum(bt)) )
    print ("error2: %s " % (np.sum(out2-bt)/np.sum(bt)) )
    
    outrec2 = np.rec.fromarrays((ids,out2),names=(BASINID,"2010"))
    mlab.rec2csv(outrec2,OUTCSV)
    


if __name__ == '__main__':
    VERBOSE=True
    test()
    
