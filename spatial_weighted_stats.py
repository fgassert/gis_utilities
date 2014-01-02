"""
spatial_weighted_stats.py

Generates spatially weighted summary statistics on input polygons over each output polygon.

Requires ArcGIS 10.1+, numpy, matplotlib.mlab

by Francis Gassert

License: The MIT License (MIT)

Copyright (c) 2013 Francis Gassert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import arcpy as ap
import string
import numpy as np
import os, sys
import matplotlib.mlab as mlab
#import time, datetime

# Overwrite default read/write encoding hack
reload(sys)
sys.setdefaultencoding( "latin-1" )

# Temporary/intermediate geoprocessing files
TMPTABLE = "in_memory/tmptable"
TMPINTERSECT = "in_memory/tmpintersect"
TMPPRJ = "tmpprj"
TMPSELECT = "tmpselect"
TMPZONES = "in_memory/tmpzones"
TMPDOWNSAMPLE = "tmpds_"

VERBOSE = True

def _string_as_list(s):
    "If s is a string, returns a list [s], or a list of the comma separated substrings in s"
    if type (s) == str:
        s = [string.strip(x) for x in string.split(s,",")]
    assert type (s) in (list, tuple)
    return s

def _resample_rasters(rasters, prefix, resolution=None):
    """
    Takes a raster or list of rasters and resamples them to a given resolution (decimal degrees)
    Appends a prefix on the raster names
    """
    rasters = _string_as_list(rasters)
    out_rasters = []
    
    # resample raster
    if resolution is not None:
        for i in range(len(weight_rasters)):
            r = weight_rasters[i]
            d = ap.Describe(r)
            print "resampling from (%s,%s) to %s dd" % (d.meanCellHeight,
                                                        d.meanCellWidth,
                                                        downsample)
            rdir,rname = os.path.split(r)
            r_ds = os.path.join(rdir,"%s%s") % (TMPDOWNSAMPLE,rname)
            ap.Resample_management(r, r_ds, downsample)
            out_rasters.append(r_ds)
    else:
        out_rasters = rasters[:]
            
    return out_rasters

def dprint(s):
    if VERBOSE:
        print s

def spatial_weighted_stats(working_directory, scratch_gdb, input_polys, input_fields,
                        weight_rasters, out_file, output_polys = None,
                        stats={"mean":1}, downsample=None, null_value=-32767):
    """
    Generates spatially weighted summary statistics on input polygons over each output polygon.
    
    Parameters:
        working_directory:      set working directory
        scratch_gdb:            set scratch geodatabase for intermediate processing
        input_polys:            polygon features or list of features containing the attributes to be summarized
        input_fields:           names of fields containing values to be aggregated
        weight_rasters:         raster or list of rasters to be used as weights
        out_file:               prefix to apply to output csv
        output_polys:           [optional] polygons to calculate summary stats for
                                leave as <None> to calculate global statistics
        stats:                  dictionary defining the statistics to be calculated (see below)
        downsample:             [optional] resolution to resample rasters to before computation (decimal degrees)
        null_value:             integer value to represent no data
        
    Returns:
        None
    
    Saves:
        Summary statistics in <out_file>_<output_polys>.csv
        One file for each file in output_polys
        
    Stats options:
        Determines the statistics to be calculated
        Default: {"mean":1}
        
        {
            "mean":boolean      weighted mean value
            "v":boolean         weighted variance (population)
            "sd":boolean        weighted standard deviation 
            "mode":boolean      mode
            "min":boolean       lowest value
            "max":boolean       highest value
            "hist":int or list of values in ascending order
                                Calculate a histogram (sum of weights within each bin) with <int> equal width bins
                                or with bins bounded by a list of values (bins include the lower but not the upper bound)
            "sum":boolean       Sum of weights
            "nan":boolean       Sum of weights where value is null
            "p":list            List of weighted percentiles (0-1) to calculate, median is 0.5 percentile
        }
    
    """
    
    os.chdir(working_directory)
    dprint( "initializing" )
    
    if not ap.Exists(scratch_gdb):
        ap.CreateFileGDB_management(".",scratch_gdb)
    
    ap.env.overwriteOutput = True
    ap.env.workspace = scratch_gdb
    ap.env.outputCoordinateSystem = ap.SpatialReference("WGS 1984")
    ap.CheckOutExtension("Spatial")
    
    weight_rasters = _resample_rasters(weight_rasters,TMPDOWNSAMPLE,downsample)
    input_fields = _string_as_list(input_fields)
    
    if output_polys is not None:
        output_polys = _string_as_list(output_polys)
    else:
        output_polys = [None]
        
    for poly in output_polys:
        arr = stats_over_polys(weight_rasters, input_polys, input_fields, poly, stats, null_value)
        try:
            csv = "%s_%s.csv" % (out_file,repr(os.path.basename(poly)))
            mlab.rec2csv(arr,csv)
            dprint( "saved %s" % csv)
        except:
            csv = "%s_1.csv" % (out_file)
            mlab.rec2csv(arr,csv)
            dprint( "saved %s" % csv)
    
    dprint( "cleaning" )
    try:
        ap.Delete_management(TMPTABLE)
        ap.Delete_management(TMPINTERSECT)
        ap.Delete_management(TMPPRJ)
        ap.Delete_management(TMPZONES)
        if downsample is not None:
            for r in weight_rasters:
                ap.Delete_management(r)
    except Exception as e:
        dprint( "ERROR: %s" % e)
    dprint( "complete" )
    

def stats_over_polys(weight_rasters, input_polys, input_fields, output_poly=None,
                     stats={"mean":1}, null_value=-32767, dice=True):
    """generate summary statistics for given polygons using a raster weight"""
    # parse stat dict  
    stat_labels = get_stat_labels(stats)
    
    # initialize column labels
    labels = ["OID"]
    labels.extend(["%s%s_%s" % (os.path.basename(r),f,s) for r in weight_rasters for f in input_fields for s in stat_labels])
    
    if output_poly is None:
        # summarize for globe
        arr = np.zeros((len(labels),1))
        arr[1:,0] = _spatial_stats(weight_rasters, input_polys, input_fields, stats, null_value)
        
    else:
        # initialize output array
        oid_name = ap.Describe(output_poly).OIDFieldName
        ids = ap.da.FeatureClassToNumPyArray(output_poly,oid_name)[oid_name]
        obs = len(ids)
        arr = np.empty((len(labels),obs))
    
        if dice:
            ap.MakeFeatureLayer_management(output_poly,TMPSELECT)
            # iterate features and summarize
            for i in range(obs):
                dprint( "feature %s of %s" % (i+1, obs))
                ap.SelectLayerByAttribute_management(TMPSELECT, None,'"%s" = %s' % (oid_name, ids[i]))
                
                dprint( "intersecting")
                ap.Intersect_analysis([input_polys, TMPSELECT], TMPINTERSECT, "NO_FID")
                arr[0,i] = ids[i]
                try:
                    arr[1:,i] = _spatial_stats(weight_rasters, TMPINTERSECT, input_fields, stats, null_value)
                except Exception as e:
                    dprint( "ERROR: %s" % e)
        else:
            dprint( "intersecting")
            ap.Intersect_analysis([input_polys, output_poly], TMPINTERSECT, "NO_FID")
            ap.MakeFeatureLayer_management(TMPINTERSECT,TMPSELECT)
            # iterate features and summarize
            for i in range(obs):
                dprint( "feature %s of %s" % (i+1, obs))
                try:
                    ap.SelectLayerByAttribute_management(TMPSELECT, None,'"%s_%s" = %s' % (output_poly, oid_name, ids[i]))
                except:
                    ap.SelectLayerByAttribute_management(TMPSELECT, None,'"%s" = %s' % (oid_name, ids[i]))
                
                arr[0,i] = ids[i]
                try:
                    arr[1:,i] = _spatial_stats(weight_rasters, TMPSELECT, input_fields, stats, null_value)
                except Exception as e:
                    dprint( "ERROR: %s" % e)
        
    #return recarray
    arr = np.rec.fromarrays(arr, names=labels)
    return arr


def _spatial_stats(weight_rasters, input_polys, input_fields, stats, null_value):
    """generate summary statistics for a masked region using a raster weight"""
    # initialize output vector to no data
    x = get_stat_len(stats)
    arr = np.ones(len(input_fields)*x*len(weight_rasters))*null_value
    
    # ensure some overlap
    if ap.GetCount_management(input_polys).getOutput(0) == '0':
        dprint( "warning: no overlap")
        return arr
    
    oid_name = ap.Describe(input_polys).OIDFieldName
    export_fields=input_fields[:]
    export_fields.append(oid_name)
    dprint( "getting field values")
    polyarr = ap.da.FeatureClassToNumPyArray(input_polys,export_fields,null_value=null_value)
    
    dprint( "generating zones")
    ap.PolygonToRaster_conversion(input_polys,oid_name,TMPZONES,"CELL_CENTER",cellsize=weight_rasters[0])
    
    # ensure raster not all null
    if ap.sa.Raster(TMPZONES).minimum is None:
        dprint( "warning: no data (no overlap or cell size too large)")
        return arr
    
    dprint( "generating weights")
    for j in range(len(weight_rasters)):
        ap.sa.ZonalStatisticsAsTable(TMPZONES, "value", weight_rasters[j], TMPTABLE, "DATA", "SUM")
        freps = ap.da.TableToNumPyArray(TMPTABLE,["value","SUM"])
        weights = freps["SUM"]
        
        # ensure some overlap
        if np.sum(weights)==0:
            dprint( "warning: no meaningful overlap" )
            continue
        
        dprint( "summarizing statistics" )
        sort_idx = np.argsort(polyarr[oid_name])
        value_idx = np.searchsorted(polyarr[oid_name][sort_idx],freps["value"])
        for i in range(len(input_fields)):
            values = polyarr[input_fields[i]][sort_idx][value_idx]
            if values.dtype.kind in "iuf":
                values = values.astype(np.double)
                values[values==null_value] = np.nan
            ras_idx = j*x*len(input_fields)
            retarr = _weighted_stats(values, weights, stats)
            arr[ras_idx+i*x:ras_idx+(i+1)*x] = retarr
                   
    return arr

def _weighted_stats(values, weights, stats):
    """generates stats in order"""
    
    f = np.isfinite(values)
    valuesf = values[f]
    weightsf = weights[f]
    nnan = ~np.isnan(values)
    valuesn = values[nnan]
    weightsn = weights[nnan]
    
    arr = np.ones(get_stat_len(stats))*null_value
    if np.sum(weightsf) == 0:
        return arr
    
    i = 0
    
    avg = np.average(valuesf,weights=weightsf)
    v = np.dot(weightsf, (valuesf-avg)**2)/weightsf.sum()
    if "mean" in stats and stats["mean"]:   #mean
        arr[i]=avg; i+=1
    if "v" in stats and stats["v"]:   #var
        arr[i]=v; i+=1
    if "sd" in stats and stats["sd"]:   #sd
        sd = np.sqrt(v)
        arr[i]=sd; i+=1
    if "mode" in stats and stats["mode"]:   #mode
        arr[i]=values[np.argmax(weights)]; i+=1
    if "min" in stats and stats["min"]:   #min
        arr[i]=np.min(valuesn); i+=1
    if "max" in stats and stats["max"]:   #max
        arr[i]=np.max(valuesn); i+=1
    if "hist" in stats and stats["hist"]:
        if type(stats["hist"]) in (list, tuple):
            x = len(stats["hist"]); y = x-1
        else:
            x = stats["hist"]; y = x
        arr[i:i+y],_ = np.histogram(valuesn, bins=stats["hist"], weights=weightsn)
        i += x
    if "sum" in stats and stats["sum"]:
        arr[i] = np.sum(weights); i+=1
    if "nan" in stats and stats["nan"]:
        arr[i] = np.sum(weights[~nnan]); i+=1
    if "p" in stats:
        x = len(stats["p"])
        if x:
            arr[i:i+x]=weighted_percentiles(valuesn, weightsn, stats["p"]);
            i+=x
            
    return arr

def get_stat_labels(stats):
    """parses and sorts stat labels"""
    assert type (stats) is dict
    labels = []
    
    for k in ("mean","v","sd","mode","min","max","hist","sum","nan","p"):
        if k in stats.iterkeys():
            if type (stats[k]) in (list,tuple):
                for i in stats[k]:
                    labels.append("%s%s" % (k,i))
            elif stats[k]>1:
                for i in range(stats[k]):
                    labels.append("%s%s" % (k,i))
            elif stats[k]:
                labels.append(k)
    return labels

def get_stat_len(stats):
    assert type (stats) is dict
    x = 0
    for k in ("mean","v","sd","mode","min","max","hist","sum","nan","p"):
        if k in stats:
            if type (stats[k]) in (list,tuple):
                x += len(stats[k])
            else:
                x += int(stats[k])
    return x
   
def weighted_percentiles(values, weights, ps):
    """gets the value at a percentile given a set of values and frequencies"""
    assert len(values) == len(weights)
    n = len(values)
    idx = np.argsort(values)
    values = values[idx]
    weight = weights[idx]
    partialsum = np.empty(n)
    partialsum[0] = weight[0]
    for i in range(1,n):
        partialsum[i] = partialsum[i-1] + weight[i]
    percentiles = (partialsum - weight/2) / np.sum(weight)
    
    if type(ps) not in (list, tuple):
        ps = (ps,)
        
    x = len(ps)
    ps = np.empty(x)
    for i in range(x):
        if ps[i] < percentiles[0]:
            ps[i] = values[0]
        elif ps[i] > percentiles[-1]:
            ps[i] = values[-1]
        else:
            # identify first value greater than or equal to q
            j = np.searchsorted(percentiles, ps[i])
            # return closer value above or below
            if percentiles[j] - ps[i] <= ps[i] - percentiles[j-1]:
                ps[i] = values[j]
            else:
                ps[i] = values[j-1]
    return ps

if __name__ == '__main__':
    
    
    # ArcGIS scratch workspace
    working_directory = "C:/Users/francis.gassert/Documents/ArcGIS/GISSync/"
    scratch_gdb = "scratch.gdb"
    
       # Original polygon data
    input_polys = "global_maps/aqueduct_global_dl_20130123.gdb/global_master_20130123"
    # Attribute fields to aggregate (must be numeric)        
    input_fields = ["BWS_s","WSV_s","SV_s","HFO_s","DRO_s","DEFAULT"]
    
    # Weights for aggregation
    weight_rasters = ["use/U_2010.gdb/Ua","use/U_2010.gdb/Ud","use/U_2010.gdb/Ui","use/U_2010.gdb/Ut"]
    
    # Polygons to aggregate data into
    output_polys = ["global_maps/global_inputs.gdb/ne_10m_countries_20120902"]
    
    # Csv to save statistics in
    out_file = "countries_01"
    
    # Resolution to which weight rasters should be resampled before aggregation in decimal degrees,
    #  <None> to keep original resolution
    downsample = None # 1.0/24 # in dd
    
    # No data value
    null_value = -32767
    
    # Statistics to calculate
    stats = {"mean":1,"v":0,"sd":1,"mode":0,"hist":0,"min":0,"max":0,"sum":0,"nan":1,"hist":[0,1,2,3,4,5],"p":[0.5]}

    spatial_weighted_stats(working_directory, scratch_gdb, input_polys, input_fields,
                        weight_rasters, out_file, output_polys,
                        stats, downsample, null_value)