"""
areagrid.py 

Generates a global grid in geograpic coordinates with the area of each gridcell as the value of the cell

Requires: ArcGIS 10.1+, Numpy

The MIT License (MIT)

Copyright (c) 2011 Francis Gassert

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

import numpy as np
import arcpy as ap

def areagrid(outraster, c = 0.083333001, r = 6371007.2, minx = -180, miny = -90, w = 360, h = 180):
    """
    Generates a global grid in geograpic coordinates
    with the area of each gridcell as the value of the cell
    
    Parameters:
        outraster = path to output raster file
        c = cellsize in decimal degrees
        r = radius of earth in desired units (e.g. 6371007.2m for sq. meters)
        
    Returns:
    """
    
    X = np.ones(round(w/c))
    Y = np.arange(round(h/c))
    
    radn = (Y*c + miny)*np.pi/180.0
    radc = c * np.pi/180.0
    
    print radn
    print radc
    
    A = (np.sin(radn+radc/2)-np.sin(radn-radc/2)) * radc * r**2
    
    print A
    print A[len(A)//2]
    
    M = np.outer(A,X)
    print M
    print M.shape


    r = ap.NumPyArrayToRaster(M,arcpy.Point(minx,miny),c,c)
    r.save(outraster)
    