"""
areagrid.py 

Generates a global grid in geograpic coordinates with the area of each gridcell as the value of the cell

Requires: Rasterio, Numpy, GDAL

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
import sys
import rasterio
from rasterio import Affine


def areagrid(outraster, c = 0.083333001, r = 6371007.2, minx = -180, miny = -90, w = 360, h = 180):
    """
    Generates a global grid in geograpic coordinates
    with the area of each gridcell as the value of the cell
    
    Parameters:
        outraster = path to output raster file
        c = cellsize in decimal degrees
        r = radius of earth in desired units (e.g. 6371007.2m for sq. meters)
        minx, miny = the south west coordinate of the raster in degrees
        w, h = the width and height of the raster in degrees
        
    Returns:
        None
    """
    c = float(c)
    r = float(r)

    # make a vector of ones [1,1,1, ... 1] of length equal to the number of cells from west to east.
    X = np.ones(round(w/c))
    # make a vector counting from 0 to the number of cells from south to north. e.g. [0,1,2,...,179] for 1 deg cells.
    Y = np.arange(round(h/c))
    
    # multiply all the numbers in the Y vector by the cell size, 
    # so it extends from 0 to <180 (if the cell size is different than 1 deg)
    # then add the southernmost coordinate (-90deg). This makes a vector of -90 to +90 degrees North
    degN = Y*c + miny
    # convert degrees vector to radians
    radN = degN*np.pi/180.0
    # convert the cell size to radians
    radc = c * np.pi/180.0
    
    # calculate the area of the cell
    # there's some implicit geometry that's been done here already'
    # but basically it averages the width of the top of a cell and the bottom of a cell
    # and mutiplies it by the height of the cell (which is constant no matter how far or south north you are)
    # then by the square of the radius
    # since the angles are in radians it works out area correctly
    # you end up with a vector of cell area from south to north
    A = (np.sin(radN+radc/2)-np.sin(radN-radc/2)) * radc * r**2
    
    # the outer product of any vector and a vector of ones just duplicates the first vector into columns in a matrix
    # basically we just copy the latitude vector across from east to west
    M = np.outer(A,X)
    
    # save the matrix as a raster
    with rasterio.open(outraster,'w',
                       'GTiff',
                       width=w,
                       height=h,
                       dtype=M.dtype,
                       crs={'init': 'EPSG:4326'},
                       transform=Affine.translation(-w*c/2, h*c/2) * Affine.scale(c, -c),
                       count=1) as dst:
        dst.write_band(1,M)

    
if __name__ == "__main__":
    areagrid(*sys.argv[1:])
