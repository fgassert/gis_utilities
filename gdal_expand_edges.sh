#!/bin/bash

set -e

TMP=__tmp.tif
TMP2=__tmp2.tif
TMP3=__tmp3.tif

if [ -z $1 ] || [ -z $2 ]; then
    echo "Usage: gdal_expand_edges <src_file> <dst_file> [attribute] [processing_resolution] [expand_cells]"
    exit
fi

SRC=$1
DST=$2

FID=OBJECT_ID
if [ ! -z $3 ]; then
    FID=$3
fi
echo "Attribute: $FID"
RES=0.00833333333
if [ ! -z $4 ]; then
    RES=$4
fi
echo "Resolution: $RES"
EXT=0.08333333333
if [ ! -z $5 ]; then
    EXT=$5
fi
echo "Expand distance: $EXT"

W=$(calc 360/$RES)
H=$(calc 180/$RES)

gdal_grid -ot UInt16 -l `basename $SRC .shp` -zfield $FID -a nearest:radius1=$EXT:radius2=$EXT -txe -180 180 -tye -90 90 -outsize $W $H $SRC $TMP
gdal_rasterize -a $FID -te -180 90 180 -90 -ts $W $H -ot UInt16 $SRC $TMP2

gdal_calc.py -A $TMP -B $TMP2 --outfile=$TMP3 --calc="A+B*(A==0)"

rm -f $DST
gdal_polygonize.py -mask $TMP3 $TMP3 -f "ESRI Shapefile" $DST $FID $FID

#rm -f $TMP $TMP2 $TMP3
