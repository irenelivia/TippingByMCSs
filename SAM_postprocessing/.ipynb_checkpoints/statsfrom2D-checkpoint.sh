#!/bin/bash -l

echo "Enter file name"
read filename
file=$filename'.nc'


cdo -P 8 daymean $file $filename'_daymean.nc'
cdo -P 8 fldvar $file $filename'_fldvar.nc'
cdo -P 8 fldmean $file $filename'_fldmean.nc'
#cdo -P 8 mermean $file $filename'_mermean.nc'
