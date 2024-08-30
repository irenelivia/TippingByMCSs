for f in *2Dcom
do echo "Converting $f file.."
~/Codes/SAM-NBI/UTIL/2Dcom2nc $f
done

echo "merging files"
cdo mergetime *.nc all_data.nc

rm *2Dcom_1.nc
