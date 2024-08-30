#!/bin/bash
set -ex



FSTEM=step_5K_2K
#FSTEM=step

PFAD=`pwd`



TMIN=289
TMAX=384
TMIN=385
TMAX=480

FSTEM=$1
TMIN=$2
TMAX=$3
LEV0=$4
LEV1=$5
#TMIN=193
#TMAX=288
#TMAX=2499

INPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/
OUTPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/paste_${TMIN}-${TMAX}/
mkdir -p $OUTPFAD

cd $INPFAD


#$(pwd)
NN=$( ls -l ${FSTEM}.????????.nc | wc -l )
NX=$( ls -l ${FSTEM}.0000????.nc | wc -l )
NY=$( ls -l ${FSTEM}.????0000.nc | wc -l )

NX=$(( $NX - 1 ))
NY=$(( $NY - 1 ))

for N in $(seq 0 $NY); do
#for N in $(seq 12 12); do

NSTRING=$(printf %04d $N)

for M in $(seq 0 $NX); do
#for M in $(seq 2 2); do

echo "selecting m= $M" 
MSTRING=$(printf %04d $M)


#ncea -F -d time,2500,2988 ${FSTEM}.${NSTRING}${MSTRING}.nc $OUTPFAD/$FSTEM.${NSTRING}${MSTRING}.nc &

cdo -P 4 seltimestep,${TMIN}/${TMAX} -sellevel,${LEV0} $FSTEM.${NSTRING}${MSTRING}.nc $OUTPFAD/$FSTEM.${NSTRING}${MSTRING}.${LEV0}.nc &
cdo -P 4 seltimestep,${TMIN}/${TMAX} -sellevel,${LEV1} $FSTEM.${NSTRING}${MSTRING}.nc $OUTPFAD/$FSTEM.${NSTRING}${MSTRING}.${LEV1}.nc &

done 
wait
done 
###### end of selecting time steps

