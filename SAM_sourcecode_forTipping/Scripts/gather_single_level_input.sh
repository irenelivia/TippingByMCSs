#!/bin/bash
set -ex

BASE=`pwd`

#module unload cdo/1.9.8-magicsxx-gcc64
#module load cdo/1.9.2-magicsxx-intel18

FSTEM=step_5K_2K
FSTEM=step
TMIN=193
TMAX=288
TMIN=193
TMAX=288

FSTEM=$1
TMIN=$2
TMAX=$3
LEV0=$4
LEV1=$5

INPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/paste_${TMIN}-${TMAX} 
OUTPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/paste_${TMIN}-${TMAX} 

cd $INPFAD
NY=$( ls -l ${FSTEM}.????0000.$LEV0.nc | wc -l )
NY=$(( $NY - 1 ))

cd $BASE
for N in $(seq 0 $NY); do
for VAR in  t q 
do    
    sh gather_RCE_byRow_level.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV0 & 
done
done
wait

for N in $(seq 0 $NY); do
for VAR in  n l r
do    
    sh gather_RCE_byRow_level.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV0 & 
done
done
wait

for N in $(seq 0 $NY); do
for VAR in  u v rgrp
do    
    sh gather_RCE_byRow_level.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV0 & 
done
done
wait

for N in $(seq 0 $NY); do
for VAR in  w
do    
    sh gather_RCE_byRow_level.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV1 & 
done
done
wait

for N in $(seq 0 $NY); do
for VAR in  t q n l r u v rgrp
do
    sh gather_only_single_level_aug.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV0 &
done
wait
done



for N in $(seq 0 $NY); do
for VAR in w
do
    sh gather_only_single_level_aug.sh $VAR $N $FSTEM $INPFAD $OUTPFAD $LEV1 &
done
wait
done

for VAR in  r u v t q n l rgrp 
do
   sh gather_RCE_byCol_level.sh $VAR $FSTEM $INPFAD $OUTPFAD $LEV0 &
done


for VAR in  w
do
   sh gather_RCE_byCol_level.sh $VAR $FSTEM $INPFAD $OUTPFAD $LEV1 &
done

