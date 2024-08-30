#!/bin/bash
set -ex

FSTEM=$1
NX=$2
NY=$3
TMIN=$4
TMAX=$5
LEV0=$6

./select_timesteps_single_level_index.sh ${FSTEM} ${TMIN} ${TMAX} ${LEV0}
./gather_single_level_index.sh ${FSTEM} ${TMIN} ${TMAX} ${LEV0}
#./reduce_stat.sh ${FSTEM}
#python precip_level.py ${FSTEM} ${NX} ${NY} ${TMIN} ${TMAX} ${LEV0}

