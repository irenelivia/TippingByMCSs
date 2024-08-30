B#!/bin/bash
set -ex

#module unload cdo/1.9.8-magicsxx-gcc64;module load cdo/1.9.2-magicsxx-intel18

FSTEM=channel_7K
FSTEM=channel_7K_diurnal_200m
#FSTEM=channel_0K_200m
#FSTEM=step_7K
FSTEM=agg_2km_decay_5d
FSTEM=agg_2km
FSTEM=agg_2km_300K
#FSTEM=agg_2km_5K_to_0K
#FSTEM=agg_2km_tanh_5d
FSTEM=agg_2km_300K_ampl0
FSTEM=agg_500m_incr_rad

#FSTEM=higher_vertical_resolution
#FSTEM=agg_2km

#FSTEM=channel_0K_1km
#FSTEM=channel_0K_diurnal_200m_narrow
FSTEM=T0_300K_ampl_0_2km_no_evap_constrad
FSTEM=no_evap_300K
FSTEM=no_evap_300K
#FSTEM=full_evap_300K
FSTEM=no_evap_300K_200m_v1
FSTEM=full_evap_300K_200m
FSTEM=01_evap_300K_200m
FSTEM=removeEvap
#FSTEM=T0_300K_ampl_10_1km_large_Gorm
#FSTEM=agg_2km_incr_rad
#FSTEM=removeEvap

FSTEM=$1
NX=$2
NY=$3
TMIN=$4
TMAX=$5
LEV0=$6
LEV1=$7

#NX=4320
#NY=48
#NX=480
#NY=480
#NX=12
#NX=960
#NX=480
#NY=480
#NX=240
#NY=240
#NX=240
#NY=240
#TMIN=193
#TMAX=288
#TMIN=97
#TMAX=192
#TMIN=289

#LEV0=50
#LEV1=100


./select_timesteps_single_level_input.sh ${FSTEM} ${TMIN} ${TMAX} ${LEV0} ${LEV1}
./gather_single_level_input.sh ${FSTEM} ${TMIN} ${TMAX} ${LEV0} ${LEV1}
#./reduce_stat.sh ${FSTEM}
#python precip_level.py ${FSTEM} ${NX} ${NY} ${TMIN} ${TMAX} ${LEV0}
