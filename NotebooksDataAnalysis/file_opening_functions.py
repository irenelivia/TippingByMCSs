#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cdo import Cdo
cdo = Cdo()

# open STAT file
def open_STAT(STAT_FILE):
    path = '/work/bb1166/Irene/DataAnalysis/Outputs/Stats'
    STAT_file = STAT_FILE
    STAT_filepath = os.path.join(path,STAT_file)
    STAT_DATA = xr.open_dataset(STAT_filepath)
    return STAT_DATA

# open 2D file
def open_2D(TWODIM_FILE):
    path = '/work/bb1166/Irene/DataAnalysis/Outputs/2D'
    TWODIM_file = TWODIM_FILE
    TWODIM_filepath = os.path.join(path,TWODIM_file)
    TWODIM_DATA = xr.open_dataset(TWODIM_filepath)
    return TWODIM_DATA


# open 3D file
def open_3D(TWODIM_FILE):
    path = '/work/bb1166/Irene/Projects/Land_sea_lagrangian/OUT_3D'
    TWODIM_file = TWODIM_FILE
    TWODIM_filepath = os.path.join(path,TWODIM_file)
    TWODIM_DATA = xr.open_dataset(TWODIM_filepath)
    return TWODIM_DATA

# define all variables to look at from stat files
def variables(STAT_DATA):
    sst = STAT_DATA.SST
    lhf = STAT_DATA.LHF
    shf = STAT_DATA.SHF
    pw = STAT_DATA.PW
    precip = STAT_DATA.PREC
    theta = STAT_DATA.THETA
    theta_25 = theta.isel(z=0)
    theta_1300 = theta.sel(z='1300', method='nearest')
    theta_3000 = theta.sel(z='3000', method='nearest')
    tabs = STAT_DATA.TABS
    tabs_25 = tabs.isel(z=0)
    tabs_1300 = tabs.sel(z='1300', method='nearest')
    tabs_3000 = tabs.sel(z='3000', method='nearest')
    qv = STAT_DATA.QV
    qv_25 = qv.isel(z=0)
    olr = STAT_DATA.RADLWUP
    olr_levmax = olr.isel(z=-1)
    return sst, lhf, shf, pw, precip, theta, theta_25, theta_1300, theta_3000, tabs, tabs_25, tabs_1300, tabs_3000, qv, qv_25, olr, olr_levmax