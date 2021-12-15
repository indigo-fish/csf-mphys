#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 2021

@author: wadh5699
"""

import os, sys
import numpy as np

Code_dir = 'home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'

sys.path.append(Code_dir)
import read_in_grib_files as rgrb
import save_file as svf

#reads in surface level pressure at 40S and 65S and saves as netcdf
years = np.arange(1981, 2001, 1)
ensemble = np.arange(0, 25, 1)
annual_surface_pressures = []
for year in years:
	ensemble_surface_pressures = []
	for member in ensemble:
		file_name = '/network/group/aopp/predict/AWH002_BEFORT_SEASONAL/CSF-20C/ecmf/guh4/sfc/fcmean/' + str(year) + '1101/var151/ecmf-guh4_' + str(member) + '_' + str(year) + '1101_fcmean_sfc.grb'
		msl, lats, lons, levs = rgrb.read_in_grib(file_name,'msl',lead=1)
		
		msl_40S = msl[192,:] #index 192 in latitudes corresponds to 40S
		msl_65S = msl[220,:] #index 220 in latitudes corresponds to 65S
		ensemble_surface_pressures.append([msl_40S, msl_65S])
	annual_surface_pressures.append(ensemble_surface_pressures)

file_name = 'msl_data.nc'
description = 'mean surface level pressure at 40S and 65S for 25 ensemble members'
dim1 = np.array(ensemble)
dim2 = np.array([-40, -65])
dim3 = np.array(lons)
times = np.array(years)

save = svf.save_file(file_name, description)
save.add_dimension(dim1, 'ensemble member')
save.add_dimension(dim2, 'latitude')
save.add_dimension(dim3, 'longitude')
save.add_times(times, 'standard', 'Gregorian_year', time_name='time')
save.add_variable(annual_surface_pressures, 'mslp for SAM', ('time', 'ensemble member', 'latitude', 'longitude'))
save.close_file()
