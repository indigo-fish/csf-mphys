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

def conversion(dataset='CSF-20C', season='DJF'):
	#sets aspects of filename/directory based on dataset and season used
	if dataset == 'CSF-20C':
		years = np.arange(1901, 2011, 1)
		code_str = 'guh4'
		if season == 'DJF':
			start_month = '11'
			ens_len = 51
		elif season == 'MAM':
			start_month = '02'
			ens_len = 25
		elif season == 'JJA':
			start_month = '05'
			ens_len = 51
		elif season == 'SON':
			start_month = '08'
			ens_len = 25
	elif dataset == 'ASF-20C':
		years = np.arange(1900, 2010, 1)
		ens_len = 51
		if season == 'DJF':
			code_str = 'b1gl'
			start_month = '11'
		elif season == 'MAM':
			code_str = 'gfi6'
			start_month = '02'
		elif season == 'JJA':
			code_str = 'gh6j'
			start_month = '05'
		elif season == 'SON':
			code_str = 'gir6'
			start_month = '08'
	
	#reads in surface level pressure at 40S and 65S and saves as netcdf
	ens_len = 10 #this is just to test the code more quickly
	years = np.arange(1950, 1980, 1) #this is just to test the code more quickly
	ensemble = np.arange(0, ens_len, 1)
	annual_surface_pressures = []
	for year in years:
		ensemble_surface_pressures = []
		for member in ensemble:
			file_name = '/network/group/aopp/predict/AWH002_BEFORT_SEASONAL/' + dataset + '/ecmf/' + code_str + '/sfc/fcmean/' + str(year) + start_month + '01/var151/ecmf-' + code_str + '_' + str(member) + '_' + str(year) + start_month + '01_fcmean_sfc.grb'
			msl, lats, lons, levs = rgrb.read_in_grib(file_name,'msl',lead=1)
			
			msl_40S = msl[192,:] #index 192 in latitudes corresponds to 40S
			msl_65S = msl[220,:] #index 220 in latitudes corresponds to 65S
			ensemble_surface_pressures.append([msl_40S, msl_65S])
		annual_surface_pressures.append(ensemble_surface_pressures)
	
	file_name = dataset + '_' + season + '_msl_data.nc'
	description = 'mean surface level pressure at 40S and 65S from ' + dataset + ' for ' + str(ens_len) + ' ensemble members during season ' + season
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

conversion()
