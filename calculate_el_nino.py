#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 2021

@author: wadh5699
"""

import os, sys
import numpy as np
from netCDF4 import num2date, Dataset
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d
import matplotlib.pyplot as plt


Code_dir = '/home/w/wadh5699/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)
#sets directories for accessing netcdf files and figures
Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
Figure_dir = Code_dir + 'Figures/'
Data_dir = Code_dir + 'Data/'

import reading_in_data_functions as rd_data
import save_file as sf
import significance_testing as sig_test

def save_nino34(dataset='CSF-20C', season='DJF'):
	#reads in ssts to be used for calculating El Nino index
	#latitude and longitude boundaries have previously been determined manually
	if dataset=='CSF-20C':
		file = '/home/p/pattersonm/sst_ecmf-guh4_1101_fcmean_sfc_1901_2010_' + season + '.nc'
		sst, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file, 'sst')
		lat_lower, lat_upper, lon_lower, lon_upper = 120, 135, 270, 341
	elif dataset=='SEAS5':
		file = '/network/aopp/hera/mad/patterson/MPhys/SEAS5_ensemble_means/SEAS5_sst_ensemble_mean_25members_' + season + '_init_November_1982_2017.nc'
		sst, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file, 'sst')
		times = np.arange(1982, 2017+1)
		csf_file = '/home/p/pattersonm/sst_ecmf-guh4_1101_fcmean_sfc_1901_2010_' + season + '.nc'
		csf_times, calendar, t_units = rd_data.read_time_dimension(csf_file, time_name='time')
		lat_lower, lat_upper, lon_lower, lon_upper = 85, 95, 190, 240
	elif dataset=='Marshall' or dataset=='ASF-20C' or dataset=='ERA5':
		file = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_2.5x2.5_195001-197812.nc'
		sst, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file, 'tos')
		file2 = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_2.5x2.5_197901-202012.nc'
		sst2, _, _, _, times2, _, _ = rd_data.read_in_variable(file2, 'tos')
		sst = np.append(sst, sst2, axis=0)
		times = np.append(times, times2)
		sst, times = rd_data.calculate_annual_mean(sst, times, calendar, t_units, season=season)
		lat_lower, lat_upper, lon_lower, lon_upper = 34, 38, 76, 96
	
	#creates timeseries of mean SST in chosen region
	unnormalized_nino = []
	for data_year in sst:
		nino_ssts = data_year[lat_lower:lat_upper, lon_lower:lon_upper]
		unnormalized_nino.append(np.mean(nino_ssts))
	
	#subtracts mean and divides by standard deviation to result in normalised index
	nino_mean = np.mean(unnormalized_nino)
	nino_std = np.std(unnormalized_nino)
	nino34 = np.array((unnormalized_nino - nino_mean) / nino_std)
	
	#saves Nino 3.4 index as netcdf
	destination = Data_dir + dataset + '_' + season + '_nino34_data.nc'
	nino_description = 'Nino 3.4 index from ' + dataset + ' SSTs during ' + season
	save = sf.save_file(destination, nino_description)
	save.add_times(times, calendar, t_units, time_name='time')
	save.add_variable(nino34, 'Nino 3.4 index', ('time'))
	save.close_file()

def read_nino34(dataset='CSF-20C', season='DJF'):
	#reads and returns Nino 3.4 index
	source = Data_dir + dataset + '_' + season + '_nino34_data.nc'
	read_source = Dataset(source)
	data = read_source.variables['Nino 3.4 index'][:]
	
	times, calendar, units = rd_data.read_time_dimension(source, time_name = 'time')
	
	return data, times, calendar, units

#go to run_sam_analysis to run code because don't want code to execute when imported into multiseason_analaysis
datasets = ['CSF-20C', 'SEAS5', 'ERA5']
for dataset in datasets:
	save_nino34(dataset=dataset, season='DJF')
	data, times, calendar, units = read_nino34(dataset=dataset, season='DJF')
	print(data)
