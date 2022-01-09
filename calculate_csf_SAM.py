#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 2021

@author: wadh5699
"""

import os, sys
import numpy as np
from netCDF4 import num2date

Code_dir = 'home/wadh5699/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)

import reading_in_data_functions as rd_data
import save_file as sf


def analysis(dataset, season):
	#reads in netcdf file of relevant surface pressures
	file_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_msl_data.nc'
	mslp_data, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file_name, 'mslp for SAM', lat_name='latitude', lon_name='longitude', time_name='time')
	
	#takes zonal mean of surface pressures and subtracts to find an unnormalized SAM index
	yearly_SAM_indices = []
	mean_SAM_indices = []
	SAM_stdevs = []
	for year in mslp_data:
		ens_SAM_indices = []
		for ensemble in year:
			msl_40S = ensemble[0] #mslp at 40S is stored at netcdf index 0
			msl_65S = ensemble[1] #mslp at 65S is stored at netcdf index 1
			
			zm_40S = np.mean(msl_40S) #takes zonal mean of mslp at 40S
			zm_65S = np.mean(msl_65S) #takes zonal mean of mslp at 65S
			
			SAM_index = zm_40S - zm_65S #subtracts surface pressures to find unnormalized SAM index
			ens_SAM_indices.append(SAM_index)
		yearly_SAM_indices.append(ens_SAM_indices) #stores vectors of SAM indices from all ensemble members for each year
		mean_SAM_indices.append(np.mean(ens_SAM_indices)) #stores ensemble mean SAM index for each year
		SAM_stdevs.append(np.std(ens_SAM_indices)) #stores standard deviation of SAM index for each year
	
	
	#normalizes SAM indices
	mean_norm = np.mean(mean_SAM_indices)
	std_norm = np.std(mean_SAM_indices)
	
	mean_SAM_indices -= mean_norm
	mean_SAM_indices /= std_norm
	yearly_SAM_indices -= mean_norm
	yearly_SAM_indices /= std_norm
	
	
	#saves normalized SAM indices as netcdf files
	mean_destination = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_mean_data.nc'
	mean_description = 'mean Marshall SAM index from ' + dataset + ' during ' + season
	save = sf.save_file(mean_destination, mean_description)
	save.add_times(times, calendar, t_units, time_name='time')
	save.add_variable(np.array(mean_SAM_indices), 'SAM index', ('time'))
	save.close_file()
	
	ensemble_destination = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_ensemble_data.nc'
	ensemble_description = 'Marshall SAM index from ' + dataset + ' ensemble during ' + season
	dim1 = np.arange(0, 25, 1)
	save2 = sf.save_file(ensemble_destination, ensemble_description)
	save2.add_dimension(dim1, 'ensemble member')
	save2.add_times(times, calendar, t_units, time_name='time')
	save2.add_variable(np.array(yearly_SAM_indices), 'SAM index', ('time', 'ensemble member'))
	save2.close_file()
	
	variation_destination = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_variation_data.nc'
	variation_description = 'standard deviation of Marshall SAM index from ' + dataset + ' during ' + season
	save3 = sf.save_file(variation_destination, variation_description)
	save3.add_times(times, calendar, t_units, time_name='time')
	save3.add_variable(np.array(SAM_stdevs), 'SAM index', ('time'))
	save3.close_file()
	
	#reads in offical Marshall SAM index data from text file
	Marshall_SAM_index, years_SAM = rd_data.read_Marshall_SAM_idx(season)
	
	#does same mean surface level pressure calculations for ERA5 data
	era5_file = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_195001-197812.nc'
	era_slp, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = rd_data.read_in_variable(era5_file, 'psl')
	era5_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_197901-202012.nc'
	era_slp2, _, _, _, era_times2, _, _ = rd_data.read_in_variable(era5_file2, 'psl')
	era_slp = np.append(era_slp, era_slp2, axis=0)
	era_times = np.append(era_times, era_times2)
	try: era_dates = num2date(era_times, calendar=era_calendar, units=era_t_units)
	except: era_dates = num2date(era_times, units=era_t_units)
	era_years = np.array([])
	for day in era_dates:
		era_years = np.append(era_years, day.timetuple()[0])
	
	era_SAM_indices = []
	for era_year_slp in era_slp:
		era_slp_65S = era_year_slp[10] #index 10 corresponds to latitude 65S
		era_slp_40S = era_year_slp[18] #index 18 corresponds to latitude 40S
		
		zm_40S = np.mean(era_slp_40S) #takes zonal mean of mslp at 40S
		zm_65S = np.mean(era_slp_65S) #takes zonal mean of mslp at 65S
		
		SAM_index = zm_40S - zm_65S #subtracts surface pressures to find unnormalized SAM index
		era_SAM_indices.append(SAM_index) #stores ERA5 sam index for each year
	
	#normalizes ERA SAM indices
	era_mean_norm = np.mean(era_SAM_indices)
	era_std_norm = np.std(era_SAM_indices)
	
	era_SAM_indices -= era_mean_norm
	era_SAM_indices /= era_std_norm
	
	
	#displays plots of ensemble and mean SAM indices
	import matplotlib.pyplot as plt
	
	plt.figure(1)
	plt.plot(times, yearly_SAM_indices, color='gray')
	plt.plot(times, mean_SAM_indices, linewidth = 2, color = 'black', label = 'mean')
	plt.plot(years_SAM, Marshall_SAM_index, linewidth = 2, color = 'red', label = 'Marshall data')
	#plt.plot(era_years, era_SAM_indices, linewidth = 2, color= 'blue', label = 'ERA5 data')
	#plt.errorbar(times, mean_SAM_indices, yerr=SAM_stdevs, color='black', label='standard deviation')
	plt.title('Unnormalized SAM Index in ' + dataset + ' Ensemble')
	plt.xlabel('Year')
	plt.ylabel('Unnormalized SAM Index')
	figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_Unnormalized_SAM.png'
	plt.savefig(figure_name)
	plt.show()



dataset = 'CSF-20C' #allows switching between CSF-20C and ASF-20C data
season = 'DJF'

analysis(dataset, season)
"""
years = np.arange(1981, 2001, 1)
ensemble = np.arange(0, 25, 1)
ens_SAM_indices = []
for member in ensemble:
	SAM_indices = []
	for year in years:		
		msl_40S = msl[192,:] #index 192 in latitudes corresponds to 40S
		msl_65S = msl[220,:] #index 220 in latitudes corresponds to 65S
		
		zm_40S = sum(msl_40S) / len(msl_40S) #zonal mean of 40S mean surface level pressure
		zm_65S = sum(msl_65S) / len(msl_65S) #zonal mean of 65S mean surface level pressure
		norm_40S = 10000 #I don't know what the normalization factor for 40S should be
		norm_65S = 10000 #I don't know what the normalization factor for 65S should be
		SAM_index = zm_40S / norm_40S - zm_65S / norm_65S #SAM index defined as difference of normalized zonal mean surface pressures
		SAM_indices.append(SAM_index)
	ens_SAM_indices.append(SAM_indices)

mean_SAM_indices = []
indices = np.arange(0, len(years), 1)
for index in indices:
	all_SAM_indices = []
	for member in ens_SAM_indices:
		all_SAM_indices.append(member[index])
	mean = sum(all_SAM_indices) / len(all_SAM_indices)
	mean_SAM_indices.append(mean)

import matplotlib.pyplot as plt

#for SAM_indices in ens_SAM_indices:
#	plt.plot(years, SAM_indices)
plt.plot(years, mean_SAM_indices)
plt.show()
figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/my_mean_SAM_indices'
plt.savefig(figure_name)
"""
