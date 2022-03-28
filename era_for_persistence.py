#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:31:50 2022

@author: wadh5699
"""

import sys
import numpy as np
from netCDF4 import num2date, Dataset

Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)

import reading_in_data_functions as rd_data

def mask_era_persistence(all_SAM_indices, times, calendar, units, season='DJF'):
	#reduces ERA SAM index data to desired season in same way as reading_in_data_functions.calculate_annual_mean
	try: dates = num2date(times, calendar=calendar, units=units)
	except: dates = num2date(times, units=units)
	months = np.array([])
	years = np.array([])
	for day in dates:
		months = np.append(months, day.timetuple()[1])
		years = np.append(years, day.timetuple()[0])
	
	return_values = np.array([])
	for i,yr in enumerate(np.unique(years)):
		if season == 'DJF': mask = (years==yr-1)&(months==11)
		elif season == 'MAM': mask = (years==yr)&(months==2)
		elif season == 'JJA': mask = (years==yr)&(months==5)
		elif season == 'SON': mask = (years==yr)&(months==8)
		elif (season==None)|(season=='ANN'): mask = (years==yr)
		else:
			print('Season is not valid')
			raise NameError
		return_values = np.append(return_values, np.nanmean(all_SAM_indices[mask]))
	
	return_values = return_values[1:]
	years = np.unique(years)
	years = years[1:]
	return return_values, years


def get_era_SAM_indices(season='DJF'):
	
	#reads in ERA5 sea level pressure data and uses it to produce same type of SAM index
	era5_file = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_195001-197812.nc'
	era_slp, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = rd_data.read_in_variable(era5_file, 'psl')
	era5_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_197901-202012.nc'
	era_slp2, _, _, _, era_times2, _, _ = rd_data.read_in_variable(era5_file2, 'psl')
	era_slp = np.append(era_slp, era_slp2, axis=0)
	era_times = np.append(era_times, era_times2)
	
	#stores arrays of zonal mean slp from ERA5 at both latitudes
	mean_pressures_40S = np.array([])
	mean_pressures_65S = np.array([])
	for era_year_slp in era_slp:
		era_slp_65S = era_year_slp[10] #index 10 corresponds to latitude 65S
		era_slp_40S = era_year_slp[18] #index 18 corresponds to latitude 40S
		zm_40S = np.mean(era_slp_40S) #takes zonal mean of mslp at 40S
		zm_65S = np.mean(era_slp_65S) #takes zonal mean of mslp at 65S
		mean_pressures_40S = np.append(mean_pressures_40S, zm_40S)
		mean_pressures_65S = np.append(mean_pressures_65S, zm_65S)
	
	#eliminates irrelevant seasons from each dataset
	mean_pressures_40S, era_years = mask_era_persistence(mean_pressures_40S, era_times, era_calendar, era_t_units, season=season)
	mean_pressures_65S, era_years = mask_era_persistence(mean_pressures_65S, era_times, era_calendar, era_t_units, season=season)
	
	#takes seasonal average of zonal mean surface pressures at both latitudes
	norm_40S = np.nanmean(mean_pressures_40S)
	norm_65S = np.nanmean(mean_pressures_65S)
	
	era_SAM_indices = np.array([])
	for pressure_40S, pressure_65S in zip(mean_pressures_40S, mean_pressures_65S):
		era_SAM_index = pressure_40S / norm_40S - pressure_65S / norm_65S #subtracts normalized surface pressures to find unnormalized SAM index
		era_SAM_indices = np.append(era_SAM_indices, era_SAM_index)
	
	#normalizes ERA SAM indices
	era_mean_norm = np.nanmean(era_SAM_indices)
	era_std_norm = np.std(era_SAM_indices)
	
	era_SAM_indices -= era_mean_norm
	era_SAM_indices /= era_std_norm
	
	return era_SAM_indices, era_years