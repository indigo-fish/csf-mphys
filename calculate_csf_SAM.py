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


Code_dir = 'home/wadh5699/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)

import reading_in_data_functions as rd_data
import save_file as sf


def get_SAM_indices(dataset='CSF-20C', season='DJF'):
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
	
	return yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units


def save_SAM_indices(dataset='CSF-20C', season='DJF'):
	
	yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units = get_SAM_indices(dataset, season)
	
	#saves normalized SAM indices as netcdf files
	mean_destination = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_mean_data.nc'
	mean_description = 'mean Marshall SAM index from ' + dataset + ' during ' + season
	save = sf.save_file(mean_destination, mean_description)
	save.add_times(times, calendar, t_units, time_name='time')
	save.add_variable(np.array(mean_SAM_indices), 'SAM index', ('time'))
	save.close_file()
	
	ensemble_destination = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_ensemble_data.nc'
	ensemble_description = 'Marshall SAM index from ' + dataset + ' ensemble during ' + season
	ens_len = len(yearly_SAM_indices[0])
	dim1 = np.arange(0, ens_len, 1)
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


def read_SAM_indices(dataset='CSF-20C', season='DJF'):
	
	mean_source = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_mean_data.nc'
	mean_read = Dataset(mean_source)
	mean_data = mean_read.variables['SAM index'][:]
	
	ensemble_source = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_ensemble_data.nc'
	ensemble_read = Dataset(ensemble_source)
	ensemble_data = mean_read.variables['SAM index'][:]
	
	variation_source = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/' + dataset + '_' + season + '_sam_variation_data.nc'
	variation_read = Dataset(variation_source)
	variation_data = variation_read.variables['SAM index'][:]
	
	times, calendar, units = rd_data.read_time_dimension(mean_source, time_name = 'time')
	
	return mean_data, times, calendar, units


def mask_era(all_SAM_indices, times, calendar, units, season='DJF'):
	
	#reduces data to appropriate season in same way as reading_in_data_functions.calculate_annual_mean
	try: dates = num2date(times, calendar=calendar, units=units)
	except: dates = num2date(times, units=units)
	months = np.array([])
	years = np.array([])
	for day in dates:
		months = np.append(months, day.timetuple()[1])
		years = np.append(years, day.timetuple()[0])
	
	return_values = np.array([])
	for i,yr in enumerate(np.unique(years)):
		if season == 'DJF': mask = ((years==(yr-1))&(months==12)) | (years==yr) & ((months==1)|(months==2))
		elif season == 'MAM': mask = (years==yr) & ((months == 3)|(months==4)|(months==5))
		elif season == 'JJA': mask = (years==yr) & ((months==6)|(months==7)|(months==8))
		elif season == 'SON': mask = (years==yr) & ((months==9)|(months==10)|(months==11))
		elif (season==None)|(season=='ANN'): mask = (years==yr)
		else:
			print('Season is not valid')
			raise NameError
		return_values = np.append(return_values, np.mean(all_SAM_indices[mask]))
	
	return return_values, np.unique(years)


def get_era_SAM_indices(season='DJF'):
	
	#reads in ERA5 sea level pressure data and uses it to produce same type of SAM index
	era5_file = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_195001-197812.nc'
	era_slp, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = rd_data.read_in_variable(era5_file, 'psl')
	era5_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_197901-202012.nc'
	era_slp2, _, _, _, era_times2, _, _ = rd_data.read_in_variable(era5_file2, 'psl')
	era_slp = np.append(era_slp, era_slp2, axis=0)
	era_times = np.append(era_times, era_times2)
	
	all_year_era_SAM_indices = np.array([])
	for era_year_slp in era_slp:
		era_slp_65S = era_year_slp[10] #index 10 corresponds to latitude 65S
		era_slp_40S = era_year_slp[18] #index 18 corresponds to latitude 40S
		
		zm_40S = np.mean(era_slp_40S) #takes zonal mean of mslp at 40S
		zm_65S = np.mean(era_slp_65S) #takes zonal mean of mslp at 65S
		
		era_SAM_index = zm_40S - zm_65S #subtracts surface pressures to find unnormalized SAM index
		all_year_era_SAM_indices= np.append(all_year_era_SAM_indices, era_SAM_index) #stores ERA5 sam index for each year
	
	#eliminates irrelevant seasons
	era_SAM_indices, era_years = mask_era(all_year_era_SAM_indices, era_times, era_calendar, era_t_units, season)
	
	#normalizes ERA SAM indices
	era_mean_norm = np.mean(era_SAM_indices)
	era_std_norm = np.std(era_SAM_indices)
	
	era_SAM_indices -= era_mean_norm
	era_SAM_indices /= era_std_norm
	
	return era_SAM_indices, era_years


def graph_SAM_indices(dataset='CSF-20C', season='DJF'):
	
	yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units = get_SAM_indices(dataset, season)
	
	#reads in offical Marshall SAM index data from text file
	Marshall_SAM_index, years_SAM = rd_data.read_Marshall_SAM_idx(season)
	
	#does same mean surface level pressure calculations for ERA5 data
	era_SAM_indices, era_years = get_era_SAM_indices(season)
	
	#displays plots of ensemble and mean SAM indices
	plt.figure(1)
	plt.plot(times, yearly_SAM_indices, color='gray')
	plt.plot(times, mean_SAM_indices, linewidth = 2, color = 'black', label = 'mean')
	plt.plot(years_SAM, Marshall_SAM_index, linewidth = 2, color = 'red', label = 'Marshall data')
	plt.plot(era_years, era_SAM_indices, linewidth = 2, color= 'blue', label = 'ERA5 data')
	#plt.errorbar(times, mean_SAM_indices, yerr=SAM_stdevs, color='black', label='standard deviation')
	plt.title('Normalized SAM Index in ' + dataset + ' Ensemble During ' + season)
	plt.xlabel('Year')
	plt.ylabel('Normalized SAM Index')
	plt.legend()
	figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_' + season + '_Normalized_SAM.png'
	plt.savefig(figure_name)
	plt.show()


def truncate_to_pairs(times1, data1, times2, data2):
	#creates arrays containing only data for years in both datasets provided
	paired1 = np.array([])
	paired2 = np.array([])
	paired_times = np.array([])
	
	for time1, point1 in zip(times1, data1):
		for time2, point2 in zip(times2, data2):
			if time1 == time2:
				paired1 = np.append(paired1, point1)
				paired2 = np.append(paired2, point2)
				paired_times = np.append(paired_times, time1)
	
	return paired1, paired2, paired_times


def correlate_pairs(data1, data2, label1, label2, season='DJF', smoothing=None):
	#produces scatter plots of SAM indices in the same years in two different datasets
	plt.plot(data1, data2, 'o', label='original data')
	res = linregress(data1, data2)
	
	lineplot = []
	for point in data1:
		lineplot.append(point * res.slope + res.intercept)
	
	plt.plot(data1, lineplot, 'r', label='fitted line')
	titlestr = 'Relationship between ' + label1 + ' and ' + label2 + ' SAM index during ' + season
	if not smoothing==None: titlestr += ' with ' + str(smoothing) + '-year average'
	plt.title(titlestr)
	plt.xlabel(label1 + ' SAM')
	plt.ylabel(label2 + ' SAM')
	plt.legend()
	plt.annotate(f"R squared: {res.rvalue**2:.6f}", (0.2, 0.2), xycoords='axes fraction')
	figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + label1 + '_' + label2 + '_' + season + '_SAM_correlation'
	if smoothing==None: figure_name += '.png'
	else: figure_name += '_' + str(smoothing) + '_average.png'
	plt.savefig(figure_name)
	plt.show()


def running_mean(data, years, timescale):
	#produces a cropped dataset of ten year averages of a longer dataset
	averaged_data = uniform_filter1d(np.array(data), timescale)
	return averaged_data[int(timescale/2):len(averaged_data) - int(timescale/2)], years[int(timescale/2):len(averaged_data) - int(timescale/2)]


def smooth_and_plot(data1, data2, label1, label2, times, timescale, season):
	#takes running mean of a pair of datasets and plots the comparison between them
	smoothed_data1, _ = running_mean(data1, times, timescale)
	smoothed_data2, _ = running_mean(data2, times, timescale)
	correlate_pairs(smoothed_data1, smoothed_data2, label1, label2, season=season, smoothing=timescale)


def compare_smoothings(dataset='CSF-20C', season='DJF'):
	#reads in seasonal forecast SAM indices
	mean_SAM_indices, times, calendar, t_units = read_SAM_indices(dataset, season)
	#reads in official Marshall SAM index data from text file
	Marshall_SAM_index, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	#reads in ERA SAM indices
	era_SAM_indices, era_years = get_era_SAM_indices(season)
	#we want to do same analysis first for Marshall index, then for ERA5 index
	pairs = [[Marshall_SAM_index, Marshall_years, 'Marshall'], [era_SAM_indices, era_years, 'ERA5']]
	
	smoothings = np.arange(1, 40, 1)
	for pair in pairs:
		r_squares = []
		paired_my_index, paired_other_index, paired_times = truncate_to_pairs(times, mean_SAM_indices, pair[1], pair[0])
		for smoothing in smoothings:
			smoothed_my_index, _ = running_mean(paired_my_index, paired_times, smoothing)
			smoothed_other_index, _ = running_mean(paired_other_index, paired_times, smoothing)
			res = linregress(smoothed_my_index, smoothed_other_index)
			r_squares.append(res.rvalue**2)
		plt.plot(smoothings, r_squares)
		plt.xlabel('years averaged')
		plt.ylabel('R squared value')
		plt.title('Correlation strength between ' + dataset + ' and ' + pair[2] + ' for different time scales')
		figure_name = figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_' + pair[2] + '_' + season + '_SAM_correlation_depending_on_averaging.png'
		plt.savefig(figure_name)
		plt.show()


def stat_analysis(dataset='CSF-20C', season='DJF'):
	#reads in seasonal forecast SAM indices
	mean_SAM_indices, times, calendar, t_units = read_SAM_indices(dataset, season)
	
	#reads in offical Marshall SAM index data from text file
	Marshall_SAM_index, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	
	#reads in ERA SAM indices
	era_SAM_indices, era_years = get_era_SAM_indices(season)
	
	#we want to do same analysis first for Marshall index, then for ERA5 index
	pairs = [[Marshall_SAM_index, Marshall_years, 'Marshall'], [era_SAM_indices, era_years, 'ERA5']]
	
	for pair in pairs:
		#creates arrays containing only SAM indices for years in both forecast dataset and Marshall data
		paired_my_index, paired_other_index, paired_times = truncate_to_pairs(times, mean_SAM_indices, pair[1], pair[0])
		#plots linear regression comparing forecast to other SAM indices
		correlate_pairs(paired_my_index, paired_other_index, dataset, pair[2], season=season)
		
		#repeats same process for different lengths of running means
		timescales = [2, 3, 5, 10, 15, 20, 30]
		for timescale in timescales:
			smooth_and_plot(paired_my_index, paired_other_index, dataset, pair[2], paired_times, timescale, season)
	
	#compares ERA5 data and Marshall data
	paired_era_index2, paired_Marshall_index2, Marshall_and_era_times = truncate_to_pairs(era_years, era_SAM_indices, Marshall_years, Marshall_SAM_index)
	correlate_pairs(paired_era_index2, paired_Marshall_index2, 'ERA5', 'Marshall', season=season)


def full_analysis(dataset='CSF-20C', season='DJF'):
	save_SAM_indices(dataset, season)
	graph_SAM_indices(dataset, season)
	stat_analysis(dataset, season)


#runs code
datasets = ['CSF-20C', 'ASF-20C']
seasons = ['DJF', 'JJA']
for dataset in datasets:
	for season in seasons:
		compare_smoothings(dataset=dataset, season=season)
"""
dataset = 'ASF-20C'
seasons = 'MAM', 'JJA', 'SON'
for season in seasons:
	full_analysis(dataset, season)

dataset = 'CSF-20C'
season = 'SON'
full_analysis(dataset, season)
"""
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
