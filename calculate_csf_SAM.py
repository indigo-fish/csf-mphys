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

def get_zm_pressures(dataset='CSF-20C', season='DJF'):
	#returns zonally averaged pressures at 40S and 65S from netcdf files, both ensemble mean and individual ensemble members
	
	#reads in netcdf file of relevant surface pressures
	file_name = Data_dir + dataset + '_' + season + '_msl_data.nc'
	mslp_data, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file_name, 'mslp for SAM', lat_name='latitude', lon_name='longitude', time_name='time')
	
	#takes zonal mean of surface pressures to be plotted
	ensemble_pressures_40S = []
	ensemble_pressures_65S = []
	mean_pressures_40S = []
	mean_pressures_65S = []
	for year in mslp_data:
		pressures_40S = []
		pressures_65S = []
		for ensemble in year:
			msl_40S = ensemble[0] #mslp at 40S is stored at netcdf index 0
			msl_65S = ensemble[1] #mslp at 65S is stored at netcdf index 1
			
			zm_40S = np.mean(msl_40S) #takes zonal mean of mslp at 40S
			zm_65S = np.mean(msl_65S) #takes zonal mean of mslp at 65S
			
			pressures_40S.append(zm_40S)
			pressures_65S.append(zm_65S)
		mean_pressures_40S.append(np.mean(pressures_40S))
		mean_pressures_65S.append(np.mean(pressures_65S))
		ensemble_pressures_40S.append(pressures_40S)
		ensemble_pressures_65S.append(pressures_65S)
	
	return mean_pressures_40S, mean_pressures_65S, ensemble_pressures_40S, ensemble_pressures_65S

def get_SAM_indices(dataset='CSF-20C', season='DJF'):
	#reads in pressure data from netcdf files, converts it to netcdf files, and returns it
	
	#reads in netcdf file of surface pressures to access times, calendar, and time units
	file_name = Data_dir + dataset + '_' + season + '_msl_data.nc'
	mslp_data, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file_name, 'mslp for SAM', lat_name='latitude', lon_name='longitude', time_name='time')
	
	#gets ensemble mean of zonal mean surface pressures at both latitudes
	mean_pressures_40S, mean_pressures_65S, _, _ = get_zm_pressures(dataset=dataset, season=season)
	
	#takes seasonal average of zonal mean surface pressures at both latitudes
	norm_40S = np.mean(mean_pressures_40S)
	norm_65S = np.mean(mean_pressures_65S)
	
	#gets arrays of full ensemble and ensemble mean surface pressures, and each year's standard deviation
	ensemble_SAM_indices = []
	mean_SAM_indices = []
	SAM_stdevs = []
	for year in mslp_data:
		ens_SAM_indices = []
		for ensemble in year:
			msl_40S = ensemble[0] #mslp at 40S is stored at netcdf index 0
			msl_65S = ensemble[1] #mslp at 65S is stored at netcdf index 1
			
			zm_40S = np.mean(msl_40S) #takes zonal mean of mslp at 40S
			zm_65S = np.mean(msl_65S) #takes zonal mean of mslp at 65S
			SAM_index = zm_40S / norm_40S - zm_65S / norm_65S #subtracts normalized surface pressures to find unnormalized SAM index
			ens_SAM_indices.append(SAM_index)
		ensemble_SAM_indices.append(ens_SAM_indices) #stores vectors of SAM indices from all ensemble members for each year
		mean_SAM_indices.append(np.mean(ens_SAM_indices)) #stores ensemble mean SAM index for each year
		SAM_stdevs.append(np.std(ens_SAM_indices)) #stores standard deviation of SAM index for each year
	
	
	#normalizes SAM indices
	mean_norm = np.mean(mean_SAM_indices)
	std_norm = np.std(mean_SAM_indices)
	
	mean_SAM_indices -= mean_norm
	mean_SAM_indices /= std_norm
	ensemble_SAM_indices -= mean_norm
	ensemble_SAM_indices /= std_norm
	SAM_stdevs /= std_norm
	
	return ensemble_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units


def save_SAM_indices(dataset='CSF-20C', season='DJF'):
	#reads in arrays of SAM indices and stores as netcdf
	yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units = get_SAM_indices(dataset, season)
	
	#saves normalized SAM indices as netcdf files
	mean_destination = Data_dir + dataset + '_' + season + '_sam_mean_data.nc'
	mean_description = 'mean Marshall SAM index from ' + dataset + ' during ' + season
	save = sf.save_file(mean_destination, mean_description)
	save.add_times(times, calendar, t_units, time_name='time')
	save.add_variable(np.array(mean_SAM_indices), 'SAM index', ('time'))
	save.close_file()
	
	ensemble_destination = Data_dir + dataset + '_' + season + '_sam_ensemble_data.nc'
	ensemble_description = 'Marshall SAM index from ' + dataset + ' ensemble during ' + season
	ens_len = len(yearly_SAM_indices[0])
	dim1 = np.arange(0, ens_len, 1)
	save2 = sf.save_file(ensemble_destination, ensemble_description)
	save2.add_dimension(dim1, 'ensemble member')
	save2.add_times(times, calendar, t_units, time_name='time')
	save2.add_variable(np.array(yearly_SAM_indices), 'SAM index', ('time', 'ensemble member'))
	save2.close_file()
	
	variation_destination = Data_dir + dataset + '_' + season + '_sam_variation_data.nc'
	variation_description = 'standard deviation of Marshall SAM index from ' + dataset + ' during ' + season
	save3 = sf.save_file(variation_destination, variation_description)
	save3.add_times(times, calendar, t_units, time_name='time')
	save3.add_variable(np.array(SAM_stdevs), 'SAM index', ('time'))
	save3.close_file()


def read_SAM_indices(dataset='CSF-20C', season='DJF'):
	#reads and returns ensemble mean SAM index from netcdf
	mean_source = Data_dir + dataset + '_' + season + '_sam_mean_data.nc'
	mean_read = Dataset(mean_source)
	mean_data = mean_read.variables['SAM index'][:]
	
	#I think because this one has different dimensions from the others I need to change the way it's read in
	#ensemble_source = Data_dir + dataset + '_' + season + '_sam_ensemble_data.nc'
	#ensemble_read = Dataset(ensemble_source)
	#ensemble_data = mean_read.variables['SAM index'][:]
	
	#variation_source = Data_dir + dataset + '_' + season + '_sam_variation_data.nc'
	#variation_read = Dataset(variation_source)
	#variation_data = variation_read.variables['SAM index'][:]
	
	times, calendar, units = rd_data.read_time_dimension(mean_source, time_name = 'time')
	
	return mean_data, times, calendar, units


def read_SAM_ensemble(dataset='CSF-20C', season='DJF'):
    ensemble_source = Data_dir + dataset + '_' + season + '_sam_ensemble_data.nc'
    ensemble_read = Dataset(ensemble_source)
    ensemble_data = ensemble_read.variables['SAM index'][:]
    
    times, calendar, units = rd_data.read_time_dimension(ensemble_source, time_name = 'time')
    
    return ensemble_data, times, calendar, units


def mask_era(all_SAM_indices, times, calendar, units, season='DJF'):
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
	mean_pressures_40S, era_years = mask_era(mean_pressures_40S, era_times, era_calendar, era_t_units, season=season)
	mean_pressures_65S, era_years = mask_era(mean_pressures_65S, era_times, era_calendar, era_t_units, season=season)
	
	#takes seasonal average of zonal mean surface pressures at both latitudes
	norm_40S = np.mean(mean_pressures_40S)
	norm_65S = np.mean(mean_pressures_65S)
	
	era_SAM_indices = np.array([])
	for pressure_40S, pressure_65S in zip(mean_pressures_40S, mean_pressures_65S):
		era_SAM_index = pressure_40S / norm_40S - pressure_65S / norm_65S #subtracts normalized surface pressures to find unnormalized SAM index
		era_SAM_indices = np.append(era_SAM_indices, era_SAM_index)
	
	#normalizes ERA SAM indices
	era_mean_norm = np.mean(era_SAM_indices)
	era_std_norm = np.std(era_SAM_indices)
	
	era_SAM_indices -= era_mean_norm
	era_SAM_indices /= era_std_norm
	
	return era_SAM_indices, era_years


def get_seas5_SAM_indices(season='DJF'):
	
	#reads in SEAS5 sea level pressure data and uses it to produce same type of SAM index
	seas5_file = '/network/aopp/hera/mad/patterson/MPhys/SEAS5_ensemble_means/SEAS5_msl_ensemble_mean_25members_DJF_init_November_1982_2017.nc'
	seas_slp, seas_lats, seas_lons, seas_levs, seas_times, seas_calendar, seas_t_units = rd_data.read_in_variable(seas5_file, 'msl', time_name='years')
	
	#stores arrays of zonal mean slp from SEAS5 at both latitudes
	mean_pressures_40S = np.array([])
	mean_pressures_65S = np.array([])
	for seas_year_slp in seas_slp:
		seas_slp_65S = seas_year_slp[10] #index 10 corresponds to latitude 65S
		seas_slp_40S = seas_year_slp[18] #index 18 corresponds to latitude 40S
		zm_40S = np.mean(seas_slp_40S) #takes zonal mean of mslp at 40S
		zm_65S = np.mean(seas_slp_65S) #takes zonal mean of mslp at 65S
		mean_pressures_40S = np.append(mean_pressures_40S, zm_40S)
		mean_pressures_65S = np.append(mean_pressures_65S, zm_65S)
	
	#takes seasonal average of zonal mean surface pressures at both latitudes
	norm_40S = np.mean(mean_pressures_40S)
	norm_65S = np.mean(mean_pressures_65S)
	
	seas_SAM_indices = np.array([])
	for pressure_40S, pressure_65S in zip(mean_pressures_40S, mean_pressures_65S):
		seas_SAM_index = pressure_40S / norm_40S - pressure_65S / norm_65S #subtracts normalized surface pressures to find unnormalized SAM index
		seas_SAM_indices = np.append(seas_SAM_indices, seas_SAM_index)
	
	#normalizes SEAS5 SAM indices
	seas_mean_norm = np.mean(seas_SAM_indices)
	seas_std_norm = np.std(seas_SAM_indices)
	
	seas_SAM_indices -= seas_mean_norm
	seas_SAM_indices /= seas_std_norm
	
	return seas_SAM_indices, seas_times, 'standard', 'Gregorian_year'


def save_seas5_SAM_indices(season='DJF'):
	#reads in arrays of SAM indices and stores as netcdf
	mean_SAM_indices, times, calendar, t_units = get_seas5_SAM_indices(season)
	
	#saves normalized SAM indices as netcdf files
	mean_destination = Data_dir + 'SEAS5_' + season + '_sam_mean_data.nc'
	mean_description = 'mean Marshall SAM index from SEAS5 during ' + season
	save = sf.save_file(mean_destination, mean_description)
	save.add_times(times, calendar, t_units, time_name='time')
	save.add_variable(np.array(mean_SAM_indices), 'SAM index', ('time'))
	save.close_file()


def graph_SAM_indices(dataset='CSF-20C', season='DJF', variance=True, trend=True, cut_years=False):
	#plots CSF SAM index, Marshall SAM index, and ERA SAM index on same graph
	#should be adjusted for optional trendline
	if dataset=='SEAS5': variance=False
	if variance:
		yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units = get_SAM_indices(dataset, season)
	else:
		mean_SAM_indices, times, calendar, t_units = read_SAM_indices(dataset, season)
	if (dataset=='CSF-20C' or dataset=='ASF-20C') and season=='DJF': times += 1
	#reads in offical Marshall SAM index data from text file
	Marshall_SAM_index, years_SAM = rd_data.read_Marshall_SAM_idx(season)
	
	#does same mean surface level pressure calculations for ERA5 data
	era_SAM_indices, era_years = get_era_SAM_indices(season)
	
	if cut_years:
		if dataset=='SEAS5': start_year, end_year = 1982, 2010
		else: start_year, end_year = 1958, 2010
		year_mask = (times >= start_year) & (times <= end_year)
		mean_SAM_indices = mean_SAM_indices[year_mask]
		times = times[year_mask]
		if variance: yearly_SAM_indices = yearly_SAM_indices[year_mask]
		Marshall_mask = (years_SAM >= start_year) & (years_SAM <= end_year)
		Marshall_SAM_index = Marshall_SAM_index[Marshall_mask]
		years_SAM = years_SAM[Marshall_mask]
		era_mask = (era_years >= start_year) & (era_years <= end_year)
		era_SAM_indices = era_SAM_indices[era_mask]
		era_years = era_years[era_mask]
	
	#displays plots of ensemble and mean SAM indices
	plt.clf()
	if variance: plt.plot(times, yearly_SAM_indices, color='gray')
	plt.plot(times, mean_SAM_indices, linewidth = 2, color = 'black', label = 'mean')
	plt.plot(years_SAM, Marshall_SAM_index, linewidth = 2, color = 'red', label = 'Marshall data')
	plt.plot(era_years, era_SAM_indices, linewidth = 2, color= 'blue', label = 'ERA5 data')
	
	if trend:
		#does linear regression to find trend for each dataset
		csf_res = linregress(times, mean_SAM_indices)
		Marshall_res = linregress(years_SAM, Marshall_SAM_index)
		era_res = linregress(era_years, era_SAM_indices)
		
		#generates line to plot for each dataset
		csf_line = []
		for point in times:
			csf_line.append(point * csf_res.slope + csf_res.intercept)
		Marshall_line = []
		for point in years_SAM:
			Marshall_line.append(point * Marshall_res.slope + Marshall_res.intercept)
		era_line = []
		for point in era_years:
			era_line.append(point * era_res.slope + era_res.intercept)
	
		#plots trend lines
		plt.plot(times, csf_line, color='black', linestyle='dashed')
		plt.plot(years_SAM, Marshall_line, color='red', linestyle='dashed')
		plt.plot(era_years, era_line, color='blue', linestyle='dashed')
	
	title = 'SAM Index in ' + dataset + ' Ensemble During ' + season
	if cut_years: title += ': ' + str(start_year) + '-' + str(end_year)
	plt.title('Normalized SAM Index in ' + dataset + ' Ensemble During ' + season)
	plt.xlabel('Year')
	plt.ylabel('Normalized SAM Index')
	plt.legend()
	
	figure_name = Figure_dir + 'Final_' + dataset + '_' + season
	if trend: figure_name += '_trend'
	if cut_years: figure_name += '_cut_years'
	figure_name += '_Normalized_SAM.png'
	print('saving figure to ' + figure_name)
	plt.savefig(figure_name)
	#plt.show()


def graph_stdev(dataset='CSF-20C', season='DJF'):
	#reads in and plots by year standard deviations from netcdf
	variation_source = Data_dir + dataset + '_' + season + '_sam_variation_data.nc'
	variation_read = Dataset(variation_source)
	variation_data = variation_read.variables['SAM index'][:]
	
	times, calendar, units = rd_data.read_time_dimension(variation_source, time_name = 'time')
	
	plt.clf()
	plt.plot(times, variation_data)
	plt.xlabel('Year')
	plt.ylabel('Normalized SAM Index')
	plt.title('Standard Deviation of SAM index in ' + dataset + ' During ' + season)
	#plt.legend()
	figure_name = Figure_dir + dataset + '_' + season + '_Standard_Deviation.png'
	plt.savefig(figure_name)
	plt.show()


def truncate_to_pairs(times1, data1, times2, data2):
	#creates arrays containing only data for years which are included in both datasets
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


def correlate_pairs(data1, data2, label1, label2, season='DJF', smoothing=None, shift_years=False):
	#produces scatter plots of SAM indices in the same years in two different datasets
	plt.clf()
	plt.plot(data1, data2, 'o', label='original data')
	res = linregress(data1, data2)
	
	#produces trendline
	lineplot = []
	for point in data1:
		lineplot.append(point * res.slope + res.intercept)
	
	plt.plot(data1, lineplot, 'r', label='fitted line')
	titlestr = 'Relationship between ' + label1 + ' and ' + label2 + ' SAM index during ' + season
	if smoothing != None: titlestr += ' with ' + str(smoothing) + '-year average'
	plt.title(titlestr)
	plt.xlabel(label1 + ' SAM')
	plt.ylabel(label2 + ' SAM')
	plt.legend()
	if smoothing == None: significance = sig_test.significance(data1, data2, 1)
	else: significance = sig_test.significance(data1, data2, smoothing)
	plt.annotate(f"p-value: {significance:.6f}", (0.2, 0.2), xycoords='axes fraction')
	figure_name = Figure_dir + label1 + '_' + label2 + '_' + season + '_SAM_correlation'
	if smoothing!=None: figure_name += '_' + str(smoothing) + '_average'
	if shift_years: figure_name += '_offset1'
	figure_name += '.png'
	print('saving figure to ' + figure_name)
	plt.savefig(figure_name)
	#plt.show()


def running_mean(data, years, timescale):
	#produces a cropped dataset of averages of a longer dataset on a particular timescale
	averaged_data = uniform_filter1d(np.array(data), timescale)
	return averaged_data[int(timescale/2):len(averaged_data) - int(timescale/2)], years[int(timescale/2):len(averaged_data) - int(timescale/2)]


def smooth_and_plot(data1, data2, label1, label2, times, timescale, season, shift_years = False):
	#takes running mean of a pair of datasets and plots the comparison between them
	smoothed_data1, _ = running_mean(data1, times, timescale)
	smoothed_data2, _ = running_mean(data2, times, timescale)
	correlate_pairs(smoothed_data1, smoothed_data2, label1, label2, season=season, smoothing=timescale, shift_years = shift_years)


def compare_smoothings(dataset='CSF-20C', season='DJF'):
	#plots R squared value as a function of number of years of averaging when comparing two datasets
	
	#reads in seasonal forecast SAM indices
	mean_SAM_indices, times, calendar, t_units = read_SAM_indices(dataset, season)
	#reads in official Marshall SAM index data from text file
	Marshall_SAM_index, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	#reads in ERA SAM indices
	era_SAM_indices, era_years = get_era_SAM_indices(season)
	#we want to do same analysis first for Marshall index, then for ERA5 index
	pairs = [[Marshall_SAM_index, Marshall_years, 'Marshall'], [era_SAM_indices, era_years, 'ERA5']]
	
	#tests correlation index for each averaging timescale between 1 and 40 years
	smoothings = np.arange(1, 40, 1)
	plt.clf()
	for pair in pairs:
		r_squares = []
		paired_my_index, paired_other_index, paired_times = truncate_to_pairs(times, mean_SAM_indices, pair[1], pair[0])
		for smoothing in smoothings:
			smoothed_my_index, _ = running_mean(paired_my_index, paired_times, smoothing)
			smoothed_other_index, _ = running_mean(paired_other_index, paired_times, smoothing)
			res = linregress(smoothed_my_index, smoothed_other_index)
			r_squares.append(res.rvalue**2)
		plt.plot(smoothings, r_squares, label=pair[2]) #plots comparison with ERA and Marshall SAM on same graph
	plt.xlabel('years averaged')
	plt.ylabel('R squared value')
	plt.legend()
	plt.title('Correlation strength between ' + dataset + ' during ' + season + ' and ERA/Marshall')
	figure_name = Figure_dir + dataset + '_' + season + '_SAM_correlation_depending_on_averaging.png'
	plt.savefig(figure_name)
	#plt.show()


def separate_pressures(dataset='CSF-20C', season='DJF'):
	#plots pressure at each of 40S and 65S throughout 20th century
	#reads in netcdf file of relevant surface pressures
	file_name = Data_dir + dataset + '_' + season + '_msl_data.nc'
	mslp_data, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file_name, 'mslp for SAM', lat_name='latitude', lon_name='longitude', time_name='time')
	
	mean_pressures_40S, mean_pressures_65S, ensemble_pressures_40S, ensemble_pressures_65S = get_zm_pressures(dataset=dataset, season=season)
	
	#reads in ERA5 sea level pressure data and does same thing with it
	era5_file = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_195001-197812.nc'
	era_slp, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = rd_data.read_in_variable(era5_file, 'psl')
	era5_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/psl/mon/psl_mon_ERA5_2.5x2.5_197901-202012.nc'
	era_slp2, _, _, _, era_times2, _, _ = rd_data.read_in_variable(era5_file2, 'psl')
	era_slp = np.append(era_slp, era_slp2, axis=0)
	era_times = np.append(era_times, era_times2)
	
	era_pressures_40S = np.array([])
	era_pressures_65S = np.array([])
	for era_year_slp in era_slp:
		era_slp_65S = era_year_slp[10] #index 10 corresponds to latitude 65S
		era_slp_40S = era_year_slp[18] #index 18 corresponds to latitude 40S
		
		zm_40S = np.mean(era_slp_40S) #takes zonal mean of mslp at 40S
		zm_65S = np.mean(era_slp_65S) #takes zonal mean of mslp at 65S
		
		era_pressures_40S = np.append(era_pressures_40S, zm_40S)
		era_pressures_65S = np.append(era_pressures_65S, zm_65S)
	
	#eliminates irrelevant seasons
	era_pressures_40S, era_years = mask_era(era_pressures_40S, era_times, era_calendar, era_t_units, season)
	era_pressures_65S, era_years = mask_era(era_pressures_65S, era_times, era_calendar, era_t_units, season)
	
	#plots 40S pressure and 65S pressure on same axes
	things_to_plot = [[mean_pressures_40S, ensemble_pressures_40S, 'zonal mean 40S', 'blue'], [mean_pressures_65S, ensemble_pressures_65S, 'zonal mean 65S', 'red']]
	plt.clf()
	for row in things_to_plot:
		plt.plot(times, row[1], color='grey')
		plt.plot(times, row[0], label=row[2], color=row[3])
	plt.plot(era_years, era_pressures_40S, color='black')
	plt.plot(era_years, era_pressures_65S, color='black')
	plt.legend()
	plt.xlabel('years')
	plt.ylabel('mean surface level pressure')
	plt.title('Zonal Mean MSLP at 40S and 65S during ' + season + ' in ' + dataset)
	figure_name = Figure_dir + dataset + '_' + season + '_40S_and_65S_Pressures.png'
	plt.savefig(figure_name)
	#plt.show()


def stat_analysis(dataset='CSF-20C', season='DJF', shift_years = False):
	#runs many of the above plotting functions one after the other
	#reads in seasonal forecast SAM indices
	mean_SAM_indices, times, calendar, t_units = read_SAM_indices(dataset, season)
	if shift_years:
		times += 1
	
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
		correlate_pairs(paired_my_index, paired_other_index, dataset, pair[2], season=season, shift_years = shift_years)
		
		#repeats same process for different lengths of running means
		timescales = [2, 3, 5, 10, 15, 20, 30]
		for timescale in timescales:
			smooth_and_plot(paired_my_index, paired_other_index, dataset, pair[2], paired_times, timescale, season, shift_years = shift_years)
	
	#compares ERA5 data and Marshall data
	paired_era_index2, paired_Marshall_index2, Marshall_and_era_times = truncate_to_pairs(era_years, era_SAM_indices, Marshall_years, Marshall_SAM_index)
	correlate_pairs(paired_era_index2, paired_Marshall_index2, 'ERA5', 'Marshall', season=season, shift_years = shift_years)


def full_analysis(dataset='CSF-20C', season='DJF'):
	save_SAM_indices(dataset=dataset, season=season)
	graph_SAM_indices(dataset=dataset, season=season)
	#stat_analysis(dataset, season)
	graph_stdev(dataset=dataset, season=season)

#go to run_sam_analysis to run code because don't want code to execute when imported into multiseason_analaysis
