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

import reading_in_data_functions as rd_data
import save_file as sf
import calculate_csf_SAM as calc1
import significance_testing as sig_test

Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
Data_dir = Code_dir + 'Data/'
Figure_dir = Code_dir + 'Figures/'

def detrend_data(dataset='CSF-20C', season='DJF', compare_SEAS5 = True, full_period = False, choose_Marshall=0):
	#produces detrended datasets for use in various analyses
	
	#reads in SAM index data
	if dataset == 'Marshall':
		mean_SAM_indices, times = rd_data.read_Marshall_SAM_idx(season)
		_, _, calendar, t_units = calc1.read_SAM_indices('CSF-20C', season)
	else:
		mean_SAM_indices, times, calendar, t_units = calc1.read_SAM_indices(dataset, season)
	if (dataset == 'CSF-20C' or dataset == 'ASF-20C') and season=='DJF':
		times += 1
	
	#truncates data series to desired common period
	if full_period:
		if dataset=='CSF-20C' and season=='DJF': start_year, end_year = 1958, 2011
		elif dataset=='CSF-20C': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C' and season=='DJF': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C': start_year, end_year = 1958, 2009
		elif dataset=='SEAS5': start_year, end_year = 1982, 2017
		elif dataset=='Marshall':
			if choose_Marshall==1: start_year, end_year = 1958, 2011
			elif choose_Marshall==2: start_year, end_year = 1958, 2010
			elif choose_Marshall==3: start_year, end_year = 1982, 2017
	else:
		if compare_SEAS5 and season=='DJF':
			start_year, end_year = 1982, 2010 #common period is 1982 (SEAS5) to December 2009 (ASF-20C)
		elif compare_SEAS5:
			start_year, end_year = 1982, 2009
		else:
			start_year, end_year = 1958, 1986
	year_mask = (times >= start_year) & (times <= end_year)
	masked_SAM = mean_SAM_indices[year_mask]
	masked_times = times[year_mask]
	
	#calculates trend
	csf_res = linregress(masked_times, masked_SAM)
	
	#constructs new datasets with trend subtracted
	mid_year = (start_year + end_year) / 2
	detrended_SAM = []
	for SAM, year in zip(masked_SAM, masked_times):
		delta_year = year - mid_year
		detrended_SAM.append(SAM - delta_year * csf_res.slope)
	
	#saves netcdf files
	if compare_SEAS5:
		dataset_destination = Data_dir + dataset + '_detrended_1982-2010_' + season + '_sam_mean_data.nc'
	else:
		dataset_destination = Data_dir + dataset + '_detrended_1958-1986_' + season + '_sam_mean_data.nc'
	if full_period:
		dataset_destination = Data_dir + dataset + '_detrended_' + str(start_year) + '-' + str(end_year) + '_' + season + '_sam_mean_data.nc'
	dataset_description = 'detrended Marshall SAM index from ' + dataset + ' during ' + season
	save = sf.save_file(dataset_destination, dataset_description)
	save.add_times(masked_times, calendar, t_units, time_name = 'time')
	save.add_variable(np.array(detrended_SAM), 'SAM index', ('time'))
	save.close_file()


def corr_without_trend(axes, dataset='CSF-20C', season='DJF', compare_SEAS5 = True, full_period=False):
	#produces plots of interannual variability with SAM subtracted
	if compare_SEAS5:
		years = '1982-2010'
	else:
		years = '1958-1986'
	if full_period:
		if dataset=='CSF-20C' and season=='DJF': start_year, end_year = 1958, 2011
		elif dataset=='CSF-20C': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C' and season=='DJF': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C': start_year, end_year = 1958, 2009
		elif dataset=='SEAS5': start_year, end_year = 1982, 2017
		years = str(start_year) + '-' + str(end_year)
	detrended_SAM, _, _, _ = calc1.read_SAM_indices(dataset + '_detrended_' + years, season)
	detrended_Marshall, _, _, _ = calc1.read_SAM_indices('Marshall' + '_detrended_' + years, season)
	
	#plots interannual variability
	axes.plot(detrended_SAM, detrended_Marshall, 'o', label='original data')
	
	res = linregress(detrended_SAM, detrended_Marshall)
	xpts = np.arange(-2.2, 2.2, 0.1)
	lineplot = []
	for point in xpts:
		lineplot.append(point * res.slope + res.intercept)
	
	axes.plot(xpts, lineplot, 'r', label='fitted line')
	axes.set(xlabel=dataset + ' SAM', ylabel='Marshall SAM')
	significance = sig_test.significance(detrended_SAM, detrended_Marshall, 1)
	axes.annotate(f"p-value: {significance:.6f}", (0.1, 0.85), xycoords='axes fraction')
	axes.annotate(f"r-value: {res.rvalue:.6f}", (0.1, 0.65), xycoords='axes fraction')


def corr_with_trend(axes, dataset='CSF-20C', season='DJF', compare_SEAS5 = True, full_period=False):
	#produces plots of interannual variability with SAM subtracted
	if compare_SEAS5:
		years = '1982-2010'
	else:
		years = '1958-1986'
	if compare_SEAS5 and season=='DJF':
		start_year, end_year = 1982, 2010 #common period is 1982 (SEAS5) to December 2009 (ASF-20C)
	elif compare_SEAS5:
		start_year, end_year = 1982, 2009
	else:
		start_year, end_year = 1958, 1986
	
	if full_period:
		if dataset=='CSF-20C' and season=='DJF': start_year, end_year = 1958, 2011
		elif dataset=='CSF-20C': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C' and season=='DJF': start_year, end_year = 1958, 2010
		elif dataset=='ASF-20C': start_year, end_year = 1958, 2009
		elif dataset=='SEAS5': start_year, end_year = 1982, 2017
		years = str(start_year) + '-' + str(end_year)
	mean_SAM_indices, times, _, _ = calc1.read_SAM_indices(dataset, season)
	if (dataset=='CSF-20C' or dataset == 'ASF-20C') and season=='DJF':
		times += 1
	Marshall_SAM_indices, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	
	year_mask = (times >= start_year) & (times <= end_year)
	masked_SAM = mean_SAM_indices[year_mask]
	masked_times = times[year_mask]
	year_mask_Marshall = (Marshall_years >= start_year) & (Marshall_years <= end_year)
	masked_Marshall = Marshall_SAM_indices[year_mask_Marshall]
	masked_Marshall_times = Marshall_years[year_mask_Marshall]
	
	#plots interannual variability
	axes.plot(masked_SAM, masked_Marshall, 'o', label='original data')
	
	res = linregress(masked_SAM, masked_Marshall)
	xpts = np.arange(-2.2, 2.2, 0.1)
	lineplot = []
	for point in xpts:
		lineplot.append(point * res.slope + res.intercept)
	
	axes.plot(xpts, lineplot, 'r', label='fitted line')
	axes.set(xlabel=dataset + ' SAM', ylabel='Marshall SAM')
	significance = sig_test.significance(masked_SAM, masked_Marshall, 1)
	axes.annotate(f"p-value: {significance:.6f}", (0.1, 0.85), xycoords='axes fraction')
	axes.annotate(f"r-value: {res.rvalue:.6f}", (0.1, 0.65), xycoords='axes fraction')


def graph_all(season='DJF', compare_SEAS5=True, detrended=False, full_period=False):
	if season=='DJF' and compare_SEAS5: datasets=['CSF-20C', 'ASF-20C', 'SEAS5']
	elif season=='DJF' and full_period: datasets=['CSF-20C', 'ASF-20C', 'SEAS5']
	else: datasets=['CSF-20C', 'ASF-20C']
	
	if compare_SEAS5: year_range = '1982-2010'
	else: year_range = '1958-1986'
	
	fig, axes = plt.subplots(len(datasets), sharex=True, sharey=True)
	for dataset, axis in zip (datasets, axes):
		if detrended: corr_without_trend(axis, dataset=dataset, season=season, compare_SEAS5=compare_SEAS5, full_period=full_period)
		else: corr_with_trend(axis, dataset=dataset, season=season, compare_SEAS5=compare_SEAS5, full_period=full_period)
	if full_period:
		if detrended: title = 'Detrended SAM Skill in ' + season + ': Full Period'
		else: title = 'SAM Skill in ' + season + ': Full Period'
	else:
		if detrended: title = 'Detrended SAM Skill in ' + season + ': ' + year_range
		else: title = 'SAM Skill in ' + season + ': ' + year_range
	fig.suptitle(title)
	
	figure_name = Figure_dir + 'Final_multiplot_SAM_correlations_'
	if detrended: figure_name += 'detrended_'
	else: figure_name += 'with_trend_'
	if full_period: figure_name += season + '_full_periods.png'
	else: figure_name += season + '_' + year_range + '.png'
	print('saving figure to ' + figure_name)
	fig.savefig(figure_name)

def confidence_interval(period, r):
	z_r = np.log((1+r)/(1-r))/2
	L = z_r - (1.96/np.sqrt(period - 3))
	U = z_r + (1.96/np.sqrt(period - 3))
	lower_bound = (np.exp(2*L) - 1) / (np.exp(2*L) + 1)
	upper_bound = (np.exp(2*U) - 1) / (np.exp(2*U) + 1)
	return lower_bound, upper_bound

def diff_period_correlations(axes, period, dataset='CSF-20C', season='DJF'):
	#produces plots of interannual variability with SAM subtracted
	if dataset=='SEAS5': start_years = np.arange(1982, 2017+1-period, 1)
	elif dataset=='ERA5': start_years = np.arange(1958, 2020+1-period, 1)
	elif dataset=='CSF-20C' and season=='DJF': start_years = np.arange(1958, 2011+1-period, 1)
	else: start_years = np.arange(1958, 2010+1-period, 1)
	
	if dataset=='ERA5': mean_SAM_indices, times = calc1.get_era_SAM_indices(season)
	else: mean_SAM_indices, times, _, _ = calc1.read_SAM_indices(dataset, season)
	if (dataset=='CSF-20C' or dataset == 'ASF-20C') and season=='DJF':
		times += 1
	Marshall_SAM_indices, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	
	r_values = []
	lower_bounds = []
	upper_bounds = []
	#p_values = []
	for start_year in start_years:
		end_year = start_year + period
		year_mask = (times >= start_year) & (times <= end_year)
		masked_SAM = mean_SAM_indices[year_mask]
		masked_times = times[year_mask]
		year_mask_Marshall = (Marshall_years >= start_year) & (Marshall_years <= end_year)
		masked_Marshall = Marshall_SAM_indices[year_mask_Marshall]
		masked_Marshall_times = Marshall_years[year_mask_Marshall]
		res = linregress(masked_SAM, masked_Marshall)
		r_values.append(res.rvalue)
		lower_bound, upper_bound = confidence_interval(period, res.rvalue)
		lower_bounds.append(lower_bound)
		upper_bounds.append(upper_bound)
		#p_values.append(sig_test.significance(masked_SAM, masked_Marshall, 1))
	
	axes.plot(start_years + int(period / 2), r_values, label=dataset)
	axes.fill_between(start_years + int(period / 2), lower_bounds, upper_bounds, alpha=.1)
	#plt.plot(start_years + int(period / 2), p_values)
	#plt.show()

def run_multi_correlations(period, season='DJF'):
	if season=='DJF': datasets=['CSF-20C', 'ASF-20C', 'SEAS5']
	else: datasets=['CSF-20C', 'ASF-20C']
	fig, axes = plt.subplots()
	
	for dataset in datasets:
		diff_period_correlations(axes, period, dataset=dataset, season=season)
	fig.suptitle('Correlation During ' + str(period) + ' year periods')
	axes.set(xlabel='Central Year',ylabel='r-value')
	axes.legend()
	figure_name = Figure_dir + 'Interval_period_correlations_' + season + '_' + str(period) + '_years.png'
	print('saving figure to ' + figure_name)
	fig.savefig(figure_name)

def diff_period_trends(axes, period, dataset='CSF-20C', season='DJF'):
	#produces plots of interannual variability with SAM subtracted
	if dataset=='SEAS5': start_years = np.arange(1982, 2017+1-period, 1)
	elif dataset=='Marshall': start_years = np.arange(1958, 2020+1-period, 1)
	elif dataset=='ERA5': start_years = np.arange(1950, 2020+1-period, 1)
	elif dataset=='CSF-20C' and season=='DJF': start_years = np.arange(1958, 2011+1-period, 1)
	else: start_years = np.arange(1958, 2010+1-period, 1)
	
	if dataset=='Marshall': mean_SAM_indices, times = rd_data.read_Marshall_SAM_idx(season)
	elif dataset=='ERA5': mean_SAM_indices, times = calc1.get_era_SAM_indices(season)
	else: mean_SAM_indices, times, _, _ = calc1.read_SAM_indices(dataset, season)
	if (dataset=='CSF-20C' or dataset == 'ASF-20C') and season=='DJF': times += 1
	
	if dataset=='Marshall': mean_SAM_indices /= np.std(mean_SAM_indices)
	trends = []
	lower_bounds = []
	upper_bounds = []
	#p_values = []
	for start_year in start_years:
		end_year = start_year + period
		year_mask = (times >= start_year) & (times <= end_year)
		masked_SAM = mean_SAM_indices[year_mask]
		masked_times = times[year_mask]
		res = linregress(masked_times, masked_SAM)
		trends.append(res.slope)
		lower_bounds.append(res.slope - 1.96 * res.stderr)
		upper_bounds.append(res.slope + 1.96 * res.stderr)
		#p_values.append(sig_test.significance(masked_SAM, masked_Marshall, 1))
	
	axes.plot(start_years + int(period / 2), trends, label=dataset)
	axes.fill_between(start_years + int(period / 2), lower_bounds, upper_bounds, alpha=.1)
	#plt.plot(start_years, p_values)
	#plt.show()

def run_multi_trends(period, season='DJF'):
	if season=='DJF': datasets=['Marshall', 'CSF-20C', 'ASF-20C', 'SEAS5']
	else: datasets=['Marshall', 'CSF-20C', 'ASF-20C']
	fig, axes = plt.subplots()
	
	for dataset in datasets:
		diff_period_trends(axes, period, dataset=dataset, season=season)
	fig.suptitle('Trend During ' + str(period) + ' year periods')
	axes.set(xlabel='Central Year',ylabel='$\Delta$ SAM/year')
	axes.legend()
	figure_name = Figure_dir + 'Final_period_trends_' + season + '_' + str(period) + '_years.png'
	print('saving figure to ' + figure_name)
	fig.savefig(figure_name)


array = [True, False]
for i in array:
	graph_all(season='DJF', compare_SEAS5=False, detrended=i, full_period=True)

"""
intervals = np.arange(10, 35, 5)
intervals = [30]
for period in intervals:
	#run_multi_correlations(period, season='DJF')
	run_multi_trends(period, season='DJF')
"""
"""
run_multi_trends(20, season='DJF')
run_multi_trends(15, season='DJF')
run_multi_trends(10, season='DJF')
run_multi_trends(25, season='DJF')
run_multi_trends(30, season='DJF')
"""
"""
datasets = ['CSF-20C', 'ASF-20C', 'SEAS5']
for dataset in datasets:
	detrend_data(dataset=dataset, season='DJF', compare_SEAS5=False, full_period=True)

Marshall_choices = [1, 2, 3]
for choice in Marshall_choices:
	detrend_data(dataset='Marshall', season='DJF', compare_SEAS5=False, full_period=True, choose_Marshall=choice)
"""
