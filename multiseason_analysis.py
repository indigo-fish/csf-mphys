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

def graph_all_SAM_indices(dataset='CSF-20C'):
	
	seasons = ['DJF', 'MAM', 'JJA', 'SON']
	
	SAM_ensemble_seasonal = []
	SAM_mean_seasonal = []
	Marshall_seasonal = []
	ERA_seasonal = []
	for season in seasons:
		yearly_SAM_indices, mean_SAM_indices, SAM_stdevs, times, calendar, t_units = calc1.get_SAM_indices(dataset=dataset, season=season)
		SAM_ensemble_seasonal.append(yearly_SAM_indices)
		SAM_mean_seasonal.append(mean_SAM_indices)
		Marshall_SAM_index, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
		Marshall_seasonal.append(Marshall_SAM_index)
		ERA_SAM_index, ERA_years = calc1.get_era_SAM_indices(season)
		ERA_seasonal.append(ERA_SAM_index)
	
	#displays plots of ensemble and mean SAM indices
	plt.figure(1)
	for seasonal_ensemble, seasonal_mean, season in zip(SAM_ensemble_seasonal, SAM_mean_seasonal, seasons):
		plt.plot(times, seasonal_ensemble, color='gray', alpha=.1)
		plt.plot(times, seasonal_mean, linewidth=2, label=season)
	#plt.errorbar(times, mean_SAM_indices, yerr=SAM_stdevs, color='black', label='standard deviation')
	plt.title('Normalized SAM Index in ' + dataset + ' Ensemble')
	plt.xlabel('Year')
	plt.ylabel('Normalized SAM Index')
	plt.legend()
	figure_name = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_All_Seasons_Normalized_SAM.png'
	plt.savefig(figure_name)
	plt.show()

def short_correlations(dataset='CSF-20C', season='DJF', timescale = 20):
	#does correlation analysis on 20 year subsets of longer datasets
	
	#reads in the 3 SAM index datasets as usual
	mean_SAM_indices, times, calendar, t_units = calc1.read_SAM_indices(dataset, season)
	Marshall_SAM_indices, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	era_SAM_indices, era_years = calc1.get_era_SAM_indices(season)
	
	pairs = [[Marshall_SAM_indices, Marshall_years, 'Marshall'], [era_SAM_indices, era_years, 'ERA5']]
	for pair in pairs:
		paired_my_index, paired_other_index, paired_times = calc1.truncate_to_pairs(times, mean_SAM_indices, pair[1], pair[0])
		start_index = 0
		plt.clf()
		while start_index < len(paired_times) - timescale:
			x = paired_my_index[start_index:start_index+timescale]
			y = paired_other_index[start_index:start_index+timescale]
			label = str(paired_times[start_index]) + '-' + str(paired_times[start_index + timescale])
			plt.scatter(x, y, label=label)
			res = linregress(x, y)
			line_y = []
			for x_pt in x:
				line_y.append(res.slope * x_pt + res.intercept)
			plt.plot(x, line_y)
			start_index += 10
		plt.xlabel(dataset + ' SAM index')
		plt.ylabel(pair[2] + ' SAM index')
		plt.title('Correlation Between ' + season + ' ' + dataset + ' and ' + pair[2] + ' In Different ' + str(timescale) + '-year Periods')
		plt.legend()
		figure_name = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_' + season + '_' + pair[2] + '_'  + str(timescale) + 'Period_Correlations.png'
		plt.savefig(figure_name)
		plt.show()

#use run_sam_analysis to run code
