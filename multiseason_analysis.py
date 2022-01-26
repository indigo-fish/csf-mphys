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
	figure_name = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/' + dataset + '_All_Seasons_Normalized_SAM.png'
	plt.savefig(figure_name)
	plt.show()

#runs code
datasets = ['CSF-20C', 'ASF-20C']
for dataset in datasets:
	graph_all_SAM_indices(dataset)

