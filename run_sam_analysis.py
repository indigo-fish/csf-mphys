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
import multiseason_analysis as mult
import detrending as det

"""
datasets = ['CSF-20C', 'ASF-20C', 'Marshall']
seasons = ['DJF', 'MAM', 'JJA', 'SON']
for dataset in datasets:
	for season in seasons:
		det.detrend_data(dataset=dataset, season=season)
		det.detrend_data(dataset=dataset, season=season, compare_SEAS5 = False)

det.detrend_data(dataset='SEAS5', season='DJF')
"""

datasets = ['CSF-20C', 'ASF-20C', 'SEAS5']
seasons = ['DJF']
for dataset in datasets:
	for season in seasons:
		"""
		det.corr_without_trend(dataset=dataset, season=season)
		det.corr_without_trend(dataset=dataset, season=season, compare_SEAS5 = False)
		det.corr_with_trend(dataset=dataset, season=season)
		det.corr_with_trend(dataset=dataset, season=season, compare_SEAS5 = False)
		"""
		calc1.graph_SAM_indices(dataset=dataset, season=season, cut_years=False, trend=False)
		calc1.graph_SAM_indices(dataset=dataset, cut_years=True, trend=True)

#det.corr_without_trend(dataset='SEAS5', season='DJF')
#det.corr_with_trend(dataset='SEAS5', season='DJF')

#calc1.stat_analysis(dataset='SEAS5', season='DJF', shift_years = False)
#calc1.stat_analysis(dataset='SEAS5', season='DJF', shift_years = True)

"""
datasets = ['CSF-20C', 'ASF-20C']
seasons = ['DJF', 'MAM', 'JJA', 'SON']
for dataset in datasets:
	for season in seasons:
		calc1.stat_analysis(dataset=dataset, season=season, shift_years = False)
		calc1.stat_analysis(dataset=dataset, season=season, shift_years = True)
"""
"""
file_name = 'Data/CSF-20C_DJF_msl_data.nc'
data, lats, lons, levs, times, calendar, t_units = rd_data.read_in_variable(file_name, 'mslp for SAM', lat_name='latitude', lon_name='longitude', time_name='time')
#calc1.save_seas5_SAM_indices('DJF')
#calc1.graph_SAM_indices('SEAS5', 'DJF', variance=False)
calc1.stat_analysis('SEAS5', 'DJF')
"""
