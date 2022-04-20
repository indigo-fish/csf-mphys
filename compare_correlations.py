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

def get_corr_difference(dataset1, detrending1, dataset2, detrending2, period, season='DJF'):
	if period == 'early':
		start_year, end_year = 1958, 1986
	elif period == 'late':
		start_year, end_year = 1982, 2010
	all_SAMs = []
	all_years = []
	print(dataset1 + ' ' + str(detrending1) + ' ' + dataset2 + ' ' + str(detrending2) + ' ' + period)
	
	if detrending1:
		dataset_arr = [dataset1, 'Marshall']
		for dataset in dataset_arr:
			if period == 'early': access_str = dataset + '_detrended_1958-1986'
			elif period == 'late': access_str = dataset + '_detrended_1982-2010'
			SAM, years, _, _ = calc1.read_SAM_indices(access_str, season)
			#if dataset == 'CSF-20C' or dataset == 'ASF-20C': years += 1
			all_SAMs.append(SAM)
			all_years.append(years)
	else:
		SAM, years, _, _ = calc1.read_SAM_indices(dataset1, season)
		print(str(years[0]) + ' ' + str(years[len(years) - 1]))
		if dataset1 == 'CSF-20C' or dataset1 == 'ASF-20C':
			years += 1
			print('shifted years')
		all_SAMs.append(SAM)
		all_years.append(years)
		SAM, years = rd_data.read_Marshall_SAM_idx(season)
		all_SAMs.append(SAM)
		all_years.append(years)
	
	if detrending2:
		dataset_arr = [dataset2, 'Marshall']
		for dataset in dataset_arr:
			if period == 'early': access_str = dataset + '_detrended_1958-1986'
			elif period == 'late': access_str = dataset + '_detrended_1982-2010'
			SAM, years, _, _ = calc1.read_SAM_indices(access_str, season)
			#if dataset == 'CSF-20C' or dataset == 'ASF-20C': years += 1
			all_SAMs.append(SAM)
			all_years.append(years)
	else:
		SAM, years, _, _ = calc1.read_SAM_indices(dataset2, season)
		all_SAMs.append(SAM)
		all_years.append(years)
		SAM, years = rd_data.read_Marshall_SAM_idx(season)
		if dataset2 == 'CSF-20C' or dataset2 == 'ASF-20C':
			years += 1
			print('shifted years')
		all_SAMs.append(SAM)
		all_years.append(years)
	
	
	for i, (SAM_array, year_array) in enumerate(zip(all_SAMs, all_years)):
		mask = (year_array >= start_year) & (year_array <= end_year)
		print(mask)
		all_SAMs[i] = SAM_array[mask]
		all_years[i] = year_array[mask]
	
	for year_array in all_years:
		print(str(year_array[0]) + ' ' + str(year_array[len(year_array) - 1]))
	
	ay = linregress(all_SAMs[0], all_SAMs[1])
	r_ay = ay.rvalue
	by = linregress(all_SAMs[2], all_SAMs[3])
	r_by = by.rvalue
	ab = linregress(all_SAMs[0], all_SAMs[2])
	r_ab = ab.rvalue
	
	c_ab = ((r_ab - .5 * r_ay * r_by) * (1 - r_ay ** 2 - r_by ** 2 - r_ab ** 2) + r_ab ** 3) / (1 - r_ay ** 2) / (1 - r_by ** 2)
	
	n = len(all_SAMs[0])
	
	z_ay = np.log((1+r_ay) / (1-r_ay)) / 2
	l1_a = z_ay - (1.96 / np.sqrt(n - 3))
	u1_a = z_ay + (1.96 / np.sqrt(n - 3))
	L_a = (np.exp(2 * l1_a) - 1) / (np.exp(2 * l1_a) + 1)
	U_a = (np.exp(2 * u1_a) - 1) / (np.exp(2 * u1_a) + 1)
	
	z_by = np.log((1+r_by) / (1-r_by)) / 2
	l1_b = z_by - (1.96 / np.sqrt(n - 3))
	u1_b = z_by + (1.96 / np.sqrt(n - 3))
	L_b = (np.exp(2 * l1_b) - 1) / (np.exp(2 * l1_b) + 1)
	U_b = (np.exp(2 * u1_b) - 1) / (np.exp(2 * u1_b) + 1)
	
	det_l = (r_by - L_b) ** 2 + (U_a - r_ay) ** 2 - 2 * c_ab * (r_by - L_b) * (U_a - r_ay)
	det_u = (U_b - r_by) ** 2 + (r_ay - L_a) ** 2 - 2 * c_ab * (U_b - r_by) * (r_ay - L_a)
	L = (r_by - r_ay) - np.sqrt(det_l)
	U = (r_by - r_ay) + np.sqrt(det_u)
	print_string = dataset1
	if detrending1: print_string += ' detrended'
	print_string += '; ' + dataset2
	if detrending2: print_string += ' detrended'
	print(print_string)
	print('[' + str(L) + ', ' + str(U) + ']')

datasets = ['CSF-20C', 'ASF-20C']
periods = ['early', 'late']
for dataset in datasets:
	for period in periods:
		get_corr_difference(dataset, False, dataset, True, period)

get_corr_difference('SEAS5', False, 'SEAS5', True, 'late')

datasets = ['CSF-20C', 'ASF-20C', 'SEAS5']
for d1 in datasets:
	for d2 in datasets:
		if d1 != d2: get_corr_difference(d1, False, d2, False, 'late')
