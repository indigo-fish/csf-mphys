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

def persistence(axes, dataset='CSF-20C', season='DJF', compare_SEAS5 = True):
	#produces plots of interannual variability with SAM subtracted
	if compare_SEAS5:
		years = '1982-2010'
	else:
		years = '1958-1986'
	mean_SAM_indices, times, _, _ = calc1.read_SAM_indices(dataset, season)
	if (dataset=='CSF-20C' or dataset == 'ASF-20C') and season=='DJF': times += 1
	#uses Marshall SAM index from SON, which is hopefully close to ERA in November
	Marshall_SAM_indices, Marshall_years = calc1.get_era_SAM_indices(season)
	Marshall_years += 1
	
	if compare_SEAS5 and season=='DJF':
		start_year, end_year = 1982, 2010 #common period is 1982 (SEAS5) to December 2009 (ASF-20C)
	elif compare_SEAS5:
		start_year, end_year = 1982, 2009
	else:
		start_year, end_year = 1958, 1986
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


def graph_all(season='DJF', compare_SEAS5=True):
	if season=='DJF' and compare_SEAS5: datasets=['CSF-20C', 'ASF-20C', 'SEAS5']
	else: datasets=['CSF-20C', 'ASF-20C']
	
	if compare_SEAS5: year_range = '1982-2010'
	else: year_range = '1958-1986'
	
	fig, axes = plt.subplots(len(datasets), sharex=True, sharey=True)
	for axis, dataset in zip(axes, datasets):
		persistence(axis, dataset=dataset, season=season, compare_SEAS5=compare_SEAS5)
	title = 'SAM Persistence in ' + season + ': ' + year_range
	fig.suptitle(title)
	
	figure_name = Figure_dir + 'Final_multiplot_SAM_persistence_'
	figure_name += season + '_' + year_range + '.png'
	print('saving figure to ' + figure_name)
	fig.savefig(figure_name)


array = [True, False]
for i in array:
	graph_all(season='DJF', compare_SEAS5=i)
