#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 2022

@author: wadh5699
"""

import os, sys
import numpy as np
from netCDF4 import num2date, Dataset
from scipy import stats
import matplotlib.pyplot as plt

Code_dir = '/home/w/wadh5699/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)

import reading_in_data_functions as rd_data
import save_file as sf
import calculate_csf_SAM as calc1
import multiseason_analysis as mult

def calc_t_test_correlation(corr, n_eff):
	#if passed correlation coefficient and effective number of datapoints, returns significance
	df = n_eff - 1
	test_statistic = np.sqrt(n_eff - 2) * corr / np.sqrt(1 - corr**2)
	pvals = stats.t.sf(np.abs(test_statistic), df) * 2
	return pvals #p < 0.05 generally taken as significant

def significance(my_SAM_indices, Marshall_SAM_indices, timescale):
	#calculates correlation coefficient and n_eff to pass to calc_t_test_correlation, and returns significance value
	res = stats.linregress(my_SAM_indices, Marshall_SAM_indices)
	autocorrelation = stats.linregress(my_SAM_indices[:len(my_SAM_indices) - 1], my_SAM_indices[1:])
	rho = autocorrelation.rvalue
	if timescale != 1:
		n_eff = len(my_SAM_indices) * (1 - rho) / (1 + rho)
	else:
		n_eff = len(my_SAM_indices)
	corr = res.rvalue
	return calc_t_test_correlation(corr, n_eff)

def plot_signif_with_averaging(dataset, season):
	#reads in SAM indices to be analysed
	mean_SAM_indices, times, calendar, t_units = calc1.read_SAM_indices(dataset, season)
	Marshall_SAM_index, Marshall_years = rd_data.read_Marshall_SAM_idx(season)
	
	plt.clf()
	paired_my_index, paired_Marshall_index, paired_times = calc1.truncate_to_pairs(times, mean_SAM_indices, Marshall_years, Marshall_SAM_index)
	timescales = np.arange(1, 10, 1)
	significances = []
	for timescale in timescales:
		smoothed_my_index, _ = calc1.running_mean(paired_my_index, paired_times, timescale)
		smoothed_Marshall_index, _ = calc1.running_mean(paired_Marshall_index, paired_times, timescale)
		significances.append(significance(smoothed_my_index, smoothed_Marshall_index, timescale))
	plt.plot(timescales, significances)
	plt.xlabel('number of years averaged together')
	plt.ylabel('p-value (p<.05 is significant)')
	plt.title('Significance of Correlation Between ' + dataset + ' and Marshall SAM Indices during ' + season)
	plt.savefig('Figures/' + dataset + '_Marshall_' + season + '_significance.png')
	plt.show()

plot_signif_with_averaging('CSF-20C', 'DJF')
plot_signif_with_averaging('SEAS5', 'DJF')
