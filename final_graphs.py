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
from calculate_csf_SAM import *

def graph_one(axes, dataset='CSF-20C', season='DJF', variance=True, trend=True, cut_years=False):
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
	if variance: axes.plot(times, yearly_SAM_indices, color='gray')
	axes.plot(times, mean_SAM_indices, linewidth = 2, color = 'black', label = dataset)
	axes.plot(years_SAM, Marshall_SAM_index, linewidth = 2, color = 'red', label = 'Marshall')
	axes.plot(era_years, era_SAM_indices, linewidth = 2, color= 'blue', label = 'ERA5')
	
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
		axes.plot(times, csf_line, color='black', linestyle='dashed')
		axes.plot(years_SAM, Marshall_line, color='red', linestyle='dashed')
		axes.plot(era_years, era_line, color='blue', linestyle='dashed')
	
	axes.set(xlabel='Year', ylabel='SAM')
	axes.label_outer()
	axes.legend(loc='center left')
	
	figure_name = Figure_dir + 'Final_' + dataset + '_' + season
	if trend: figure_name += '_trend'
	if cut_years: figure_name += '_cut_years'
	#plt.show()

def graph_all(season='DJF', variance=True, trend=True, cut_years=False):
	if season=='DJF': datasets=['CSF-20C', 'ASF-20C', 'SEAS5']
	else: datasets=['CSF-20C', 'ASF-20C']
	fig, axes = plt.subplots(len(datasets), sharex=True, sharey=True)
	for dataset, axis in zip(datasets, axes):
		graph_one(axis, dataset=dataset, season=season, variance=variance, trend=trend, cut_years=cut_years)
	title = 'SAM Indices in ' + season
	fig.suptitle(title)
	
	figure_name = Figure_dir + 'Final_Multiseries_' + season
	if trend: figure_name += '_trend'
	if cut_years: figure_name += '_cut_years'
	figure_name += '_Normalized_SAM.png'
	print('saving figure to ' + figure_name)
	fig.savefig(figure_name)
	#fig.show()

graph_all(season='DJF', variance=True, trend=False, cut_years=False)
graph_all(season='DJF', variance=False, trend=True, cut_years=True)
