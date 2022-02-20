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

datasets = ['CSF-20C', 'ASF-20C']
seasons = ['DJF', 'MAM', 'JJA', 'SON']
for dataset in datasets:
	for season in seasons:
		mult.short_correlations(dataset=dataset, season=season)
