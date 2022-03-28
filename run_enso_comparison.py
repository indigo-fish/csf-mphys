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
import calculate_csf_SAM as cal_sam
import multiseason_analysis as mult
import detrending as det
import calculate_el_nino as cal_nino

def compare_sam_enso(dataset='CSF-20C', season='DJF', compare_SEAS5 = True, full_period = False):
	sam, sam_times, _, _ = cal_sam.read_SAM_indices(dataset=dataset, season=season)
	nino, nino_times, _, _ = cal_nino.read_nino34(dataset=dataset, season=season)
	
	
	sam, nino, times = cal_sam.truncate_to_pairs(sam_times, sam, nino_times, nino)
	
