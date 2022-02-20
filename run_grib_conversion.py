#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 2021

@author: wadh5699
"""

import os, sys
import numpy as np

Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'

sys.path.append(Code_dir)
import read_in_grib_files as rgrb
import save_file as svf
import msl_grib_to_netcdf as grb2nc



#use run_grib_conversion() to run so no problems with importing
datasets = ['CSF-20C', 'ASF-20C']
seasons = ['DJF', 'MAM', 'JJA', 'SON']
for dataset in datasets:
	for season in seasons:
		grb2nc.convert_ssts(dataset=dataset, season=season)
