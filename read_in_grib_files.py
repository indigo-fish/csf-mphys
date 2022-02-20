#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 12:12:05 2021

@author: pattersonm
"""

import xarray as xr
import numpy as np

def read_in_grib(file_name,variable_name,lead=1):
    """ Read in grib file
    
    Args:
        file_name = name of file
        variable_name = name of variable in file e.g. 't2m'
        lead = number of months ahead to return, e.g. for runs initialised in May, lead=1 returns 
    the mean of June, July and August.
    
    Returns:
        Mean of three months
        latitudes,longitudes,levels
    """
    ds = xr.open_dataset(file_name, engine="cfgrib")
    data = ds[variable_name].data
    t = np.arange(data.shape[0])
    
    if lead==0: mask = (t==0)|(t==1)|(t==2)
    elif lead==1: mask = (t==1)|(t==2)|(t==3)
    elif lead==2: mask = (t==2)|(t==3)|(t==4)
    elif lead==3: mask = (t==3)|(t==4)|(t==5)

    lats = ds.variables['latitude'].data
    lons = ds.variables['longitude'].data
    try: levs = ds.variables['isobaricInhPa'].data
    except: levs = None
    
    return np.mean(data[mask],axis=0),lats,lons,levs

"""
# example to read in mean sea level pressure from CSF-20C for December, january February following initialisation in November 1981
file_name = '/network/group/aopp/predict/AWH002_BEFORT_SEASONAL/CSF-20C/ecmf/guh4/sfc/fcmean/19811101/var151/ecmf-guh4_0_19811101_fcmean_sfc.grb'

msl, lats, lons, levs = read_in_grib(file_name,'msl',lead=1)

# plot this quickly
import matplotlib.pyplot as plt
cs = plt.contourf(lons,lats,msl,cmap='RdBu_r')
plt.colorbar(cs)
plt.show()
figure_name = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/mslp_from_grib_19811101'
plt.savefig(figure_name)
"""

"""
#calculates mean zonal surface level pressure at 40S and 65S
years = np.arange(1981, 2001, 1)
ensemble = np.arange(0, 25, 1)
ens_SAM_indices = []
for member in ensemble:
	SAM_indices = []
	for year in years:
		file_name = '/network/group/aopp/predict/AWH002_BEFORT_SEASONAL/CSF-20C/ecmf/guh4/sfc/fcmean/' + str(year) + '1101/var151/ecmf-guh4_' + str(member) + '_' + str(year) + '1101_fcmean_sfc.grb'
		msl, lats, lons, levs = read_in_grib(file_name,'msl',lead=1)
		
		msl_40S = msl[192,:] #index 192 in latitudes corresponds to 40S
		msl_65S = msl[220,:] #index 220 in latitudes corresponds to 65S
		
		zm_40S = sum(msl_40S) / len(msl_40S) #zonal mean of 40S mean surface level pressure
		zm_65S = sum(msl_65S) / len(msl_65S) #zonal mean of 65S mean surface level pressure
		norm_40S = 10000 #I don't know what the normalization factor for 40S should be
		norm_65S = 10000 #I don't know what the normalization factor for 65S should be
		SAM_index = zm_40S / norm_40S - zm_65S / norm_65S #SAM index defined as difference of normalized zonal mean surface pressures
		SAM_indices.append(SAM_index)
	ens_SAM_indices.append(SAM_indices)

mean_SAM_indices = []
indices = np.arange(0, len(years), 1)
for index in indices:
	all_SAM_indices = []
	for member in ens_SAM_indices:
		all_SAM_indices.append(member[index])
	mean = sum(all_SAM_indices) / len(all_SAM_indices)
	mean_SAM_indices.append(mean)

import matplotlib.pyplot as plt

#for SAM_indices in ens_SAM_indices:
#	plt.plot(years, SAM_indices)
plt.plot(years, mean_SAM_indices)
plt.show()
figure_name = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/Figures/my_mean_SAM_indices'
plt.savefig(figure_name)
"""
