""" Script to read in the Marshall SAM index and calculate the associated geopotential height pattern
via linear regression analysis.

Plots this using two different maps (usually only one is necessary)"""

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import cartopy.crs as ccrs
import os, sys

Code_dir = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
start_year,end_year = 1958,2010 # choose the period over which to calculate regression pattern

sys.path.append(Code_dir)
import plotting_functions
import reading_in_data_functions
import analysis_functions
import calculate_csf_SAM

Code_dir = '/home/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
Figure_dir = Code_dir + 'Figures/' # directory to save figures to 
Data_dir = Code_dir + 'Data/' #directory to access data in

def sst_regression(season, dataset='CSF-20C'):
    variable = 'sst'
    
    #reads in previously calculated netcdf of SAM indices from given dataset
    ensemble_SAM_index, my_times, calendar, units = calculate_csf_SAM.read_SAM_ensemble(dataset=dataset, season=season)
    # reads in Marshall SAM index
    Marshall_SAM_idx, Marshall_years = reading_in_data_functions.read_Marshall_SAM_idx(season=season)
    
    # Read in seasonal sea surface temperature data (previously converted from grib files)
    file_str = Data_dir + dataset + '_' + season + '_sst_data.nc'
    var, lats, lons, levs, times, calendar, t_units = reading_in_data_functions.read_in_variable(file_str, 'sst')
    
    # calculate the seasonal mean SST
    # needed for ERA5 data, not for CSF/ASF data
    
    # restrict geopotential height and SAM index to chosen years
    year_mask_var = (times>=start_year)&(times<=end_year)
    var_am = var[year_mask_var,:,:]
    year_mask_Marshall = (Marshall_years>=start_year)&(Marshall_years<=end_year)
    Marshall_SAM_idx = Marshall_SAM_idx[year_mask_Marshall]
    my_year_mask = (my_times>=start_year)&(my_times<=end_year)
    ensemble_SAM_index = ensemble_SAM_index[my_year_mask]
    mean_sst = []
    mean_sam = []
    for (sst_yr, sam_yr) in zip(var_am, ensemble_SAM_index):
        mean_sst.append(np.mean(sst_yr))
        mean_sam.append(np.mean(sam_yr))
    
    plt.plot(mean_sst, mean_sam)
    print(mean_sst)
    print(mean_sam)
    
    """
    #regression map is expecting 2 time series, so pretend ensembles are time series
    unwrapped_sst = []
    unwrapped_sam = []
    for (sst_yr, sam_yr) in zip(var_am, ensemble_SAM_index):
        for (sst_ensemble, sam_ensemble) in zip(sst_yr, sam_yr):
            unwrapped_sst.append(sst_ensemble)
            unwrapped_sam.append(sam_ensemble)
    """
    """
    unwrapped_sst = []
    unwrapped_sam = []
    for (sst_yr, sam_yr) in zip(var_am, ensemble_SAM_index):
        unwrapped_sst.append(sst_yr[0])
        unwrapped_sam.append(sam_yr[0])
    unwrapped_sst = np.array(unwrapped_sst)
    unwrapped_sam = np.array(unwrapped_sam)
    
    # calculate regression pattern
    #regress_coeff,corr,pvals = analysis_functions.regress_map(my_SAM_index,var_am)
    regress_coeff,corr,pvals = analysis_functions.regress_map(unwrapped_sam, unwrapped_sst)
    print(corr)
    
    # make some plots
    plt.figure(figsize=(15,15))
    gs = gridspec.GridSpec(3,1,height_ratios=[5,10,0.5])
    
    clevs = np.linspace(np.min(regress_coeff), np.max(regress_coeff), 10)
    #clevs = np.append(-a[::-1],a) # contour level
    subplot_labels = ['a)','b)']
    
    for i, projection in enumerate([ccrs.PlateCarree(central_longitude=0.),ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0) ]):
    
        ax = plt.subplot(gs[i,0],projection=projection)
        cs = plotting_functions.plot_filled_contours(regress_coeff,lons,lats,clevs,ax,title='')
        plotting_functions.add_significance(pvals,lons,lats,clevs=np.array([0,0.05]))  # plot hatching to show where p values are less than 0.05, i.e. stat significant
        ax.text(-0.05,1, subplot_labels[i], transform=ax.transAxes,fontsize=25, va='top', ha='right')
        if i == 0:
            title_str = variable + ' SAM pattern in ' + dataset + ' during ' + season
            plt.title(title_str,fontsize=30) 
            plotting_functions.add_latlon_labels(ax,xticks=np.arange(-180,181,60),yticks=np.arange(-80,81,20)) # add latitude longitude labels
        ax.set_extent([-180,179,-90,20],crs=ccrs.PlateCarree())
    
    # colour bar
    ax = plt.subplot(gs[2,:])
    plotting_functions.colorbar(ax,cs)
    
    plt.subplots_adjust(hspace=0.1,wspace=0.1) # force subplots to be close together
    
    # save figure
    figure_name = Figure_dir + 'SAM_pattern_' + dataset + '_' + season + '_' + variable + str(start_year)+'-'+str(end_year)+'.png'
    print('saving to %s' % (figure_name))
    plt.savefig(figure_name,bbox_inches='tight')
    """
    plt.show()


seasons = ['DJF']
for season in seasons:
    sst_regression(season=season)

