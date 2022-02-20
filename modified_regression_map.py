""" Script to read in the Marshall SAM index and calculate the associated geopotential height pattern
via linear regression analysis.

Plots this using two different maps (usually only one is necessary)"""

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import cartopy.crs as ccrs
import os, sys

Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
Figure_dir = Code_dir + 'Figures/' # directory to save figures to
start_year,end_year = 1982,2009 # choose the period over which to calculate regression pattern

sys.path.append(Code_dir)
import plotting_functions
import reading_in_data_functions
import analysis_functions
import calculate_csf_SAM

def runregression(variable, season, level_to_choose=None, dataset=None):
    
    # read in the SAM index
    if dataset == None: SAM_idx, years_SAM = reading_in_data_functions.read_Marshall_SAM_idx(season=season)
    else: SAM_idx, years_SAM, calendar, units = calculate_csf_SAM.read_SAM_indices(dataset=dataset, season=season)
    
    # Read in monthly geopotential height data on the 500hPa level
    # this is split across two files so need to read in both and combine them
    file_str_ERA5 = '/network/group/aopp/met_data/MET001_ERA5/data/' + variable + '/mon/' + variable + '_mon_ERA5_2.5x2.5_195001-197812.nc'
    if variable == 'zg':
        var,lats,lons,levs,times,calendar,t_units = reading_in_data_functions.read_in_variable(file_str_ERA5,'geopotential_height',chosen_level=level_to_choose)
    else:
        var,lats,lons,levs,times,calendar,t_units = reading_in_data_functions.read_in_variable(file_str_ERA5, variable, chosen_level=level_to_choose)
    file_str_ERA5_2 = '/network/group/aopp/met_data/MET001_ERA5/data/' + variable + '/mon/' + variable + '_mon_ERA5_2.5x2.5_197901-202012.nc'
    var2,_,_,_,times2,_,_ = reading_in_data_functions.read_in_variable(file_str_ERA5_2,variable,chosen_level=level_to_choose) # weirdly, these two files use different variable names zg vs geopotential_height
    var = np.append(var,var2,axis=0)
    var = var * 100
    times = np.append(times,times2)
        
    
    # calculate the seasonal mean geopotential height
    var_am, years_var = reading_in_data_functions.calculate_annual_mean(var,times,calendar,t_units,season=season)
    
    
    # restrict geopotential height and SAM index to chosen years
    year_mask_var = (years_var>=start_year)&(years_var<=end_year)
    var_am = var_am[year_mask_var,:,:]
    year_mask_SAM = (years_SAM>=start_year)&(years_SAM<=end_year)
    SAM_idx = SAM_idx[year_mask_SAM]
    
    
    # calculate regression pattern
    regress_coeff,corr,pvals = analysis_functions.regress_map(SAM_idx,var_am)
    
    # make some plots
    plt.figure(figsize=(15,15))
    gs = gridspec.GridSpec(3,1,height_ratios=[5,10,0.5])
    
    if variable=='tos': clevs = np.linspace(-50, 50, 10)
    else: clevs = np.linspace(min(corr), max(corr), 10)
    #clevs = np.append(-a[::-1],a) # contour level
    subplot_labels = ['a)','b)']
    
    for i, projection in enumerate([ccrs.PlateCarree(central_longitude=0.),ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0) ]):
    
        ax = plt.subplot(gs[i,0],projection=projection)
        cs = plotting_functions.plot_filled_contours(regress_coeff,lons,lats,clevs,ax,title='')
        plotting_functions.add_significance(pvals,lons,lats,clevs=np.array([0,0.05]))  # plot hatching to show where p values are less than 0.05, i.e. stat significant
        ax.text(-0.05,1, subplot_labels[i], transform=ax.transAxes,fontsize=25, va='top', ha='right')
        if i == 0:
            if dataset==None: title_str = variable + str(level_to_choose) + ' SAM pattern in ' + season + ' from ERA5 (m)'
            else: title_str = variable + str(level_to_choose) + ' SAM pattern in ' + season + ' from ERA5 (m) and ' + dataset
            plt.title(title_str,fontsize=30) 
            plotting_functions.add_latlon_labels(ax,xticks=np.arange(-180,181,60),yticks=np.arange(-80,81,20)) # add latitude longitude labels
        ax.set_extent([-180,179,-90,20],crs=ccrs.PlateCarree())
    
    # colour bar
    ax = plt.subplot(gs[2,:])
    plotting_functions.colorbar(ax,cs)
    
    plt.subplots_adjust(hspace=0.1,wspace=0.1) # force subplots to be close together
    
    # save figure
    figure_name = Figure_dir + 'SAM_pattern_' + variable + str(level_to_choose) + '_'
    if dataset==None: figure_name += 'ERA5_' + str(start_year) + '-' + str(end_year) + '_' + season + '.png'
    else: figure_name += dataset + '_ERA5_' + str(start_year) + '-' + str(end_year) + '_' + season + '.png' 
    print('saving to %s' % (figure_name))
    plt.savefig(figure_name,bbox_inches='tight')
    #plt.show()

levels = [0] #[250, 500, 750, 1000]
seasons = ['DJF']
datasets = [None, 'CSF-20C', 'ASF-20C']
for level in levels:
    for season in seasons:
        for dataset in datasets:
            runregression('tos', season, dataset=dataset)

