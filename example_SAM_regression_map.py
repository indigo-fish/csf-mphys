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
season = 'DJF' # choose a season from DJF/MAM/JJA/SON
start_year,end_year = 1958,1978 # choose the period over which to calculate regression pattern

sys.path.append(Code_dir)
import plotting_functions
import reading_in_data_functions
import analysis_functions

# read in the SAM index
SAM_idx, years_SAM = reading_in_data_functions.read_Marshall_SAM_idx(season=season)


# Read in monthly geopotential height data on the 500hPa level
# this is split across two files so need to read in both and combine them
file_str_ERA5 = '/network/group/aopp/met_data/MET001_ERA5/data/zg/mon/zg_mon_ERA5_2.5x2.5_195001-197812.nc'
zg,lats,lons,levs,times,calendar,t_units = reading_in_data_functions.read_in_variable(file_str_ERA5,'geopotential_height',chosen_level=500)
file_str_ERA5_2 = '/network/group/aopp/met_data/MET001_ERA5/data/zg/mon/zg_mon_ERA5_2.5x2.5_197901-202012.nc'
zg2,_,_,_,times2,_,_ = reading_in_data_functions.read_in_variable(file_str_ERA5_2,'geopotential_height',chosen_level=500) # weirdly, these two files use different variable names zg vs geopotential_height
zg = np.append(zg,zg2,axis=0)
times = np.append(times,times2)

# calculate the seasonal mean geopotential height
zg_am, years_zg = reading_in_data_functions.calculate_annual_mean(zg,times,calendar,t_units,season=season)


# restrict geopotential height and SAM index to chosen years
year_mask_zg = (years_zg>=start_year)&(years_zg<=end_year)
zg_am = zg_am[year_mask_zg,:,:]
year_mask_SAM = (years_SAM>=start_year)&(years_SAM<=end_year)
SAM_idx = SAM_idx[year_mask_SAM]


# calculate regression pattern
regress_coeff,corr,pvals = analysis_functions.regress_map(SAM_idx,zg_am)


# make some plots
plt.figure(figsize=(15,7))
gs = gridspec.GridSpec(2,1,height_ratios=[10,0.5])

#a = np.arange(5,31,5)
#clevs = np.append(-a[::-1],a) # contour level
clevs = np.linspace(-20, 15, 10)
#subplot_labels = ['a)','b)']

#for i, projection in enumerate([ccrs.PlateCarree(central_longitude=0.),ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0) ]):
for i, projection in enumerate([ccrs.PlateCarree(central_longitude=0.)]):
    ax = plt.subplot(gs[i,0],projection=projection)
    cs = plotting_functions.plot_filled_contours(regress_coeff,lons,lats,clevs,ax,title='')
    plotting_functions.add_significance(pvals,lons,lats,clevs=np.array([0,0.05]))  # plot hatching to show where p values are less than 0.05, i.e. stat significant
    #ax.text(-0.05,1, subplot_labels[i], transform=ax.transAxes,fontsize=25, va='top', ha='right')
    if i == 0:
        plt.title('Z500 SAM pattern in '+season+' from ERA5 (m)',fontsize=30) 
        plotting_functions.add_latlon_labels(ax,xticks=np.arange(-180,181,60),yticks=np.arange(-80,81,20)) # add latitude longitude labels
    ax.set_extent([-180,179,-90,-1],crs=ccrs.PlateCarree())

# colour bar
ax = plt.subplot(gs[1,:])
plotting_functions.colorbar(ax,cs)

plt.subplots_adjust(hspace=0.1,wspace=0.1) # force subplots to be close together

# save figure
figure_name = Figure_dir + 'SAM_pattern_Z500_ERA5_'+str(start_year)+'_'+str(end_year)+'_'+season+'.png'
print('saving to %s' % (figure_name))
plt.savefig(figure_name,bbox_inches='tight')
#plt.show()
