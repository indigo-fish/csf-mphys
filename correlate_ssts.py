import numpy as np
from scipy.stats import linregress
import os,sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs

Code_dir = '/home/w/wadh5699/Desktop/Example_Scripts/Amelia_example_scripts/'
sys.path.append(Code_dir)
import plotting_functions
import reading_in_data_functions
from interpolate_grid2 import interpolate_grid
from cal_grid_point_correlations import cal_grid_point_correlations

Data_dir = Code_dir + 'Data/'
Figure_dir = Code_dir + 'Figures/'

def make_corr_map(dataset, season):
	dataset_file = '/home/p/pattersonm/sst_ecmf-guh4_1101_fcmean_sfc_1901_2010_DJF.nc'
	#read in both sst maps
	ds_var, ds_lats, ds_lons, ds_levs, ds_times, ds_calendar, ds_t_units = reading_in_data_functions.read_in_variable(dataset_file, 'sst')
	era_file = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_1x1_195001-197812.nc'
	era_var, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = reading_in_data_functions.read_in_variable(era_file, 'tos')
	print(np.shape(era_var))
	era_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_1x1_197901-202012.nc'
	era_var2,_,_,_,era_times2,_,_= reading_in_data_functions.read_in_variable(era_file2, 'tos')
	era_var = np.append(era_var, era_var2,axis=0)
	era_times = np.append(era_times,era_times2)
	
	
	start_year, end_year = 1958, 2010 #choose the period over which to calculate regression pattern
	year_mask_ds = (ds_times >= start_year)&(ds_times <= end_year)
	ds_am = ds_var[year_mask_ds,:,:]
	year_mask_era = (era_times >= start_year)&(era_times <= end_year)
	era_am = era_var[year_mask_era,:,:]
	
	interpolated_ds = interpolate_grid(ds_am, ds_lons, ds_lats, era_lons, era_lats)
	
	corrs_map, pvals_map = cal_grid_point_correlations(interpolated_ds, era_am, return_pvals=True)
	
	#make plots
	plt.figure(figsize=(15,15))
	gs = gridspec.GridSpec(3,1,height_ratios=[5,10,0.5])
	
	clevs = np.linspace(np.nanmin(corrs_map), np.nanmax(corrs_map), 10)
	subplot_labels = ['a)', 'b)']
	
	for i, projection in enumerate([ccrs.PlateCarree(central_longitude=0.),ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0) ]):
		
		ax = plt.subplot(gs[i,0], projection=projection)
		cs = plotting_functions.plot_filled_contours(corrs_map, ds_lons, ds_lats, clevs, ax, title='')
		plotting_functions.add_significance(pvals_map, ds_lons, ds_lats, clevs=np.array([0,0.05])) #plot hatching to show where p values are less than 0.05, ie statistically significant
		ax.text(-0.05,1, subplot_labels[i], transform = ax.transAxes, fontsize=25, va='top', ha='right')
		if i == 0:
			title_str = 'SST correlation between ' + dataset + ' and ERA5 during ' + season
			plt.title(title_str, fontsize=30)
			plotting_functions.add_latlon_labels(ax,xticks=np.arange(-180,181,60),yticks=np.arange(-80,81,20)) #add latitude longitude labels
		ax.set_extent([-180,179,-90,20],crs=ccrs.PlateCarree())
	
	#colour bar
	ax = plt.subplot(gs[2,:])
	plotting_functions.colorbar(ax,cs)
	
	plt.subplots_adjust(hspace=0.1, wspace=0.1) #force subplots to be close together
	
	#save figure
	figure_name = Figure_dir + 'SST_correlations_CSF-20C_ERA_DJF' + str(start_year) + '-' + str(end_year) + '.png'
	print('saving to %s' % (figure_name))
	plt.savefig(figure_name, bbox_inches='tight')
	
	plt.show()

dataset = 'CSF-20C'
season = 'DJF'
make_corr_map(dataset, season)
