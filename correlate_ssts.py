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

Data_dir = Code_dir + 'Data/'
Figure_dir = Code_dir + 'Figures/'

def cal_grid_point_correlations(X,Y,return_pvals=False):
	""" Calculate correlations between two fields at each grid-point. Assumes
	fields have time as first dimension then two spatial dimensions."""
	
	# check that time series and field have the same time dimension
	for i in np.arange(X.ndim):
		if X.shape[i] != Y.shape[i]:
			print("Fields 1 and 2 do not have the same dimensions")
			print("Field 1 size:", X.shape, "Field 2 size:", Y.shape)
			raise ValueError
	
	corrs_map = np.zeros_like(X[0,:,:])
	pvals_map = np.zeros_like(X[0,:,:])
	
	# calculate anomalies
	X_anom = X - np.mean(X, axis=0)
	Y_anom = Y - np.mean(Y, axis=0)
	
	print("Calculating grid-point correlations...")
	milestone = 20
	for i in np.arange(X.shape[1]):
		if 100*float(i+1)/float(X.shape[1]) >= milestone:
			print("%d %% complete" % (milestone))
			milestone = milestone + 20
		
		for j in np.arange(X.shape[2]):
			slope, intercept, r_value, p_value, std_err = linregress(X_anom[:,i,j], Y_anom[:,i,j])
			corrs_map[i,j] = r_value
			pvals_map[i,j] = p_value
	
	if return_pvals == True: return corrs_map, pvals_map
	else: return corrs_map

def make_corr_map(dataset, season):
	dataset_file = Data_dir + 'CSF-20C_' + season + '_sst_data.nc'
	#read in both sst maps
	ds_var, ds_lats, ds_lons, ds_levs, ds_times, ds_calendar, ds_t_units = reading_in_data_functions.read_in_variable(dataset_file, 'sst')
	era_file = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_1x1_195001-197812.nc'
	era_var, era_lats, era_lons, era_levs, era_times, era_calendar, era_t_units = reading_in_data_functions.read_in_variable(era_file, 'tos')
	print(np.shape(era_var))
	era_file2 = '/network/group/aopp/met_data/MET001_ERA5/data/tos/mon/tos_mon_ERA5_1x1_197901-202012.nc'
	era_var2,_,_,_,era_times2,_,_= reading_in_data_functions.read_in_variable(era_file2, 'tos')
	era_var = np.append(era_var, era_var2,axis=0)
	era_times = np.append(era_times,era_times2)
	print(np.shape(era_var))
	
	start_year, end_year = 1958, 2009 #choose the period over which to calculate regression pattern
	year_mask_ds = (ds_times >= start_year)&(ds_times <= end_year)
	ds_am = ds_var[year_mask_ds,:,:]
	year_mask_era = (era_times >= start_year)&(era_times <= end_year)
	era_am = era_var[year_mask_era,:,:]
	print('year masked' + str(np.shape(era_am)))
	
	#takes ensemble mean of CSF SST at each grid point
	ens_mean_ds = []
	for ds_yr in ds_am:
		start_dim = np.shape(ds_yr)
		flattened = ds_yr.flatten()
		temp_array = [np.array(x) for x in flattened] #creates arrays of same gridpoint from all ensembles
		mean_array = [np.mean(k) for k in zip(*temp_array)] #averages same gridpoint in all ensembles
		ens_mean_ds.append(np.reshape(mean_array, start_dim)) #stores result in new array able to be analysed
		print('one year' + str(np.shape(mean_array)))
	ens_mean_ds = np.array(ens_mean_ds)
	print(np.shape(ens_mean_ds))
	
	corrs_map, pvals_map = cal_grid_point_correlations(ens_mean_ds, era_am, return_pvals=True)
	
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
	figure_name = Figure_dir + 'SST_pattern_CSF-20C_ASF-20C_DJF' + str(start_year) + '-' + str(end_year) + '.png'
	print('saving to %s' % (figure_name))
	plt.savefig(figure_name, bbox_inches='tight')
	
	plt.show()

dataset = 'CSF-20C'
season = 'DJF'
make_corr_map(dataset, season)
