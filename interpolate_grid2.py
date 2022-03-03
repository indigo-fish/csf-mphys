from scipy import interpolate
import numpy as np

def interpolate_grid(old_grid_data,lons_old,lats_old,lons_new,lats_new,method='cubic'):
    """ Assumes that data is of the form [time,latitude,longitude] or [latitude,longitude] """
    print('Interpolating data onto new grid...')

    if old_grid_data.ndim == 2:

        X, Y = np.meshgrid(lons_old, lats_old)
        XI, YI = np.meshgrid(lons_new,lats_new)

        new_grid =  interpolate.griddata((X.flatten(),Y.flatten()),old_grid_data.flatten() , (XI,YI),method=method)

    elif old_grid_data.ndim == 3:

        new_grid = np.zeros([0,lats_new.shape[0],lons_new.shape[0]])
        X, Y = np.meshgrid(lons_old, lats_old)
        XI, YI = np.meshgrid(lons_new,lats_new)

        for i in np.arange(old_grid_data.shape[0]):

            new_grid_temp = interpolate.griddata((X.flatten(),Y.flatten()),old_grid_data.flatten() , (XI,YI),method=method)
            new_grid = np.append(new_grid,new_grid_temp.reshape(1,lats_new.shape[0],lons_new.shape[0]),axis=0)

    return new_grid
