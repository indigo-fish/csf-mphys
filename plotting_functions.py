""" A set of functions to plot maps """

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point


def add_latlon_labels(ax,xticks=np.arange(-80,81,10),yticks=np.arange(0,360,60),fontsize=20):
    """Add latitude and longitude labels to map plot. This function only works for the PlateCarree projection.
    
    Args:
        ax: the subplot axes which must be predefined (e.g. ax = plt.subplot(gs[0,0],projection=ccrs.PlateCarree(central_longitude=0.)))
        xticks, yticks: arrays with the longitude and latitude labels to plot.
        fontsize: font size for the longitude / latitude labels  """

    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    lon_formatter = LongitudeFormatter(number_format='.0f',degree_symbol='',dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.0f',degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)


def plot_contours(data,lons,lats,clevs,ax,title='', colors='k', transform=ccrs.PlateCarree(),
                 clabel=None, title_fontsize=25, clabel_fontsize=20,aspect=1):
    """ Add unfilled contours to a map
    
    Args:
        data: the 2d grid to be plotted
        lons, lats: corresponding longitude and latitude arrays
        clevs: contour level array
        title: title for plot
        colors: color of contours, defaults to black"""

    plt.title(title,fontsize=title_fontsize) 

    # make data loop round latitude circles. Otherwise a line of white appears at the prime meridian
    data, lons_new = add_cyclic_point(data, coord=lons)

    # add the contours
    cs = plt.contour(lons_new,lats,data,clevs,colors=colors,transform=transform)

    # add contour labels 
    if clabel is not None:
        plt.clabel(cs, fmt='%1.0f', colors='k', fontsize=clabel_fontsize)

    ax.coastlines()
    ax.set_aspect(aspect) # adjust aspect ratio of plot, defaults to 1


def plot_filled_contours(data,lons,lats,clevs,ax,title='', cmap='RdBu_r', transform=ccrs.PlateCarree(), 
        extend='both', title_fontsize=25, aspect=1):
    """ Add filled contours to a map 

    Args:
        data: the 2d grid to be plotted
        lons, lats: corresponding longitude and latitude arrays
        clevs: contour level array
        title: title for plot
        colors: color of filled contours, defaults to red-blue
    Returns:
        cs: contour levels, which can be used for colour bar plotting"""

    plt.title(title,fontsize=title_fontsize)

    # make data loop round latitude circles. Otherwise a line of white appears at the prime meridian
    #data, lons_new = add_cyclic_point(data, coord=lons)
    data = np.append(data, data[:,0].reshape(lats.shape[0], 1), axis=1)
    lons = np.append(lons, 359.9)

    # add contours
    cs = plt.contourf(lons,lats,data,clevs,cmap=cmap,extend=extend,transform=transform)

    ax.coastlines()
    ax.set_aspect(aspect)

    return cs


def add_significance(sig,lons,lats,clevs=np.array([0,0.05]),display='hatching'):
    """ Plot p values or some other measure of statistical significance.

    Args:
        sig: map of p-values (i.e. in the form [latitude,longitude])
        lons,lats: corresponding longitude, latitude arrays
        clevs: contour levels. Use np.array([0,0.05]) to plot regions with p-values less than 0.05, or np.array(0.05,1.01])
            to plot regions with p-values greater than 0.05.
        display: type of plotting, can be 'hatching' or 'contour' with the latter giving a thick black contour.
    
    """

    # make data loop round latitude circles. Otherwise a line of white appears at the prime meridian
    sig, lons_new = add_cyclic_point(sig, coord=lons)

    if display == 'hatching':
        plt.contourf(lons_new, lats, sig, clevs,transform=ccrs.PlateCarree(),colors='none',hatches=['\\\\'])
    elif display == 'contour':
        plt.contour(lons_new, lats, sig, clevs,transform=ccrs.PlateCarree(),colors='k',linewidths=2)



def colorbar(ax,cs,labelsize=20,orientation='horizontal'):
    """ Make a colour bar.

    Args: 
        ax: predefined subplot axes (e.g. ax = plt.subplot(gs[0]))
        cs: contour levels to plot, can get this from contour functions
    Returns 
        cb: colorbar axis """

    cb = plt.colorbar(cs,cax=ax,orientation=orientation)
    cb.ax.tick_params(labelsize=labelsize)
    return cb

