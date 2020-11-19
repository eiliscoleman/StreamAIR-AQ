from __future__ import print_function, unicode_literals

import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import json
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import seaborn as sns
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
#import mercantile as mti
from copy import copy

import argparse
from argparse import RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description = "read in args", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('filein', help = 'file to plot')
parser.add_argument('iproot', help = 'Path for input')
parser.add_argument('oproot', help = 'Path for output')
parser.add_argument('datestr', help = 'datestr used to label output files')
args = parser.parse_args()

tilesize = 768
dpi = 288


dir_path=args.iproot
datestr=args.datestr
print(dir_path)

png_out_dir=args.oproot
fid=args.filein

f_path=(dir_path+'/'+fid)
print(f_path)
print('this is the file to be plotted', args.filein, 'will be sent to', png_out_dir, 'with date', datestr)
data = xr.open_dataset(f_path)
LON, LAT=np.meshgrid(data.lon.values, data.lat.values)



# Wind
print ('Wind' )
wind_out_path=(png_out_dir+'/'+datestr+':00-wind-10m.png')
fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
ax=fig.add_subplot(111, projection=ccrs.Mercator())
lev_range=np.arange(0,30,0.5)
levels=lev_range
wind10=xu.sqrt(data.u10.values[0,:,:]**2+data.v10.values[0,:,:]**2)
cs=ax.contourf(LON, LAT, wind10, levels, cmap=plt.cm.jet, extend="min", transform=ccrs.PlateCarree())
plt.box(on=None)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
#ax.coastlines('10m')
plt.axis('off')
ax.figsize=(tilesize/dpi, tilesize/dpi)
ax.dpi=dpi
ax.outline_patch.set_visible(False)
#ax.add_feature(cfeature.BORDERS,linewidth=0.25)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(wind_out_path, dpi=288, tilesize=768, transparent=True)
plt.close()
 
 
  ##Plot White
#print('Prs isobars - bk')
##png_out_path=(png_out_dir+'/'+datestr+':00-psl_white.png')
#png_out_path=(png_out_dir+'/outline_mercator.png')
#levels_n=np.arange(1000,1040,5)
#fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#ax=fig.add_subplot(111, projection=ccrs.Mercator(), facecolor='red')
#fig.patch.set_alpha(0)
#cs=ax.contour(LON, LAT,data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='k', alpha=1, transform=ccrs.PlateCarree())
#plt.clabel(cs, inline=1, fmt='%.0f', fontsize=2, color='k')
#plt.box(on=None)
## we need this line!
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
#ax.coastlines('10m')
#plt.axis('off')
##ax.figsize=(tilesize/dpi, tilesize/dpi)
#ax.dpi=dpi
##the line below actually controls outilne
### to see how to set the cartopy object trasparent:
##https://stackoverflow.com/questions/30076560/make-a-transparent-plot-with-cartopy
#ax.outline_patch.set_visible(False)
#ax.background_patch.set_visible(False)
#ax.patch.set_alpha(0)
##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)
##plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
#fig.savefig(png_out_path, transparent=True )
#plt.close()

#print('SO2')
#png_out_path=(png_out_dir+'/'+datestr+':00-so2.png')
#lev_range=np.arange(0,8,.002)
#levels=np.linspace(0, 100, num=11)

  ##plt.contourf(LON, LAT, data.so2_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min")
#fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#ax=fig.add_subplot(111, projection=ccrs.Mercator())  
#cs=ax.contourf(LON, LAT,data.so2_concentration.values[0,0,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min", transform=ccrs.PlateCarree())
#plt.box(on=None)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
##ax.coastlines('10m')
#plt.axis('off')
#ax.figsize=(tilesize/dpi, tilesize/dpi)
#ax.dpi=dpi
#ax.outline_patch.set_visible(False)
##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)
#plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
#plt.close()

#Plot Outline
print('Outline')
#png_out_path=(png_out_dir+'/'+datestr+':00-psl_white.png')
png_out_path=(png_out_dir+'/outline_mercator_0.25line.png')
levels_n=np.arange(1000,1040,5)
fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
ax=fig.add_subplot(111, projection=ccrs.Mercator())
cs=ax.contour(LON, LAT,data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='w', alpha=0, transform=ccrs.PlateCarree())
fig.patch.set_alpha(0)
plt.box(on=None)
# we need this line!
plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
ax.coastlines('10m', linewidth=4)
plt.axis('off')
#ax.figsize=(tilesize/dpi, tilesize/dpi)
ax.dpi=dpi
#the line below actually controls outilne
## to see how to set the cartopy object trasparent:
#https://stackoverflow.com/questions/30076560/make-a-transparent-plot-with-cartopy
#ax.add_feature(cfeature.coastlines('10m'), linewidth=0.25)
ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)
ax.patch.set_alpha(0)
#ax.add_feature(cfeature.BORDERS,linewidth=0.25)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
#plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
fig.savefig(png_out_path, transparent=True )
plt.close()

 
