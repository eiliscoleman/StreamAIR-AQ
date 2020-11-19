from __future__ import print_function, unicode_literals

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
from mpl_toolkits.basemap import Basemap
#import mercantile as mti

import argparse
from argparse import RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description = "read in args", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('filein', help = 'file to plot')
parser.add_argument('oproot', help = 'Path for output')
args = parser.parse_args()

tilesize = 768
dpi = 288


dir_path=args.oproot
print(dir_path)

png_out_dir=args.oproot
fid=args.filein

f_path=(png_out_dir+'/'+fid)
print(f_path)
print('this is the file to be plotted', args.filein)
data = xr.open_dataset(f_path)
LON, LAT=np.meshgrid(data.lon.values, data.lat.values)

print ('T2', data.t2.shape )
t2_out_path=(png_out_dir+'/test_t2.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=range(0,40)  
levels=lev_range
plt.contourf(LON, LAT, data.t2.values[0,:,:]-273.15, levels, cmap=plt.cm.get_cmap('hot_r'), extend="both")
plt.savefig(t2_out_path, dpi=288, tilesize=768)
plt.close()
# Wind
print ('Wind' )
wind_out_path=(png_out_dir+'/test_wind.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=range(0,30)
levels=lev_range
wind10=xu.sqrt(data.u10.values[0,:,:]**2+data.v10.values[0,:,:]**2)
plt.contourf(LON, LAT, wind10, levels, cmap=plt.cm.GnBu, extend="both")
plt.savefig(wind_out_path, dpi=288, tilesize=768)
plt.close()
 

print ('RAIN and pressure')
out_path=(png_out_dir+'/test_rain_psl.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,50,1)  
levels=lev_range
cs=plt.contourf(LON, LAT, data.rain.values[0,:,:], levels, cmap=plt.cm.Blues, extend="both")
cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n=20, linewidths=0.5, colors='k')
levels_n=np.arange(1000,1040,1)
cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n, linewidths=0.2, colors='k')
plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
plt.savefig(out_path, dpi=1152, tilesize=768)
plt.close()



#                                     Plot Rh
print('Rh')
rh_png_out_path=(png_out_dir+'/test_rh.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,101,1)
levels=lev_range
plt.contourf(LON, LAT, data.rh.values[0,:,:]  , levels, cmap=plt.cm.Spectral, extend="both")
plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
plt.close()

#                                     Plot O3
print('O3', data.o3_concentration.shape )
rh_png_out_path=(png_out_dir+'/test_o3.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0, 80, 2)
levels=lev_range
plt.contourf(LON, LAT, data.o3_concentration.values[0,0,:,:]*1000., levels, cmap=plt.cm.jet, extend="both")
plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
plt.close()
 
#                                    Plot NOx
print('NOx')
rh_png_out_path=(png_out_dir+'/test_nox.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(1, 60, 1)
levels=lev_range
plt.contourf(LON, LAT, data.nox_concentration.values[0,0,:,:]*1000., levels, cmap=plt.cm.rainbow, extend="both")
plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
plt.close()
#                 
#                                    Plot SO2
print('SO2')
png_out_path=(png_out_dir+'/test_so2.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.logspace(0,2,num=200)
levels=lev_range
plt.contourf(LON, LAT, data.so2_concentration.values[0,0,:,:]*1000000., levels, cmap=plt.cm.gnuplot2, extend="both")
plt.savefig(png_out_path, dpi=288, tilesize=768)
plt.close()

#                                    Plot PM10
print('PM10')
png_out_path=(png_out_dir+'/test_pm10.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,20,.2)
levels=lev_range
plt.contourf(LON, LAT, data.pm10.values[0,0,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both")
plt.savefig(png_out_path, dpi=288, tilesize=768)
plt.close()

#                                   Plot PM25 
print('PM25')
pm25_png_out_path=(png_out_dir+'/test_pm25.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,20,.2)
levels=lev_range
plt.contourf(LON, LAT, data.pm25.values[0,0,:,:],  levels, cmap=plt.cm.get_cmap('viridis_r'), extend="both")
plt.savefig(pm25_png_out_path, dpi=288, tilesize=768)
plt.close()
