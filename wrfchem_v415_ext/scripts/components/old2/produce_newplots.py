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

print ('T2', data.t2.shape )
t2_out_path=(png_out_dir+'/'+datestr+':00-t2.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=range(-20,40)  
levels=lev_range
plt.contourf(LON, LAT, data.t2.values[0,:,:]-273.15, levels, cmap=plt.cm.get_cmap('jet'), extend="both")
plt.savefig(t2_out_path, dpi=288, tilesize=768)
plt.close()
# Wind
print ('Wind' )
wind_out_path=(png_out_dir+'/'+datestr+':00-wind-10m.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,30,0.05)
levels=lev_range
wind10=xu.sqrt(data.u10.values[0,:,:]**2+data.v10.values[0,:,:]**2)
plt.contourf(LON, LAT, wind10, levels, cmap=plt.cm.jet, extend="min")
plt.savefig(wind_out_path, dpi=288, tilesize=768)
plt.close()
 

#print ('RAIN and pressure')
#out_path=(png_out_dir+'/test_rain_psl.png')
#plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
#plt.axis('off')
##bounds=[0, .01, .1, .3, .5, .7, 1, 2, 3, 6, 12, 20, 50, 100]
#bounds=[0, .1, .3, .5, .7, 1, 2, 3,  5, 6, 12]
##lev_range=np.logspace(0,2, num=8)  
#lev_range=np.array(bounds)
#levels=lev_range
#palette=plt.cm.get_cmap('cubehelix_r')
#palette.set_over('r')
#cs=plt.contourf(LON, LAT, data.rain.values[0,:,:], levels, cmap=palette, extend="max")
#cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n=20, linewidths=0.5, colors='k')
#levels_n=np.arange(1000,1040,2)
##cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n, linewidths=0.2, colors='k')
#plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
#plt.savefig(out_path, dpi=1152, tilesize=768)
#plt.close()

###                                     Plot AQI
#print('aqi')
#rh_png_out_path=(png_out_dir+'/'+datestr+':00-aqi.png')
#plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
#plt.axis('off')
##bounds=np.arange[0, 11, 1]
#lev_range=np.arange(0,11,1)
#levels=lev_range
#aqiraw=data.rh.values[0,:,:]
#aqiscaled=aqiraw[:,:]*0.1
#aqiround=np.around([aqiscaled], decimals=0)
#print(aqiround.shape)
#plt.contourf(LON, LAT, aqiround[0,:,:], levels, cmap=plt.cm.viridis, extend="both")
#plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
#plt.close()


##                                     Plot Rh
print('Rh')
rh_png_out_path=(png_out_dir+'/'+datestr+':00-rh.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(30,101,1)
levels=lev_range
plt.contourf(LON, LAT, data.rh.values[0,:,:], levels, cmap=plt.cm.get_cmap('gist_rainbow_r'), extend="both")
plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
plt.close()

#The following plots need to be plotted for all levels
#                                     Plot O3
#for lev in range(0,29):
for lev in range(0,1):

  levstr=str(lev)
  levpadded=levstr.zfill(2)
  print(levstr)
  print('O3', data.o3_concentration.shape, lev )
#  rh_png_out_path=(png_out_dir+'/'+datestr+'-o3-lev'+levpadded+'.png')
  rh_png_out_path=(png_out_dir+'/'+datestr+':00-o3.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  lev_range=np.arange(30, 50, .2)
  levels=lev_range
  plt.contourf(LON, LAT, data.o3_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.jet, extend="min")
  plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
  plt.close()
 
#                                    Plot NOx
  print('NOx')
  rh_png_out_path=(png_out_dir+'/'+datestr+':00-nox.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  lev_range=np.arange(1, 60, 1)
  levels=lev_range
  plt.contourf(LON, LAT, data.nox_concentration.values[0,lev,:,:]*1000., levels  , cmap=plt.cm.gnuplot, extend="min")
  plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
  plt.close()
#                 
#                                    Plot SO2

  print('SO2')
  png_out_path=(png_out_dir+'/'+datestr+':00-so2.png')
  #png_out_path=('/mnt/raid/rong-ming/wrfchem/'+datestr+':00-so2.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  #bounds=[ 1.e-5, 2.e-5, 3.e-5, 4.e-5. 5.e-5, 6.e-5, 7.e-5, 8.e-5, 9.e-5, 1.e-4 ]
  lev_range=np.arange(0,8,.002)
  #levels=lev_range
  ##levels=np.logspace(0,2, num=9)
  levels=np.linspace(0, 100, num=11)
  print(levels)
  print(data.so2_concentration.values[0,lev,:,:]*1000.)
  plt.contourf(LON, LAT, data.so2_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min")
  plt.savefig(png_out_path, dpi=288, tilesize=768)
  plt.close()
  #print('SO2')
  #png_out_path=(png_out_dir+'/'+datestr+':00-so2.png')
  #plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  #plt.axis('off')
  #lev_range=np.arange(0,8,.002)
  #levels=lev_range
  #plt.contourf(LON, LAT, data.so2_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min")
  #plt.savefig(png_out_path, dpi=288, tilesize=768)
  #plt.close()

#                                    Plot PM10
  print('PM10')
  png_out_path=(png_out_dir+'/'+datestr+':00-pm10.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  lev_range=np.arange(0,20,.2)
  levels=lev_range
  plt.contourf(LON, LAT, data.pm10.values[0,lev,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  plt.savefig(png_out_path, dpi=288, tilesize=768)
  plt.close()

#                                   Plot PM25 
  print('PM25')
  pm25_png_out_path=(png_out_dir+'/'+datestr+':00-pm25.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  lev_range=np.arange(0,20,.2)
  levels=lev_range
  plt.contourf(LON, LAT, data.pm25.values[0,lev,:,:],  levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  plt.savefig(pm25_png_out_path, dpi=288, tilesize=768)
  plt.close()
  
  #                                   Plot Black Pressure lines 
  print('Prs isobars - bk')
  pm25_png_out_path=(png_out_dir+'/'+datestr+':00-psl_black.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  levels_n=np.arange(1000,1040,5)
  cs=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='k')
  plt.clabel(cs, inline=1, fmt='%.0f', fontsize=2, color='k')
  plt.savefig(pm25_png_out_path, dpi=288, tilesize=768, transparent=True)
  plt.close()

#                                   Plot White Pressure Lines 
  print('Prs isobars - white')
  pm25_png_out_path=(png_out_dir+'/'+datestr+':00-psl_white.png')
  plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
  plt.axis('off')
  levels_n=np.arange(1000,1040,5)
  cs=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='w')
  plt.clabel(cs, inline=1, fmt='%.0f', fontsize=2, color='w')
  plt.savefig(pm25_png_out_path, dpi=288, tilesize=768, transparent=True)
  plt.close()

