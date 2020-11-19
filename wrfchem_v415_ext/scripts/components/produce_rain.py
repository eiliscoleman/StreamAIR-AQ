#!/usr/bin/env python3
import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import json
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as mc
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import seaborn as sns
#from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import glob
from copy import copy
#import wrf
import argparse
from argparse import RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description = "read in args", formatter_class = RawDescriptionHelpFormatter)

#parser.add_argument('inroot', help = 'Path for rain input')
#parser.add_argument('snowinroot', help = 'Path for snow input')
#parser.add_argument('oproot', help = 'Path for output')
#parser.add_argument('psl_op_dir', help = 'Path for psl output')

#args = parser.parse_args()

#import mercantile as mti
#from wrf import to_np, getvar, smooth2d, get_cartopy
#tilesize = 2304
#dpi = 864
#dpi3 = 864
tilesize = 768
dpi = 288
#pngjson_out_dir=args.oproot
#dir_op=args.inroot
#dir_snow_op=args.snowinroot
##dir_op_back=('/mnt/raid/rong-ming/wrfchem/data_back/output-wrf')
#dir_dfile_op=args.psl_op_dir

dfiledir=('/mnt/raid/wrf-chem/wrfchem_v415_ext/data/output-wrf-remapped')
pngjson_out_dir=('/mnt/raid/wrf-chem/wrfchem_v415_ext/web')

rain_prev = {}
rainrate = {}
rainmm = {}
# initialise the snow variables
snow_prev = {}
snowrate = {}
snowratehr = {}
snowmm = {}
no_files=len(os.listdir(dfiledir))
rain_dict={}
snow_dict={}
coords_dict={}
 
for  i, file in enumerate(sorted(os.listdir(dfiledir))):
   print('file name', i, file)
   #fpath=os.path.join(dfiledir, file)
   newfpath=os.path.join(dfiledir, file)
   data=xr.open_dataset(newfpath)
   LON, LAT=np.meshgrid(data.lon.values, data.lat.values)
   
   rainmm[i]=data.rain.values
   snowmm[i]=data.snow.values
   datet=data.time.values
   dd=' '.join(map(str, datet))
   dates=str(dd).split(".")
   dt=dates[0].split('[')
   strr=('0.'+(dates[1]))
   hours=(float(strr)*24)
   hr=round(hours, 0)
   hint=int(hr)
   dt=pd.to_datetime(str(dates[0]))
   dttf=dt.replace(hour=hint)

   pngfilename = (dttf.strftime("%Y-%m-%d-%H")+':00:00-rain_hires.png')
   isopngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-psl_hires.png')
   snowpngfilename = (dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_hires.png')
   snowconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow.png')
   convfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-rain.png')
   isoconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-psl.png')
   #snowpngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl_hires.png')
   #snowconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl.png')
   print(pngfilename)
   if i == 0:
     rain_prev[i] = np.empty(data.rain.shape)
     rainrate[i]=rainmm[i]
     snow_prev[i] = np.empty(data.snow.shape)
     snowrate[i]=snowmm[i]
   else:   
     rain_prev[i]=rainmm[i-1] 
     rainrate[i]=rainmm[i]-rain_prev[i]
     snow_prev[i]=snowmm[i-1] 
     snowrate[i]=snowmm[i]-snow_prev[i]

   print ('plotting RAIN and pressure')
   out_path=(pngjson_out_dir+'/'+pngfilename)
   conv_out_path=(pngjson_out_dir+'/'+convfilename)
   snow_out_path=(pngjson_out_dir+'/'+snowpngfilename)
   snow_conv_out_path=(pngjson_out_dir+'/'+snowconvfilename)
   
   fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   ax=fig.add_subplot(111, projection=ccrs.Mercator())
   #bounds=[0, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 6.0, 12.0]
   bounds=[0, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 6.0, 12.0]
   lev_range=np.array(bounds)
   levels=lev_range
   palette=['w', 'lightcyan', 'paleturquoise', 'turquoise', 'lightskyblue', 'royalblue', 'orange', 'orangered', 'firebrick', 'darkmagenta']
   cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
   cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
   cm.set_over('purple')
   
   #norm = mc.BoundaryNorm(bounds, cm.N)
   #norm = mc.BoundaryNorm(bounds, plt.cm.gnuplot2)
   norm = mc.BoundaryNorm(boundaries=bounds, ncolors=256)
   rainfall=rainrate[i]

   #print('rainshape', rainfall.shape)
   #cs=ax.contourf(LON, LAT, rainfall[0,:,:], levels, cmap=plt.cm.BuPu, norm=norm, transform=ccrs.PlateCarree())
   cs=ax.contourf(LON, LAT, rainfall[0, :,:], norm=norm, levels=levels, cmap=plt.cm.gnuplot2, transform=ccrs.PlateCarree(), extend='max')
   
   cs.cmap.set_over('w')
   plt.box(on=None)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
   #ax.coastlines('10m', linewidth=0.15)
   #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
   #cbar=plt.colorbar(cs, ax=ax)
   #mpl.colorbar.ColorbarBase(ax, cmap=plt.cm.gnuplot2, orientation='horizontal', extend='max')
   plt.axis('off')
   ax.figsize=(tilesize/dpi, tilesize/dpi)
   ax.dpi=dpi
   ax.outline_patch.set_visible(False)
   ax.background_patch.set_visible(False)
   ax.patch.set_alpha(0)
   #ax.add_feature(cfeature.BORDERS,linewidth=0.25)
   ax.axes.get_xaxis().set_visible(False)
   ax.axes.get_yaxis().set_visible(False)
   plt.savefig(out_path, dpi=dpi, tilesize=tilesize, transparent=True)
   plt.close()
   #os.system('convert -resize 768x576 ' + out_path+' '+conv_out_path)
   os.system('convert +dither -colors 68 ' + out_path+' '+conv_out_path)
   #os.system('rm ' + out_path)

   
   print ('plotting just pressure isobars')
   iso_out_path=(pngjson_out_dir+'/'+isopngfilename)
   isoconv_out_path=(pngjson_out_dir+'/'+isoconvfilename)
   plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
   plt.axis('off')
   levels_n=np.arange(900,1240,5)
   cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='k')
   plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
   plt.savefig(iso_out_path, dpi=dpi, tilesize=tilesize)
   plt.close()
   
   os.system('convert -resize 768x576 ' + out_path+' '+conv_out_path)
   os.system('convert +dither -colors 68 ' + iso_out_path+' '+isoconv_out_path)
   os.system('rm ' + iso_out_path)
   if i == 0:
       os.system('rm ' + conv_out_path)
       os.system('rm ' + out_path)
   print ('plotting snow filled contour plot')
   fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   ax=fig.add_subplot(111, projection=ccrs.Mercator())
   
   bounds=[0, 0.1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0, 400]
   lev_range=np.array(bounds)
   levels=lev_range
   
   palette=['w', 'lightcyan', 'paleturquoise', 'turquoise', 'lightskyblue', 'royalblue', 'orange', 'orangered', 'firebrick', 'darkmagenta']
   #palette.set_under('w')
   #palette.set_over('#5B2C6F')
   cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
   cm.set_over('purple')
   norm = mc.BoundaryNorm(bounds, cm.N)
   snowfall=snowrate[i]
   
   cs=ax.contourf(LON, LAT, snowfall[0,:,:], levels, cmap=cm, norm=norm, transform=ccrs.PlateCarree())
   plt.box(on=None)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
   #ax.coastlines('10m', linewidth=0.15)
   plt.axis('off')
   ax.figsize=(tilesize/dpi, tilesize/dpi)
   ax.dpi=dpi
   ax.outline_patch.set_visible(False)
   ax.background_patch.set_visible(False)
   ax.patch.set_alpha(0)
   #ax.add_feature(cfeature.BORDERS,linewidth=0.25)
   ax.axes.get_xaxis().set_visible(False)
   ax.axes.get_yaxis().set_visible(False)
   plt.savefig(snow_conv_out_path, dpi=dpi, tilesize=tilesize, transparent=True)
   plt.close()
   
   #Make the json files, going to pngjson_out_dir
   
   
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   
   rain_json_path = os.path.join(pngjson_out_dir, '%s-rain.json' % (timestr,))
   snow_json_path = os.path.join(pngjson_out_dir, '%s-snow.json' % (timestr,))
   
   for i, lonv in enumerate(data.lon.values):
       for j, latv in enumerate(data.lat.values):
           coords_string = (str(latv)+','+str(lonv))
           coords_dict[j,i] = coords_string
           rain_dict[j,i]=rainfall[0,j,i]
           snow_dict[j,i]=snowfall[0,j,i]
   coords_df = pd.DataFrame.from_dict(coords_dict, orient='index')
   rain_df = pd.DataFrame.from_dict(rain_dict, orient='index')
   snow_df = pd.DataFrame.from_dict(snow_dict, orient='index')
   rain_df['rain']=rain_df.values
   rain_df = rain_df.round({'rain': 2})
   snow_df['snow']=snow_df.values
   snow_df = snow_df.round({'snow': 2})
   #print(rain_df.shape)
   rain_df['coords'] = coords_df
   snow_df['coords'] = coords_df
   text_output = rain_df.set_index('coords').to_dict()['rain']
   snow_text_output=snow_df.set_index('coords').to_dict()['snow']
   with open(rain_json_path, 'w') as fp:
       json.dump(text_output, fp)
   with open(snow_json_path, 'w') as fp:
       json.dump(snow_text_output, fp)

   if i == 0:
       os.system('rm ' + rain_json_path)
       os.system('rm ' + snow_json_path)
       print('i should have removed', rain_json_path)
     