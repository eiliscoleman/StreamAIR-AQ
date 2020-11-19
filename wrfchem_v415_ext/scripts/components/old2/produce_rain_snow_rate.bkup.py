#!/usr/bin/env python3
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

parser.add_argument('inroot', help = 'Path for rain input')
parser.add_argument('snowinroot', help = 'Path for snow input')
parser.add_argument('oproot', help = 'Path for output')
parser.add_argument('psl_op_dir', help = 'Path for psl output')

args = parser.parse_args()

#import mercantile as mti
#from wrf import to_np, getvar, smooth2d, get_cartopy
tilesize = 2304
dpi = 864
#dpi3 = 864

pngjson_out_dir=args.oproot
dir_op=args.inroot
dir_snow_op=args.snowinroot
#dir_op_back=('/mnt/raid/rong-ming/wrfchem/data_back/output-wrf')
dir_dfile_op=args.psl_op_dir
rain_prev = {}
rainrate = {}
rainmm = {}
# initialise the snow variables
snow_prev = {}
snowrate = {}
snowratehr = {}
snowmm = {}
no_files=len(os.listdir(dir_op))
 
for  i, file in enumerate(sorted(os.listdir(dir_op))):
   print(i, file)
   fpath=os.path.join(dir_op, file)
   newfpath=os.path.join(dir_op, file)
   data=xr.open_dataset(newfpath)
   LON, LAT=np.meshgrid(data.lon.values, data.lat.values)
   rainmm[i]=data.rain.values
   datet=data.Times.values[0]
   
   #Define the date. as read from filenames (format of "Times" variable in the files (%Y%m%d.%f) is problematic to read in )
   datestamp = fpath.split('_')[2]
   yst=datestamp[:4]
   mst=datestamp[4:-4]
   dst=datestamp[6:-2]
   hrst=datestamp[8:]
   print(datestamp, datestamp[:-2])
   dttf=datetime.strptime(datestamp[:-2], "%Y%m%d%H")
   print(dttf)
   
   pngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-rain_psl_hires.png')
   isopngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-psl_hires.png')
   convfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-rain_psl.png')
   isoconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-psl.png')
   #snowpngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl_hires.png')
   #snowconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl.png')
   print(pngfilename)
   if i == 0:
     rain_prev[i] = np.empty(data.rain.shape)
     rainrate[i]=rainmm[i]
   else:   
     rain_prev[i]=rainmm[i-1] 
     rainrate[i]=rainmm[i]-rain_prev[i]
   print(rainrate[i].shape)
   # Read in original d files in order to plot pressure isobars 
   dfileloc=(dir_dfile_op+'/d01_'+dttf.strftime("%Y%m%d%H")+'00_psl')
   data2=xr.open_dataset(dfileloc)
   # plot rainrate
   print ('plotting RAIN and pressure')
   out_path=(pngjson_out_dir+'/'+pngfilename)
   conv_out_path=(pngjson_out_dir+'/'+convfilename)
   
   
   plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   ax=fig.add_subplot(111, projection=ccrs.Mercator())
   bounds=[0, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 6.0, 12.0]
   lev_range=np.array(bounds)
   levels=lev_range
   
   
   
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
   plt.axis('off')
   #bounds=[0, .01, .1, .25, .5, 1, 2.5, 5, 10, 20, 50, 100]
   #bounds=[0, .01, .1, .2, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
   
   print(levels)
   #palette=plt.cm.get_cmap('cubehelix_r')
   #palette=copy(mc.ListedColormap(['#EAFAF1', '#ABEBC6', '#52BE80', '#5DADE2',
                               #'#2980B9', '#F39C12', '#BA4A00', '#A93226']))
   palette=['w', 'lightcyan', 'paleturquoise', 'turquoise', 'lightskyblue', 'royalblue', 'orange', 'orangered', 'firebrick', 'darkmagenta']
   cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
   cm.set_over('purple')
   norm = mc.BoundaryNorm(bounds, cm.N)
   #palette.set_under('w') #anything under 0.1 will be white
   #palette.set_over('#5B2C6F')
   rainfall=rainrate[i]
   print(rainfall.shape)
   cs=plt.contourf(LON, LAT, rainfall[0,:,:], levels, cmap=cm, norm=norm)
   print(data2.p_sl.shape)
   levels_n=np.arange(960,1040,2)
   #changed linewidth from 0.2 to 0.1
   # commented out the 2 lines underneath to prevent isobars  on same plot
   #cs2=plt.contour(LON, LAT, data2.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.1, colors='k')
   #levels_n=np.arange(1000,1040,1)

   #plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
   plt.savefig(out_path, dpi=dpi, tilesize=tilesize)
   plt.close()
   #os.system('convert -resize 768x576 ' + out_path+' '+conv_out_path)
   os.system('convert +dither -colors 68 ' + out_path+' '+conv_out_path)
   os.system('rm ' + out_path)

   
   print ('plotting just pressure isobars')
   iso_out_path=(pngjson_out_dir+'/'+isopngfilename)
   isoconv_out_path=(pngjson_out_dir+'/'+isoconvfilename)
   plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
   plt.axis('off')
   levels_n=np.arange(1000,1040,5)
   cs2=plt.contour(LON, LAT, data2.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='k')
   plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
   plt.savefig(iso_out_path, dpi=dpi, tilesize=tilesize)
   plt.close()
   #os.system('convert -resize 768x576 ' + out_path+' '+conv_out_path)
   os.system('convert +dither -colors 68 ' + iso_out_path+' '+isoconv_out_path)
   os.system('rm ' + iso_out_path)
   
   
   #Make the json files, going to pngjson_out_dir
   
   nx = 162
   ny = 114
   dx = 0.25
   dy = 0.25
   la1 = 34.5
   la2 = 62.75
   lo1 = -21.75
   lo2 = 18.5
   gridsize=162*114
   longs=np.arange(lo1,lo2+dx,dx)
   latit=np.arange(la1, la2+dy,dy)
   #times = newf.variables['Times']
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   dirpath = args.oproot
   
   outpath = os.path.join(dirpath, '%s-rain-wrf-chem-3.5.json' % (timestr,))
   gc=0
   outf = open(outpath, 'w')
   outf.write('{ ')

   rawdatavals=np.zeros([ny,nx])
   for xways in range(0, nx):
     for yways in range(0,ny):
       strlat=str(latit[yways])
       strlon=str(longs[xways])
       rawdatavals[yways, xways]=rainfall[0,yways,xways]
       rounddatavals=round(rawdatavals[yways, xways], 2)
       strval=str(rounddatavals)
       gc=gc+1
       if gc==gridsize:
         strtowrite=(' "'+strlat+', '+strlon+'":"'+strval+'" }')
       else:
         strtowrite=(' "'+strlat+','+strlon+'":"'+strval+'",')
       outf.write(strtowrite)    
   print(rawdatavals.shape, np.mean(rawdatavals), np.max(rawdatavals))
   outf.close()
     

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Snowfall


for  i, file in enumerate(sorted(os.listdir(dir_snow_op))):
   print(i, file)
   fpath=os.path.join(dir_snow_op, file)
   newfpath=os.path.join(dir_snow_op, file)
   data=xr.open_dataset(newfpath)
   LON, LAT=np.meshgrid(data.lon.values, data.lat.values)
   snowmm[i]=data.snow.values
   datet=data.Times.values[0]
   
   #Define the date. as read from filenames (format of "Times" variable in the files (%Y%m%d.%f) is problematic to read in )
   datestamp = fpath.split('_')[2]
   yst=datestamp[:4]
   mst=datestamp[4:-4]
   dst=datestamp[6:-2]
   hrst=datestamp[8:]
   print(datestamp[:], datestamp[:-2])

   dttf=datetime.strptime(datestamp[:-2], "%Y%m%d%H")
   print(dttf)
   
   pngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl_hires.png')
   convfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl.png')
   #snowpngfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl_hires.png')
   #snowconvfilename=(dttf.strftime("%Y-%m-%d-%H")+':00:00-snow_psl.png')
   print(pngfilename)
   if i == 0:
       snow_prev[i] = np.empty(data.snow.shape)
   else:   
     snow_prev[i]=snowmm[i-1] 
   #snowrate[i]=snowmm[i]-snow_prev[i]
   snowratehr[i]=snowmm[i]-snow_prev[i]
   print('shape of snowratehr', snowratehr[i].shape)
   print('type of snowratehr', type(snowratehr[i]))
   snowrate[i]=snowmm[i]
   # Read in original d files in order to plot pressure isobars 
   dfileloc=(dir_dfile_op+'/d01_'+dttf.strftime("%Y%m%d%H")+'00_psl')
   print(dfileloc)
   data2=xr.open_dataset(dfileloc)
   # plot rainrate
   print ('plotting Snow and pressure')
   out_path=(pngjson_out_dir+'/'+pngfilename)
   conv_out_path=(pngjson_out_dir+'/'+convfilename)
   plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
   plt.axis('off')
   #bounds=[0, .01, .1, .25, .5, 1, 2.5, 5, 10, 20, 50, 100]
   #bounds=[0, .01, .1, .2, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
   bounds=[0, 0.1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0, 400]
   lev_range=np.array(bounds)
   levels=lev_range
   #palette=plt.cm.get_cmap('cubehelix_r')
   #palette=copy(mc.ListedColormap(['w','#EAFAF1', '#ABEBC6', '#52BE80', '#5DADE2',
                               #'#2980B9', '#F39C12', '#BA4A00', '#A93226', '#5B2C6F']))
   palette=['w', 'lightcyan', 'paleturquoise', 'turquoise', 'lightskyblue', 'royalblue', 'orange', 'orangered', 'firebrick', 'darkmagenta']
   #palette.set_under('w')
   #palette.set_over('#5B2C6F')
   cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
   cm.set_over('purple')
   norm = mc.BoundaryNorm(bounds, cm.N)
   snowfall=snowrate[i]
   #print(snowfall.shape)
   cs=plt.contourf(LON, LAT, snowfall[0,:,:], levels, cmap=cm, extend="max", norm=norm)
   print(data2.p_sl.shape)
   
   #cpmment out lines underneath to supress isobar plotting
   #levels_n=np.arange(960,1040,2)
   #cs2=plt.contour(LON, LAT, data2.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.1, colors='k')
   #levels_n=np.arange(1000,1040,1)

   #plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
   plt.savefig(out_path, dpi=dpi, tilesize=tilesize)
   plt.close()
   #os.system('convert -resize 768x576 ' + out_path+' '+conv_out_path)
   os.system('convert +dither -colors 68 ' + out_path+' '+conv_out_path)
   os.system('rm ' + out_path)


   
   #Make the snow json files, going to pngjson_out_dir
   
   nx = 162
   ny = 114
   dx = 0.25
   dy = 0.25
   la1 = 34.5
   la2 = 62.75
   lo1 = -21.75
   lo2 = 18.5
   gridsize=162*114
   longs=np.arange(lo1,lo2+dx,dx)
   latit=np.arange(la1, la2+dy,dy)
   #times = newf.variables['Times']
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   dirpath = args.oproot
   
   outpath = os.path.join(dirpath, '%s-snow-wrf-chem-3.5.json' % (timestr,))
   outpathhr = os.path.join(dirpath, '%s-snow-hr-wrf-chem-3.5.json' % (timestr,))
   gc=0
   outf = open(outpath, 'w')
   outf.write('{ ')

   rawdatavals=np.zeros([ny,nx])
   for xways in range(0, nx):
     for yways in range(0,ny):
       strlat=str(latit[yways])
       strlon=str(longs[xways])
       rawdatavals[yways, xways]=snowfall[0,yways,xways]
       #convdatavals=rawdatavals
       rounddatavals=round(rawdatavals[yways, xways], 2)
       strval=str(rounddatavals)
       gc=gc+1
       if gc==gridsize:
         strtowrite=(' "'+strlat+', '+strlon+'":"'+strval+'" }')
       else:
         strtowrite=(' "'+strlat+','+strlon+'":"'+strval+'",')
       outf.write(strtowrite)    
   print(rawdatavals.shape, np.mean(rawdatavals), np.max(rawdatavals))
   outf.close()  
   
   gc=0
   outf = open(outpathhr, 'w')
   outf.write('{ ')

   rawdatavals=np.zeros([ny,nx])
   for xways in range(0, nx):
     for yways in range(0,ny):
       strlat=str(latit[yways])
       strlon=str(longs[xways])
       rawdatavals[yways, xways]=snowratehr[i][0,yways,xways]
       #convdatavals=rawdatavals
       rounddatavals=round(rawdatavals[yways, xways], 2)
       strval=str(rounddatavals)
       gc=gc+1
       if gc==gridsize:
         strtowrite=(' "'+strlat+', '+strlon+'":"'+strval+'" }')
       else:
         strtowrite=(' "'+strlat+','+strlon+'":"'+strval+'",')
       outf.write(strtowrite)    
   outf.close()  

