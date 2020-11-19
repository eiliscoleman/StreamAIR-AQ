#!/usr/bin/env python3
import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfea
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
import scipy.interpolate
#from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
###from nco import Nco
from cdo import Cdo
import shutil
import argparse
from argparse import RawDescriptionHelpFormatter
#import mercantile as mti


#details on EU AQI: https://www.eea.europa.eu/themes/air/air-quality-index/index

#from wrf import to_np, getvar, smooth2d, get_cartopy
parser = argparse.ArgumentParser(description = "read in args", formatter_class = RawDescriptionHelpFormatter)

parser.add_argument('basdir_path', help = 'Path for basdir')

args = parser.parse_args()

#nco=Nco()
cdo=Cdo()

tilesize = 768
dpi = 288
dpi3 = 864
basdir=args.basdir_path
dir_op_back=(basdir+'/data_back/output-wrf-remapped')
dir_op=(basdir+'/data/output-wrf-remapped')
dir_png_json_op=(basdir+'/web')
#dir_aqi=('/mnt/raid/rong-ming/wrfchem/data/output-wrf-aqi')
# Create list of AQI for various pollutants
species=['o3', 'nox', 'so2', 'pm25', 'pm10']
#bboundsugm3=(
  #[0,34, 67,101,121,141,161,188,214,241],
  #[0,68,135,201,268,335,401,468,535,601],
  #[0,30,60,90,120,150,180,237,296,355],
  #[0,12,24,36,42,48,54,59,65,71],
  #[0,17,34,51,59,67,76,84,92,101]
#)
bboundsugm3=(
  [0,50, 100,130, 240, 380],
  [0,40,90,120, 230, 340],
  [0,100,200,350,500, 750],
  [0,10,20,25,50, 75],
  [0,20,40,50,100, 150]
)
#aqi_colours=['yellow', 'yellowgreen', 'forestgreen', 'gold', 'orange', 'coral', 'orangered', 'firebrick', 'brown', 'm']
aqi_colours=['cyan', 'mediumseagreen', 'yellow', 'tomato', 'darkred', 'indigo' ]
listAQI={
  'spcs':pd.Series(species, index=np.arange(0,len(species))),
  #'desc':pd.Series(['good1', 'good2', 'good3', 'fair1', 'fair2', 'fair3', 'poor1', 'poor2', 'poor3', 'very poor'],index=np.arange(1,11,1))
  'desc':pd.Series(['Good', 'Fair', 'Moderate', 'Poor', 'Very Poor', 'Extremely Poor'],index=np.arange(1,7,1))
}

listAQI['o3']={
  'boundsugm3':pd.Series([0,50, 100,130, 240, 380], index=np.arange(1,7,1)),
  'hrspan':1,
  'varname':'o3_concentration'
}
listAQI['nox']={
  'boundsugm3':pd.Series([0,40,90,120, 230, 340],index=np.arange(1,7,1)),
  'hrspan':1,
  'varname':'nox_concentration'
}
listAQI['so2']={
  'boundsugm3':pd.Series( [0,100,200,350,500, 750],index=np.arange(1,7,1)),
  'hrspan':1,
  'varname':'so2_concentration'
}
listAQI['pm25']={
  'boundsugm3':pd.Series([0,10,20,25,50, 75],index=np.arange(1,7,1)),
  'hrspan':1,
  'varname':'pm25'
}
listAQI['pm10']={
  'boundsugm3':pd.Series([0,20,40,50,100, 150],index=np.arange(1,7,1)),
  'hrspan':1,
  'varname':'pm10'
}

#Create colormap
cm=LinearSegmentedColormap.from_list('aqi_cmap', aqi_colours, N=len(aqi_colours))

#print(listAQI['spcs']['o3']['boundsugm3'])
for file in sorted(os.listdir(dir_op)):
   fpath=os.path.join(dir_op, file)
   fpath_bk=os.path.join(dir_op_back, file)
   # file already copied from aq_who routine, called before this current routine
   #shutil.copyfile(fpath, fpath_bk)
   ##print('file to be opened', fpath)
   data_ref=xr.open_dataset(fpath)
   # Interpolate the data to regular grid
   # Step 1: apply the wrf out raw grid to the current file
   
   
   
   
   #print( data_ref.variables['lat'].shape, data_ref.variables['lon'].shape)
   
   nolats, nolons=len(data_ref['lat'][:]), len(data_ref['lon'][:])
   #print('lat,lon', data_ref['lat'].shape, nolats, nolons)
   lats = data_ref.variables['lat'][:]
   lons = data_ref.variables['lon'][:]
   #print(fpath.split("reg_",1)[1])
   datestamp=fpath.split("reg_",1)[1]
   #datestamp=fpath[-15:-5]
   yst=datestamp[:4]
   mst=datestamp[4:-4]
   dst=datestamp[6:-2]
   hrst=datestamp[8:]
   dttf=datetime.strptime(datestamp, "%Y%m%d%H")
   datestr=dttf.strftime("%Y-%m-%d-%H:00")
   #filedate=dtti.strftime("%Y%m%d%H")
   print(datestr)
   
   
   #-----------------------
#     Loop structure:
#     
#    for species in list:
#       get hrspan y
#       open relevant files  y 
#       get average y
#       convert units y
#       apply aq limits for each species
#       calucate final aqi
#       rm d file from the output-wrf directory
#   -----------------------
   #spc_aqi=np.empty([141,156,5])
   
   spc_aqi=np.empty([nolats,nolons,len(species)])
   spi=0
   for   sp in listAQI['spcs']:
     ##print('species is', sp)


     ave_period=listAQI[sp]['hrspan']
     spc_hrtmp=np.empty([nolats,nolons,ave_period])
     
     
     
     for hour in range(ave_period):
       dtti=dttf+timedelta(hours=-hour)
       filedate=dtti.strftime("%Y%m%d%H")
       filepath=(dir_op_back+'/d01reg_'+filedate)
       #print(filepath)
       datatmp=xr.open_dataset(filepath)
       spc_ave=np.empty([nolats,nolons])
       spc_ave_conv=np.empty([nolats,nolons])       
       if sp == 'o3':   
         #print(datatmp.o3_concentration.shape)
         spc_hrtmp[:,:,hour]=datatmp.o3_concentration[0,0,:,:]
         #factor to convert ppm to ug/m3
         conv= 1995.7 # 1000pp = 1 ppm, 1 ppb = 2.00 ug/m3
       elif sp == 'nox':
         spc_hrtmp[:,:,hour]=datatmp.nox_concentration[0,0,:,:]
         conv= 1912.5 # 1000pp = 1 ppm, 1 ppb = 1.9125 ug/m3         
       elif sp == 'so2':
         
         spc_hrtmp[:,:,hour]=datatmp.so2_concentration[0,0,:,:]
         conv= 2660.9 # 1000pp = 1 ppm, 1 ppb = 2.6609 ug/m3         
       elif sp == 'pm25':
              
         spc_hrtmp[:,:,hour]=datatmp.pm25[0,0,:,:]
         #print(spc_hrtmp[88,77,hour], hour)
         conv= 1 # already in ug/m3
       elif sp == 'pm10':
         spc_hrtmp[:,:,hour]=datatmp.pm10[0,0,:,:]         
         conv= 1 # already in ug/m3
     
     #get average over the timespan of each array
     spc_ave=np.mean(spc_hrtmp, axis=2)
     spc_ave_conv=spc_ave*conv
     #print(sp, spc_ave_conv.shape, spc_ave_conv[88,77])
     lims=listAQI[sp]['boundsugm3']
     
     #np.digitize returns indices of specified bins difined in boundsugm3
     
     spc_aqi[:,:,spi]=np.digitize(spc_ave_conv, lims, right=True)
     vals=spc_aqi[:,:,spi]
     #print(spc_aqi[88,77,spi])
     #searchval = 3
     #xx,yy=np.where(vals==3)
     #print('aq=',searchval, xx.shape, yy.shape,'for species number', spi, 'name', sp, spc_ave_conv[xx,yy])
     spi=spi+1
     #spc_aqi=array with aqi of each species in a dimension
   #Total AQI = max of each spc_aqi array along the last axis
   tot_aq=np.amax(spc_aqi, axis=2)
   #xxx,yyy=np.where(tot_aq==10)
   #print(xxx,yyy)
   #print('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTtotal aqi is 10', xxx.shape, yyy.shape, tot_aq[xxx,yyy])
   
   ##print(tot_aq)
   
   ##print('Plotting aqi')
   png_out_path=(dir_png_json_op+'/'+datestr+':00-aqi.png')
   
   fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   ax=fig.add_subplot(111, projection=ccrs.Mercator())
   lev_range=np.arange(0,7,1)
   levels=lev_range
   bounds=np.linspace(0,7,1)
   norm = mc.BoundaryNorm(bounds, cm.N)
   print(tot_aq.shape)
   #cs=ax.contourf(lons, lats, tot_aq, levels=levels, cmap=cm, extend="max", spacing='proportional', ticks=levels, transform=ccrs.PlateCarree())
   #won't work with transform: look at https://github.com/SciTools/cartopy/pull/885 to fix!
   cs=ax.contourf(lons, lats, tot_aq[:,:], levels=levels, cmap=cm, extend="max",  transform=ccrs.PlateCarree())
   ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
   plt.box(on=None)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
   #ax.coastlines('10m', linewidth=0.15)
   plt.axis('off')
   ax.figsize=(tilesize/dpi, tilesize/dpi)
   ax.dpi=dpi
   ax.outline_patch.set_visible(False)
   ax.background_patch.set_visible(False)
   #ax.background_patch.set_alpha(0)
   #ax.add_feature(cfeature.BORDERS,linewidth=0.25)
   ax.axes.get_xaxis().set_visible(False)
   ax.axes.get_yaxis().set_visible(False)
   cs=cm.set_over('m')
   #cbar = plt.colorbar(shrink=0.9, orientation='horizontal')
   plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
   #plt.close()
   
   
   
   #output the dataset as an xarray
   #aqi_out_path=(dir_aqi+'/'+datestr+':00-aqi.nc')
   #out=xr.Dataset()
   #out.coords['time'] = dttf
   #out.coords['x'] = range(data_ref['lon'].shape[0])
   #out.coords['y'] = range(data_ref['lat'].shape[0])
   ##prin
   #out['lat'] = (('y', 'x'), lats)
   #out['lon'] = (('y', 'x'), lons)
   #out['aqi'] = (('y', 'x'), tot_aq[:])
   #out.to_netcdf(aqi_out_path)
   
   # make jsons
   #nx = 162
   #ny = 114
   #dx = 0.25
   #dy = 0.25
   #la1 = 34.5
   #la2 = 62.75
   #lo1 = -21.75
   #lo2 = 18.5
   #gridsize=162*114
   #longs=np.arange(lo1,lo2+dx,dx)
   #latit=np.arange(la1, la2+dy,dy)
   #times = newf.variables['Times']
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   dirpath = dir_png_json_op
   
   outpath = os.path.join(dirpath, '%s-aqi.json' % (timestr,))
   #gc=0
   #outf = open(outpath, 'w')
   #outf.write('{ ')

   
   coords_dict = {}
   tot_aq_dict = {}
   #for xways in range(0, nx):
     #for yways in range(0,ny):
       #strlat=str(latit[yways])
       #strlon=str(longs[xways])
       ##print(tot_aq.shape, rawdatavals.shape)
       #rawdatavals[yways, xways]=tot_aq[yways,xways]
       #rounddatavals=int(round(rawdatavals[yways, xways], 0))
       #strval=str(rounddatavals)
       #gc=gc+1
       #if gc==gridsize:
         #strtowrite=(' "'+strlat+', '+strlon+'": "'+strval+'" }')
       #else:
         #strtowrite=(' "'+strlat+','+strlon+'": "'+strval+'",')
       #outf.write(strtowrite)    
   #outf.close()
   #fpath=os.path.join(dir_op, file)
   #fpath_bk=os.path.join(dir_op_back, file)
   #print('moving from', fpath, fpath_bk)
   #shutil.move(fpath, fpath_bk)
   
   for i, lonv in enumerate(data_ref.lon.values):
       for j, latv in enumerate(data_ref.lat.values):
           coords_string = (str(latv)+','+str(lonv))
           coords_dict[j,i] = coords_string
           tot_aq_dict[j,i]= tot_aq[j,i]
           
   coords_df =pd.DataFrame.from_dict(coords_dict, orient='index')
   tot_aq_df = pd.DataFrame.from_dict(tot_aq_dict, orient='index')
   tot_aq_df['tot_aq']=tot_aq_df.values
   tot_aq_df = tot_aq_df.round({'tot_aq': 2})
   tot_aq_df['coords'] = coords_df
   text_output = tot_aq_df.set_index('coords').to_dict()['tot_aq']
   with open(outpath, 'w') as fp:
       json.dump(text_output, fp)
   print('moving from', fpath, fpath_bk)
   shutil.move(fpath, fpath_bk)
   
     
     
     #create dynamic variable name for the array "comp" (AQI for individual species)
     
     
     
   #move the file to back directory)
   #newfpath=os.path.join(dir_op_back, file)
   #os.rename(fpath, newfpath)
