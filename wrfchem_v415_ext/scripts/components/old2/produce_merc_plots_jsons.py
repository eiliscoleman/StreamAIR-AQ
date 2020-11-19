from __future__ import print_function, unicode_literals

import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import simplejson as json
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import xarray as xr
#import xarray.ufuncs as xu
#import seaborn as sns
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
print(LON.shape)
print(data.lon.shape)
#coords=np.empty(LON.shape)
coords={}
coords_dict = {}
t2_dict={}
t2_dict2={}
t2_dat=data.t2.values[0,:,:]-273.15
#for i, lonv in enumerate(data.lon.values):
    ##print(i, lonv)
    #for j, latv in enumerate(data.lat.values):

        #stri=(str(latv)+','+str(lonv))
        ##print(j*i)
        #coords[j,i]=stri
        #t22[j,i]=data.t2.values[0,j,i]-273.15
#print(type(coords))        
#co_p=pd.DataFrame.from_dict(coords, orient='index')
#t2_0=pd.DataFrame.from_dict(t22, orient='index')
#t2_0.columns=['t2']
#print(t2_0.shape)
#t2_0['coords']=co_p
#t2_new=t2_0.set_index('coords')
#col_list=list(t2_0)
#col_list[0], col_list[1]=col_list[1], col_list[0]
#print(col_list)
#t222=t2_0.reindex(columns=col_list)
##t2_0.set_index('coords', inplace=True)
##print(t2_0.index)
#print(t222.keys)
#t222.set_index('coords', inplace=True)
json_path=(png_out_dir+'/'+datestr+':00-t2.json')
#export_df_lat=t222.to_json(json_path,orient='table')


for i, lonv in enumerate(data.lon.values):
    for j, latv in enumerate(data.lat.values):
        coords_string = (str(latv)+','+str(lonv))
        coords_dict[j,i] = coords_string
        t2_dict[j,i] = data.t2.values[0,j,i]-273.15
#for key, value in t2_dict.items() :
    #print(key, value)
    #t2_dict2
    
    
#coords_df = pd.DataFrame.from_dict(coords_dict, orient='index')
#t2_df = pd.DataFrame.from_dict(t2_dict, orient='index', columns=["t2"])
#t2_df = t2_df.round({'t2': 2})
#t2_df['coords'] = coords_df

#text_output = t2_df.set_index('coords').rename(columns={'t2':'t2'}).to_dict()['t2']

#with open(json_path, 'w') as fp:
    #json.dump(text_output, fp)

#############################
coords_df = pd.DataFrame.from_dict(coords_dict, orient='index')
#index = [0]
t2_df = pd.DataFrame.from_dict(t2_dict, orient='index')
print(t2_df.columns)
t2_df['t2']=t2_df.values
t2_df = t2_df.round({'t2': 2})
t2_df['coords'] = coords_df

text_output = t2_df.set_index('coords').to_dict()['t2']
with open(json_path, 'w') as fp:
    # works, except the whole thing is led by {0: {array we want}}
    json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
    #below commamd almost works!, except  each value is preceded by the column name {0}
    #export_df=text_output_df.to_json(json_path, orient='index')
    #export_df=text_output_df.to_json(json_path, orient='values')
##################################################
##df_t2=pd.DataFrame(t2_dat, index=data.lat.values, columns=data.lon.values)
#print(co_p.values)
#print(co_p.values.shape)
#df_t2=pd.DataFrame(t2_dat, index=co_p.values)
#print(df_t2.head())
#, index=co_p.values)
#df_t2_2=df_t2.stack()
#df_t2['coords']=pd.DataFrame(co_p)
#print(df_t2.keys)
#df_t2.set_index('coords', inplace = True)
#print(df_t2.head())
#ix=pd.DataFrame(co_p)


#export_df=df_t2.to_json(json_path)
#export_df_lon=df_lon.to_json(json_path_lon)
#export_df_lat=df_lat.to_json(json_path_lat)
#, orient='index')
#, orient='index')
#with open(json_path, 'w') as outfile:
    #json.dump(t2_dat[:], outfile)
#with open(json_lon_path, 'w') as outfile:
    #json.dump(LON, outfile)
#with open(json_lat_path, 'w') as outfile:
    #json.dump(LAT, outfile)
## Wind
#print ('Wind' )
#wind_out_path=(png_out_dir+'/'+datestr+':00-wind-10m.png')
#fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#ax=fig.add_subplot(111, projection=ccrs.Mercator())
#lev_range=np.arange(0,30,0.5)
#levels=lev_range
#wind10=xu.sqrt(data.u10.values[0,:,:]**2+data.v10.values[0,:,:]**2)
#cs=ax.contourf(LON, LAT, wind10, levels, cmap=plt.cm.jet, extend="min", transform=ccrs.PlateCarree())
#plt.box(on=None)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
##ax.coastlines('10m', linewidth=0.15)
#plt.axis('off')
#ax.figsize=(tilesize/dpi, tilesize/dpi)
#ax.dpi=dpi
#ax.outline_patch.set_visible(False)
#ax.outline_patch.set_visible(False)
#ax.background_patch.set_visible(False)
#ax.patch.set_alpha(0)
##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)
#plt.savefig(wind_out_path, dpi=288, tilesize=768, transparent=True)
#plt.close()
 

 

##print ('RAIN and pressure')
##out_path=(png_out_dir+'/test_rain_psl.png')
##plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
##plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
##plt.axis('off')
###bounds=[0, .01, .1, .3, .5, .7, 1, 2, 3, 6, 12, 20, 50, 100]
##bounds=[0, .1, .3, .5, .7, 1, 2, 3,  5, 6, 12]
###lev_range=np.logspace(0,2, num=8)  
##lev_range=np.array(bounds)
##levels=lev_range
##palette=plt.cm.get_cmap('cubehelix_r')
##palette.set_over('r')
##cs=plt.contourf(LON, LAT, data.rain.values[0,:,:], levels, cmap=palette, extend="max")
##cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n=20, linewidths=0.5, colors='k')
##levels_n=np.arange(1000,1040,2)
###cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n, linewidths=0.2, colors='k')
##plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
##plt.savefig(out_path, dpi=1152, tilesize=768)
##plt.close()

####                                     Plot AQI
##print('aqi')
##rh_png_out_path=(png_out_dir+'/'+datestr+':00-aqi.png')
##plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
##plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
##plt.axis('off')
###bounds=np.arange[0, 11, 1]
##lev_range=np.arange(0,11,1)
##levels=lev_range
##aqiraw=data.rh.values[0,:,:]
##aqiscaled=aqiraw[:,:]*0.1
##aqiround=np.around([aqiscaled], decimals=0)
##print(aqiround.shape)
##plt.contourf(LON, LAT, aqiround[0,:,:], levels, cmap=plt.cm.viridis, extend="both")
##plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
##plt.close()


####                                     Plot Rh
#print('Rh')
#rh_png_out_path=(png_out_dir+'/'+datestr+':00-rh.png')
#plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
#plt.axis('off')
#lev_range=np.arange(30,101,1)
#levels=lev_range
#fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#ax=fig.add_subplot(111, projection=ccrs.Mercator())
#cs=ax.contourf(LON, LAT,data.rh.values[0,:,:], levels, cmap=plt.cm.get_cmap('gist_rainbow_r'), extend="both", transform=ccrs.PlateCarree())
#plt.box(on=None)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
##ax.coastlines('10m', linewidth=0.15)
#plt.axis('off')
#ax.figsize=(tilesize/dpi, tilesize/dpi)
#ax.dpi=dpi
#ax.outline_patch.set_visible(False)
#ax.background_patch.set_visible(False)
#ax.patch.set_alpha(0)
##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)
##plt.contourf(LON, LAT, data.t2.values[0,:,:]-273.15, levels, cmap=plt.cm.get_cmap('jet'), extend="both")
#plt.savefig(rh_png_out_path, dpi=288, tilesize=768, transparent=True)
#plt.close()
##The following plots need to be plotted for all levels
##                                     Plot O3
##for lev in range(0,29):
#for lev in range(0,1):

  #levstr=str(lev)
  #levpadded=levstr.zfill(2)
  #print(levstr)
  #print('O3', data.o3_concentration.shape, lev )
##  rh_png_out_path=(png_out_dir+'/'+datestr+'-o3-lev'+levpadded+'.png')
  #png_out_path=(png_out_dir+'/'+datestr+':00-o3.png')
  #lev_range=np.arange(30, 70, .2)
  #levels=lev_range
  ##plt.contourf(LON, LAT, data.o3_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.jet, extend="min")
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator())
  #cs=ax.contourf(LON, LAT,data.o3_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.jet, extend="both", transform=ccrs.PlateCarree())
  #plt.box(on=None)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
  #plt.axis('off')
  #ax.figsize=(tilesize/dpi, tilesize/dpi)
  #ax.dpi=dpi
  #ax.outline_patch.set_visible(False)
  #ax.background_patch.set_visible(False)
  #ax.patch.set_alpha(0)
  ##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  #ax.axes.get_xaxis().set_visible(False)
  #ax.axes.get_yaxis().set_visible(False)
  #plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  #plt.close()
##                                    Plot NOx
  #print('NOx')
  #png_out_path=(png_out_dir+'/'+datestr+':00-nox.png')
  #lev_range=np.arange(1, 60, 1)
  #levels=lev_range
  ##plt.contourf(LON, LAT, data.nox_concentration.values[0,lev,:,:]*1000., levels  , cmap=plt.cm.gnuplot, extend="min")
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator())
  #cs=ax.contourf(LON, LAT,data.nox_concentration.values[0,lev,:,:]*1000., levels  , cmap=plt.cm.gnuplot, extend="min", transform=ccrs.PlateCarree())
  #plt.box(on=None)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
  #plt.axis('off')
  #ax.figsize=(tilesize/dpi, tilesize/dpi)
  #ax.dpi=dpi
  #ax.outline_patch.set_visible(False)
  #ax.background_patch.set_visible(False)
  #ax.patch.set_alpha(0)
  ###ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  #ax.axes.get_xaxis().set_visible(False)
  #ax.axes.get_yaxis().set_visible(False)
  #plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  #plt.close()
##                 
 ##Plot SO2

  #print('SO2')
  #png_out_path=(png_out_dir+'/'+datestr+':00-so2.png')
  #lev_range=np.arange(0,8,.002)
  #levels=np.linspace(0, 100, num=11)

  ##plt.contourf(LON, LAT, data.so2_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min")
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator())
  #cs=ax.contourf(LON, LAT,data.so2_concentration.values[0,lev,:,:]*1000., levels, cmap=plt.cm.inferno, extend="min", transform=ccrs.PlateCarree())
  #plt.box(on=None)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
  #plt.axis('off')
  #ax.figsize=(tilesize/dpi, tilesize/dpi)
  #ax.dpi=dpi
  #ax.outline_patch.set_visible(False)
  #ax.background_patch.set_visible(False)
  #ax.patch.set_alpha(0)
  ##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  #ax.axes.get_xaxis().set_visible(False)
  #ax.axes.get_yaxis().set_visible(False)
  #plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  #plt.close()

###                                    Plot PM10
  ##print('PM10')
  ##png_out_path=(png_out_dir+'/'+datestr+':00-pm10.png')
 
  ##lev_range=np.arange(0,20,.2)
  ##levels=lev_range
  ###plt.contourf(LON, LAT, data.pm10.values[0,lev,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  ##fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  ##ax=fig.add_subplot(111, projection=ccrs.Mercator())
  ##cs=ax.contourf(LON, LAT,data.pm10.values[0,lev,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both", transform=ccrs.PlateCarree())
  ##plt.box(on=None)
  ##plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ###ax.coastlines('10m', linewidth=0.15)
  ##plt.axis('off')
  ##ax.figsize=(tilesize/dpi, tilesize/dpi)
  ##ax.dpi=dpi
  ##ax.outline_patch.set_visible(False)
  ##ax.background_patch.set_visible(False)
  ##ax.patch.set_alpha(0)
  ###ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  ##ax.axes.get_xaxis().set_visible(False)
  ##ax.axes.get_yaxis().set_visible(False)
  ##plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  ##plt.close()
  ### change the extended PM plots to the main output
  ##                                    Plot PM10_30
  #print('PM10_30')
  #png_out_path=(png_out_dir+'/'+datestr+':00-pm10.png')
 
  #lev_range=np.arange(0,30.2,.2)
  #levels=lev_range
  ##plt.contourf(LON, LAT, data.pm10.values[0,lev,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator())
  #cs=ax.contourf(LON, LAT,data.pm10.values[0,lev,:,:], levels, cmap=plt.cm.get_cmap('magma'), extend="both", transform=ccrs.PlateCarree())
  #plt.box(on=None)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
  #plt.axis('off')
  #ax.figsize=(tilesize/dpi, tilesize/dpi)
  #ax.dpi=dpi
  #ax.outline_patch.set_visible(False)
  #ax.background_patch.set_visible(False)
  #ax.patch.set_alpha(0)
  ##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  #ax.axes.get_xaxis().set_visible(False)
  #ax.axes.get_yaxis().set_visible(False)
  #plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  #plt.close()

###                                   Plot PM25 
  ##print('PM25')
  ##png_out_path=(png_out_dir+'/'+datestr+':00-pm25.png')
  ##lev_range=np.arange(0,20,.2)
  ##levels=lev_range
  
  ###plt.contourf(LON, LAT, data.pm25.values[0,lev,:,:],  levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  ##fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  ##ax=fig.add_subplot(111, projection=ccrs.Mercator())
  ##cs=ax.contourf(LON, LAT,data.pm25.values[0,lev,:,:],  levels, cmap=plt.cm.get_cmap('magma'), extend="both", transform=ccrs.PlateCarree())
  ##plt.box(on=None)
  ##plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ###ax.coastlines('10m', linewidth=0.15)
  ##plt.axis('off')
  ##ax.figsize=(tilesize/dpi, tilesize/dpi)
  ##ax.dpi=dpi
  ##ax.outline_patch.set_visible(False)
  ##ax.background_patch.set_visible(False)
  ##ax.patch.set_alpha(0)
  ###ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  ##ax.axes.get_xaxis().set_visible(False)
  ##ax.axes.get_yaxis().set_visible(False)
  ###plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  ##plt.savefig(png_out_path, dpi=dpi, tilesize=tilesize, transparent=True)
  ##plt.close()
  
  ##                                   Plot PM25_30 
  #print('PM25_30')
  #png_out_path=(png_out_dir+'/'+datestr+':00-pm25.png')
  #lev_range=np.arange(0,30.2,.2)
  #levels=lev_range
  
  ##plt.contourf(LON, LAT, data.pm25.values[0,lev,:,:],  levels, cmap=plt.cm.get_cmap('magma'), extend="both")
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator())
  #cs=ax.contourf(LON, LAT,data.pm25.values[0,lev,:,:],  levels, cmap=plt.cm.get_cmap('magma'), extend="both", transform=ccrs.PlateCarree())
  #plt.box(on=None)
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
  #plt.axis('off')
  #ax.figsize=(tilesize/dpi, tilesize/dpi)
  #ax.dpi=dpi
  #ax.outline_patch.set_visible(False)
  #ax.background_patch.set_visible(False)
  #ax.patch.set_alpha(0)
  ##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
  #ax.axes.get_xaxis().set_visible(False)
  #ax.axes.get_yaxis().set_visible(False)
  ##plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
  #plt.savefig(png_out_path, dpi=dpi, tilesize=tilesize, transparent=True)
  #plt.close()
  
  ##Plot Black lines
  #print('Prs isobars - bk')
  #png_out_path=(png_out_dir+'/'+datestr+':00-psl_black.png')
  #svg_out_path=(png_out_dir+'/'+datestr+':00-psl_black.svg')
  #levels_n=np.arange(900,1240,5)
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator(), facecolor='red')
  #fig.patch.set_alpha(0)
  #cs=ax.contour(LON, LAT,data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='k', alpha=1, transform=ccrs.PlateCarree())
  #plt.clabel(cs, inline=1, fmt='%.0f', fontsize=2, color='k')
  #plt.box(on=None)
  ## we need this line!
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
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
  #fig.savefig(png_out_path, transparent=True)
  #fig.savefig(svg_out_path, format='svg', bbox_inches='tight', pad_inches=0, dpi=dpi, tilesize=tilesize, transparent=True)
  #plt.close()
  
  
  ##Plot Black lines
  #print('Prs isobars - bk - finer')
  ##png_out_path=(png_out_dir+'/'+datestr+':00-psl_black.png')
  #svg_out_path=(png_out_dir+'/'+datestr+':00-psl_black_hires.svg')
  #levels_n=np.arange(900,1240,2)
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator(), facecolor='red')
  #fig.patch.set_alpha(0)
  #cs=ax.contour(LON, LAT,data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.1, colors='k', alpha=1, transform=ccrs.PlateCarree())
  #plt.clabel(cs, inline=1, fmt='%.0f', fontsize=1, color='k')
  #plt.box(on=None)
  ## we need this line!
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
  ##ax.coastlines('10m', linewidth=0.15)
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
  ##fig.savefig(png_out_path, transparent=True )
  #fig.savefig(svg_out_path, format='svg', bbox_inches='tight', pad_inches=0, dpi=24, tilesize=64, transparent=True)
  #plt.close()
  
  

  ##Plot white lines
  #print('Prs isobars - wh')
  #png_out_path=(png_out_dir+'/'+datestr+':00-psl_white.png')
  #levels_n=np.arange(900,1240,5)
  #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
  #ax=fig.add_subplot(111, projection=ccrs.Mercator(), facecolor='red')
  #fig.patch.set_alpha(0)
  #cs=ax.contour(LON, LAT,data.p_sl[0,:,:]*1e-2, levels=levels_n, linewidths=0.2, colors='w', alpha=1, transform=ccrs.PlateCarree())
  #plt.clabel(cs, inline=1, fmt='%.0f', fontsize=2, color='w')
  #plt.box(on=None)
  ## we need this line!
  #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
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
