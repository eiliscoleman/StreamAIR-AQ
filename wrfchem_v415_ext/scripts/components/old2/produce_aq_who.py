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
import numpy.ma as ma
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
bboundsugm3=(100,200,100,25,50)
time_ave=(8,1,1,24,24)
    

listwho={
  'spcs':pd.Series(species, index=np.arange(0,len(species))),
  'bnds':pd.Series(bboundsugm3, index=np.arange(0,len(species))),
  'time_ave':pd.Series(time_ave, index=np.arange(0,len(species))),
}
listwho2={
  #'spcs':pd.Series(species, index=np.arange(0,len(species))),
  'bnds':pd.Series(bboundsugm3, index=species),
  'time_ave':pd.Series(time_ave, index=species),
}
#print(listwho2.keys())

df_who=pd.DataFrame(data=listwho)

df_who2=pd.DataFrame(data=listwho2)


##Create colormap
#cm=LinearSegmentedColormap.from_list('aqi_cmap', aqi_colours, N=len(aqi_colours))

##print(listAQI['spcs']['o3']['boundsugm3'])
for file in sorted(os.listdir(dir_op)):
   fpath=os.path.join(dir_op, file)
   fpath_bk=os.path.join(dir_op_back, file)
   shutil.copyfile(fpath, fpath_bk)
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
   exc=np.empty([nolats,nolons])
   ex_mx=np.empty([nolats,nolons,len(species)])
   spc_aq_who=np.empty([nolats,nolons,len(species)])
   spi=0
   for  spi, sp in enumerate(listwho['spcs']):
     print('species is', sp)
     print('bounds', df_who2.loc[sp]['bnds'])
     ave_period=df_who2.loc[sp]['time_ave']
     print('hours', ave_period)
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
         conv= 1995.7 # 1000ppb = 1 ppm, 1 ppb = 2.00 ug/m3
       elif sp == 'nox':
         spc_hrtmp[:,:,hour]=datatmp.nox_concentration[0,0,:,:]
         conv= 1912.5 # 1000ppb = 1 ppm, 1 ppb = 1.9125 ug/m3   
       elif sp == 'so2':
         
         spc_hrtmp[:,:,hour]=datatmp.so2_concentration[0,0,:,:]
         conv= 2660.9 # 1000ppb = 1 ppm, 1 ppb = 2.6609 ug/m3         
       elif sp == 'pm25':
              
         spc_hrtmp[:,:,hour]=datatmp.pm25[0,0,:,:]
         #print(spc_hrtmp[88,77,hour], hour)
         conv= 1 # already in ug/m3
       elif sp == 'pm10':
         spc_hrtmp[:,:,hour]=datatmp.pm10[0,0,:,:]         
         conv= 1 # already in ug/m3
     
     #get average over the timespan of each array
     spc_ave=np.mean(spc_hrtmp, axis=2)
     # below: the time averaged concentration in ug/m3
     spc_ave_conv=spc_ave*conv
     #print('spc_ave', spc_ave_conv[0:100])

     rexc, cexc=np.where(spc_ave_conv>df_who2.loc[sp]['bnds'])
     print(len(rexc), len(cexc))
     # this is the martrix that'll be output to give the exceedance values
     ex_mx[rexc,cexc,spi]=((spc_ave_conv[rexc, cexc])-df_who2.loc[sp]['bnds'])*100./df_who2.loc[sp]['bnds']
     for lines, rows in enumerate(rexc):
         #print(lines)
         #print(rexc[lines], cexc[lines])
         print(spc_ave_conv[rexc[lines], cexc[lines]], '%',ex_mx[rexc[lines],cexc[lines],spi] )
      # make jsons
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
     dirpath = dir_png_json_op
     opspfilename=timestr+'-'+sp+'-aq-who-wrf-chem-3.9.json'
     print(opspfilename)
     outpath = os.path.join(dirpath,opspfilename)
     gc=0
     outf = open(outpath, 'w')
     outf.write('{ ')
     rawdatavals=np.zeros([ny,nx])
     for xways in range(0, nx):
      for yways in range(0,ny):
        strlat=str(latit[yways])
        strlon=str(longs[xways])
        #print(tot_aq.shape, rawdatavals.shape)
        rawdatavals[yways, xways]=ex_mx[yways,xways, spi]
        #rawdatavals[yways, xways]=tot_aq[yways,xways]
        rounddatavals=round(rawdatavals[yways, xways], 1)
        strval=str(rounddatavals)
        gc=gc+1
        if gc==gridsize:
          strtowrite=(' "'+strlat+', '+strlon+'": "'+strval+'" }')
        else:
          strtowrite=(' "'+strlat+','+strlon+'": "'+strval+'",')
        outf.write(strtowrite)    
     outf.close()
     #fpath=os.path.join(dir_op, file)
     #fpath_bk=os.path.join(dir_op_back, file)
     #print('moving from', fpath, fpath_bk)
     #shutil.move(fpath, fpath_bk)
   
     
     
     #create dynamic variable name for the array "comp" (AQI for individual species)
     
     
     
   #move the file to back directory)
   #newfpath=os.path.join(dir_op_back, file)
   #os.rename(fpath, newfpath)

    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     #print(np.where(exc))
     #for point in exc:
         #print('exceedance in', sp, 'at', point)
     ##print(sp, spc_ave_conv.shape, spc_ave_conv[88,77])
     #lims=listAQI[sp]['boundsugm3']
     
     ##np.digitize returns indices of specified bins difined in boundsugm3
     
     #spc_aqi[:,:,spi]=np.digitize(spc_ave_conv, lims, right=True)
     #vals=spc_aqi[:,:,spi]
     ##print(spc_aqi[88,77,spi])
     ##searchval = 3
     ##xx,yy=np.where(vals==3)
     ##print('aq=',searchval, xx.shape, yy.shape,'for species number', spi, 'name', sp, spc_ave_conv[xx,yy])
     #spi=spi+1
     ##spc_aqi=array with aqi of each species in a dimension
   ##Total AQI = max of each spc_aqi array along the last axis
   #tot_aq=np.amax(spc_aqi, axis=2)
   ##xxx,yyy=np.where(tot_aq==10)
   ##print(xxx,yyy)
   ##print('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTtotal aqi is 10', xxx.shape, yyy.shape, tot_aq[xxx,yyy])
   
   ###print(tot_aq)
   
   ###print('Plotting aqi')
   #png_out_path=(dir_png_json_op+'/'+datestr+':00-aqi.png')
   
   #fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   #ax=fig.add_subplot(111, projection=ccrs.Mercator())
   #lev_range=np.arange(0,11,1)
   #levels=lev_range
   #bounds=np.linspace(0,11,1)
   #norm = mc.BoundaryNorm(bounds, cm.N)
   #print(tot_aq.shape)
   ##cs=ax.contourf(lons, lats, tot_aq, levels=levels, cmap=cm, extend="max", spacing='proportional', ticks=levels, transform=ccrs.PlateCarree())
   ##won't work with transform: look at https://github.com/SciTools/cartopy/pull/885 to fix!
   #cs=ax.contourf(lons, lats, tot_aq, levels=levels, cmap=cm, extend="max",  transform=ccrs.PlateCarree())
   #plt.box(on=None)
   #plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
   ##ax.coastlines('10m', linewidth=0.15)
   #plt.axis('off')
   #ax.figsize=(tilesize/dpi, tilesize/dpi)
   #ax.dpi=dpi
   #ax.outline_patch.set_visible(False)
   #ax.background_patch.set_visible(False)
   ##ax.background_patch.set_alpha(0)
   ##ax.add_feature(cfeature.BORDERS,linewidth=0.25)
   #ax.axes.get_xaxis().set_visible(False)
   #ax.axes.get_yaxis().set_visible(False)
   #cs=cm.set_over('m')
   ##cbar = plt.colorbar(shrink=0.9, orientation='horizontal')
   #plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
   #plt.close()

   
   ## make jsons
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
   ##times = newf.variables['Times']
   #timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   #dirpath = dir_png_json_op
   
   #outpath = os.path.join(dirpath, '%s-aqi-wrf-chem-3.5.json' % (timestr,))
   #gc=0
   #outf = open(outpath, 'w')
   #outf.write('{ ')

   #rawdatavals=np.zeros([ny,nx])
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
   ##fpath=os.path.join(dir_op, file)
   ##fpath_bk=os.path.join(dir_op_back, file)
   #print('moving from', fpath, fpath_bk)
   #shutil.move(fpath, fpath_bk)
   
     
     
     #create dynamic variable name for the array "comp" (AQI for individual species)
     
     
     
   #move the file to back directory)
   #newfpath=os.path.join(dir_op_back, file)
   #os.rename(fpath, newfpath)
