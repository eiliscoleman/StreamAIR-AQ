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
dir_op2=(basdir+'/data/output-wrf-remapped2')
dir_png_json_op=(basdir+'/web')
#dir_aqi=('/mnt/raid/rong-ming/wrfchem/data/output-wrf-aqi')
# Create list of AQI for various pollutants
species=['o3', 'nox', 'so2', 'pm25', 'pm10']
#bboundsugm3=(100,200,20,25,50)
# bboundsppb=(50,106,7.6,25,50)
bboundsugm3=(100,106,7.6,25,50)
time_ave=(8,1,24,24,24)
#time_ave=(1,1,1,1,1)
lst_cmaps=('Greens','gist_yarg','bone_r','pink_r','afmhot_r')
    

listwho={
  'spcs':pd.Series(species, index=np.arange(0,len(species))),
  'bnds':pd.Series(bboundsugm3, index=np.arange(0,len(species))),
  'time_ave':pd.Series(time_ave, index=np.arange(0,len(species))),
}
listwho2={
  #'spcs':pd.Series(species, index=np.arange(0,len(species))),
  'bnds':pd.Series(bboundsugm3, index=species),
  'time_ave':pd.Series(time_ave, index=species),
  'ls_cm':pd.Series(lst_cmaps, index=species),
}
df_who=pd.DataFrame(data=listwho)

df_who2=pd.DataFrame(data=listwho2)


##Create colormap
#cm=LinearSegmentedColormap.from_list('aqi_cmap', aqi_colours, N=len(aqi_colours))

##print(listAQI['spcs']['o3']['boundsugm3'])
for file in sorted(os.listdir(dir_op)):
   fpath=os.path.join(dir_op, file)
   fpath2=os.path.join(dir_op2, file)
   print(fpath)
   fpath_bk=os.path.join(dir_op_back, file)
   #print(fpath_bk)
   shutil.copyfile(fpath, fpath_bk)
   ##print('file to be opened', fpath)
   data_ref=xr.open_dataset(fpath)
   nolats, nolons=len(data_ref['lat'][:]), len(data_ref['lon'][:])
   lats = data_ref.variables['lat'][:]
   lons = data_ref.variables['lon'][:]
   datestamp=fpath.split("reg_",1)[1]
   yst=datestamp[:4]
   mst=datestamp[4:-4]
   dst=datestamp[6:-2]
   hrst=datestamp[8:]
   dttf=datetime.strptime(datestamp, "%Y%m%d%H")
   datestr=dttf.strftime("%Y-%m-%d-%H:00")
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
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   #print(datestr)
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
   ex_mx=np.empty([nolats,nolons,len(species)])
   tot_ex_mx=np.empty([nolats,nolons])
   spi=0
   for  spi, sp in enumerate(listwho['spcs']):
     ave_period=df_who2.loc[sp]['time_ave']
     spc_hrtmp=np.empty([nolats,nolons,ave_period])
     opspfilename=timestr+'-'+sp+'-ave-who-wrf-chem-3.9.json'
     dirpath = dir_png_json_op
     outpath = os.path.join(dirpath,opspfilename)
     for hour in range(ave_period):
       
       dtti=dttf+timedelta(hours=-hour)
       filedate=dtti.strftime("%Y%m%d%H")
       filepath=(dir_op_back+'/d01reg_'+filedate)
       datatmp=xr.open_dataset(filepath)       
       spc_ave_conv=np.empty([nolats,nolons])       
       if sp == 'o3':   
         #print(datatmp.o3_concentration.shape)
         spc_hrtmp[:,:,hour]=datatmp.o3_concentration[0,0,:,:]
         #factor to convert ppm to ug/m3
         conv= 1000 # 1000ppb = 1 ppm, 1 ppb = 2.00 ug/m3
       elif sp == 'nox':
         spc_hrtmp[:,:,hour]=datatmp.nox_concentration[0,0,:,:]
         conv= 1000 # 1000ppb = 1 ppm, 1 ppb = 1.9125 ug/m3   
       elif sp == 'so2':
         spc_hrtmp[:,:,hour]=datatmp.so2_concentration[0,0,:,:]
         conv= 1000 # 1000ppb = 1 ppm, 1 ppb = 2.6609 ug/m3         
       elif sp == 'pm25':
              
         spc_hrtmp[:,:,hour]=datatmp.pm25[0,0,:,:]
         #print(spc_hrtmp[88,77,hour], hour)
         conv= 1 # already in ug/m3
       elif sp == 'pm10':
         spc_hrtmp[:,:,hour]=datatmp.pm10[0,0,:,:]         
         conv= 1 # already in ug/m3
       datatmp.close()
       
     #get average over the timespan of each array
     spc_ave=np.empty([nolats,nolons])
     spc_ave=np.mean(spc_hrtmp, axis=2)
     # below: the time averaged concentration in ug/m3
     spc_ave_conv=np.empty([nolats,nolons])
     spc_ave_conv=spc_ave*conv
     #print('spc_ave', spc_ave_conv[0:100])
     #rexc and cexc are problematic!
     rexc, cexc=np.empty([0]), np.empty([0])
     rexc, cexc=np.where(spc_ave_conv>df_who2.loc[sp]['bnds'])
     ##print(len(species))
     ##print('shape, yakno', spi, ex_mx.shape)
     #ex_mx[rexc,cexc,spi]=((spc_ave_conv[rexc, cexc])-df_who2.loc[sp]['bnds'])*100./df_who2.loc[sp]['bnds']
     bndconstant=df_who2.loc[sp]['bnds']
     ex_mx[:,:,spi]=100.*(spc_ave_conv[:, :]-bndconstant)/bndconstant
     gc=0
     outf = open(outpath, 'w')
     outf.write('{ ')
     rawdatavals=np.zeros([ny,nx], dtype=np.float128)
     for xways in range(0, nx):
           for yways in range(0,ny):
               strlat=str(latit[yways])
               strlon=str(longs[xways])
               rawdatavals[yways, xways]=spc_ave_conv[yways,xways]
               #print(spc_ave[yways,xways])
               rounddatavals=round(rawdatavals[yways, xways], 1)
               strval=str(rounddatavals)
               gc=gc+1
               if gc==gridsize:
                   strtowrite=(' "'+strlat+', '+strlon+'": "'+strval+'" }')
               else:
                   strtowrite=(' "'+strlat+','+strlon+'": "'+strval+'",')
               outf.write(strtowrite) 
     outf.close()
     #shutil.move(fpath, fpath2)
     
     
     
     ##print(df_who2.loc[sp]['bnds'])
     #print(len(rexc), len(cexc))
   tot_ex_mx=np.amax(ex_mx, axis=2)
   #print(tot_ex_mx.shape, sp)
   #tot_ex_mx=ex_mx[:,:,1]
   ##print(tot_ex_mx.shape, sp)
   png_out_path=(dir_png_json_op+'/'+datestr+':00-who_tot_ex.png')
   ##print(png_out_path)
   fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
   ax=fig.add_subplot(111, projection=ccrs.Mercator())
   lev_range=np.arange(0,210,10)
   strcmap=df_who2.loc[sp]['ls_cm']
   strcmap='afmhot_r'
   ##print(strcmap)
   #extent = mtransforms.Bbox.union([col.get_datalim(self.transData) for col in result.collections] if col.get_paths())
   ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
   #cs=ax.contourf(lons, lats,  ex_mx[:,:,2], levels=lev_range,  cmap=plt.cm.get_cmap(strcmap),  transform=ccrs.PlateCarree(),extend='both')
   cs=ax.contourf(lons, lats,  tot_ex_mx[:,:], levels=lev_range,  cmap=plt.cm.get_cmap(strcmap),  transform=ccrs.PlateCarree(),extend='both')
   ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
   plt.box(on=None)
   plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
   #ax.coastlines('10m', linewidth=0.15)
   plt.axis('off')
   ax.figsize=(tilesize/dpi, tilesize/dpi)
   ax.dpi=dpi
   ax.outline_patch.set_visible(False)
   
   ax.background_patch.set_visible(False)
   ax.axes.get_xaxis().set_visible(False)
   ax.axes.get_yaxis().set_visible(False)
   plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=True)
   plt.close("all")
   plt.cla()
   plt.clf()
   plt.close(fig)
   #data_ref.close()
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
   timestr = (dttf.strftime("%Y-%m-%d-%H")+':00:00')
   dirpath = dir_png_json_op
   opspfilename=timestr+'-exc-aq-who-wrf-chem-3.9.json'
   outpath = os.path.join(dirpath,opspfilename)
   gc=0
   outf = open(outpath, 'w')
   outf.write('{ ')
   rawdatavals=np.zeros([ny,nx], dtype=np.float128)
   for xways in range(0, nx):
      for yways in range(0,ny):
        strlat=str(latit[yways])
        strlon=str(longs[xways])
        #print(tot_aq.shape, rawdatavals.shape)
        rawdatavals[yways, xways]=tot_ex_mx[yways,xways]
        #rawdatavals[yways, xways]=tot_aq[yways,xways]
        rounddatavals=round(rawdatavals[yways, xways])
        strval=str(rounddatavals)
        gc=gc+1
        if gc==gridsize:
          strtowrite=(' "'+strlat+', '+strlon+'": "'+strval+'" }')
        else:
          strtowrite=(' "'+strlat+','+strlon+'": "'+strval+'",')
        outf.write(strtowrite) 
   outf.close()
   shutil.move(fpath, fpath2)
   
