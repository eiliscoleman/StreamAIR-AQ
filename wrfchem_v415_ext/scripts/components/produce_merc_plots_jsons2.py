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
import xarray.ufuncs as xu
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


# create an array for each variables to be made into json files
coords_dict = {}
t2_dict={}
wdir_dict={}
w10_dict={}
uv_dict={}
o3_dict={}
rh_dict={}
psl_dict={}
swdown_dict={}
nox_dict={}
so2_dict={}
pm10_dict={}
pm25_dict={}
org25_dict={}
sulf25_dict={}
nitr25_dict={}
nh425_dict={}
cl25_dict={}
bc_dict={}
cldfra_dict={}
fog_dict={}
qsnow_dict={}
qgraup_dict={}
#
#
#
t2_json_path=(png_out_dir+'/'+datestr+':00-t2.json')
wdir_json_path=(png_out_dir+'/'+datestr+':00-wdir.json')
w10_json_path=(png_out_dir+'/'+datestr+':00-w10.json')
uv_json_path=(png_out_dir+'/'+datestr+':00-uv.json')
o3_json_path=(png_out_dir+'/'+datestr+':00-o3.json')
rh_json_path=(png_out_dir+'/'+datestr+':00-rh.json')
psl_json_path=(png_out_dir+'/'+datestr+':00-psl.json')
swdown_json_path=(png_out_dir+'/'+datestr+':00-swdown.json')
nox_json_path=(png_out_dir+'/'+datestr+':00-nox.json')
so2_json_path=(png_out_dir+'/'+datestr+':00-so2.json')
pm10_json_path=(png_out_dir+'/'+datestr+':00-pm10.json')
pm25_json_path=(png_out_dir+'/'+datestr+':00-pm25.json')
org25_json_path=(png_out_dir+'/'+datestr+':00-org25.json')
sulf25_json_path=(png_out_dir+'/'+datestr+':00-sulf25.json')
nitr25_json_path=(png_out_dir+'/'+datestr+':00-nitr25.json')
nh425_json_path=(png_out_dir+'/'+datestr+':00-nh425.json')
cl25_json_path=(png_out_dir+'/'+datestr+':00-cl25.json')
bc_json_path=(png_out_dir+'/'+datestr+':00-bc.json')
cldfra_json_path=(png_out_dir+'/'+datestr+':00-cldfra.json')
fog_json_path=(png_out_dir+'/'+datestr+':00-fog.json')
qsnow_json_path=(png_out_dir+'/'+datestr+':00-qsnow.json')
qgraup_json_path=(png_out_dir+'/'+datestr+':00-qgraup.json')
#
#
#
wind10=xu.sqrt(data.u10.values[0,:,:]**2+data.v10.values[0,:,:]**2)
#normalise the u and v components
uu_n=data.u10.values[0,:,:]/wind10
vv_n=data.u10.values[0,:,:]/wind10
#get wind dir
wdir=np.rad2deg(np.arctan2(vv_n,uu_n))
o3levs=data.o3_concentration.values[0,:,:,:]
uv_index_levs=((o3levs[:,:,:]*data.pb[0,:,:,:]/6950.0))

#uv_index=np.sum(uv_index_levs, axis=0)


#o3=data.o3_concentration.values[0,0,:,:]*1000

#rh=data.rh.values[0,:,:]
#psl=data.p_sl.values[0,:,:]*1e-2
## Below: the second dimension determines the levels to be summed over to determine the cloud fraction area. at present, the cloud fractions is calculated by 
#cldfra_levs=data.cldfra.values[0,5:11,:,:]
#print(cldfra_levs.shape)
#cldfra=np.amax(cldfra_levs, axis=0)
#print(cldfra.shape)

#for i, lonv in enumerate(data.lon.values):
    #for j, latv in enumerate(data.lat.values):
        #coords_string = (str(latv)+','+str(lonv))
        #coords_dict[j,i] = coords_string
        #t2_dict[j,i] = data.t2.values[0,j,i]-273.15
        #wdir_dict[j,i] = wdir[j,i]
        #w10_dict[j,i] = wind10[j,i]
        #uv_dict[j,i] = uv_index[j,i]
        #o3_dict[j,i]=o3[j,i]
        #rh_dict[j,i]=rh[j,i]
        #psl_dict[j,i]=psl[j,i]
        #swdown_dict[j,i]=data.swdown.values[0,j,i]
        #nox_dict[j,i]=data.nox_concentration.values[0,0,j,i]*1000
        #so2_dict[j,i]=data.so2_concentration.values[0,0,j,i]*1000
        #pm10_dict[j,i]=data.pm10.values[0,0,j,i]
        #pm25_dict[j,i]=data.pm25.values[0,0,j,i]
        #org25_dict[j,i]=data.org25.values[0,0,j,i]
        #sulf25_dict[j,i]=data.sulf25.values[0,0,j,i]
        #nitr25_dict[j,i]=data.nitr25.values[0,0,j,i]
        #nh425_dict[j,i]=data.nh425.values[0,0,j,i]
        #cl25_dict[j,i]=data.cl25.values[0,0,j,i]
        #bc_dict[j,i]=data.bc.values[0,0,j,i]
        #cldfra_dict[j,i]=cldfra[j,i]
        #fog_dict[j,i]=data.cldfra.values[0,0,j,i]
        #qsnow_dict[j,i]=data.qsnow.values[0,0,j,i]
        #qgraup_dict[j,i]=data.qgraup.values[0,0,j,i]
##############################
##first2pairs = {k: w10_dict[k] for k in w10_dict.keys()[:2]}
##print(w10_dict[0].type)
##print(uv_dict[0].type)
#coords_df = pd.DataFrame.from_dict(coords_dict, orient='index')

#t2_df = pd.DataFrame.from_dict(t2_dict, orient='index')
#t2_df['t2']=t2_df.values
#t2_df = t2_df.round({'t2': 2})
#t2_df['coords'] = coords_df
#text_output = t2_df.set_index('coords').to_dict()['t2']
#with open(t2_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## wind 10 json
#w10_df = pd.DataFrame.from_dict(w10_dict, orient='index')
#w10_df['w10']=w10_df.values
#w10_df = w10_df.round({'w10': 2})
#w10_df['coords'] = coords_df
#text_output = w10_df.set_index('coords').to_dict()['w10']
#with open(w10_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)

## wind dir json
#wdir_df = pd.DataFrame.from_dict(wdir_dict, orient='index')
#print(wdir_df.shape)
#wdir_df['wdir']=wdir_df.values
#print(wdir_df.shape)
#wdir_df = wdir_df.round({'wdir': 2})
#print(wdir_df.shape)
#wdir_df['coords'] = coords_df
#print(wdir_df.shape)
#text_output = wdir_df.set_index('coords').to_dict()['wdir']
#with open(wdir_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## uv dir json
#uv_df = pd.DataFrame.from_dict(uv_dict, orient='index')
#print(uv_df.shape)
#uv_df['uv']=uv_df.values
#print(uv_df.shape)
#uv_df = uv_df.round({'uv': 0})
#print(uv_df.shape)
#uv_df['coords'] = coords_df
#text_output = uv_df.set_index('coords').to_dict()['uv']
#print(uv_df.shape)
#keys_values = text_output.items()
#new_tx = {str(key): str(value) for key, value in keys_values}
##text_output = text_output.applymap(str)
#with open(uv_json_path, 'w') as fp:
    #json.dump(new_tx, fp, skipkeys=True, namedtuple_as_object=False)
## o3 dir json
#o3_df = pd.DataFrame.from_dict(o3_dict, orient='index')
#o3_df['o3']=o3_df.values
#o3_df = o3_df.round({'o3': 0})
#o3_df['coords'] = coords_df
#text_output = o3_df.set_index('coords').to_dict()['o3']
#with open(o3_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## rh dir json
#rh_df = pd.DataFrame.from_dict(rh_dict, orient='index')
#rh_df['rh']=rh_df.values
#rh_df = rh_df.round({'rh': 0})
#rh_df['coords'] = coords_df
#text_output = rh_df.set_index('coords').to_dict()['rh']
#with open(rh_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## psl dir json
#psl_df = pd.DataFrame.from_dict(psl_dict, orient='index')
#psl_df['psl']=psl_df.values
#psl_df = psl_df.round({'psl': 0})
#psl_df['coords'] = coords_df
#text_output = psl_df.set_index('coords').to_dict()['psl']
#with open(psl_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## swdown dir json
#swdown_df = pd.DataFrame.from_dict(swdown_dict, orient='index')
#swdown_df['swdown']=swdown_df.values
#swdown_df = swdown_df.round({'swdown': 0})
#swdown_df['coords'] = coords_df
#text_output = swdown_df.set_index('coords').to_dict()['swdown']
#with open(swdown_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## nox dir json
#nox_df = pd.DataFrame.from_dict(nox_dict, orient='index')
#nox_df['nox']=nox_df.values
#nox_df = nox_df.round({'nox': 1})
#nox_df['coords'] = coords_df
#text_output = nox_df.set_index('coords').to_dict()['nox']
#with open(nox_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## so2 dir json
#so2_df = pd.DataFrame.from_dict(so2_dict, orient='index')
#so2_df['so2']=so2_df.values
#so2_df = so2_df.round({'so2': 1})
#so2_df['coords'] = coords_df
#text_output = so2_df.set_index('coords').to_dict()['so2']
#with open(so2_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## pm10 dir json
#pm10_df = pd.DataFrame.from_dict(pm10_dict, orient='index')
#pm10_df['pm10']=pm10_df.values
#pm10_df = pm10_df.round({'pm10': 1})
#pm10_df['coords'] = coords_df
#text_output = pm10_df.set_index('coords').to_dict()['pm10']
#with open(pm10_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)    
## pm25 dir json
#pm25_df = pd.DataFrame.from_dict(pm25_dict, orient='index')
#pm25_df['pm25']=pm25_df.values
#pm25_df = pm25_df.round({'pm25': 1})
#pm25_df['coords'] = coords_df
#text_output = pm25_df.set_index('coords').to_dict()['pm25']
#with open(pm25_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)      
## org25 dir json
#org25_df = pd.DataFrame.from_dict(org25_dict, orient='index')
#org25_df['org25']=org25_df.values
#org25_df = org25_df.round({'org25': 4})
#org25_df['coords'] = coords_df
#text_output = org25_df.set_index('coords').to_dict()['org25']
#with open(org25_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)    
## sulf25 dir json
#sulf25_df = pd.DataFrame.from_dict(sulf25_dict, orient='index')
#sulf25_df['sulf25']=sulf25_df.values
#sulf25_df = sulf25_df.round({'sulf25': 4})
#sulf25_df['coords'] = coords_df
#text_output = sulf25_df.set_index('coords').to_dict()['sulf25']
#with open(sulf25_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False) 
## nitr25 dir json
#nitr25_df = pd.DataFrame.from_dict(nitr25_dict, orient='index')
#nitr25_df['nitr25']=nitr25_df.values
#nitr25_df = nitr25_df.round({'nitr25': 4})
#nitr25_df['coords'] = coords_df
#text_output = nitr25_df.set_index('coords').to_dict()['nitr25']
#with open(nitr25_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)  
## nh425 dir json
#nh425_df = pd.DataFrame.from_dict(nh425_dict, orient='index')
#nh425_df['nh425']=nh425_df.values
#nh425_df = nh425_df.round({'nh425': 4})
#nh425_df['coords'] = coords_df
#text_output = nh425_df.set_index('coords').to_dict()['nh425']
#with open(nh425_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)  
## cl25 dir json
#cl25_df = pd.DataFrame.from_dict(cl25_dict, orient='index')
#cl25_df['cl25']=cl25_df.values
#cl25_df = cl25_df.round({'cl25': 4})
#cl25_df['coords'] = coords_df
#text_output = cl25_df.set_index('coords').to_dict()['cl25']
#with open(cl25_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)      
## bc dir json
#bc_df = pd.DataFrame.from_dict(bc_dict, orient='index')
#bc_df['bc']=bc_df.values
#bc_df = bc_df.round({'bc': 4})
#bc_df['coords'] = coords_df
#text_output = bc_df.set_index('coords').to_dict()['bc']
#with open(bc_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)    
## cldfra dir json
#cldfra_df = pd.DataFrame.from_dict(cldfra_dict, orient='index')
#cldfra_df['cldfra']=cldfra_df.values
#cldfra_df = cldfra_df.round({'uv': 0})
#cldfra_df['coords'] = coords_df
#text_output = cldfra_df.set_index('coords').to_dict()['cldfra']
#keys_values = text_output.items()
#new_tx = {str(key): str(value) for key, value in keys_values}
#with open(cldfra_json_path, 'w') as fp:
    #json.dump(new_tx, fp, skipkeys=True, namedtuple_as_object=False)
## fog dir json
#fog_df = pd.DataFrame.from_dict(fog_dict, orient='index')
#fog_df['fog']=fog_df.values
#fog_df = fog_df.round({'uv': 0})
#fog_df['coords'] = coords_df
#text_output = fog_df.set_index('coords').to_dict()['fog']
##keys_values = text_output.items()
##new_tx = {str(key): str(value) for key, value in keys_values}
#with open(fog_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## qsnow dir json
#qsnow_df = pd.DataFrame.from_dict(qsnow_dict, orient='index')
#qsnow_df['qsnow']=qsnow_df.values
#qsnow_df = qsnow_df.round({'uv': 0})
#qsnow_df['coords'] = coords_df
#text_output = qsnow_df.set_index('coords').to_dict()['qsnow']
##keys_values = text_output.items()
##new_tx = {str(key): str(value) for key, value in keys_values}
#with open(qsnow_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
## qgraup dir json
#qgraup_df = pd.DataFrame.from_dict(qgraup_dict, orient='index')
#qgraup_df['qgraup']=qgraup_df.values
#qgraup_df = qgraup_df.round({'uv': 0})
#qgraup_df['coords'] = coords_df
#text_output = qgraup_df.set_index('coords').to_dict()['qgraup']
##keys_values = text_output.items()
##new_tx = {str(key): str(value) for key, value in keys_values}
#with open(qgraup_json_path, 'w') as fp:
    #json.dump(text_output, fp, skipkeys=True, namedtuple_as_object=False)
##---------------------------------------------------------
    
##print ('T2', data.t2.shape )
#t2_out_path=(png_out_dir+'/'+datestr+':00-t2.png')
#lev_range=range(-20,40)  
#levels=lev_range
#fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#ax=fig.add_subplot(111, projection=ccrs.Mercator())
#cs=ax.contourf(LON, LAT, data.t2.values[0,:,:]-273.15, levels, cmap=plt.cm.jet, extend="both", transform=ccrs.PlateCarree())
#plt.box(on=None)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
#ax.coastlines('10m', linewidth=0.15)
#plt.axis('off')
#ax.figsize=(tilesize/dpi, tilesize/dpi)
#ax.dpi=dpi
#ax.outline_patch.set_visible(False)
#ax.background_patch.set_visible(False)
#ax.patch.set_alpha(0)
#ax.add_feature(cfeature.BORDERS,linewidth=0.25)
#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)
##plt.contourf(LON, LAT, data.t2.values[0,:,:]-273.15, levels, cmap=plt.cm.get_cmap('jet'), extend="both")
#plt.savefig(t2_out_path, dpi=288, tilesize=768, transparent=True)
#plt.close()
##-----------------------------------------------------

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
for lev in range(0,1):

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
  
  # plot the PM constituents
  pm_species=['org25', 'sulf25', 'nh425', 'cl25', 'bc', 'nitr25']
  for pms in pm_species:
      outfile=(datestr+'_'+pms+'.png')
      print(outfile)
      png_out_path=(png_out_dir+'/'+outfile)
      lev_range=range(0,10)  
      levels=lev_range
      fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
      ax=fig.add_subplot(111, projection=ccrs.Mercator())
      #print(data[sp].shape)
      cs=ax.contourf(LON, LAT, data[pms].values[0,0,:,:], levels, cmap=plt.cm.jet, extend="both", transform=ccrs.PlateCarree())
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
      fig.colorbar(cs, ax=ax, orientation='horizontal', shrink=0.5)
      plt.savefig(png_out_path, dpi=288, tilesize=768, transparent=False)
      plt.close()
      

