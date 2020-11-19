#!/usr/bin/env python3
#from mpl_toolkits.basemap import Basemap
import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import json
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import seaborn as sns

#very useful tutorial ! http://earthpy.org/cartopy_backgroung.html
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


tilesize = 768
dpi = 288

remapdir='/mnt/raid/wrf-chem/wrfchem_v39/data_back/output-wrf-remapped'

fplot=os.path.join(remapdir, 'd01reg_2019011800')
#print(fplot)

ds=xr.open_dataset(fplot)
#print(ds.lon)
lon0=ds.lon.min()
lon1=ds.lon[1]
loninc=(lon1-lon0)
#print(lon0, lon1, 'loninc', loninc, loninc.values/2)
tagon=0
#loninc.values/2
ex_dom=([ds.lon.min()-tagon, ds.lon.max()+tagon,
         ds.lat.min()-tagon, (ds.lat.max()-1.5 +tagon)])
print(ex_dom, ds.lat.max())

png_out_path='/mnt/raid/wrf-chem/wrfchem_v39/imgs/outline.png'
lons = ds.variables['lon'][:]
lats = ds.variables['lat'][:]
terr=ds.variables['terrain']
print(ds.variables['terrain'].shape)
plt.axis('off')
plt.figure()
plt.box(on=None)
#plt.set_frame_on(False)
#(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
#margins = {  #     vvv margin in inches
    #"left"   :     .00000001 ,
    #"bottom" :     0.00000001 ,
    #"right"  : 1 - .00000001,
    #"top"    : 1 - .00000001 
#}
#plt.subplots_adjust(**margins)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
#plt.axis('off')
#plt.background_patch=None
#plt.contourf(lons, lats, terr[0,:,:])
cs3=plt.axes(projection=ccrs.PlateCarree(), frameon=False)
cs3.axis('off')
cs3.figsize=(tilesize/dpi, tilesize/dpi)
cs3.dpi=dpi
cs3.outline_patch.set_visible(False)
  
#cs3=plt.axes(projection=ccrs.PlateCarree())
cs3.add_feature(cfeature.BORDERS,linewidth=0.25)

cs3.axes.get_xaxis().set_visible(False)

cs3.axes.get_yaxis().set_visible(False)
#plt.get_xaxis().set_visible(False)

#plt.get_yaxis().set_visible(False)
##gl = cs3.gridlines(draw_labels=True, linestyle='-', color='gray',
                      #linewidth=0.25)
##gl.xlabels_top = None
#gl.ylabels_right = None
#cs3.add_feature(cfeature.LAND)


cs3.set_extent(ex_dom)
cs3.coastlines('10m')


cs3.set_axis_off()

#
levels=(25, 50,100,150,200,250,300,350,400,500,1000, 1500, 2000)
#plt.contourf(lons, lats, terr[0,:,:], levels=levels)
plt.savefig(png_out_path, bbox_inches='tight', pad_inches=0, dpi=288, tilesize=768)
plt.close()


