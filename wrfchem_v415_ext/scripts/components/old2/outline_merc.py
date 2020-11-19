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
tagon=0
#loninc.values/2
ex_dom=([ds.lon.min()-tagon, ds.lon.max()+tagon,ds.lat.min()-tagon, (ds.lat.max() +tagon-1.5)])
#ex_dom=         
print(ex_dom, ds.lat.max())

png_out_path='/mnt/raid/wrf-chem/wrfchem_v39/imgs/outline_merc_fine.svg'
lons = ds.variables['lon'][:]
lats = ds.variables['lat'][:]
terr=ds.variables['t2']
print(ds.variables['terrain'].shape)



#plotting
levels=(-5, -4, -2, 0,1, 2 , 5, 7, 10, 12, 13, 14, 15, 20)
fig=plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
ax=fig.add_subplot(111, projection=ccrs.PlateCarree())
# NB: in order to contourf onto a particular projection: must define using "transform" the projection relevant to the dataset. Eg: if the dataset is in equirectangular grid (latlon): then the contourf trasnform must be defined as ccrs.PlateCarree()
#cs=ax.contourf(lons, lats, terr[0,:,:], transform=ccrs.PlateCarree())
ax.set_extent(ex_dom)
plt.box(on=None)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)


#cs.background_patch=None



#ax.set_extent(ex_dom)
ax.coastlines('10m', linewidth=0.05)
plt.axis('off')
ax.figsize=(tilesize/dpi, tilesize/dpi)
ax.dpi=dpi
ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)
ax.background_patch.set_alpha(0)
#ax.add_feature(cfeature.BORDERS,linewidth=0.25)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
#ax.add_feature(cfeature.LAND)
#ax.set_axis_off()
#

plt.savefig(png_out_path, format='svg', bbox_inches='tight', pad_inches=0, dpi=288, tilesize=768, transparent=True)
plt.close()


