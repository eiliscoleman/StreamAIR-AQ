from __future__ import print_function, unicode_literals
#import cartopy as cpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
#from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
#import mercantile as mti
from copy import copy



tilesize = 768
dpi = 288
#TERESA: Change teh fpath to look at different files. the directory YYYYMMDD contains simulaions related to that date, and the file name refers tot he specefic hour.  
#change the op_png to the path you want: it'll have to be in your own directory becuase that's where you can write. 
fpath='/mnt/raid/wrf-chem/wrfchem_v415/data_back/20200513/output-wrf' 
fnm='d01_202005141800.nc'
op_png='/mnt/raid/wrf-chem/wrfchem_v415/bc_grid_140520_1800_2.png'
f_path=os.path.join(fpath, fnm)
    
    
data = xr.open_dataset(f_path)
LON= data.lon
LAT=data.lat

extents = {
     #'ireland': [-12, -3, 51, 55.5]
#IRE domain
     'europe':[LON.min(),LON.max(),LAT.min(),LAT.max()-1]
     }
#lev_range=range(0,10)  
#levels=lev_range

for dom in extents:


    print('loms', LON.shape)
    print('lats', LAT.shape)
    print ('terrain', data.terrain.values[:,:].shape )
    #t2_out_path=('hgt_raw_'+dom+'.png')
    #lev_range=np.linspace(0, 1.5, num=15, endpoint=False)
    #lev_range=np.logspace(0, 20, num=15, endpoint=False)
    bounds=[0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 5, 7, 10, 20, 30]
    lev_range=np.array(bounds)
    levels=lev_range
    print(levels)
    norm = mc.BoundaryNorm(boundaries=bounds, ncolors=256)
    fig=plt.figure(figsize=(10,6), dpi=dpi)
    ax=fig.add_subplot(111, projection=ccrs.Mercator())
    print('2')
    print(data.bc.shape)
    cs=ax.contourf(LON, LAT, data.bc.values[0, :,:]*1.2, norm=norm, levels=levels, cmap=plt.cm.jet, extend="both", transform=ccrs.PlateCarree())
    #cs=ax.contourf(LON, LAT, data.bc.values[0, :,:], cmap=plt.cm.jet, extend="both", transform=ccrs.PlateCarree())
    #cs=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  #linewidth=2, color='gray', alpha=0.5, linestyle='--')
    cbar = fig.colorbar(cs,ticks=lev_range, shrink=1.2)
    #cbar.ax.set_yticklabels(lev_range)
#   plt.box(on=None)
#plt.subplots_adjust(bottom=0, left=0, right=1, top=1,  hspace = 0, wspace = 0)
    #extents=[LON.min(),LON.max(),LAT.min(),LAT.max()]
    print(extents)
    ax.set_extent(extents[dom])
    ax.coastlines('10m', linewidth=0.15)
    #plt.axis('off')
    #ax.figsize=(tilesize/dpi, tilesize/dpi)
    ax.dpi=dpi
    plt.savefig(op_png, dpi=288, tilesize=768, transparent=True)
    plt.close()
