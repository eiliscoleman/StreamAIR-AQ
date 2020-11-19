from __future__ import print_function, unicode_literals

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
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
#import mercantile as mti

import argparse
from argparse import RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description = "read in args", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('filein', help = 'file to plot')
parser.add_argument('oproot', help = 'Path for output')
args = parser.parse_args()

tilesize = 768
dpi = 288


dir_path=args.oproot
print(dir_path)

png_out_dir=args.oproot
fid=args.filein

f_path=(png_out_dir+'/'+fid)
print(f_path)
print('this is the file to be plotted', args.filein)
data = xr.open_dataset(f_path)
LON, LAT=np.meshgrid(data.lon.values, data.lat.values)



print ('RAIN and pressure')
out_path=(png_out_dir+'/test_rain_psl.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(0,50,1)  
levels=lev_range
cs=plt.contourf(LON, LAT, data.rain.values[0,:,:], levels, cmap=plt.cm.Blues, extend="both")
cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n=20, linewidths=0.5, colors='k')
levels_n=np.arange(1000,1040,1)
cs2=plt.contour(LON, LAT, data.p_sl[0,:,:]*1e-2, levels_n, linewidths=0.2, colors='k')
plt.clabel(cs2, inline=1, fmt='%.0f', fontsize=2)
plt.savefig(out_path, dpi=1152, tilesize=768)
plt.close()

