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


 
#                                    Plot NOx
print('NOx')
rh_png_out_path=(png_out_dir+'/test_nox.png')
plt.figure(figsize=(tilesize/dpi, tilesize/dpi), dpi=dpi)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1)
plt.axis('off')
lev_range=np.arange(1, 60, 1)
levels=lev_range
plt.contourf(LON, LAT, data.nox_concentration.values[0,0,:,:]*1000., levels, cmap=plt.cm.rainbow, extend="both")
plt.savefig(rh_png_out_path, dpi=288, tilesize=768)
plt.close()
#                 
#                                    Plot SO2
