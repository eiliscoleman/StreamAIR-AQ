#!/usr/bin/env python3
import json
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import seaborn as sns
from matplotlib.font_manager import FontProperties
#from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
from copy import copy
#import seaborn as sns
#import mercantile as mti
#from wrf import to_np, getvar, smooth2d, get_cartopy
tilesize = 768
dpi = 288
dpi3 = 864

png_out_dir=('/mnt/raid/wrf-chem/wrfchem_v415/scripts/components/cbars/pngs')
#

print('pm25 Colourbar hh')

#cmap.set_under('w')
#cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
cmap=plt.cm.magma
font = FontProperties()
font.set_name('Ariel')
#cm.set_over('purple')
bounds=[0, .01, .1, .2, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
#cmap.set_over('#5B2C6F')
#bounds=[0, .01, .1, .25, .5, 1, 2.5, 5, 10, 20, 50, 100]
#bounds=[0, .1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
#bounds=[0.0001, 0.1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0, 400]
#norm = mc.BoundaryNorm(bounds, 256)
#norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
cbar_png_out_path=(png_out_dir+'/pm25_cbar.png')
fig=plt.figure(figsize=(10, 1))
#add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='black', linewidth='0', label='None') )
#ax2.annotate(' Rain (mm/hr)', (0.001, 0.5),  color='white', weight='bold', fontsize=18)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
#Rectange((x,y, width, height))
ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='black', linewidth='0', label='None') )
ax2.annotate('PM  $\mu$g m', (0.005, 0.45),  color='white', weight='bold', fontsize=18)
ax2.annotate('25', (0.299, 0.4),  color='white', weight='bold', fontsize=8)
ax2.annotate('-3', (0.85, 0.65),  color='white', weight='bold', fontsize=10)
ax1=fig.add_axes([0.15, 0.05, 0.85, 0.8475])
norm = mpl.colors.Normalize(vmin=0, vmax=30)
cbarticks=np.arange(0,32, 2)
cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbarticks=np.arange(0,22,2)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=palette, norm=norm, orientation='horizontal', extend='max')
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, orientation='horizontal')
cbar.ax.set_xticklabels([' '])
cbar.ax.get_yaxis().set_ticks([])
divfac=cbarticks.shape[0]-.35
for j, lab in enumerate(cbarticks):
    print(j, lab)
    if j<9:
        labcol="white"
    else:
        labcol="indigo"
            
    cbar.ax.text( ( (j + 0.075)) / (divfac), .5, lab, ha='center', va='center', color=labcol, fontsize="20", fontweight="bold")
#for j, lab in enumerate(cbarticks):
    #print(j, lab)
    #cbar.ax.text( ( (j + 0.09)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (8 + 0.09)) / (divfac), .5, "16", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (9 + 0.09)) / (divfac), .5, "18", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (10 + 0.09)) / (divfac), .5, "20", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
cbar.ax.tick_params(labelsize=6)
plt.savefig(cbar_png_out_path, transparent=True)
plt.close()


