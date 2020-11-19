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

print('wind Colourbar hh')

#cmap.set_under('w')
#cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
cmap=plt.cm.jet
font = FontProperties()
font.set_name('Ariel')
#cm.set_over('purple')
bounds=[0,10,20,30,40,50,60,70,80,90,100]
#cmap.set_over('#5B2C6F')
#bounds=[0, .01, .1, .25, .5, 1, 2.5, 5, 10, 20, 50, 100]
#bounds=[0, .1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
#bounds=[0.0001, 0.1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0, 400]
#norm = mc.BoundaryNorm(bounds, 256)
#norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
cbar_png_out_path=(png_out_dir+'/wind_cbar_kph.png')
fig=plt.figure(figsize=(10, 1))
#add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='black', linewidth='0', label='None') )
#ax2.annotate(' Rain (mm/hr)', (0.001, 0.5),  color='white', weight='bold', fontsize=18)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
ax2=fig.add_axes([0.0001, 0.005, 0.2, 0.85])
#Rectange((x,y, width, height))
ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.99999, 0.99999, color='navy', linewidth='0', label='None') )
fig.patch.set_visible(False)
#ax2.axis('off')
#plt.box(on=None)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["bottom"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.annotate('Windspeed kph', (0.01, 0.45),  color='white', weight='bold', fontsize=16)
#ax2.annotate(' o', (0.785, .6),  color='white', weight='bold', fontsize=8)
ax1=fig.add_axes([0.2, 0.008, 0.8, 0.847])
#ax1.axis('off')
#plt.box(on=None)
#ax1.patch.set_visible(False)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.spines["bottom"].set_visible(False)
ax1.spines["left"].set_visible(False)
norm = mpl.colors.Normalize(vmin=0, vmax=108)
cbarticks=np.arange(0,110,2)
cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal', ticks=bounds)
cbar.ax.set_xticklabels([' '])
cbar.ax.get_yaxis().set_ticks([])
divfac=divfac=50-0.4
for j, lab in enumerate(bounds):
    print(j, lab)
    #cbar.ax.text( 0.00, .5, "0.0", ha='left', va='center', color="white", fontsize=16, fontweight="bold")
    #print(float(lab)/108.0)
    #cbar.ax.text( float(lab)/float(108)-0.025, .5, str(lab), ha='left', va='center',   color="white", fontsize=18, fontweight="bold")
    if j>0:
        cbar.ax.text( float(lab)/float(108), .5, str(lab), ha='left', va='center',   color="white", fontsize=18, fontweight="bold")
    else:
        cbar.ax.text( 0.01, .5, '0', ha='left', va='center',   color="white", fontsize=18, fontweight="bold")
cbar.ax.tick_params(labelsize=6)   
#cbar.ax.text( ( (3 + 0.05)) / (divfac), .5, "10", ha='left', va='center', color="blue", fontsize=18, fontweight="bold")
#cbar.ax.text( ( (4 + 0.05)) / (divfac), .5, "20", ha='left', va='center', color="blue", fontsize=18, fontweight="bold")
#cbar.ax.text( ( (5 + 0.05)) / (divfac), .5, "30", ha='left', va='center', color="blue", fontsize=18, fontweight="bold")
#cbar.ax.text( ( (6 + 0.05)) / (divfac), .5, "40", ha='left', va='center', color="white", fontsize=18, fontweight="bold")
#cbar.ax.tick_params(labelsize=6)#ax2.annotate(r'$\ m s^{-1}$'
#ax1.axis('off')
#plt.box(on=None)
#ax1.patch.set_visible(False)
cbar.outline.set_visible(False)
plt.savefig(cbar_png_out_path, transparent=True)
plt.close()  


