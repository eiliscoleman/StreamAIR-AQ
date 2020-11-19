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

png_out_dir=('/mnt/raid/wrf-chem/wrfchem_v415/cbars')
#

print('rain Colourbar hh')


#New lines to overwrite previous colour scheme
palette=['w', 'lightcyan', 'paleturquoise', 'turquoise', 'lightskyblue', 'royalblue', 'orange', 'orangered', 'firebrick', 'darkmagenta']
#cmap.set_under('w')
#cm=LinearSegmentedColormap.from_list('palette', palette, N=len(palette))
cmap=plt.cm.gnuplot2
font = FontProperties()
font.set_name('Ariel')
#cm.set_over('purple')

#cmap.set_over('#5B2C6F')
#bounds=[0, .01, .1, .25, .5, 1, 2.5, 5, 10, 20, 50, 100]
#bounds=[0, .1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0]
bounds=[0.0001, 0.1, .3, .5, .7, 1.0, 2.0, 3.0, 6.0, 12.0, 400]
norm = mc.BoundaryNorm(bounds, 256)
#norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
cbar_png_out_path=(png_out_dir+'/rain_cbar_gnuplot2.png')
fig=plt.figure(figsize=(10, 1))
#add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='black', linewidth='0', label='None') )
#ax2.annotate(' Rain (mm/hr)', (0.001, 0.5),  color='white', weight='bold', fontsize=18)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
ax1=fig.add_axes([0.15, 0.05, 0.8, 0.8475])
#norm = mpl.colors.LogNorm(vmin=1, vmax=100)
#cbarticks=np.logspace(0,2,num=10)
cbarticks=np.array(bounds)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=palette, norm=norm, orientation='horizontal', extend='max')
cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
cbar.ax.set_xticklabels([' '])
cbar.ax.get_yaxis().set_ticks([])
divfac=cbarticks.shape[0]
#print divfac

cbar.ax.text( ( (10.4)) / (divfac), .5, '12.0', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (9.2 )) / (divfac), .5, '6.0', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (8.29 )) / (divfac), .5, '3.0', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (7.2 )) / (divfac), .5, '2.0', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (6  )) / (divfac), .5, '1.0', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (5 )) / (divfac), .5, '0.7', ha='center', va='center', weight='bold', color="black", fontsize="20", fontweight="bold")
cbar.ax.text( ( (3.7 )) / (divfac), .5, '0.5', ha='center', va='center', weight='bold', color="white", fontsize="20", fontweight="bold")
cbar.ax.text( ( (2.8 )) / (divfac), .5, '0.3', ha='center', va='center', weight='bold', color="white", fontsize="20", fontweight="bold")
cbar.ax.text( ( (1.5 )) / (divfac), .5, '0.1', ha='center', va='center', weight='bold', color="white", fontsize="20", fontweight="bold")
ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
#Rectange((x,y, width, height))
ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='black', linewidth='0', label='None') )
ax2.annotate(' Rain mm hr', (0.001, 0.45),  color='white', weight='bold', fontsize=20)
ax2.annotate('           -1', (0.9, 0.6),  color='white', weight='bold', fontsize=12)
#ax2.annotate(r'$\ m s^{-1}$'
plt.savefig(cbar_png_out_path, transparent=True)
plt.close() 

#print('RH Colourbar')
#cmap=plt.cm.get_cmap('gist_rainbow_r')
#cbar_png_out_path=(png_out_dir+'/extended_labelled/toserver/rh_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='deeppink', linewidth='0', label='None') )
#ax2.annotate('RH $\%$ ', (0.06, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
#norm = mpl.colors.Normalize(vmin=30, vmax=100)
#cbarticks=np.arange(30,110, 10)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]
#for j, lab in enumerate(cbarticks):
    ##print(j, lab)
    #cbar.ax.text( ( (j + 0.1)) / (divfac-1.+.35), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (5 + 0.1)) / (divfac-1.+.35), .5, "80", ha='center', va='center', color="navy", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (6 + 0.1)) / (divfac-1.+.35), .5, "90", ha='center', va='center', color="navy", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (7 + 0.1)) / (divfac-1.+.35), .5, "100", ha='center', va='center', color="navy", fontsize="20", fontweight="bold")
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close() 

#print('O3 Colourbar')
#cmap=plt.cm.jet
#cbar_png_out_path=(png_out_dir+'/extended_labelled/o3_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='navy', linewidth='0', label='None') )
#ax2.annotate('$\ O_3$ $\ ppb $ ', (0.06, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
#norm = mpl.colors.Normalize(vmin=30, vmax=50)
#cbarticks=np.arange(30,55, 5)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]-.6
#for j, lab in enumerate(cbarticks):
    #cbar.ax.text( ( (j + 0.07)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (2 + 0.07)) / (divfac), .5, "40", ha='center', va='center', color="navy", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (3 + 0.07)) / (divfac), .5, "45", ha='center', va='center', color="navy", fontsize="20", fontweight="bold")
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close()  



#print('NOx Colourbar')
#cmap=plt.cm.gnuplot
#cbar_png_out_path=(png_out_dir+'/extended_labelled/nox_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='k', linewidth='0', label='None') )
#ax2.annotate('$\ NO_x$ $\ ppb $ ', (0.06, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
##norm = mpl.colors.Normalize(vmin=0, vmax=40)
#cbarticks=np.arange(0,18, 2)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,  orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]-.45
#for j, lab in enumerate(cbarticks):
    ##print(j, lab)
    #cbar.ax.text( ( (j + 0.07)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (7 + 0.07)) / (divfac), .5, "14", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (8 + 0.07)) / (divfac), .5, "16", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")    
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close()  

#print('SO2 Colourbar')
#cmap=plt.cm.inferno
#cbar_png_out_path=(png_out_dir+'/extended_labelled/so2_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='k', linewidth='0', label='None') )
#ax2.annotate('$\ SO_2$ $\ ppb $ ', (0.06, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
#norm = mpl.colors.Normalize(vmin=0, vmax=8)
#cbarticks=np.arange(0,9, 1)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]-.45
#for j, lab in enumerate(cbarticks):
    ##print(j, lab)
    #cbar.ax.text( ( (j + 0.07)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (6 + 0.07)) / (divfac), .5, "6", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (7 + 0.07)) / (divfac), .5, "7", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")     
#cbar.ax.text( ( (8 + 0.07)) / (divfac), .5, "8", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")    
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close()  

#print('PM10 Colourbar')
#cmap=plt.cm.get_cmap('magma')
#cbar_png_out_path=(png_out_dir+'/extended_labelled/pm10_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='k', linewidth='0', label='None') )
#ax2.annotate('$\ PM_{10}$ $\mu g m^{-3} $ ', (0.03, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
#norm = mpl.colors.Normalize(vmin=0, vmax=20)
#cbarticks=np.arange(0,22, 2)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]-.35
#for j, lab in enumerate(cbarticks):
    #print(j, lab)
    #cbar.ax.text( ( (j + 0.075)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (8 + 0.075)) / (divfac), .5, "16", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (9 + 0.075)) / (divfac), .5, "18", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (10 + 0.075)) / (divfac), .5, "20", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close()   

#print('PM25 Colourbar')
#cmap=plt.cm.get_cmap('magma')
#cbar_png_out_path=(png_out_dir+'/extended_labelled/pm25_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.2, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='k', linewidth='0', label='None') )
#ax2.annotate('$\ PM_{25}$ $\mu g m^{-3} $ ', (0.03, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.20, 0.05, 0.8, 0.85])
#norm = mpl.colors.Normalize(vmin=0, vmax=20)
#cbarticks=np.arange(0,22, 2)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]-.35
#for j, lab in enumerate(cbarticks):
    #print(j, lab)
    #cbar.ax.text( ( (j + 0.075)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (8 + 0.075)) / (divfac), .5, "16", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (9 + 0.075)) / (divfac), .5, "18", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.text( ( (10 + 0.075)) / (divfac), .5, "20", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold") 
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close() 

#print('AQI Colourbar')
#cmap=plt.cm.get_cmap('viridis')
#cbar_png_out_path=(png_out_dir+'/extended_labelled/aqi_cbar.png')
#fig=plt.figure(figsize=(10, 1))
##add_axes(Left, bottom, width, height)
#ax2=fig.add_axes([0.0001, 0.05, 0.15, 0.85])
##Rectange((x,y, width, height))
#ax2.add_patch(patches.Rectangle((0.00001, 0.0001), 0.9999, 0.9999, color='indigo', linewidth='0', label='None') )
#ax2.annotate('AQI ', (0.12, 0.45),  color='white', weight='bold', fontsize=20)
##ax2.annotate(r'$\ m s^{-1}$', (0.25, 0.3),  color='w', weight='bold', fontsize=10)
#ax1=fig.add_axes([0.15, 0.05, 0.9, 0.85])
##norm = mpl.colors.Normalize(vmin=30, vmax=100)
#cbarticks=np.arange(0,11, 1)
#cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, orientation='horizontal')
#cbar.ax.set_xticklabels([' '])
#cbar.ax.get_yaxis().set_ticks([])
#divfac=cbarticks.shape[0]
#for j, lab in enumerate(cbarticks):
    ##print(j, lab)
    #cbar.ax.text( ( (j + 0.1)) / (divfac), .5, lab, ha='center', va='center', color="white", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (10 + 0.1)) / (divfac), .5, "10", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (9 + 0.1)) / (divfac), .5, "9", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")
#cbar.ax.text( ( (8 + 0.1)) / (divfac), .5, "8", ha='center', va='center', color="indigo", fontsize="20", fontweight="bold")
##cbar.ax.text( ( (7 + 0.1)) / (divfac-1.+.35), .4, "100", ha='center', va='center', color="white", fontsize="12", fontweight="bold")
#cbar.ax.tick_params(labelsize=6)
#plt.savefig(cbar_png_out_path, transparent=True)
#plt.close() 
