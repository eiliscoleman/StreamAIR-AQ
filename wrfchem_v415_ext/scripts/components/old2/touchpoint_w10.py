from __future__ import print_function, unicode_literals

import os
import json
import numpy as np
from matplotlib.mlab import griddata
from netCDF4 import Dataset
import argparse
from argparse import RawDescriptionHelpFormatter
from numpy.random import uniform, seed
import xarray.ufuncs as xu

parser = argparse.ArgumentParser(description = "", epilog = "", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('--lat', default = 'XLAT', dest = 'latitude_name', help = 'Name for latitude variable')
parser.add_argument('--lon', default = 'XLONG', dest = 'longitude_name', help = 'Name for longitude variable')
parser.add_argument('-u', '--ucomponent', default = 'u10', dest = 'ucomponent', help = 'Name for u-component of wind')
parser.add_argument('-v', '--vcomponent', default = 'v10', dest = 'vcomponent', help = 'Name for v-component of wind')

parser.add_argument('wrffile', help = 'Path to wrfout or wrfinput file that contains latitude_name, longitude_name, w10')
parser.add_argument('outroot', help = 'Path for output', default = '../public/data/weather/')
#parser.add_argument('o3_conc',default = 'o3_concentration', help = 'Name for O3 conc')
args = parser.parse_args()
outroot = args.outroot
uhdr = {}
vhdr = {}
data = [uhdr, vhdr]
newf = Dataset(args.wrffile)

nx = 162
ny = 114
dx = 0.25
dy = 0.25
la1 = 62.75
la2 = 34.5
lo1 = -21.75
lo2 = 18.5
gridsize=162*114
longs=np.arange(lo1,lo2+dx,dx)
latit=np.arange(la1, la2-dy,-dy)
times = newf.variables['Times']
timestr = 9999
dirpath = args.outroot
outpath = os.path.join(dirpath, '%s-w10-surface-level-gfs-1.0.json' % (timestr,))
gc=0
outf = open(outpath, 'w')
outf.write('{ ')
uu10=newf.variables[args.ucomponent]
vv10=newf.variables[args.vcomponent]
w10=xu.sqrt(uu10[0,:,:]**2+vv10[0,:,:]**2)
for xways in range(0, nx):
  for yways in range(0,ny):
    strlat=str(latit[yways])
    strlon=str(longs[xways])
    rawdatavals=w10[yways,xways]
    rounddatavals=round(rawdatavals,3)
    strval=str(rounddatavals)
    gc=gc+1
    if gc==gridsize:
      strtowrite=(' "'+strlat+', '+strlon+'": {"value":"'+strval+'m/s"} }')
    else:
      strtowrite=(' "'+strlat+','+strlon+'": {"value":"'+strval+'m/s"},')
    outf.write(strtowrite)    
outf.close()
