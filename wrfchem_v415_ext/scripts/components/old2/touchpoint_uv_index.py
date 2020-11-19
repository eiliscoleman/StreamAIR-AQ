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
parser.add_argument('-o3', '--o3_conc', default = 'o3_concentration', dest = 'o3_conc', help = 'Name for o3 concentration')
parser.add_argument('-pb', '--pb', default = 'pb', dest = 'pb', help = 'Name for ground state pressure')

parser.add_argument('wrffile', help = 'Path to wrfout or wrfinput file that contains latitude_name, longitude_name, o3, pb')
parser.add_argument('outroot', help = 'Path for output', default = '../public/data/weather/')
#parser.add_argument('o3_conc',default = 'o3_concentration', help = 'Name for O3 conc')
args = parser.parse_args()
outroot = args.outroot
newf = Dataset(args.wrffile)
nx = 162
ny = 114
dx = 0.25
dy = 0.25
la1 = 34.5
la2 = 62.75
lo1 = -21.75
lo2 = 18.5
gridsize=162*114
longs=np.arange(lo1,lo2+dx,dx)
latit=np.arange(la1, la2+dy,dy)
times = newf.variables['Times']
timestr = 9999
dirpath = args.outroot
outpath = os.path.join(dirpath, '%s-uv-index.json' % (timestr,))
gc=0
outf = open(outpath, 'w')
outf.write('{ ')
o3=newf.variables[args.o3_conc]
pb=newf.variables[args.pb]
uv_index_levs=((o3[0,:,:,:]*pb[0,:,:,:]/6950.0))
uv_index=np.sum(uv_index_levs, axis=0)


for xways in range(0, nx):
  for yways in range(0,ny):
    strlat=str(latit[yways])
    strlon=str(longs[xways])
    rawdatavals=uv_index[yways,xways]
    rounddatavals=round(rawdatavals, 1)
    strval=str(rounddatavals)
    gc=gc+1
    if gc==gridsize:
      strtowrite=(' "'+strlat+', '+strlon+'":"'+strval+'" }')
    else:
      strtowrite=(' "'+strlat+','+strlon+'":"'+strval+'",')
    outf.write(strtowrite)    
outf.close()
