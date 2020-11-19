from __future__ import print_function, unicode_literals

import os
import json
import numpy as np
from matplotlib.mlab import griddata
from netCDF4 import Dataset
import argparse
from argparse import RawDescriptionHelpFormatter
from numpy.random import uniform, seed

parser = argparse.ArgumentParser(description = "", epilog = "", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('--lat', default = 'XLAT', dest = 'latitude_name', help = 'Name for latitude variable')
parser.add_argument('--lon', default = 'XLONG', dest = 'longitude_name', help = 'Name for longitude variable')
parser.add_argument('-t2', '--temperature', default = 't2', dest = 'temperature', help = 'Name for surface temperaure')

parser.add_argument('wrffile', help = 'Path to wrfout or wrfinput file that contains latitude_name, longitude_name, t2')
parser.add_argument('outroot', help = 'Path for output', default = '../public/data/weather/')
#parser.add_argument('o3_conc',default = 'o3_concentration', help = 'Name for O3 conc')
args = parser.parse_args()
outroot = args.outroot
uhdr = {}
vhdr = {}
data = [uhdr, vhdr]
newf = Dataset(args.wrffile)

nx = 210
ny = 150
dx = 0.25
dy = 0.25
la1 = 72
la2 = 34.5
lo1 = -22.75
lo2 = 29.75
gridsize=210*150
longs=np.arange(lo1,lo2+dx,dx)
latit=np.arange(la1, la2-dy,-dy)
times = newf.variables['Times']
timestr = 9999
dirpath = args.outroot
outpath = os.path.join(dirpath, '%s-t2-surface-level-gfs-1.0.json' % (timestr,))
gc=0
outf = open(outpath, 'w')
outf.write('{ ')
for xways in range(0, nx):
  for yways in range(0,ny):
    strlat=str(latit[yways])
    strlon=str(longs[xways])
    rawdatavals=newf.variables[args.temperature][0,yways,xways]
    convdatavals=rawdatavals-273.15
    rounddatavals=round(convdatavals, 1)
    strval=str(rounddatavals)
    gc=gc+1
    if gc==gridsize:
      strtowrite=(' "'+strlat+', '+strlon+'":"'+strval+'" }')
    else:
      strtowrite=(' "'+strlat+','+strlon+'":"'+strval+'",')
    outf.write(strtowrite)    
outf.close()
