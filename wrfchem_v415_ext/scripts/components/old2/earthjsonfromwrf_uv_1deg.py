from __future__ import print_function, unicode_literals

import os
import json

import numpy as np

#import matplotlib.mlab as ml
from matplotlib.mlab import griddata

from netCDF4 import Dataset
import argparse
from argparse import RawDescriptionHelpFormatter
from numpy.random import uniform, seed
#from scipy.interpolate import griddata

parser = argparse.ArgumentParser(description = """earthjsonfromwrf.py
Creates a json in the format required by earth[1] from
Weather Research and Forecasting (WRF) outputs. The script
was made for and tested with Python3, but to be compatible
with Python2. See the example for a typical usage
[1]http://github.com/cambecc/earth
""", epilog = """
Example:
    Create input for earth from wrfinput_d01. This script was run
    from the root directory of the earth-master so that the outroot
    is in the standard weather location. The second line starts
    a node server with earth loaded. The third line (on a mac)
    starts the default webbrowser with the data loaded.
    
    $ python earthjsonfromwrf.py ../wrf/wrfinput_d01 public/data/weather/
    Latitude Error: 0.015507
    Latitude Sum Error: 0.00580698
    Longitude Error: 0.0334339
    Longitude Sum Error: 0.0344685
    Add "#2015/06/23/1800Z/wind/surface/level/orthographic" to url to see this time
    $ node dev-server.js 8089 &
    $ open http://localhost:8089/#2015/06/23/1800Z/wind/surface/level/orthographic
""", formatter_class = RawDescriptionHelpFormatter)
parser.add_argument('--lat', default = 'XLAT', dest = 'latitude_name', help = 'Name for latitude variable')
parser.add_argument('--lon', default = 'XLONG', dest = 'longitude_name', help = 'Name for longitude variable')
parser.add_argument('-u', '--ucomponent', default = 'u10', dest = 'ucomponent', help = 'Name for u-component of wind')
parser.add_argument('-v', '--vcomponent', default = 'v10', dest = 'vcomponent', help = 'Name for v-component of wind')

parser.add_argument('wrffile', help = 'Path to wrfout or wrfinput file that contains latitude_name, longitude_name, ucomponent, and vcomponent')
parser.add_argument('outroot', help = 'Path for output', default = '../public/data/weather/')
args = parser.parse_args()

outroot = args.outroot

# Copied a bunch of meta data. Some should (and has not)
# been updated. For example, should the earth radius be updated?
# If it is part of the display, maybe not. If it is part of the
# data meta-data, it should be set to 6,370,000.0
uhdr = {"header":{"discipline":0,"disciplineName":"Meteorological products","gribEdition":2,"gribLength":131858,"center":0,"centerName":"WRF OUTPUT","subcenter":0,"refTime":"2014-01-31T00:00:00.000Z","significanceOfRT":1,"significanceOfRTName":"Start of forecast","productStatus":0,"productStatusName":"Operational products","productType":1,"productTypeName":"Forecast products","productDefinitionTemplate":0,"productDefinitionTemplateName":"Analysis/forecast at horizontal level/layer at a point in time","parameterCategory":2,"parameterCategoryName":"Momentum","parameterNumber":2,"parameterNumberName":"U-component_of_wind","parameterUnit":"m.s-1","genProcessType":2,"genProcessTypeName":"Forecast","forecastTime":3,"surface1Type":103,"surface1TypeName":"Specified height level above ground","surface1Value":10,"surface2Type":255,"surface2TypeName":"Missing","surface2Value":0,"gridDefinitionTemplate":0,"gridDefinitionTemplateName":"Latitude_Longitude","numberPoints":65160,"shape":6,"shapeName":"Earth spherical with radius of 6,371,229.0 m","gridUnits":"degrees","resolution":48,"winds":"true","scanMode":0,"nx":360,"ny":181,"basicAngle":0,"subDivisions":0,"lo1":0,"la1":90,"lo2":359,"la2":-90,"dx":1,"dy":1}}
vhdr = {"header":{"discipline":0,"disciplineName":"Meteorological products","gribEdition":2,"gribLength":131858,"center":0,"centerName":"WRF OUTPUT","subcenter":0,"refTime":"2014-01-31T00:00:00.000Z","significanceOfRT":1,"significanceOfRTName":"Start of forecast","productStatus":0,"productStatusName":"Operational products","productType":1,"productTypeName":"Forecast products","productDefinitionTemplate":0,"productDefinitionTemplateName":"Analysis/forecast at horizontal level/layer at a point in time","parameterCategory":2,"parameterCategoryName":"Momentum","parameterNumber":3,"parameterNumberName":"V-component_of_wind","parameterUnit":"m.s-1","genProcessType":2,"genProcessTypeName":"Forecast","forecastTime":3,"surface1Type":103,"surface1TypeName":"Specified height level above ground","surface1Value":10,"surface2Type":255,"surface2TypeName":"Missing","surface2Value":0,"gridDefinitionTemplate":0,"gridDefinitionTemplateName":"Latitude_Longitude","numberPoints":65160,"shape":6,"shapeName":"Earth spherical with radius of 6,371,229.0 m","gridUnits":"degrees","resolution":48,"winds":"true","scanMode":0,"nx":360,"ny":181,"basicAngle":0,"subDivisions":0,"lo1":0,"la1":90,"lo2":359,"la2":-90,"dx":1,"dy":1}}

data = [uhdr, vhdr]
newf = Dataset(args.wrffile)
lat = newf.variables[args.latitude_name][0]
lon = newf.variables[args.longitude_name][0]
#dys = np.diff(lat, axis = 0).mean(1)
#dy = float(dys.mean())
#print('Latitude Error:', np.abs((dy / dys) - 1).max())
#print('Latitude Sum Error:', (dy / dys - 1).sum())
#dxs = np.diff(lon, axis = 1).mean(0)
#dx = float(dxs.mean())
#print('Longitude Error:', np.abs(dx / dxs - 1).max())
#print('Longitude Sum Error:', (dx / dxs - 1).sum())
#nx = float(lon.shape[1])
#ny = float(lat.shape[0])
#la1 = float(lat[-1, -1])
#la2 = float(lat[0, 0])
#lo1 = float(lon[0, 0])
#lo2 = float(lon[-1, -1])

#x  = np.arange(lo1,lo2,dx)
#print(x)
#y  = np.arange(la2,la1,dy)

nx = 41
ny = 29
dx = 1
dy = 1
la1 = 62.75
la2 = 34.75
lo1 = -21.75
lo2 = 18.25

#x1  = np.arange(lo1,lo2,dx1)
#print (x1)

#y1  = np.arange(la2,la1,dy1)
#x1, y1 = np.meshgrid(x1, y1)
#print (y1)
#print(lat) 
#u=newf.variables[args.ucomponent]
#u=newf.variables[args.ucomponent][0]
#print(u) 
#test
#npts = 200
#x = uniform(-2, 2, npts)
#y = uniform(-2, 2, npts)
#z = x*np.exp(-x**2 - y**2)
# define grid.
#xi = np.linspace(-2.1, 2.1, 100)
#yi = np.linspace(-2.1, 2.1, 200)
# grid the data.
#zi = griddata(x, y, z, xi, yi, interp='linear')
#print(z)

#exit()

# interpolation
#ui = ml.griddata(x,y,u,x1,y1,interp='linear') 
#zi = griddata(x, y, u, x1, y1, interp='linear')
#zi = griddata((x,y),u,(x1,y1),method='linear')
#print(zi)
#exit()

#times = newf.variables['Times'][:].copy().view('S19')
times = newf.variables['Times']
print(times[0])
#exit()

for ti, time in enumerate(times):
    #2012/02/07/0100Z/wind/surface/level/orthographic=-74.01,4.38,29184
#    datestr = (time[0]).decode('ascii').replace('-', '/')
#    timestr = (time[0]).decode('ascii') + '00'
#    timestr = 1603072017
#    timestr = times[0]
    timestr = 9999
#    print('Add "#' + datestr + '/' + timestr + 'Z/wind/surface/level/orthographic" to url to see this time')
#    dirpath = os.path.join(args.outroot, *datestr.split('/'))
    dirpath = args.outroot

    print (dirpath)
#    os.makedirs(dirpath, exist_ok = True)
    outpath = os.path.join(dirpath, '%s-wind-surface-level-gfs-1.0.json' % (timestr,))

    for u0_or_v1 in [0, 1]:
        # Update header data for some properties
        # that are now to affect.
        h = data[u0_or_v1]['header']
        h['la1'] = la1
        h['la2'] = la2
        h['lo1'] = lo1
        h['lo2'] = lo2
        h['nx'] = nx
        h['ny'] = ny
        h['dx'] = dx
        h['dy'] = dy
        h['forecastTime'] = 0
        #h['refTime'] = time[0].decode('ascii').replace('_', 'T') + '.000Z'
        #"2014-01-31T00:00:00.000Z"
        
        h['gribLength'] = 1538 + nx * ny * 2
        if u0_or_v1 == 0:
            data[u0_or_v1]['data'] = newf.variables[args.ucomponent][ti].ravel().tolist()
        elif u0_or_v1 == 1:
            data[u0_or_v1]['data'] = newf.variables[args.vcomponent][ti].ravel().tolist()
    if ti == 0:
        outf = open(outpath, 'w')
        json.dump(data, outf)
        outf.close()
        
    outf = open(outpath, 'w')
    json.dump(data, outf)
    outf.close()
