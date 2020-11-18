#!/usr/bin/env python3
import datetime as dt
from subprocess import Popen, PIPE
from datetime import timedelta
from shutil import copyfile
import glob as g
import json
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import os
import re
import fileinput
from cdo import *
cdo=Cdo()
#Find the input data from which we deduce ampunt of fuel burned
#_____________________________________________________


cam_dir='/mnt/raid/wrf-chem/cam-chem'
i=dt.datetime.now()
# loop to get for 4 day forecast:
for fc in range(0,5):
    fci=i+timedelta(days=fc)
    #https://www.acom.ucar.edu/waccm/DATA/f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.001 uses FINN emissions
    #https://www.acom.ucar.edu/waccm/DATA/f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.002 uses QFED emissions 
    
    #define path name for the cam files
    cam_fname='f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.002.cam.h3.'+str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)+'-00000.nc'
    
    lcam_fname=('rawcam'+str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)+'.nc')
    lcam_fpath=os.path.join(cam_dir, lcam_fname)
    longlcam_fpath=os.path.join(cam_dir, cam_fname)
    wgetcmd=('wget -P '+cam_dir+' https://www.acom.ucar.edu/waccm/DATA/'+cam_fname)
    mvcmd=('mv '+longlcam_fpath+' '+lcam_fpath)
    if not os.path.isfile(lcam_fpath):
        print('wgetting:', cam_fname)
        os.system(wgetcmd)
        os.system(mvcmd)
    else:
        print('already have data for', str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)) 


#merge all the raw cams files
mergecmd=('ncrcat '+cam_dir+'/rawcam* '+cam_dir+'/cam_cat001.nc')
print('merging the rawcam files (this will take a few minutes)')
os.system(mergecmd)

# remove the camsfile for the first timestep as this won't be needed for the next simulation
rm_cmd=('rm '+cam_dir+'/rawcam*')
