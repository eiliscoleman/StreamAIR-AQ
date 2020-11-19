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

work_root_dir='/mnt/raid/wrf-chem/wrfchem_v39_cri'
scripts_dir=os.path.join(work_root_dir,'scripts/components')

cam_dir='/mnt/raid/wrf-chem/cams'
moz_dir='/mnt/raid/wrf-chem/mozbc'
wrfrundir=os.path.join(work_root_dir,'WRFV3/run')

# link the operational wrfinput file to the emission preprocessor directory. This is necessary to generate the hourly emissions. 
wrfinfile=os.path.join(wrfrundir,'wrfinput_d01')
#link_dest=os.path.join(cam_dir, 'wrfinput_d01')

#dwrf=xr.open_dataset(wrfinfile)

#if os.path.exists(link_dest):
        #os.remove(link_dest)
        #print('removed the input file')
#os.symlink(wrfinfile, link_dest)


#ii=cdo.showdate(input=os.path.join(work_root_dir,'WRFV3/run/wrfinput_d01'))
#i=dt.datetime.strptime(ii[0], '%Y-%m-%d')
i=dt.datetime.now()
# loop to get for 4 day forecast:
for fc in range(0,1):
    fci=i+timedelta(days=fc)
    
    unamepwd='liz.coleman:fj0yVQEV'
    ftpcmd=('ftp://'+unamepwd+'dissemination.ecmwf.int/DATA/CAMS_GLOBAL/'+str(fci.year)+str(fci.month).zfill(2)+str(fci.day).zfill(2)+'00/*fc_pl*aermr10.nc')
    print(ftpcmd)

    #wget  ftp://liz.coleman:fj0yVQEV@dissemination.ecmwf.int/DATA/CAMS_GLOBAL/2019072500/*fc_pl*aermr10.nc
    
    ##define path name for the cam files
    #cam_fname='f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.002.cam.h3.'+str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)+'-00000.nc'
    
    #lcam_fname=('rawcam'+str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)+'.nc')
    #lcam_fpath=os.path.join(cam_dir, lcam_fname)
    #longlcam_fpath=os.path.join(cam_dir, cam_fname)
    #wgetcmd=('wget -P '+cam_dir+' https://www.acom.ucar.edu/waccm/DATA/'+cam_fname)
    #mvcmd=('mv '+longlcam_fpath+' '+lcam_fpath)
    #if not os.path.isfile(lcam_fpath):
        #print('wgetting:', cam_fname)
        #os.system(wgetcmd)
        
        ##exec("mv %s %s " % (longlcam_fpath, lcam_fpath))
        #os.system(mvcmd)
    #else:
        #print('already have data for', str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)) 




##merge all the raw cams files
#mergecmd=('ncrcat '+cam_dir+'/rawcam* '+cam_dir+'/cam_cat001.nc')
#print('merging the rawcam files (this will take a few minutes)')
#os.system(mergecmd)

## remove the camsfile for the first timestep as this won't be needed for the next simulation
#rm_cmd=('rm '+cam_dir+'/rawcam*')


        
