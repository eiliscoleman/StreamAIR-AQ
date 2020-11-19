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

cam_dir='/mnt/raid/wrf-chem/cam-chem'
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
        
        #exec("mv %s %s " % (longlcam_fpath, lcam_fpath))
        os.system(mvcmd)
    else:
        print('already have data for', str(fci.year)+'-'+str(fci.month).zfill(2)+'-'+str(fci.day).zfill(2)) 


#merge all the raw cams files
mergecmd=('ncrcat '+cam_dir+'/rawcam* '+cam_dir+'/cam_cat001.nc')
print('merging the rawcam files (this will take a few minutes)')
os.system(mergecmd)

# remove the camsfile for the first timestep as this won't be needed for the next simulation
rm_cmd=('rm '+cam_dir+'/rawcam*')
#os.system(rm_cmd)











#moz_namelist=os.path.join(moz_dir, 'CRIMECH-mosaic8bin-fromWMACC-noPM.inp')

#def replace_line(file_name, line_num, text):
        #lines = open(file_name, 'r').readlines()
        #lines[line_num] = text
        #out = open(file_name, 'w')
        #out.writelines(lines)
        #out.close()
##replace_line(upd_file, 0, 'Mage\n')
#f=open(upd_file ,'r')
#namelist_data=f.readlines()

#for ll, line in enumerate(namelist_data):
    #if 'year' in line:
        #nl_year = line.split('=')[-1]
        #new_nl_year=i.year
        #new_nl_year_line=line.split('=')[0]+'= '+str(new_nl_year)+',\n'
        #replace_line(upd_file, ll, new_nl_year_line)
    #if 'month' in line:
        #print('month line', line, 'll', ll)
        #nl_month = line.split('=')[-1]
        #new_nl_month=i.month
        #new_nl_month_line=line.split('=')[0]+'= '+str(new_nl_month)+',\n'
        #replace_line(upd_file, ll, new_nl_month_line)
        
#f.close()

## The namelist file now updated, run the preprocessor
## First: change to directory 
#os.chdir(emis_preproc_dir)
### comment 15 lines below to run without emission preprocessor
###-------------------------------------------------------------------
## if wrfchemi exists in the emission preprocessor directory, or enhacend emission directory clean them out
#for iii in os.listdir(emis_preproc_dir):
    #if re.match('wrfchemi_d01', iii):
        ##print('clean these out', os.path.join(emis_preproc_dir, iii))
        #os.remove(os.path.join(emis_preproc_dir, iii))
#for iv in os.listdir(enh_emis_dir):
    #if re.match('wrfchemi_d01', iv):
        ##print('clean these out too', os.path.join(enh_emis_dir, iv))
        #os.remove(os.path.join(enh_emis_dir, iv))
#prep_exec='./main'
#print(prep_exec)

#p = subprocess.Popen(prep_exec, stdout=PIPE)
#output, err = p.communicate(b"input data that is passed to subprocess' stdin")
#rc=p.returncode
###-------------------------------------------------------------------------
#os.chdir(scripts_dir)
#rundirwrfchemi_nodate=os.path.join(wrfrundir, 'wrfchemi_d01_')
#rundirwrfchemi_i=(rundirwrfchemi_nodate+str(i.year)+'-'+str(i.month).zfill(2)+'-'+str(i.day).zfill(2)+'_00:00:00')
#rundirwrfchemis=os.path.join(wrfrundir, 'wrfchemi_d01_')
#wrfchemis=os.path.join(emis_preproc_dir, 'wrfchemi_d*')


## Alter the emissions according to temperature for each day of the simulation, evaluated at 18:00
## revert range below to 0,5
#for n in range(0,5):
    #j=dt.datetime.now()+timedelta(days=n)
    #dayemisname=('wrfchemi_d01_'+str(j.year)+'-'+str(j.month).zfill(2)+'-'+str(j.day).zfill(2)+'_00:00:00')
    #orig_emis_fname=os.path.join(emis_preproc_dir, dayemisname)
    #enh_emis_fname=os.path.join(enh_emis_dir, dayemisname)
    #print('adding solid fuel to file', orig_emis_fname)
    #copyfile(orig_emis_fname, enh_emis_fname)
    #if n<4:
        #print('Dealing with date', j, '+days from sim start', n)
        #exec("f_day%s='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'0000'+str(j.month).zfill(2)+str(j.day).zfill(2)+'18001'" % (n))
        ## convert the file to a netcdf so can be opened by python (this can be skipped if twe get pygrib installed)
        #tempout=('IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'0000'+str(j.month).zfill(2)+str(j.day).zfill(2)+'18001.nc')
        #exec("cdo.copy(input=os.path.join(ecmwf_path, f_day%s), output=os.path.join(scripts_dir, tempout), options = '-f nc',)" % (n))
        #data=xr.open_dataset(tempout)
        ##lat/lon: data.lat/lon[:].values)
        ## Note: we use a different formulation here because ECMWF data is a regular grid
        #arraylon=np.asarray(data.lon[:])
        #arraylat=np.asarray(data.lat[:])
        #lat_ix=np.argmin((arraylat-dub_lat)**2)
        #lon_ix=np.argmin((arraylon-dub_lon)**2)
        #print(lon_ix, lat_ix) #verified as corresponding with the location of Dublin
        #idx = {
            #'x': lon_ix,
            #'y': lat_ix
            #}
        #dub_temp=data.var235[:,lat_ix, lon_ix].values-273.15
        #print('Todays temp in Dublin', dub_temp)
        #if dub_temp<12.0:
                    ##Assume linear relationship between temperature and proportion of houses burning solid fuel as secondary heat source: 
                    #enh_fac=(12.- dub_temp)/12
                    #print('temp below threshold: enhancement factor', enh_fac)
                    #ds=Dataset(enh_emis_fname, 'r+')
                    #for demis in radm_spcs_list:
                        ##print(demis, ds.variables[demis].shape)
                        ##exec("print(%s_ugperm2pers.shape)"% (demis))
                        #exec("ds.variables['%s'][:,0,%i,%i]=ds.variables['%s'][:,0,%i,%i]+(%s_ugperm2pers[:].values)*enh_fac" % (demis, idxx['y'], idxx['x'], demis, idxx['y'], idxx['x'], demis)) 
            ## get rid of wrfchemi files residing in run directory if they exist
        #os.remove(tempout)
    #rundirwrfchemi_i=os.path.join(wrfrundir, dayemisname)
    #print('unlink these and link new ones', rundirwrfchemi_i)
    #if os.path.islink(rundirwrfchemi_i):
            #print('remove the existing wrfchemi files in run dir ')
            #os.unlink(rundirwrfchemi_i)
    #linkdesti=os.path.join(wrfrundir, dayemisname)
    #os.symlink(enh_emis_fname, linkdesti)

        
