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
#Stats to determine fuel burned
#houses per km2 in Dublin area
house_density=3500
#proportion of houses burning solid fuel
burning_factor=0.4
burners_per_m2=house_density*burning_factor*1e-6

#solid fuel national consumtion
sf_types=['bit_coal', 'anthracite', 'sod_peat', 'briquettes', 'coke', 'biomass']
sf_annual_ktoe=pd.Series([101., 69., 128.,69., 6., 33.])
sf_tot_ktoe=sf_annual_ktoe.sum()
#proportion of solid fuel burned in homes
sf_prop=sf_annual_ktoe/sf_tot_ktoe
# houses per m2 burning the fuel sf_types
sf_burned_per_m2=burners_per_m2*sf_prop
#kg of fuel burned nightly per m2
kg_fuel_consumed=pd.Series([10., 10., 10,10, 10, 10])

#from International Energy Agency: 1 toe = 41.868 gigajoules (GJ)
sf_annual_TJ=sf_annual_ktoe[:]*41.868
#calorific vals MJ/kg
sf_calorific_vals=pd.Series([27.84, 27.84, 13.1,18.55, 32.1, 16.])
# following emission factors (kgPM2.5 / TJ) are estimated based on literature: large uncertainty source as it depends on state of the fuel, use of stove/fireplace etc. (see ozgen et al 14) 
sf_ef=pd.Series([30, 30, 60,60, 30, 270])

#resultant PM emitted
sf_pm_ugperm2=kg_fuel_consumed*sf_burned_per_m2*sf_calorific_vals*1e3*sf_ef


# Emission time series over 24 hours
emis_hr_weight=pd.Series([0.5,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.2,.5,1.3,1.75,1.25,1,.75,.5,])/8

#Total PM emitted per s (/3600) from combustion of all fuels over the hours
E_PM_10_ugperm2pers=emis_hr_weight*sf_pm_ugperm2.sum()/3600
# calculate proportions of PM that are in different species. For the future, this  should be done according to the proportions being burned and emitted species, but for the moment we use proportions as determined by Lin during dublin pollution events. 
#PM2.5 = 75% of all PM
E_PM25_ugperm2pers=E_PM_10_ugperm2pers*.75
E_BC_1_ugperm2pers=E_PM_10_ugperm2pers*.13
E_OC_DOM_ugperm2pers=E_PM25_ugperm2pers*.65
E_OC_25_10_ugperm2pers=E_PM_10_ugperm2pers*.65-E_OC_DOM_ugperm2pers
E_OIN_10_ugperm2pers=E_PM_10_ugperm2pers*.22
E_OIN_25_ugperm2pers=E_PM25_ugperm2pers*.22

#RADM species to be adjusted. Note:RADM2 doesn't have BC
#radm_spcs_list={'E_PM_10', "E_PM25I","E_PM25J", "E_ORGI", "E_ORGJ"}
cri_spcs_list={'E_BC_1', 'E_OC_25_10', 'E_OC_DOM', 'E_OIN_10', 'E_OIN_25', 'E_PM_10', 'E_PM25'}
#_____________________________________________________




work_root_dir='/mnt/raid/wrf-chem/wrfchem_v39_cams_cri'
scripts_dir=os.path.join(work_root_dir,'scripts/components')
#ecmwf_path=os.path.join(work_root_dir, 'data/input-wps')
ecmwf_path='/mnt/raid/wrf-chem/ECMWF-op-VOLCEX'
emis_preproc_dir='/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA_cams_crimech'
enh_emis_dir=os.path.join(emis_preproc_dir, 'solid_fuel_emis') 
wrfrundir=os.path.join(work_root_dir,'WRFV3/run')

# link the operational wrfinput file to the emission preprocessor directory. This is necessary to generate the hourly emissions. 
wrfinfile=os.path.join(wrfrundir,'wrfinput_d01')
link_dest=os.path.join(emis_preproc_dir, 'wrfinput_d01')
dub_lon, dub_lat= -6.249910,  53.426448
dwrf=xr.open_dataset(wrfinfile)
# Get the index of dublin grid cell from the wrfinput file as the produced wrfchemi files don't contain xlat, xlon
lat=np.asarray(dwrf.XLAT[0,:,:])
longg=np.asarray(dwrf.XLONG[0,:,:])

diffarray_lats=lat[:]-[dub_lat]
diffarray_lons=longg[:]-[dub_lon]
print(diffarray_lats.shape, diffarray_lons.shape)
#diffarraylon1=longg-dub_lon
diffarrayabs=xu.sqrt(xu.square(diffarray_lons)+xu.square(diffarray_lats))
ixlat,ixlon=np.where(diffarrayabs == np.min(diffarrayabs))
idxx = {
            'x': ixlon,
            'y': ixlat
        }
print(idxx)

#if os.path.exists(link_dest):
        #os.remove(link_dest)
        #print('removed the input file')
#os.symlink(wrfinfile, link_dest)
#ii='datestr'
ii=cdo.showdate(input=os.path.join(work_root_dir,'WRFV3/run/wrfinput_d01'))
i=dt.datetime.strptime(ii[0], '%Y-%m-%d')
##i=dt.datetime.strptime('2018-11-02', '%Y-%m-%d')
##2. Update the namelist with current date
#upd_file=os.path.join(emis_preproc_dir, 'namelist')

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

# The namelist file now updated, run the preprocessor
# First: change to directory 
os.chdir(emis_preproc_dir)
## comment 15 lines below to run without emission preprocessor
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
##-------------------------------------------------------------------------

####-------------------------------------------------------------------------
##New section  of code:
##if startdate month ne enddate month:
## update the namelist with the new month
## re run main

#i=dt.datetime.strptime(ii[0], '%Y-%m-%d')
#iend=i+timedelta(days=4)
#print(iend, i.month, iend.month)
#if i.month != iend.month:
    #print("my good god...it's the turn of the month")
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
            #new_nl_month=iend.month
            #new_nl_month_line=line.split('=')[0]+'= '+str(new_nl_month)+',\n'
            #replace_line(upd_file, ll, new_nl_month_line)
        
    #f.close()

    #prep_exec='./main'
    #print(prep_exec, 'for the new month')
    #p = subprocess.Popen(prep_exec, stdout=PIPE)
    #output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    #rc=p.returncode
    #os.chdir(scripts_dir)
    #rundirwrfchemi_nodate=os.path.join(wrfrundir, 'wrfchemi_d01_')
    #rundirwrfchemi_i=(rundirwrfchemi_nodate+str(i.year)+'-'+str(i.month).zfill(2)+'-'+str(i.day).zfill(2)+'_00:00:00')
    #rundirwrfchemis=os.path.join(wrfrundir, 'wrfchemi_d01_')
    #wrfchemis=os.path.join(emis_preproc_dir, 'wrfchemi_d*')


# Alter the emissions according to temperature for each day of the simulation, evaluated at 18:00
# revert range below to 0,5


for n in range(0,4):
    j=dt.datetime.now()+timedelta(days=n)
    dayemisname=('wrfchemi_d01_'+str(j.year)+'-'+str(j.month).zfill(2)+'-'+str(j.day).zfill(2)+'_00:00:00')
    orig_emis_fname=os.path.join(emis_preproc_dir, dayemisname)
    enh_emis_fname=os.path.join(enh_emis_dir, dayemisname)
    print('adding solid fuel to file', orig_emis_fname)
    copyfile(orig_emis_fname, enh_emis_fname)
    if n<3:
         IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'1200'+str(j.month).zfill(2)+str(j.day).zfill(2)+'18001'
    else:
        print('the day increment', n)
        IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'1200'+str(j.month).zfill(2)+str(j.day).zfill(2)+'12001'
    print('Dealing with date', j, '+days from sim start', n)
    exec("f_day%s=IKD_f" % (n))
    # convert the file to a netcdf so can be opened by python (this can be skipped if twe get pygrib installed)
    tempout=(IKD_f+'.nc')
    print(tempout)
    exec("print(os.path.join(ecmwf_path, f_day%s))" % (n))
    print(os.path.join(scripts_dir, tempout))
    exec("cdo.copy(input=os.path.join(ecmwf_path, f_day%s), output=os.path.join(scripts_dir, tempout), options = '-f nc',)" % (n))
    print(tempout)
    #data=xr.open_dataset('/mnt/raid/wrf-chem/wrfchem_v39_cri/scripts/components/IKD09250000092518001.nc')
    data=xr.open_dataset(os.path.join(scripts_dir, tempout))
        #lat/lon: data.lat/lon[:].values)
        # Note: we use a different formulation here because ECMWF data is a regular grid
    arraylon=np.asarray(data.lon[:])
    arraylat=np.asarray(data.lat[:])
    lat_ix=np.argmin((arraylat-dub_lat)**2)
    lon_ix=np.argmin((arraylon-dub_lon)**2)
    print(lon_ix, lat_ix) #verified as corresponding with the location of Dublin
    idx = {
            'x': lon_ix,
            'y': lat_ix
            }
    dub_temp=data.var235[:,lat_ix, lon_ix].values-273.15
    print('Todays temp in Dublin', dub_temp)
    ds=Dataset(enh_emis_fname, 'r+')
    # firstly, apply fix to NO, NO2 ratios
    ds.variables['E_NO2'][:,0,:,:]=ds.variables['E_NO'][:,0,:,:]*0.22
    ds.variables['E_NO'][:,0,:,:]=ds.variables['E_NO'][:,0,:,:]*0.88
    if dub_temp<12.0:
                    #Assume linear relationship between temperature and proportion of houses burning solid fuel as secondary heat source: 
                    enh_fac=(12.- dub_temp)/12
                    print('temp below threshold: enhancement factor', enh_fac)
                    
                    
                    
                    for demis in cri_spcs_list:
                        #print(demis, ds.variables[demis].shape)
                        #exec("print(%s_ugperm2pers.shape)"% (demis))
                        exec("ds.variables['%s'][:,0,%i,%i]=ds.variables['%s'][:,0,%i,%i]+(%s_ugperm2pers[:].values)*enh_fac" % (demis, idxx['y'], idxx['x'], demis, idxx['y'], idxx['x'], demis)) 
    # get rid of wrfchemi files residing in run directory if they exist
    #os.remove('/mnt/raid/wrf-chem/wrfchem_v39_cri/scripts/components/IKD09250000092518001.nc')
#os.remove(os.path.join(scripts_dir, tempout))
#rundirwrfchemi_i=os.path.join(wrfrundir, dayemisname)
#print('unlink these and link new ones', rundirwrfchemi_i)
#if os.path.islink(rundirwrfchemi_i):
            #print('remove the existing wrfchemi files in run dir ')
            #os.unlink(rundirwrfchemi_i)
#linkdesti=os.path.join(wrfrundir, dayemisname)
#os.symlink(enh_emis_fname, linkdesti)

        
