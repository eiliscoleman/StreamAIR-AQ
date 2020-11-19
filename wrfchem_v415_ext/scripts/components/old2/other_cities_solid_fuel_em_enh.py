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


#work_root_dir='../'
work_root_dir='/mnt/raid/wrf-chem/wrfchem_v415_ext'
scripts_dir=os.path.join(work_root_dir,'scripts/components')
#ecmwf_path=os.path.join(work_root_dir, 'data/input-wps')
ecmwf_path='/mnt/raid/wrf-chem/ECMWF-op-VOLCEX'
emis_preproc_dir='/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA_v415_ext'
enh_emis_dir=os.path.join(emis_preproc_dir, 'solid_fuel_emis2') 
wrfrundir=os.path.join(work_root_dir,'WRFV3/run')



sf_types=['bit_coal', 'anthracite', 'sod_peat', 'briquettes', 'coke', 'biomass']
sf_annual_ktoe=pd.Series([101., 69., 128.,69., 6., 33.])
sf_tot_ktoe=sf_annual_ktoe.sum()
#proportion of solid fuel burned in homes
sf_prop=sf_annual_ktoe/sf_tot_ktoe
#kg of fuel burned nightly per m2
kg_fuel_consumed=pd.Series([10., 10., 10,10, 10, 10])
#Stats to determine fuel burned
#houses per km2 in Dublin area
##from International Energy Agency: 1 toe = 41.868 gigajoules (GJ)
sf_annual_TJ=sf_annual_ktoe[:]*41.868
sf_calorific_vals=pd.Series([27.84, 27.84, 13.1,18.55, 32.1, 16.])
# following emission factors (kgPM2.5 / TJ) are estimated based on literature: large uncertainty source as it depends on state of the fuel, use of stove/fireplace etc. (see ozgen et al 14) 
sf_ef=pd.Series([30, 30, 60,60, 30, 270])
## Emission time series over 24 hours
emis_hr_weight=pd.Series([0.5,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.2,.5,1.3,1.75,1.25,1,.75,.5,])/8


#RADM species to be adjusted. Note:RADM2 doesn't have BC
#radm_spcs_list={'E_PM_10', "E_PM25I","E_PM25J", "E_ORGI", "E_ORGJ"}
cri_spcs_list={'E_BC_1', 'E_OC_25_10', 'E_OC_DOM', 'E_OIN_10', 'E_OIN_25', 'E_PM_10', 'E_PM25'}
#_____________________________________________________
#in km2
gridcellarea=25*25
#galway
g_city_pop_dens=1475.2
g_county_pop_dens=42
g_city_area=54.2

#gpop_dens=g_city_area/gridcellarea*g_city_pop_dens+(gridcellarea-g_city_area)/gridcellarea*g_county_pop_dens
gpop_dens=g_city_pop_dens
print(gpop_dens)


gpphouse=2.8
ghouse_density=gpop_dens/gpphouse
#cork
gpop_dens=1475.2
gpphouse=2.8
ghouse_density=gpop_dens/gpphouse

cork_city_popdens=1123
ck_house_density=cork_city_popdens/gpphouse
#cork_county_popdens=72


lk_popdens=1591
lk_house_density=lk_popdens/gpphouse

city_list={"Dublin":{"house_density":3500, "lat": 53.426448,"lon":-6.249910},
                   "Galway":{"house_density":ghouse_density, "lat":53.27, "lon":-9.087},
                   "Cork":{"house_density":ck_house_density, "lat":51.89, "lon":-8.4756},
                   "Limerick":{"house_density":lk_house_density, "lat":52.6638, "lon":-8.627}
                   }
# deal with date and time
ii=cdo.showdate(input=os.path.join(work_root_dir,'WRFV3/run/wrfinput_d01'))
i=dt.datetime.strptime(ii[0], '%Y-%m-%d')
# Case: if i=29/02/2020 (a leap year) : the emission preprocessor won't be able to handle this, so set the date a day forward, 
lydate=dt.datetime.strptime('2020-02-26', '%Y-%m-%d')
print(lydate)
if i==lydate:
    print('its dday')
    i=i+timedelta(days=1)
    print('newi', i)
    st,en=(-1,4)
else:
    st,en=(0,5)
    
    
    # copy the file to be readable in python
for n in range(st,en):
    j=i+timedelta(days=n)
    print(j)
    dayemisname=('wrfchemi_d01_'+str(j.year)+'-'+str(j.month).zfill(2)+'-'+str(j.day).zfill(2)+'_00:00:00')
    dayemisname12=('wrfchemi_d01_'+str(j.year)+'-'+str(j.month).zfill(2)+'-'+str(j.day).zfill(2)+'_12:00:00')
    orig_emis_fname=os.path.join(emis_preproc_dir, dayemisname)
    enh_emis_fname=os.path.join(enh_emis_dir, dayemisname)
    print('adding solid fuel to file', orig_emis_fname)
    print('the day increment', n)
    copyfile(orig_emis_fname, enh_emis_fname)
    if n<4:
        IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'1200'+str(j.month).zfill(2)+str(j.day).zfill(2)+'18001'
            #IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'0000'+str(j.month).zfill(2)+str(j.day).zfill(2)+'18001'
    else:
        IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'1200'+str(j.month).zfill(2)+str(j.day).zfill(2)+'12001'
            #IKD_f='IKD'+str(i.month).zfill(2)+str(i.day).zfill(2)+'0000'+str(j.month).zfill(2)+str(j.day).zfill(2)+'12001'
    print('Dealing with date', j, '+days from sim start', n)
    exec("f_day%s=IKD_f" % (n))
    ## convert the file to a netcdf so can be opened by python (this can be skipped if twe get pygrib installed)
    tempout=(IKD_f+'.nc')
    print(tempout)
    exec("print(os.path.join(ecmwf_path, f_day%s))" % (n))
    print(os.path.join(scripts_dir, tempout))
    exec("cdo.copy(input=os.path.join(ecmwf_path, f_day%s), output=os.path.join(scripts_dir, tempout), options = '-f nc',)" % (n))
    print(tempout)
    data=xr.open_dataset(os.path.join(scripts_dir, tempout))
    #lat/lon: data.lat/lon[:].values)
    # Note: we use a different formulation here because ECMWF data is a regular grid
    arraylon=np.asarray(data.lon[:])
    arraylat=np.asarray(data.lat[:])
    wrfinfile=os.path.join(wrfrundir,'wrfinput_d01')
    dwrf=xr.open_dataset(wrfinfile)
    ds=Dataset(enh_emis_fname, 'r+')
    for ss in city_list:
        h_d=city_list[ss]['house_density']
        burning_factor=0.4
        burners_per_m2=h_d*burning_factor*1e-6
        sf_burned_per_m2=burners_per_m2*sf_prop
        #resultant PM emitted
        sf_pm_ugperm2=kg_fuel_consumed*sf_burned_per_m2*sf_calorific_vals*1e3*sf_ef
        #Total PM emitted per s (/3600) from combustion of all fuels over the hours
        # calculate proportions of PM that are in different species. For the future, this  should be done according to the proportions being burned and emitted species, but for the moment we use proportions as determined by Lin during dublin pollution events. 
        E_PM_10_ugperm2pers=emis_hr_weight*sf_pm_ugperm2.sum()/3600
        #PM2.5 = 75% of all PM
        E_PM25_ugperm2pers=E_PM_10_ugperm2pers*.75
        E_BC_1_ugperm2pers=E_PM_10_ugperm2pers*.13
        E_OC_DOM_ugperm2pers=E_PM25_ugperm2pers*.65
        E_OC_25_10_ugperm2pers=E_PM_10_ugperm2pers*.65-E_OC_DOM_ugperm2pers
        E_OIN_10_ugperm2pers=E_PM_10_ugperm2pers*.22
        E_OIN_25_ugperm2pers=E_PM25_ugperm2pers*.22
        
        #link_dest=os.path.join(emis_preproc_dir, 'wrfinput_d01')
        loc_lon, loc_lat=city_list[ss]['lon'],  city_list[ss]['lat']
        ## Get the index of dublin grid cell from the wrfinput file as the produced wrfchemi files don't contain xlat, xlon
        lat=np.asarray(dwrf.XLAT[0,:,:])
        longg=np.asarray(dwrf.XLONG[0,:,:])
        diffarray_lats=lat[:]-[loc_lat]
        diffarray_lons=longg[:]-[loc_lon]
        diffarrayabs=xu.sqrt(xu.square(diffarray_lons)+xu.square(diffarray_lats))
        ixlat,ixlon=np.where(diffarrayabs == np.min(diffarrayabs))
        idxx = {
            'x': ixlon,
            'y': ixlat
            }
        lat_ix=np.argmin((arraylat-loc_lat)**2)
        lon_ix=np.argmin((arraylon-loc_lon)**2)
        print(lon_ix, lat_ix) #verified as corresponding with the location 
        idx = {
            'x': lon_ix,
            'y': lat_ix
            }
        dub_temp=data.var235[:,lat_ix, lon_ix].values-273.15
        print('Todays temp in ',ss, '  ', dub_temp)
        ## firstly, apply fix to NO, NO2 ratios
        ds.variables['E_NO2'][:,0,:,:]=ds.variables['E_NO'][:,0,:,:]*0.22
        ds.variables['E_NO'][:,0,:,:]=ds.variables['E_NO'][:,0,:,:]*0.88
        if dub_temp<12.0:
                    ##Assume linear relationship between temperature and proportion of houses burning solid fuel as secondary heat source: 
                    enh_fac=(12.- dub_temp)/12
                    print('temp below threshold: enhancement factor', enh_fac)
                    #ds=Dataset(enh_emis_fname, 'r+')
                    
                    
                    for demis in cri_spcs_list:
                        ##print(demis, ds.variables[demis].shape)
                        ##exec("print(%s_ugperm2pers.shape)"% (demis))
                        exec("ds.variables['%s'][:,0,%i,%i]=ds.variables['%s'][:,0,%i,%i]+(%s_ugperm2pers[:].values)*enh_fac" % (demis, idxx['y'], idxx['x'], demis, idxx['y'], idxx['x'], demis)) 
            ## get rid of wrfchemi files residing in run directory if they exist
    os.remove(os.path.join(scripts_dir, tempout))
    rundirwrfchemi_i=os.path.join(wrfrundir, dayemisname)
    rundirwrfchemi_i12=os.path.join(wrfrundir, dayemisname12)
    print('unlink these and link new ones', rundirwrfchemi_i)
    if os.path.islink(rundirwrfchemi_i):
            print('remove the existing wrfchemi files in run dir ')
            os.unlink(rundirwrfchemi_i)
    if os.path.islink(rundirwrfchemi_i12):
            print('remove the existing wrfchemi files in run dir ')
            os.unlink(rundirwrfchemi_i12)
    #os.unlink(rundirwrfchemi_i12)
    linkdesti=os.path.join(wrfrundir, dayemisname)
    linkdesti12=os.path.join(wrfrundir, dayemisname12)
    os.symlink(enh_emis_fname, linkdesti)
    # copy emission file, renamed to 12hourly interval
    os.symlink(enh_emis_fname, linkdesti12)

        
