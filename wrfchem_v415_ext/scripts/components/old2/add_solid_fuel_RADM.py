#!/usr/bin/env python3
import datetime as dt
import glob as g
import json
import numpy as np
import os
from netCDF4 import Dataset
#import xarray as xr
from datetime import datetime as dt
from datetime import timedelta
from shutil import copyfile
from numpy import loadtxt
import pandas as pd

#Input file directory
dirpath='/ichec/work/ngear011c/liz/wrfchem/aq_exps/emit/RADM'


#datain=reference file that will be altered
datain='wrfchemi_d01_2017-01-20_00:00:00'
ndate='2017-01-20'
ddt=dt.strptime(ndate, '%Y-%m-%d')
reffpath=os.path.join(dirpath, datain)
ndays=5
#houses per km2 in Dublin area
house_density=3500
#proportion of houses burning solid fuel
burning_factor=0.4

dirpath2=('/ichec/work/ngear011c/liz/wrfchem/aq_exps/emit/RADM/perturbed'+str(burning_factor))

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
kg_fuel_consumed=pd.Series([15., 15., 15,15, 15, 15])*2

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
pm25_ugperm2pers=E_PM_10_ugperm2pers*.75
E_PM25I_ugperm2pers=pm25_ugperm2pers*.15
E_PM25J_ugperm2pers=pm25_ugperm2pers*.85
E_ORGI_ugperm2pers=E_PM25I_ugperm2pers*.65
E_ORGJ_ugperm2pers=E_PM25J_ugperm2pers*.65
print('orgj', E_ORGJ_ugperm2pers[:])

#RADM species to be adjusted. Note:RADM2 doesn't have BC
radm_spcs_list={'E_PM_10', "E_PM25I","E_PM25J", "E_ORGI", "E_ORGJ"}

print(ddt)

for ii in range(0,ndays):
    currentday= ddt + timedelta(days=ii)
    datestr=dt.strftime(currentday, '%Y-%m-%d')
    #%H:%M:%S')
    newfname=('wrfchemi_d01_' + datestr +'_00:00:00')
    nfpath=os.path.join(dirpath2,newfname )
    print(nfpath)
    copyfile(reffpath, nfpath)
    ds=Dataset(nfpath, 'r+')
    times=ds.variables['Times']
    #print(ds.variables['E_CO'].shape)
    print('******************************************************************************')
    #print(ds.variables.keys())
    print('PM10-shape', ds.variables['E_PM_10'].shape)
    ##overwrite
    for i in range(0,24):
        #alter the year:
        ds.variables['Times'][i,0:4]=[datestr[0], datestr[1], datestr[2], datestr[3]]
        #alter the month:
        ds.variables['Times'][i,5:7]= [datestr[5],datestr[6]]
        #alter the day:
        ds.variables['Times'][i,8:10]= [datestr[8],datestr[9]]
    for demis in radm_spcs_list:
        print(demis, ds.variables[demis].shape)
        exec("print(%s_ugperm2pers.shape)"% (demis))
        #exec("ds.variables['%s'][:,0,62,40]=ds.variables['%s'][:,0,62,40]+%s_ugperm2pers[:]" % (demis, demis, demis))
        exec("ds.variables['%s'][:,0,62,40]=ds.variables['%s'][:,0,62,40]+%s_ugperm2pers[:].values" % (demis, demis, demis))    
    ds.close()


