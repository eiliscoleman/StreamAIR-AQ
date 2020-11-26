# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 12:55:11 2016

@author: fqu12abu
"""
#This script is to create NetCDF files, with emissions from the latest version of the TNO_MACC 
#inventory, which can be read by the UoM TNO emissions preprocessor. 
#It produces individual netcdf files per pollutant (NOx, CO, NMVOC, SO2, NH3, PM10, PM25)  grouped by sectors 
#(pow, res, inc, pei, exf, sol, tra1-2-3-4-5, nrt, was, and agr). 
#The script does not create all the files at the same time, so the pollutant name needs to be changed in order 
#to create a new file. 
#It is important to note that the TNO emissions preprocessor requests individual files from ship emissions also. 
#Ship emissions are no longer compiled in a supplementary inventory but incorporated in the new version of TNO-MACC. 
#Therefore, currently, I produce empty files for ship emissions for each pollutant in order to force the program 
#to run. The ideal thing would be to modify the routine in the program 

##Libraries used
import netCDF4 as nc
from netCDF4 import Dataset as ds
import numpy as np
import xray
from monthdelta import monthdelta
from datetime import datetime
from netCDF4 import num2date, date2num
dates = []


############(A) Open original TNO emissions inventory################
#tno_2011 = 'path to the TNO_MACCII'  #input file
tno_2011 = '/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/TNO_MACC/TNO_MACC_III_emissions_2011.nc'
tno = nc.Dataset(tno_2011,'r')

#######################(B)Create new netcdf file#################################
#1) Create an empty NetCDF file 
#dataset = ds('Path_directory_to_place_the_file','w','r', format='NETCDF4_CLASSIC') #output file
dataset = ds('/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/OUT/tno_ship_co.nc','w','r', format='NETCDF4_CLASSIC') #output file

#2)Create the dimensions 

lat = dataset.createDimension('lat', 672) #number of latitudes
lon = dataset.createDimension('lon', 720) #number of longitudes
time = dataset.createDimension('time', None)

#3)Create variables
######1D variables first######

lat = dataset.createVariable('lat',np.float32, ('lat'), fill_value=False)
lon = dataset.createVariable('lon',np.float32, ('lon'), fill_value=False)
time = dataset.createVariable('time',np.float32, ('time'), fill_value=False)

######3D variables (SNAP sectors) ########

ships = dataset.createVariable('ships',np.float64, ('time','lat','lon'), fill_value=False)


#4) Add attributes units to 1D and 3D variables
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time.units = 'days since 1900-01-01 00:00'
time.calendar ='gregorian'
time.long_name = 'Time'

ships.units = 'Tg?'
ships.long_name = 'International shipping'
 
                          
############# (C) Append data from the original emissions files to the new NetCDF file, example for NOx #######

#1)Extract and add latitude  and longitude values 

data = xray.open_dataset(tno_2011)  #Open NetCDF file (TNO-MACCII 2011) with x-ray library
#print data
print(data)

latitude_tno = tno['latitude'][:]  # Extract latitude from the original file
longitude_tno = tno['longitude'][:] 

lat [:] = latitude_tno  #Append latitude to the new file
lon [:] = longitude_tno   

#2)Extract other variables (NOx, co, etc) using the NETCDF4 library########
##### i) Create a 3D array with zeros
co_pow = np.zeros(shape=(12,672,720))


####iii)Append values to the variable time (12 months)
#dates = [datetime(2000,01,01)+n*monthdelta(1) for n in range(co_pow.shape[0])]     
dates = [datetime(2000,1,1)+n*monthdelta(1) for n in range(co_pow.shape[0])]     
#time[:] = date2num(dates,units=time.units,calendar=time.calendar)
time[:] = date2num(dates,units=time.units,calendar=time.calendar)
#print "time values (in units %s): " % time.units+"\n",time[:]
print ("time values (in units %s): " % time.units+"\n",time[:])
dates = num2date(time[:],units=time.units,calendar=time.calendar)
print ("dates corresponding to time values:\n",dates)
        
###iv)Append the emissions values in the new NetCDF file  #############
ships [:] = 0

dataset.close()

#Check netCDF version
#print dataset.file_format # At UEA the TNO emissions preprocessor program compiles/runs, only, with the netCDF3.6.3 library. I am not able to produce a netcd3 file directly since is an old library. Therefore, I have applied the following conversion

#############(D)Convert from netCDF4 to netCDF3_Classic  #######
#input file
dsin = ds("Path_directory_to_newly-created_file") 

#output file
dsout = ds("Path_directory_to_place_the_file", "w", format="NETCDF3_CLASSIC")

#Copy dimensions
for dname, the_dim in dsin.dimensions.iteritems():
#    print dname, len(the_dim)
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy variables
for v_name, varin in dsin.variables.iteritems():
    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
  #  print varin.datatype
    
# Copy variable attributes
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
    outVar[:] = varin[:]
# close the output file
dsout.close()

#Check netCDF version
print (dsout.file_format)
#print dsout.file_format




















