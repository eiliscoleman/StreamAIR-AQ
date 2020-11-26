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
import pandas as pd
import xray
from monthdelta import monthdelta
from datetime import datetime
from netCDF4 import num2date, date2num
dates = []


############(A) Open original TNO emissions inventory################
#tno_2011 = 'path to the TNO_MACCII'  #input file
tno_2011 = '/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/TNO_MACC/TNO_MACC_III_emissions_2011.nc'
tno = nc.Dataset(tno_2011,'r')

############ Read TNO EC and OC Fractions of PM emissions inventory################

# Assign spreadsheet filename to `file`
file = '/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/TNO_MACC/PM_split_for_TNO_MACC_III_2011.xlsx'

# Load spreadsheet
xl = pd.ExcelFile(file)

# Print the sheet names
print(xl.sheet_names)

# Load a sheet into a DataFrame by name: df1
#df1 = xl.parse('Sheet1')
#df1 = xl.parse('coarse')
df1 = xl.parse('fine')
x=[]
x.append(df1['EC_fine'])
print(x)
#######################(B)Create new netcdf file#################################
#1) Create an empty NetCDF file 
#dataset = ds('Path_directory_to_place_the_file','w','r', format='NETCDF4_CLASSIC') #output file
#dataset = ds('/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/OUT/tno_ec_1_25.nc','w','r', format='NETCDF4_CLASSIC') #output file
#dataset = ds('/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/OUT/tno_oc_1_25.nc','w','r', format='NETCDF4_CLASSIC') #output file
#dataset = ds('/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/OUT/tno_oc_25_10.nc','w','r', format='NETCDF4_CLASSIC') #output file
dataset = ds('/mnt/raid/wrf-chem/emis/WRF_EMIS_UoM_UEA/preprocessor_script/OUT/tno_ec_25_10.nc','w','r', format='NETCDF4_CLASSIC') #output file

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

pow = dataset.createVariable('pow',np.float64, ('time','lat','lon'), fill_value=False)
res = dataset.createVariable('res',np.float32, ('time','lat','lon'), fill_value=False)
inc = dataset.createVariable('inc',np.float32, ('time','lat','lon'), fill_value=False)
pei = dataset.createVariable('pei',np.float32, ('time','lat','lon'), fill_value=False)
exf = dataset.createVariable('exf',np.float32, ('time','lat','lon'), fill_value=False)
sol = dataset.createVariable('sol',np.float32, ('time','lat','lon'), fill_value=False)
tra1 = dataset.createVariable('tra1',np.float32, ('time','lat','lon'), fill_value=False)
tra2 = dataset.createVariable('tra2',np.float32, ('time','lat','lon'), fill_value=False)
tra3 = dataset.createVariable('tra3',np.float32, ('time','lat','lon'), fill_value=False)
tra4 = dataset.createVariable('tra4',np.float32, ('time','lat','lon'), fill_value=False)
tra5 = dataset.createVariable('tra5',np.float32, ('time','lat','lon'), fill_value=False)
nrt = dataset.createVariable('nrt',np.float32, ('time','lat','lon'), fill_value=False)
was = dataset.createVariable('was',np.float32, ('time','lat','lon'), fill_value=False)
agr = dataset.createVariable('agr',np.float32, ('time','lat','lon'), fill_value=False)


#4) Add attributes units to 1D and 3D variables
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time.units = 'days since 1900-01-01 00:00'
time.calendar ='gregorian'
time.long_name = 'Time'

pow.units = 'Kg yr-1'
pow.long_name = 'Power generation'
res.units = 'Kg yr-1'
res.long_name= 'Residential, comercial and other combustion'
inc.units = 'Kg yr-1'
inc.long_name = 'Industrial combustion'
pei.units = 'Kg yr-1'
pei.long_name= 'Processed emission industrial'
exf.units = 'Kg yr-1'
exf.long_name = 'Extraction and distribution of fossil fuels'
sol.units = 'Kg yr-1'
sol.long_name = 'Solvent use'
tra1.units = 'Kg yr-1'
tra1.long_name = 'Road transport, gasoline'
tra2.units = 'Kg yr-1'
tra2.long_name = 'Road transport, diesel' 
tra3.units = 'Kg yr-1'
tra3.long_name = 'Road trasnport, LPG'
tra4.units = 'Kg yr-1'
tra4.long_name = 'Road trasnport, non-exhaust, volatilisation'
tra5.units = 'Kg yr-1'
tra5.long_name = 'Road transport, non-exhaust, wear'
nrt.units = 'Kg yr-1'
nrt.long_name = 'Non-road transport'
was.units = 'Kg yr-1'
was.long_name = 'Waste tratment and disposal'
agr.units = 'Kg yr-1'
agr.long_name = 'Agriculture'
 
                          
############# (C) Append data from the original emissions files to the new NetCDF file, example for NOx #######

#1)Extract and add latitude  and longitude values 

data = xray.open_dataset(tno_2011)  #Open NetCDF file (TNO-MACCII 2011) with x-ray library
#print data
print(data)

latitude_tno = tno['latitude'][:]  # Extract latitude from the original file
longitude_tno = tno['longitude'][:] 

lat [:] = latitude_tno  #Append latitude to the new file
lon [:] = longitude_tno   

#2)Extract other variables (NOx, pm10, etc) using the NETCDF4 library########
tno_net = ds(tno_2011) 
emis_cat_index = tno_net.variables['emission_category_index'][:]
lat_index =tno_net.variables['latitude_index'][:]
lon_index =tno_net.variables['longitude_index'][:]
pm10_data = tno_net.variables['pm10'][:]  #Array with the emissions per pm10mpound 

#3)Loop through every emission category (emiss_cat) and pick the emissions values by sector, latitude, and longitude 
##### i) Create a 3D array with zeros 
pm10_pow = np.zeros(shape=(12,672,720))
pm10_res = np.zeros(shape=(12,672,720))
pm10_inc = np.zeros(shape=(12,672,720))
pm10_pei = np.zeros(shape=(12,672,720)) #Since the new inventory have sectors three and four merged as 'sector 34', I have assigned emissions from sector 34 to the category 3 and left category four empty
pm10_exf = np.zeros(shape=(12,672,720))
pm10_sol = np.zeros(shape=(12,672,720))
pm10_tra1 = np.zeros(shape=(12,672,720))
pm10_tra2 = np.zeros(shape=(12,672,720))
pm10_tra3 = np.zeros(shape=(12,672,720))
pm10_tra4 = np.zeros(shape=(12,672,720))
pm10_tra5 = np.zeros(shape=(12,672,720))
pm10_nrt = np.zeros(shape=(12,672,720))
pm10_was = np.zeros(shape=(12,672,720))
pm10_agr = np.zeros(shape=(12,672,720))

####ii)Pick values from each sector and fill the 3D arrays with the new emission values
for i in range(len(lat_index)):
    sector = emis_cat_index[i]
    if sector == 1:
        for k in range(12):
            pm10_pow[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]
    if sector == 2:
        for k in range(12):
            pm10_res[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]
    if sector == 3:
        for k in range(12):
            pm10_inc[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]
    if sector == 4:
        for k in range(12):
            pm10_exf[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 5:
        for k in range(12):
            pm10_sol[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 6:
        for k in range(12):
            pm10_tra1[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]        
    if sector == 7:
       for k in range(12):
           pm10_tra2[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]  
    if sector == 8:
        for k in range(12):
            pm10_tra3[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 9:
        for k in range(12):
            pm10_tra4[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 10:
        for k in range(12):
            pm10_tra5[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 11:
        for k in range(12):
            pm10_nrt[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i] 
    if sector == 12:
        for k in range(12):
            pm10_was[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]        
    if sector == 13:
        for k in range(12):
            pm10_agr[k,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]        
                        
print (pm10_nrt[5,105,0]) 
print (pm10_nrt.shape) #To check if I am picking the right values 
#print pm10_nrt[5,105,0] 
#print pm10_nrt.shape #To check if I am picking the right values 

####iii)Append values to the variable time (12 months)
#dates = [datetime(2000,01,01)+n*monthdelta(1) for n in range(pm10_pow.shape[0])]     
dates = [datetime(2000,1,1)+n*monthdelta(1) for n in range(pm10_pow.shape[0])]     
time[:] = date2num(dates,units=time.units,calendar=time.calendar)
#time[:] = date2num(dates,units=time.units,calendar=time.calendar)
#print "time values (in units %s): " % time.units+"\n",time[:]
print ("time values (in units %s): " % time.units+"\n",time[:])
dates = num2date(time[:],units=time.units,calendar=time.calendar)
print ("dates pm10rresponding to time values:\n",dates)
        
###iv)Append the emissions values in the new NetCDF file  #############
#for ec
pow [:] = pm10_pow * 0.10147832
res [:] = pm10_res * 0.11758079
inc [:] = pm10_inc * 0.05799639 * 0.5
pei [:] = pm10_pei * 0.05799639 * 0.5
exf [:] = pm10_exf * 0.45873689
sol [:] = pm10_sol * 0.0
tra1 [:] = pm10_tra1 * 0.29486546
tra2 [:] = pm10_tra2 * 0.68190075
tra3 [:] = pm10_tra3 * 0.3
tra4 [:] = pm10_tra4 * 0.0
tra5 [:] = pm10_tra5 * 0.04866798
nrt [:] = pm10_nrt * 0.44384309
was [:] = pm10_was * 0.2
agr [:] = pm10_agr * 0.00278432
#for oc
#pow [:] = pm10_pow * 0.1597653
#res [:] = pm10_res * 0.00965305
#inc [:] = pm10_inc * 0.08713204 * 0.5
#pei [:] = pm10_pei * 0.08713204 * 0.5
#exf [:] = pm10_exf * 0.0
#sol [:] = pm10_sol * 0.0
#tra1 [:] = pm10_tra1 * 0.0
#tra2 [:] = pm10_tra2 * 0.0
#tra3 [:] = pm10_tra3 * 0.0
#tra4 [:] = pm10_tra4 * 0.0
#tra5 [:] = pm10_tra5 * 0.15966201
#nrt [:] = pm10_nrt * 0.0
#was [:] = pm10_was * 0.0
#agr [:] = pm10_agr * 0.62946225
dataset.close()

#Check netCDF version
#print dataset.file_format # At UEA the TNO emissions preprocessor program pm25mpiles/runs, only, with the netCDF3.6.3 library. I am not able to produce a netcd3 file directly since is an old library. Therefore, I have applied the following pm25nversion

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




















