
#!/usr/bin/env python3
import datetime as dt
from subprocess import Popen, PIPE
from datetime import timedelta
import shutil
import glob 
import json
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import os
import re
import fileinput
import csv
from cdo import *
import multiprocessing as mp
cdo=Cdo()
# Written by Liz Coleman

# list no of processors available for code
print("Number of processors: ", mp.cpu_count())
#pool = mp.Pool(mp.cpu_count())

# List the species required in the download of the CAMS data
#spcs_list=['no2']
spcs_list=['aermr01','aermr02', 'aermr03', 'aermr04', 'aermr05', 'aermr06', 'aermr07', 'aermr08', 'aermr09', 'aermr10', 'aermr11', 'aermr12', 'aermr13', 'aermr14', 'aermr15', 'aermr16', 'aermr17', 'aermr18', 'c2h6', 'c3h8', 'c5h8', 'ch4_c', 'co', 'den', 'go3', 'h2o2', 'hcho', 'hno3', 'lnsp', 'no2', 'no', 'oh', 'pan', 'so2']
#gas_list=spcs_list
gas_list=['c2h6', 'c3h8', 'c5h8', 'ch4_c', 'co','go3', 'h2o2', 'hcho', 'hno3', 'no2', 'no','oh', 'pan', 'so2']
gas_mw=[30., 44.1, 68.1, 28. , 16,48., 34., 30., 63., 46., 30., 17., 121., 64.1]
aer_list=['aermr01','aermr02', 'aermr03', 'aermr04', 'aermr05', 'aermr06', 'aermr07', 'aermr08', 'aermr09', 'aermr10', 'aermr11', 'aermr16', 'aermr17', 'aermr18']
air_mw=28.97
domain_red=('-46,46,30,70')

i=dt.datetime.now()
datetag=(str(i.year)+str(i.month).zfill(2)+str(i.day).zfill(2)+'00')
# This overrides the automated function and specifices dates of the CAMS data to be 24th September, 2019
#datetag=('2019092400')

cur_dir=os.getcwd()
dl_dir=('op_data/dissemination.ecmwf.int/DATA/CAMS_NREALTIME/'+str(datetag))
box_dir=('op_data/dissemination.ecmwf.int/DATA/CAMS_NREALTIME/'+str(datetag)+'/box')
dl_data_path=os.path.join(cur_dir, dl_dir)
box_dir_path=os.path.join(cur_dir, box_dir)


#dl_path=
#call_path=('/DATA/CAMS_NREALTIME/'+datetag+'/z_cams_c_ecmf_'+datetag+'0000_prod_fc_ml_*_'+'species here'+'.nc')
#wget -m ftp://liz.coleman:fj0yVQEV@dissemination.ecmwf.int:/DATA/CAMS_GLOBAL/2019052200/z_cams_c_ecmf_20190522000000_prod_fc_ml_*_no2.nc
#f(x): function to download the data
#-------------------------------------------------------
#FUNCTION 1: f(x): DOWNLOAD ALL THE FILES
#-------------------------------------------------------
def f(x):
    print(x)
    os.system('wget -m ftp://liz.coleman:fj0yVQEV@dissemination.ecmwf.int:/DATA/CAMS_NREALTIME/'+datetag+'/z_cams_c_ecmf_'+datetag+'0000_prod_fc_ml_*_'+str(x)+'.nc -P op_data/')

    
#The command under executes the function to grab the cams data, sorted by species parallel
#if __name__ == '__main__':
    #p = mp.Pool(mp.cpu_count())
    #p.map(f,spcs_list)
    ##print(call_cmd.split(",")[0])
#----------------------------------------------    
#-------------------------------------------------------
#FUNCTION 2:  ff(x) =  take a sample of the CAMS data rather than the whole global domain. Extend of the reduced domain is define in the domain_red parameter.
#-------------------------------------------------------
if not os.path.exists(box_dir_path):
    print('we must make the box directory')
    os.makedirs(box_dir_path)
 
def ff(x):
    #make the box directory if it doesn't exist
    #if not os.path.exists(box_dir_path):
        #print('we must make the box directory')
        #os.makedirs(box_dir_path)
    for rawf in sorted(glob.glob(dl_data_path+'/z_cams_c_ecmf_'+datetag+'0000_prod_fc_ml_*_'+str(x)+'.nc')):
        #print('rawfile', rawf)
        box_file=(rawf.split('.nc')[0]+'_box.nc')
        box_fname=os.path.split(box_file)[1]
        box_dest=os.path.join(box_dir_path, box_fname)
        print(box_dest)
        if not os.path.exists(box_dest):
            print('we must make the boxed file')
            cdo.sellonlatbox(domain_red, input=rawf, output=box_dest)
        else:
            print('boxed file already exists')
            print(box_dest)
#The command under reduces map domain and moves to a specific direcotry 
#if __name__ == '__main__':
    #p = mp.Pool(mp.cpu_count())
    #p.map(ff,spcs_list)
#-------------------------------------------------
#-------------------------------------------------------
#FUNCTION 3: g(fx)
# using "no2" files to sample the dataset, obtain the filenames for each time stamp, merge all species into one file per timestep
#-------------------------------------------------------


f_list=sorted(glob.glob(box_dir_path+'/*no2*nc'))
print(f_list)
ff1={}
for ii, ff in enumerate(f_list):
    ff1[ii]=ff.split('no2')

#ff1[:][0] is an array of all the timestep roots. ncat this array, combining all species in one timestep taskfarm by spreading out each timestep to a different processor. 

ffstr={}
merged_dir_path=os.path.join(dl_dir, 'merged')
if not os.path.exists(merged_dir_path):
    print('we must make the merged directory')
    os.makedirs(merged_dir_path)


for fff in ff1:
    ffstr[fff]=str(ff1[fff][0])
    print('these are the files to merge', ffstr[fff])


#print(ffstr.keys())
# note the merged files are merged from the reduced domain cams
def g(fx):
    
    #fx_key=fx(10:20)
    #os.system('ncrcat' +fxy+'*'
    os.system('ls  ' +str(fx)+'*')
    num=(fx.split('ml_')[1]).split('_')[0]
    #print('number', (fx.split('ml_')[1]).split('_')[0])
    mergcm=('cdo merge '+str(fx)+'*  '+ dl_dir+'/merged/cams_'+str(num)+'.nc')
    print('thhhhhhhhhhiis', mergcm)
    os.system(mergcm)
    print('done', str(fx))

# THe section to merge the files is below: comment out to make script run faster    
#if __name__ == '__main__':
    #p = mp.Pool(mp.cpu_count())
    #p.map(g,ffstr.values())
    

#-------------------------------------------------------
#FUNCTION 4: q(fs)
# Create new MOZART-type file from the existing merged CAMS file: one for each timestep
#-------------------------------------------------------  


# CAMS datafile with raw data (all be merged into one file) before this process)
cams_dir=os.path.join(dl_dir, 'merged')
#for fs in (sorted(os.listdir(os.path.join(dl_dir, 'merged/cams*nc')))):
#os.chdir(cams_dir)    

hyams=pd.read_csv(os.path.join(cur_dir,'hyam'), header=None)
hybms=pd.read_csv(os.path.join(cur_dir,'hybm'), header=None)

#hyams=pd.read_csv(os.path.join(cur_dir,'a'), header=None)
#hybms=pd.read_csv(os.path.join(cur_dir,'b'), header=None)
cams_fpool=glob.glob(os.path.join(cams_dir, "*.nc"))



def q(fs):
    #print(fs)
    o_data=Dataset(fs, 'r')
    #print(o_data.variables)
    #Below: example of how ot print out from variables in the script
    #print(o_data['time'][:])
    
    #o_date=model forecast date in YYYYMMDD
    o_date=datetag[0:8]
    #o_datesec=seconds of forecast since simulation time (this should be retrievable from the data, instead of reading from filename)
    cams_hr=fs.split('cams_')[1][:3]
    cams_hr_int=int(cams_hr)
    cams_hr_inc=int(cams_hr_int/3)
    cams_inc=str(cams_hr_inc).zfill(4)
    #print(cams_inc)
    o_datesec=int(fs.split('cams_')[1][:3])*3600
    lnsp_data=o_data['lnsp']
    o_ps=np.exp(lnsp_data[0::])
    #New filename  
    n_dir=os.path.join(cur_dir, 'for_mozbc')
    #n_fname=('moz_0'+cams_hr+'.nc')
    n_fname=('moz_'+cams_inc+'.nc')
    n_fpath=os.path.join(n_dir, n_fname)
    #print(n_fpath)
    ds=Dataset(n_fpath, 'w')
    # create the new dimensions
    #lat = ds.createDimension('lat', 451)
    lat = ds.createDimension('lat', 101)
    #lon = ds.createDimension('lon', 900)
    lon = ds.createDimension('lon', 231)
    time = ds.createDimension('time', None)
    lev = ds.createDimension('lev', 137)
    ilev = ds.createDimension('ilev', 138)
    ref = ds.createDimension('ref', 0)
    # Create the new variables
    #lat
    lat=ds.createVariable('lat', np.float32, ('lat'))
    #print(o_data.variables['latitude'].shape)
    lat[:]=o_data.variables['latitude'][:]
    #lon
    lon=ds.createVariable('lon', np.float32, ('lon'))
    #ds.variables['lon'][:]=o_data.variables['longitude'][:]
    print(o_datesec)
    lon[:]=(o_data.variables['longitude'][:]+360)%360
    print(lon[:])
    latrev=np.flip(lat, axis= 0)
    lat[:]=latrev
    # hyam
    hyam=ds.createVariable('hyam', np.float32, ('lev'))
    #hyam[:]=hyams[:]
    #try1
    hyam[:]=hyams[:]/100000.0
    ## hybm
    hybm=ds.createVariable('hybm', np.float32, ('lev'))
    #hybm=ds.createVariable('hybm', np.float32, ('ilev'))
    hybm[:]=hybms[:]
    # date
    date=ds.createVariable('date', np.int32, ('time'))
    date[:]=o_date[:]
    # datesec
    datesec=ds.createVariable('datesec', np.int32, ('time'))
    datesec[:]=o_datesec
    # PS
    PS=ds.createVariable('PS', np.float32, ('time', 'lat', 'lon'))
    PS[:]=o_ps
    PS.units='Pa'
    # P0 Ref pressure
    P0=ds.createVariable('P0', np.float32, ('ref'))
    P0[:]=1
    #P0[:]=100000
    P0.units='Pa'
    P0.long_name='reference pressure'
    for ss in gas_list:
        print('index', gas_list.index(ss))
        ss_mw=gas_mw[gas_list.index(ss)]
        fact=(air_mw/ss_mw)
        print('include this species:', ss)
        print('fac', fact)
        sfact=str(fact)
        exec("%s=ds.createVariable(ss, np.float32, ('time', 'lev', 'lat', 'lon'))" % (ss))
        exec("%s[:]=o_data.variables[ss][:,:,:,:]*%s" % (ss, sfact))
        #no2[:]=o_data.variables['no2'][:,:,:,:]*0.63
        exec("%s.units='mol/mol'" % (ss))
        exec("%s.mixing_ratio='dry'" % (ss))
        exec("%srev=np.flip(%s, axis= 2)" % (ss, ss))
        exec("%s[:]=%srev" % (ss, ss))  
    for ass in aer_list:
        print('include this species:', ass)
        exec("%s=ds.createVariable(ass, np.float32, ('time', 'lev', 'lat', 'lon'))" % (ass))
        exec("%s[:]=o_data.variables[ass][:,:,:,:]" % (ass))
        exec("%s.units='kg/kg'" % (ass))
        exec("%s.mixing_ratio='dry'" % (ass))
        exec("%srev=np.flip(%s, axis= 2)" % (ass, ass))
        exec("%s[:]=%srev" % (ss, ss))
    ds.close()
    
    # create the mozart-type files
if __name__ == '__main__':
    p = mp.Pool(mp.cpu_count())
    p.map(q,cams_fpool)

#shutil.rmtree(box_dir_path)
