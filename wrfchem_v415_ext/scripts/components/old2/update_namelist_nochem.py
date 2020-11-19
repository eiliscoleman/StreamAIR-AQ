#!/usr/bin/env python3
import datetime as dt
from subprocess import Popen, PIPE
from datetime import timedelta
from shutil import copyfile
import glob as g
import json
#import numpy as np
#from netCDF4 import Dataset
#import pandas as pd
#import xarray as xr
#import xarray.ufuncs as xu
import os
import re
import fileinput
#from cdo import *
#cdo=Cdo()

pathbase=os.getcwd()

#work_root_dir='/mnt/raid/wrf-chem/wrfchem_v39_cams_cri'
work_root_dir=os.path.dirname(pathbase)
scripts_dir=os.path.join(work_root_dir,'scripts/components')
wrfrundir=os.path.join(work_root_dir,'WRFV3/run')




#2. Update the namelist with chemical scheme
nl_file=os.path.join(wrfrundir, 'namelist.input')

def replace_line(file_name, line_num, text):
        lines = open(file_name, 'r').readlines()
        lines[line_num] = text
        out = open(file_name, 'w')
        out.writelines(lines)
        out.close()

f=open(nl_file ,'r')
namelist_data=f.readlines()

for ll, line in enumerate(namelist_data):
    if 'chem_opt' in line:
        
        new_nl_chem_opt=' chem_opt                            = 0,\n'
        replace_line(nl_file, ll, new_nl_chem_opt)
        
f.close()
