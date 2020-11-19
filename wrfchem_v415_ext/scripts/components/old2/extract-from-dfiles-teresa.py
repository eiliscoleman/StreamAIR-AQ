#!/usr/bin/env python3
import bottleneck as bn
import datetime as dt
import glob as g
import json
import os
import glob
import numpy as np
import pandas as pd
import wrf
import xarray as xa
import xarray.ufuncs as xau


VARS = ['lat', 'lon', 'u10', 'v10', 'pb', 'p_sl', 't2', 'so2_concentration', 'o3_concentration', 'nox_concentration', 'pm25', 'pm10', 'rh', 'swdown', 'cldfra', 'bc', 'rain']

#rainmm= {}
#def pressure_rh(d):
    #h_agl_staggered = (d.PHB + d.PH)/wrf.G0
    #h = bn.move_mean(h_agl_staggered.values, 2, axis=0)[1:]
    #t = d['T'].values
    #p = (d.P + d.PB).values
    #q = d.QVAPOR.values
    #return (wrf.slp(h, p, t, q) * 1e-2, wrf.rh(p[0], t[0], q[0]))

def to_dayofyear(d):
    return d.dayofyear + d.hour/24 + d.minute/60/24

coords={"Mace Head": {"lon": -9.9, "lat": 53.33},
             "Shannon": {"lon": -8.914729, "lat": 52.699406},
             "Dublin": {"lon": -6.26, "lat": 53.35},
             "Catania": {"lon": 15.065998, "lat": 37.468496},
             "Galway": {"lon": -9.06, "lat": 53.27},
             "Cork": {"lon": -8.4863, "lat": 51.8969},
             "London": {"lon": -0.1278, "lat": 51.5074},
             "Paris": {"lon": 2.3522, "lat": 48.8566},
             "Carnsore":{"lon": -6.3522, "lat": 52.1706},
             "Malinhead":{"lon": -7.3321, "lat": 55.3528}}
print(coords['Mace Head'])
#def main(config_path):
    #config = {}
    #with open(config_path) as f_config:
        #config = json.load(f_config)

    #fmt_date = 'd01_{}'.format('%Y%m%d%H00.nc')
path_in = ('/mnt/raid/wrf-chem/wrfchem_v415/data_back')
path_out = ('../data/output-csv-archive')
    # for each directory in the data_back folder 
#for ddd in (0,0):
        #print('loop1')
        #next(os.walk(path_in))[1]:
        #print(ddd)

sr_path=(path_in+'/20*')
        # for each directory in this path starting with 20*
for d2 in glob.glob(sr_path):
            print('loop2')
            print(d2)
            #loc_op={'Mace Head':'mh'}
            loc_op={'Malinhead':'mal', 'Carnsore':'cp'}
            # fpaths of all the files 
            f_paths = sorted(g.glob(os.path.join(d2, 'output-wrf/d01_*')))
            datestamp=os.path.basename(os.path.normpath(d2))
            print(datestamp)
            # for every location in this 
            for kk in loc_op.keys():
                print('loop3')
                loc_fname=(datestamp+'_'+kk+'.csv')
                op_fname=os.path.join(path_out,loc_fname)
                out=pd.DataFrame()
                print(op_fname)
                rainmm = {}
                rainaccum = {}
                #rainmmt = {}
                snowmm = {}
                snowmmt = {}
                ddstr = {}
                
                for i, f_path in enumerate(sorted(f_paths)):
                
                    tmp= {}
                    print('di', i,f_path)
                    d = xa.open_dataset(f_path)[VARS]
                    index = [pd.Timestamp(dt.datetime.strptime(os.path.basename(f_path), 'd01_%Y%m%d%H00.nc'
                    ))]
                    if i == 0:
                        index0 = index[0]
                        dayofyear = to_dayofyear(index[0])
                    
                    mh = coords[kk]
                    diffarraylon=np.asarray(d['lon'][:]-mh['lon'])
                    diffarraylat=np.asarray(d['lat'][:]-mh['lat'])
                    diffarrayabs=xau.sqrt(xau.square(diffarraylon)+xau.square(diffarraylat))
                    #note: the array is stored in array dims [lat, lon]
                    ixlat,ixlon=np.where(diffarrayabs == np.min(diffarrayabs))
                    idx = {
                        'x': ixlon,
                        'y': ixlat
                        }
                    d = d.isel(x=idx['x'], y=idx['y'])
                    
                    tmp['dayofyear'] = dayofyear
                    for vv in VARS:
                         if d[vv].ndim==3:
                                tmp[vv]=d[vv][0,0,0].values
                         else:
                                tmp[vv]=d[vv][0,0].values
                        #tmp[vv] = 0.0
                        # first: put in the rain loop
                         #if vv=='rain':
                            #print('look out! its the', vv, 'vbl')
                            #cmmd=('stat -c %y ' + str(f_path))
                            ##print(cmmd)
                            #dstr=str(os.system(cmmd))
                            #print('dstr', dstr)
                            
                            #ddstr[i]=os.popen(cmmd).read()[:10]
                            #ddstr[i]=os.popen(cmmd).read()[50:-6]
                            #print('dstr[i]', dstr[i])
                    if i==0:
                        
                        #rainaccum[i]=np.empty(d['rain'].shape)
                        rainaccum[i]=d['rain'][0,0].values
                        print(rainaccum[i])
                        rainmm[i]=np.nan
                        #tmp['rain']=np.nan   
                    else:
                        rainaccum[i]=d['rain'][0,0].values
                        #rainmmt[i]=rainmm[i-1]
                        rainmm[i]=rainaccum[i]-rainaccum[i-1]
                        print('accum index', i-1)
                        #rainmm[i]=rainaccum[i-1]
                        print('inex', i, 'rain today', rainaccum[i], 'yday', rainaccum[i-1], 'new', rainmm[i])
                    tmp['rain']=rainmm[i]
                                 #if ddstr[i-1]==ddstr[i]:
                                     #print('same dataset')
                                     #print(d[vv].shape)
                                     
                                     #print('rain', tmp[vv])
                                 #else:
                                    #print('new dataset alert!!!!!! must disregard current rain rate value')
                                    #rainmm[i]=d['rain'][0,0].values
                                    #tmp['rain']=np.nan
                                    #print('fuuuuuuuuuuuuuuucl', tmp[vv])
                                    ##i=i+1
                                #print('shape', d['rain'].shape, 'i', i)
                                #print(i-1, rainmm[0])
                                #rainmm[i]=d['rain'][0,0].values - rainmm[i-1]
                                
                                #-rainmm[i-1]
                                #tmp['rain']=rainmm[i]
                                #else:
                                    #print('new dataset alert!!!!!! must disregard current rain rate value')
                                    #rainmm[i]=d['rain'][0,0].values
                                    #tmp[vv]=np.nan
           
                        #else:
                    #tmp['bc']  =d['bc'][0,0,0].values
                    tmp['winddirection_deg']  = np.asscalar(270 - xau.rad2deg(xau.arctan2(d.v10, d.u10)).values) % 360
                    tmp['windspeed_mPs'] = np.asscalar(xau.sqrt(d.u10**2 + d.v10**2).values)
                    #print(tmp.shape)
                    tmp = pd.DataFrame(tmp, index=index).sort_index(axis=1)
                    #tmp.drop_duplicates(subset ="dayofyear", keep = 'last', inplace = True)
                    out = pd.concat((out, tmp))
                #print(out)
                #out.drop_duplicates(subset ="dayofyear", keep = 'last', inplace = True)
                out.to_csv(op_fname)
                print('done ', op_fname)

#if __name__ == '__main__':
    #main('./config.json')

