#!/usr/bin/env python3
import bottleneck as bn
import datetime as dt
import glob as g
import json
import os
import numpy as np
import pandas as pd
import wrf
import xarray as xa
import xarray.ufuncs as xau


VARS = ['lat', 'lon', 'u10', 'v10', 'pb', 'p_sl', 't2','rain', 'snow', 'so2_concentration', 'o3_concentration', 'nox_concentration', 'pm25', 'pm10', 'rh', 'swdown', 'cldfra', 'org25']

#def pressure_rh(d):
    #h_agl_staggered = (d.PHB + d.PH)/wrf.G0
    #h = bn.move_mean(h_agl_staggered.values, 2, axis=0)[1:]
    #t = d['T'].values
    #p = (d.P + d.PB).values
    #q = d.QVAPOR.values
    #return (wrf.slp(h, p, t, q) * 1e-2, wrf.rh(p[0], t[0], q[0]))

def to_dayofyear(d):
    return d.dayofyear + d.hour/24 + d.minute/60/24



def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    #fmt_date = 'd01_{}'.format('%Y%m%d%H00.nc')
    path_in = ('/mnt/raid/wrf-chem/wrfchem_v39_cams_cri/data/output-wrf')
    
    #loc_op={'Mace Head':'mh', 'Dublin':'db','Galway':'gw', 'Cork':'ck', 'London':'ln', 'Paris':'pr', 'Carnsore':'cp', 'Malinhead':'mal'}    
    loc_op={'Dublin':'db'}
    out = pd.DataFrame()
    f_paths = sorted(g.glob(os.path.join(path_in, 'd01_*')))
    
    
    rainmm = {}
    rainmmt = {}
    snowmm = {}
    snowmmt = {}
    ddstr = {}
    
    #print(loc_op.keys())
    for kk in loc_op.keys():
        print('kk', kk)
        print('content', loc_op[kk])
        opdir=('output-'+loc_op[kk])
        path_out = config[opdir]
        mh = config['coords'][kk]
        
    for i, f_path in enumerate(f_paths):
        index = [pd.Timestamp(dt.datetime.strptime(os.path.basename(f_path), 'd01_%Y%m%d%H00.nc'))]
        print('innnndex', index)
        print(i, f_path)
        if i == 0:
            index0 = index[0]
            print('1st inx', index0)
        dayofyear = to_dayofyear(index[0])
        
        #d = xa.open_dataset(f_path).sel(time=0)[VARS]
        d = xa.open_dataset(f_path)[VARS]
        #print(loc_op.keys())
        for kk in loc_op.keys():
            tmp = {}
            print('kk', kk)
            print('content', loc_op[kk])
            opdir=('output-'+loc_op[kk])
            path_out = config[opdir]
            mh = config['coords'][kk]
    
            # Code to pick out the defined gridcell
            diffarraylon=np.asarray(d['lon'][:]-mh['lon'])
            diffarraylat=np.asarray(d['lat'][:]-mh['lat'])
            diffarrayabs=xau.sqrt(xau.square(diffarraylon)+xau.square(diffarraylat))
            #note: the array is stored in array dims [lat, lon]
            ixlat,ixlon=np.where(diffarrayabs == np.min(diffarrayabs))
            idx = {
                'x': ixlon,
                'y': ixlat
            }
                
            #tmp_p, tmp_rh = pressure_rh(d)
            d = d.isel(x=idx['x'], y=idx['y'])
            # windspeed
            tmp['dayofyear'] = dayofyear
            for vv in VARS:
           
                # first: put in the rain loop
                if vv=='rain':
                
                    print('look out! its the', vv, 'vbl')
                    cmmd=('stat -c %y ' + str(f_path))
                    dstr=str(os.system(cmmd))
                    ddstr[i]=os.popen(cmmd).read()[:10]
                    if i==0:
                        rainmm[i]=d['rain'][0,0].values
                        tmp[vv]=np.nan                    
                    else:
                        if ddstr[i-1]==ddstr[i]:
                            rainmm[i]=d['rain'][0,0]-rainmm[i-1]
                            tmp[vv]=rainmm[i].values
                        else:
                            print('new dataset alert!!!!!! must disregard current rain rate value')
                            rainmm[i]=d['rain'][0,0].values
                            tmp[vv]=np.nan
                            
                
                    #print('uuup', ddstr)
                
                    ## if i=0: then disregard the first value, and subtract all subsequent. if new date (find out by executing stat -c '%y' d01_201909110000.nc), then deal with it the same way etc. 
           
           
                else:
                    if d[vv].ndim==3:
                        #print(d[vv].shape)
                        tmp[vv]=d[vv][0,0,0].values
                    else:
                        #print(vv, d[vv].shape)
                        tmp[vv]=d[vv][0,0].values
            tmp['winddirection_deg']  = np.asscalar(
            270 - xau.rad2deg(xau.arctan2(d.v10, d.u10)).values) % 360
            tmp['windspeed_mPs'] = np.asscalar(xau.sqrt(d.u10**2 + d.v10**2).values)
            tmp = pd.DataFrame(tmp, index=index).sort_index(axis=1)
            
            #print(tmp.keys())
            #for kk in tmp.keys():
                #print(kk, tmp[kk].values)
            out = pd.concat((out, tmp))
            out.drop_duplicates(subset ="dayofyear", 
                                keep = 'first', inplace = True)
            #print('done {}'.format(f_path))
     #'BC1', 'BC2', 'OC1', 'OC2']
    #        out = out.rename({'pressure_sea_hPa': 'pressure_hPa'})
#        (out[['dayofyear', 'pressure_hPa', 'relativehumidity_percent',
#          'temperature2m_C', 'winddirection_deg', 'windspeed_mPs']]
#        .to_dataframe()
#        .to_csv(os.path.join(path_out,
#                          '{}.csv'.format(index0.strftime('%Y%m%d%H%M')))))
            print('are you printing')
            print(os.path.join(path_out))
            out.to_csv(os.path.join(path_out,
            '{}.csv'.format(index0.strftime('%Y%m%d%H%M'))))

if __name__ == '__main__':
    main('./config.json')