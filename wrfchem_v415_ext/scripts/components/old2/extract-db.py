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


VARS = ['XLAT', 'XLONG', 'U10', 'V10', 'PB', 'P', 'T', 'T2', 'PHB',
        'PH', 'QVAPOR', 'Q2', 'PSFC', 'QCLOUD', 'QRAIN', 'QICE', 
        'QSNOW', 'QGRAUP', 'CLDFRA', 'RAINNC', 'RAINC', 'SWDOWN', 
        'SNOWNC', 'so2', 'o3', 'no', 'no2', 'PM2_5_DRY', 'PM10', 
        'BC1', 'BC2', 'OC1', 'OC2', 'PBLH']



def pressure_rh(d):
    h_agl_staggered = (d.PHB + d.PH)/wrf.G0
    h = bn.move_mean(h_agl_staggered.values, 2, axis=0)[1:]
    t = d['T'].values
    p = (d.P + d.PB).values
    q = d.QVAPOR.values
    return (wrf.slp(h, p, t, q) * 1e-2, wrf.rh(p[0], t[0], q[0]))


def to_dayofyear(d):
    return d.dayofyear + d.hour/24 + d.minute/60/24


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    fmt_date = 'wrfout_d01_{}'.format(config['fmt']['date'])
    path_in = config['output-wrf-raw']
    path_out = config['output-db']
    mh = config['coords']['Dublin']

    out = pd.DataFrame()
    f_paths = sorted(g.glob(os.path.join(path_in, 'wrfout_d01*')))
    rainmm = {}
    rainmmt = {}
    snowmm = {}
    snowmmt = {}
    
    for i, f_path in enumerate(f_paths):
        index = [pd.Timestamp(dt.datetime.strptime(
            os.path.basename(f_path), fmt_date
        ))]
        if i == 0:
            index0 = index[0]
        dayofyear = to_dayofyear(index[0])
        tmp = {}
        d = xa.open_dataset(f_path).sel(Time=0)[VARS]
        # Code to pick out the defined gridcell
        diffarraylon=np.asarray(d.XLONG[:]-mh['lon'])
        diffarraylat=np.asarray(d.XLAT[:]-mh['lat'])
        diffarrayabs=xau.sqrt(xau.square(diffarraylon)+xau.square(diffarraylat))
        #note: the array is stored in array dims [lat, lon]
        ixlat,ixlon=np.where(diffarrayabs == np.min(diffarrayabs))
        idx = {
            'x': ixlon,
            'y': ixlat
        }
        
        print(idx['x'], idx['y'])
        tmp_p, tmp_rh = pressure_rh(d)
        d = d.isel(west_east=idx['x'], south_north=idx['y'])
        # windspeed
        tmp['dayofyear'] = dayofyear
#        tmp['winddirection_deg'] = np.asscalar(xau.rad2deg(
#            xau.arctan2(-d.U10, -d.V10)
#        ).values)
        tmp['pressure_hPa'] = tmp_p[idx['y'], idx['x']]
        tmp['qcloud'] = d.isel(bottom_top=0).QCLOUD[0,:].values
        tmp['qgraup'] = d.isel(bottom_top=0).QGRAUP[0,:].values
        tmp['qicd'] = d.isel(bottom_top=0).QICE[0,:].values
        tmp['qrain'] = d.isel(bottom_top=0).QRAIN[0,:].values
        tmp['qsnow'] = d.isel(bottom_top=0).QSNOW[0,:].values
        if i == 0:
          rainmm[i] = np.asscalar((d.RAINNC.values + d.RAINC.values))
          rainmmt[i] = np.asscalar((d.RAINNC.values + d.RAINC.values))
          snowmm[i] = np.asscalar(d.SNOWNC.values)
          snowmmt[i] = np.asscalar(d.SNOWNC.values)
        else:
          rainmmt[i] = np.asscalar((d.RAINNC.values + d.RAINC.values))
          rainmm[i] = np.asscalar(d.RAINNC.values + d.RAINC.values)-rainmmt[i-1]
          snowmmt[i] = np.asscalar(d.SNOWNC.values)
          snowmm[i] = np.asscalar(d.SNOWNC.values)-snowmmt[i-1]
        tmp['rain_mm'] = rainmm[i]  
        tmp['relativehumidity_percent'] = tmp_rh[idx['y'], idx['x']]
        tmp['temperature2m_C'] = np.asscalar((d.T2 - 273.15).values)
        tmp['winddirection_deg'] = np.asscalar(
            270 - xau.rad2deg(xau.arctan2(d.V10, d.U10)).values) % 360
        tmp['windspeed_mPs'] = np.asscalar(xau.sqrt(d.U10**2 + d.V10**2).values)
        tmp['zcldfra'] = d.isel(bottom_top=0).CLDFRA[0,:].values
        tmp['znox'] = d.isel(bottom_top=0).no2[0,:].values + d.isel(bottom_top=0).no[0,:].values
        tmp['zo3'] = d.isel(bottom_top=0).o3[0,:].values*1000.0
        tmp['zpm25'] = d.isel(bottom_top=0).PM2_5_DRY[0,:].values
        tmp['zpm10'] = d.isel(bottom_top=0).PM10[0,:].values
        tmp['zbc1'] = d.isel(bottom_top=0).BC1[0,:].values
        tmp['zbc2'] = d.isel(bottom_top=0).BC2[0,:].values
        tmp['zoc1'] = d.isel(bottom_top=0).OC1[0,:].values
        tmp['zoc2'] = d.isel(bottom_top=0).OC2[0,:].values
        tmp['zso2'] = d.isel(bottom_top=0).so2[0,:].values
        tmp['zpblh'] = np.asscalar(d.PBLH.values)
        tmp['zuv_index'] = np.sum((d.o3*d.PB/6950.0).values)
        tmp = pd.DataFrame(tmp, index=index).sort_index(axis=1)
        out = pd.concat((out, tmp))
        print('done {}'.format(f_path))
 #'BC1', 'BC2', 'OC1', 'OC2']
#        out = out.rename({'pressure_sea_hPa': 'pressure_hPa'})
#        (out[['dayofyear', 'pressure_hPa', 'relativehumidity_percent',
#          'temperature2m_C', 'winddirection_deg', 'windspeed_mPs']]
#        .to_dataframe()
#        .to_csv(os.path.join(path_out,
#                          '{}.csv'.format(index0.strftime('%Y%m%d%H%M')))))

    out.to_csv(os.path.join(
        path_out,
        '{}.csv'.format(index0.strftime('%Y%m%d%H%M'))
    ))

if __name__ == '__main__':
    main('./config.json')
