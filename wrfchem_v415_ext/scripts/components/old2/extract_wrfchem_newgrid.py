#!/usr/bin/env python3
import bottleneck as bn
import datetime as dt
import glob as g
import json
import numpy as np
import os
import wrf
import xarray as xr


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

#    out, times, h_agl = xr.Dataset(), [], []
    #for f_path in sorted(g.glob(os.path.join(
        #config['output-wrf-raw'], 'wrfout_*'
    #))):
    #dir_data=('/mnt/raid/rong-ming/wrfchem/data_in')
    #dir_out=('/mnt/raid/rong-ming/wrfchem/data_in_dfiles')
    #no_files=len(os.listdir(dir_data))
    for i, f_path in enumerate(sorted(g.glob(os.path.join(
        config['output-wrf-raw'], 'wrfout_*')
    ))):
        print(f_path)
        f_basename = os.path.basename(f_path)
        print(f_basename)
        domain = f_basename.split('_')[1]
        #fpath=os.path.join(dir_data, f_path)
        d = xr.open_dataset(f_path).isel(Time=0)
        #d = xr.open_dataset(fpath)
        time = dt.datetime.strptime(
            f_basename,
            'wrfout_{}_%Y-%m-%d_%H:%M:%S'.format(domain)
        )
        print(time)
        h_agl_staggered = (d.PHB + d.PH)/wrf.G0
        h_agl = bn.move_mean(h_agl_staggered
                              .mean(dim='south_north')
                              .mean(dim='west_east'), 2)[1:]
        h = bn.move_mean(h_agl_staggered.values, 2, axis=0)[1:]
        t = d['T'].values
        p = (d.P + d.PB).values
        q = d.QVAPOR.values

        out = xr.Dataset()
        out.coords['time'] = time
##        out.coords['h_agl'] = np.mean(h_agl, axis=0)
        out.coords['x'] = range(d.XLONG.shape[1])
        out.coords['y'] = range(d.XLAT.shape[0])
        out['lat'] = (('y', 'x'), d.XLAT.values[:])
        out['lon'] = (('y', 'x'), d.XLONG.values[:])
        out['terrain'] = (('y', 'x'), d.HGT.values)
        out['u10'] = (('y', 'x'), d.U10.values)
        out['v10'] = (('y', 'x'), d.V10.values)
        out['rain'] = (('y', 'x'), d.RAINC + d.RAINNC)
        out['t2'] = (('y', 'x'), d.T2.values)
        out['so2_concentration'] = (
            ('h_agl', 'y', 'x'), d.so2.values) 
        out['o3_concentration'] = (
            ('h_agl', 'y', 'x'), d.o3.values) 
        out['nox_concentration'] = (
            ('h_agl', 'y', 'x'), (d.no2 + d.no).values) 
        out['pm25'] = (
            ('h_agl', 'y', 'x'), d.PM2_5_DRY.values) 
        out['pm10'] = (
            ('h_agl', 'y', 'x'), d.PM10.values) 
##            wrf.x_to_yOm3(d.so2.values, (d.PB + d.P).values,
##                          d['T'].values, mm=64) 
##        )
        out['p_sl'] = (('y', 'x'), wrf.slp(h, p, t, q))
        out['rh'] = (('y', 'x'), wrf.rh(p, t, q)[0])
        out.to_netcdf(
            os.path.join(config['output-wrf'],
                         '{domain}_{date}.nc'
                         .format(domain=domain,
                                 date=(time.strftime('%Y%m%d%H%M'))))
        )


if __name__ == '__main__':
    main('./config.json')