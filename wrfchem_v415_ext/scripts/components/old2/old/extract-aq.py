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
    for f_path in sorted(g.glob(os.path.join(
        config['output-wrf-raw'], 'wrfout_*'
    ))):
        print(f_path)
        f_basename = os.path.basename(f_path)
        domain = f_basename.split('_')[1]
        d = xr.open_dataset(f_path).isel(Time=0)
        time = dt.datetime.strptime(
            f_basename,
            'wrfout_{}_%Y-%m-%d_%H:%M:%S'.format(domain)
        )
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
#        out.coords['h_agl'] = np.mean(h_agl, axis=0)
        out.coords['lat'] = d.XLAT.values[:, 0]
        out.coords['lon'] = d.XLONG.values[0, :]
        out['terrain'] = (('lat', 'lon'), d.HGT.values)
        out['u10'] = (('lat', 'lon'), d.U10.values)
        out['v10'] = (('lat', 'lon'), d.V10.values)
        out['rain'] = (('lat', 'lon'), d.RAINC + d.RAINNC)
        out['t2'] = (('lat', 'lon'), d.T2.values)
        out['so2_concentration'] = (
            ('h_agl', 'lat', 'lon'), d.so2.values) 
        out['o3_concentration'] = (
            ('h_agl', 'lat', 'lon'), d.o3.values) 
        out['nox_concentration'] = (
            ('h_agl', 'lat', 'lon'), (d.no2 + d.no).values) 
        out['pm25'] = (
            ('h_agl', 'lat', 'lon'), d.PM2_5_DRY.values) 
        out['pm10'] = (
            ('h_agl', 'lat', 'lon'), d.PM10.values) 
#            wrf.x_to_yOm3(d.so2.values, (d.PB + d.P).values,
#                          d['T'].values, mm=64) 
#        )
        out['p_sl'] = (('lat', 'lon'), wrf.slp(h, p, t, q))
        out['rh'] = (('lat', 'lon'), wrf.rh(p, t, q)[0])
        out.to_netcdf(
            os.path.join(config['output-wrf'],
                         '{domain}_{date}.nc'
                         .format(domain=domain,
                                 date=(time.strftime('%Y%m%d%H%M'))))
        )


if __name__ == '__main__':
    main('./config.json')
