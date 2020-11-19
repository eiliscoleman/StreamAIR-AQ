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
        'PH', 'QVAPOR', 'Q2', 'PSFC']


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
    path_out = config['output-mh']
    mh = config['coords']['Mace Head']

    out = pd.DataFrame()
    f_paths = sorted(g.glob(os.path.join(path_in, 'wrfout_d01*')))
    for i, f_path in enumerate(f_paths):
        index = [pd.Timestamp(dt.datetime.strptime(
            os.path.basename(f_path), fmt_date
        ))]
        if i == 0:
            index0 = index[0]
        dayofyear = to_dayofyear(index[0])
        tmp = {}
        d = xa.open_dataset(f_path).sel(Time=0)[VARS]
        idx = {
            'x': np.asscalar(abs(d.XLONG[0, :] - mh['lon']).argmin()),
            'y': np.asscalar(abs(d.XLAT[:, 0] - mh['lat']).argmin())
        }
        tmp_p, tmp_rh = pressure_rh(d)
        d = d.isel(west_east=idx['x'], south_north=idx['y'])

        # windspeed
        tmp['dayofyear'] = dayofyear
        tmp['windspeed_mPs'] = np.asscalar(xau.sqrt(d.U10**2 + d.V10**2).values)
        tmp['winddirection_deg'] = np.asscalar(xau.rad2deg(
            xau.arctan2(-d.U10, -d.V10)
        ).values)
        tmp['pressure_hPa'] = tmp_p[idx['y'], idx['x']]
        tmp['temperature2m_C'] = np.asscalar((d.T2 - 273.15).values)
        tmp['relativehumidity_percent'] = tmp_rh[idx['y'], idx['x']]
        tmp = pd.DataFrame(tmp, index=index).sort_index(axis=1)
        out = pd.concat((out, tmp))
        print('done {}'.format(f_path))

    out.to_csv(os.path.join(
        path_out,
        '{}.csv'.format(index0.strftime('%Y%m%d%H%M'))
    ))


if __name__ == '__main__':
    main('./config.json')
