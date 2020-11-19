#!/usr/bin/env python3
import cartopy as cpy
import cartopy.crs as ccrs
import json
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import seaborn as sns


class MidpointNormalize(mc.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        left = None if clip or self.clip else self.vmin - 1
        right = None if clip or self.clip else self.vmax + 1
        return np.ma.masked_array(
            np.interp(value, x, y, left=left, right=right)
        )


LEVELS = {
    'aml-e': np.linspace(0, 10 + 0.1, 0.5),
    'aml-i': np.linspace(0, 500 + 1, 50),
    'conc': [1e-2, 2, 4],
    't2-ir': np.linspace(-6, 16 + 1, 2),
    't2-e': np.linspace(-20, 38 + 1, 2),
    'rain': [0.1, 0.3, 0.5, 0.7, 1, 2, 3, 6, 12],
    'wind': np.linspace(0, 31, 3),
    'rh': np.linspace(0, 101, 10),
    'so2': np.linspace(0, 8, 0.5),
    'o3': np.linspace(0, 80, 5),
    'nox': np.linspace(0, 70, 5),
    'pm25': np.linspace(0, 30, 2), 
    'pm10': np.linspace(0, 30, 2)
}
PROJ = ccrs.Mercator()
PROJ_T = ccrs.PlateCarree()
CMAP = {
    'a': sns.cubehelix_palette(start=0.7, light=0.95, rot=0.75, as_cmap=True),
    'c': sns.cubehelix_palette(start=0.4, light=0.95, as_cmap=True),
    't2': sns.diverging_palette(240, 10, l=40, n=len(LEVELS['t2-e']),
                                as_cmap=True),
    'rain': sns.cubehelix_palette(start=0, light=0.95, rot=1, as_cmap=True),
    'wind': sns.cubehelix_palette(start=2, light=0.95, rot=1, as_cmap=True),
    'rh': sns.cubehelix_palette(start=0, dark=0.3, light=0.8,
                                rot=0.75, as_cmap=True)
}
CMAP['a'].set_over('firebrick')
CMAP['c'].set_over('indigo')
CMAP['t2'].set_under('indigo')
CMAP['t2'].set_over('firebrick')
CMAP['rain'].set_under('w', alpha=0)


def round_to(x, b=10):
    return b * round(x/b)


def set_map(extent=None):
    ax = plt.axes(projection=PROJ)
    ax.coastlines('50m')
    ax.add_feature(cpy.feature.BORDERS, linewidth=0.25)
    gl = ax.gridlines(draw_labels=True, linestyle='-', color='gray',
                      linewidth=0.25)
    gl.xlabels_top = None
    gl.ylabels_right = None
    if extent:
        ax.set_extent(extent)
    return ax


def plot1d(d, t, savename, title='', config={}):
    def axvline(**kwargs):
        for x, c in zip(kwargs['xs'], kwargs['colors']):
            plt.axvline(x=x, color=c, linestyle='--')

    df = pd.DataFrame()
    for airport, info in config['coords'].items():
        if info['airport'] != True:
            continue

        tmp = d.sel(lat=info['lat'], lon=info['lon'], method='nearest')
        df = df.append(pd.DataFrame(
            {'airport': airport,
             'conc': tmp.ash_concentration.values * 1e-3,
             'h_agl': tmp.h_agl.values * 1e-3}
        ))

    g = sns.FacetGrid(df, col='airport', size=5, aspect=0.7)
    (g.map(plt.fill_betweenx, 'h_agl', 'conc')
     .map(axvline, xs=[2, 4], colors=[sns.color_palette()[4],
                                      sns.color_palette()[2]])
     .set_titles('{col_name} airport')
     .set_xlabels('Ash concentration (mg/m$^3$)')
     .set_ylabels('Height AGL (km)')
     .set(xlim=(0, 10)))
    for ax in g.axes.flat:
        plt.text(
            0.01, 0.01,
            t.tz_localize(tz='Europe/Dublin').strftime('%Y.%m.%d %a %H:%M'),
            bbox=dict(facecolor='0.75', alpha=0.75),
            transform=ax.transAxes
        )
    g.fig.suptitle(title)
    g.fig.subplots_adjust(left=0.085, bottom=0.11, top=0.9)
    plt.savefig(os.path.join(config['imgs'], '{}'.format(savename)))
    plt.close()


def plot2d(lon, lat, d, fig, newfig=True, t=None, levels=None, levels_n=None,
           norm=None, cmap=None, label='', title='', extent=None,
           what='contourf', colorbar=True, extend='neither', format=None,
           config={}):
    if newfig:
        plt.figure(num=fig, figsize=(6, 5))
        ax = set_map(extent)
    else:
        plt.figure(num=fig)
        ax = plt.gca()

    if what == 'quiver':
        l = xu.sqrt(d.u10**2 + d.v10**2)
        u10 = d.u10/l
        v10 = d.v10/l
        plt.quiver(lon.values, lat.values, u10.values, v10.values,
                   scale=40, color='k', alpha=0.35, zorder=10,
                   transform=PROJ_T)
    elif what == 'contour':
        if levels_n:
            cs = plt.contour(lon.values, lat.values, d.values, levels_n,
                             linewidths=0.75, colors='k', transform=PROJ_T)
        else:
            cs = plt.contour(lon.values, lat.values, d.values,
                             linewidths=0.75, colors='k', transform=PROJ_T)
        plt.clabel(cs, inline=1, fmt='%.0f', fontsize=8)
    elif what == 'contourf':
#        plt.contourf(lon.values, lat.values, d.values, 9, 
        plt.contourf(lon.values, lat.values, d.values,  
                     norm=norm, cmap=cmap,
                     extend=extend, transform=PROJ_T)
#        plt.contourf(lon.values, lat.values, d.values,
#                     levels=levels, norm=norm, cmap=cmap,
#                     extend=extend, transform=PROJ_T)
        if colorbar:
            if format:
                cbar = plt.colorbar(shrink=0.9, format=format)
            else:
                cbar = plt.colorbar(shrink=0.9)
            cbar.set_label(label)

    if newfig:
        ax.set_xlabel('lon')
        ax.set_ylabel('lat')
        plt.title(title)
        plt.text(
            0.01, 0.01,
            t.tz_localize(tz='Europe/Dublin').strftime('%Y.%m.%d %a %H:%M'),
            bbox=dict(facecolor='0.75', alpha=0.75),
            transform=ax.transAxes
        )
        plt.tight_layout()
        plt.subplots_adjust(left=0.075, bottom=0.075)


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    doms = sorted(set(map(lambda x: x.split('_')[0],
                         os.listdir(os.path.join(config['output-wrf'])))))
    ds = [xr.open_mfdataset(os.path.join(config['output-wrf'],
                                         '{}*.nc'.format(dom)),
                            concat_dim='time') for dom in doms]
    extents = {
        'ireland': [-12, -3, 51, 55.5],
        'europe': [ds[0].lon.min(), ds[0].lon.max(),
                   ds[0].lat.min(), ds[0].lat.max() - 1]
    }

    for z in zip(*map(lambda x: list(x.groupby('time')), ds)):
        for i, (t, d) in enumerate(z):
            is_fst = i == 0
            is_lst = i == len(doms) - 1

            t = pd.to_datetime(t)
            t_save = t.strftime('%Y%m%d%H%M')

            if is_fst:
                print(t)


            print('\t{} - temperature & pressure (Ireland)'.format(i))
            figs = ['t2-p-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.t2 - 273.15, fig=figs[-1], newfig=is_fst,
                   t=t, levels=LEVELS['t2-ir'],
                   norm=MidpointNormalize(midpoint=0), cmap=CMAP['t2'],
                   extent=extents['ireland'], extend='both', label='$^o$C',
                   title='Temperature and Pressure',
                   colorbar=is_lst, config=config)
            if is_fst:
                plot2d(
                    d.lon.loc[extents['ireland'][0]:extents['ireland'][1]],
                    d.lat.loc[extents['ireland'][2]:extents['ireland'][3]],
                    d.p_sl.sel(
                        lon=slice(extents['ireland'][0],
                                  extents['ireland'][1]),
                        lat=slice(extents['ireland'][2],
                                  extents['ireland'][3])
                    ) * 1e-2, fig=figs[-1], newfig=False,
                    levels_n=10, what='contour', config=config
                )
            print('\t{} - temperature & pressure (Europe)'.format(i))
            figs += ['t2-p-e_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.t2 - 273.15, fig=figs[-1], newfig=is_fst,
                   t=t, levels=LEVELS['t2-e'],
                   norm=MidpointNormalize(midpoint=0), cmap=CMAP['t2'],
                   extent=extents['europe'], extend='both', label='$^o$C',
                   title='Temperature and Pressure',
                   colorbar=is_lst, config=config)
            if is_fst:
                plot2d(d.lon, d.lat, d.p_sl * 1e-2, fig=figs[-1],
                       newfig=False, levels_n=20, what='contour',
                       config=config)

            print('\t{} - rain (Ireland)'.format(i))
            figs += ['rain-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.rain, fig=figs[-1], newfig=is_fst,
                   t=t, levels=LEVELS['rain'], cmap=CMAP['rain'],
                   extent=extents['ireland'], label='mm/h', format='%.1f',
                   title='Precipitation', colorbar=is_lst, config=config)
            print('\t{} - rain (Europe)'.format(i))
            figs += ['rain-e_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.rain, fig=figs[-1], newfig=is_fst,
                   t=t, levels=LEVELS['rain'], cmap=CMAP['rain'],
                   extent=extents['europe'], label='mm/h', format='%.1f',
                   title='Precipitation', colorbar=is_lst, config=config)

            print('\t{} - wind (Ireland)'.format(i))
            step = 2
            figs += ['wind-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, xu.sqrt(d.u10**2 + d.v10**2),
                   fig=figs[-1], newfig=is_fst, t=t, levels=LEVELS['wind'],
                   extent=extents['ireland'], cmap=CMAP['wind'],
                   label='m/s', format='%.0f', title='Wind speed & direction',
                   colorbar=is_lst, config=config)
            if is_fst:
                plot2d(d.lon[::step], d.lat[::step],
                       d.isel(lat=slice(None, None, step),
                              lon=slice(None, None, step)),
                       fig=figs[-1], newfig=False, t=t,
                       what='quiver', config=config)
            print('\t{} - wind (Europe)'.format(i))
            step = 4
            figs += ['wind-e_{}'.format(t_save)]
            plot2d(d.lon, d.lat, xu.sqrt(d.u10**2 + d.v10**2),
                   fig=figs[-1], newfig=is_fst, t=t, levels=LEVELS['wind'],
                   extent=extents['europe'], cmap=CMAP['wind'],
                   label='m/s', format='%.0f', title='Wind speed & direction',
                   colorbar=is_lst, config=config)
            if is_fst:
                plot2d(d.lon[::step], d.lat[::step],
                       d.isel(lat=slice(None, None, step),
                              lon=slice(None, None, step)),
                       fig=figs[-1], newfig=False, t=t,
                       what='quiver', config=config)

            print('\t -PM2.5 (Ireland)')
            figs += ['pm25-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.pm25[0,:,:], fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'],
                   levels=LEVELS['pm25'], cmap=CMAP['a'], 
                   label='PM2.5 (ug/m$^3$)',
                   format='%.1f', title='PM2.5', colorbar=is_lst,
                   config=config)
            print('\t -PM2.5 (Europe)')
            figs += ['pm25-eu_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.pm25[0,:,:], fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'],
                   levels=LEVELS['pm25'], cmap=CMAP['a'], 
                   label='PM2.5 (ug/m$^3$)',
                   format='%.1f', title='PM2.5', colorbar=is_lst,
                   config=config)

            print('\t -PM10 (Ireland)')
            figs += ['pm10-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.pm10[0,:,:], fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'],
                   levels=LEVELS['pm10'], cmap=CMAP['a'], 
                   label='PM10 (ug/m$^3$)',
                   format='%.1f', title='PM10', colorbar=is_lst,
                   config=config)
            print('\t -PM10 (Europe)')
            figs += ['pm10-eu_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.pm10[0,:,:], fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'],
                   levels=LEVELS['pm10'], cmap=CMAP['a'], 
                   label='PM10 (ug/m$^3$)',
                   format='%.1f', title='PM10', colorbar=is_lst,
                   config=config)

            print('\t -SO2 (Ireland)')
            figs += ['so2-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.so2_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'], extend='both',
                   levels=LEVELS['so2'], cmap=CMAP['a'], 
                   label='SO2 (ppbv)',
                   format='%.1f', title='SO2', colorbar=is_lst,
                   config=config)
            print('\t -SO2 (Europe)')
            figs += ['so2-eu_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.so2_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'], extend='both',
                   levels=LEVELS['so2'], cmap=CMAP['a'], 
                   label='SO2 (ppbv)',
                   format='%.1f', title='SO2', colorbar=is_lst,
                   config=config)

            print('\t -O3 (Ireland)')
            figs += ['o3-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.o3_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'],
                   levels=LEVELS['o3'], cmap=CMAP['a'], 
                   label='O3 (ppbv)',
                   format='%.1f', title='O3', colorbar=is_lst,
                   config=config)
            print('\t -O3 (Europe)')
            figs += ['o3-eu_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.o3_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'],
                   levels=LEVELS['o3'], cmap=CMAP['a'], 
                   label='O3 (ppbv)',
                   format='%.1f', title='O3', colorbar=is_lst,
                   config=config)

            print('\t -NOx (Ireland)')
            figs += ['nox-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.nox_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'], extend='both',
                   levels=LEVELS['nox'], cmap=CMAP['a'], 
                   label='NOx (ppbv)',
                   format='%.1f', title='NOx', colorbar=is_lst,
                   config=config)
            print('\t -NOx (Europe)')
            figs += ['nox-eu_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.nox_concentration[0,:,:] * 1e3, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'], extend='both',
                   levels=LEVELS['nox'], cmap=CMAP['a'], 
                   label='NOx (ppbv)',
                   format='%.1f', title='NOx', colorbar=is_lst,
                   config=config)

#            plot2d(
#            d.lon, d.lat, d.rh, 
#            fig=figs[-1], newfig=is_fst, t=t,
#            extend=extents['europe'],
#            levels=LEVELS['pm25'], cmap=CMAP['a'], 
#            label='PM2.5 (ug/m$^2$)',
#            title='PM2.5',
#            config=config)

            print('\t{} - relative humidity (Ireland)'.format(i))
            figs += ['rh-ir_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.rh, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['ireland'],
                   levels=LEVELS['rh'], cmap=CMAP['rh'], label='%',
                   format='%.1f', title='Relative Humidity', colorbar=is_lst,
                   config=config)
            print('\t{} - relative humidity (Europe)'.format(i))
            figs += ['rh-e_{}'.format(t_save)]
            plot2d(d.lon, d.lat, d.rh, fig=figs[-1],
                   newfig=is_fst, t=t, extent=extents['europe'],
                   levels=LEVELS['rh'], cmap=CMAP['rh'], label='%',
                   format='%.1f', title='Relative Humidity', colorbar=is_lst,
                   config=config)

        for fig in figs:
            plt.figure(fig)
            plt.savefig(os.path.join(config['imgs'], fig))
            plt.close(fig)


if __name__ == '__main__':
    main('./config.json')
