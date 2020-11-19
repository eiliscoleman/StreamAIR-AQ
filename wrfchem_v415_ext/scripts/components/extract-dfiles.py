#!/usr/bin/env python3
import bottleneck as bn
import datetime as dt
import glob as g
import json
import numpy as np
import os
import wrf
import xarray as xr
import xarray.ufuncs as xu


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    #    out, times, h_agl = xr.Dataset(), [], []
    for i, f_path in enumerate(sorted(g.glob(os.path.join(
        config['output-wrf-raw'], 'wrfout_*')
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
        #cosalpha and sinalpha to account for earth's rotation relative to the grid
        cosa=d.COSALPHA.values
        sina=d.SINALPHA.values
        
        out = xr.Dataset()
        out.coords['time'] = time
##        out.coords['h_agl'] = np.mean(h_agl, axis=0)
        out.coords['x'] = range(d.XLONG.shape[1])
        out.coords['y'] = range(d.XLAT.shape[0])
        out['lat'] = (('y', 'x'), d.XLAT.values[:])
        out['lon'] = (('y', 'x'), d.XLONG.values[:])
        out['terrain'] = (('y', 'x'), d.HGT.values)
        out['u10'] = (('y', 'x'), d.U10.values*cosa - d.V10.values * sina)
        out['v10'] = (('y', 'x'), d.V10.values * cosa + d.U10.values*sina)
        out['rain'] = (('y', 'x'), d.RAINC + d.RAINNC)
        out['snow'] = (('y', 'x'), d.SNOWNC)
        out['t2'] = (('y', 'x'), d.T2.values)
        out['so2_concentration'] = (
            ('h_agl', 'y', 'x'), d.so2.values) 
        out['o3_concentration'] = (
            ('h_agl', 'y', 'x'), d.o3.values) 
        out['nox_concentration'] = (
            ('h_agl', 'y', 'x'), (d.no2 + d.no).values) 
        # This was changed today to dimish pm2.5 by a factor of 4
        wind10=xu.sqrt(d.U10.values**2+d.V10.values**2)
        divfac=(wind10/5).clip(1,20)
        #.clip(1,10)
        #cappedpm25=d.PM2_5_DRY.values/divfac
        #out['pm25'] = (('h_agl', 'y', 'x'), cappedpm25)
        #newpm10=d.PM10.values-d.PM2_5_DRY.values+cappedpm25
        shp=d.so2.shape
        #print(shp)
        ss_01=(d.na_a01+d.cl_a01).values
        ss_02=(d.na_a02+d.cl_a02).values
        ss_03=(d.na_a03+d.cl_a03).values
        ss_04=(d.na_a04+d.cl_a04).values
        ss_05=(d.na_a05+d.cl_a05).values
        ss_06=(d.na_a06+d.cl_a06).values
        ss_25=ss_01+ss_02+ss_03+ss_04+ss_05+ss_06
        ss_new=ss_25/divfac
        cappedpm25=d.PM2_5_DRY.values-ss_25+ss_new
        out['pm25'] = (('h_agl', 'y', 'x'), cappedpm25)
        newpm10=d.PM10.values-d.PM2_5_DRY.values+cappedpm25
        out['pm10'] = (('h_agl', 'y', 'x'), newpm10)
        #out['ss_25'] = (('h_agl', 'y', 'x'), ss_25)
        #out['newss'] = (('h_agl', 'y', 'x'), ss_new)
        #out['ss_25'] = (('h_agl', 'y', 'x'), ss_25)
        #out['wind10'] = (('y', 'x'), wind10)
        #out['divfac']=(('y', 'x'), divfac)
        #organics
        org_01=d.oc_a01.values
        org_02=d.oc_a02.values
        org_03=d.oc_a03.values
        org_04=d.oc_a04.values
        org_05=d.oc_a05.values
        org_06=d.oc_a06.values
        org_07=d.oc_a07.values
        org_08=d.oc_a08.values
        org25=org_01+org_02+org_03+org_04+org_05+org_06+org_07+org_08
        out['org25'] = (('h_agl', 'y', 'x'), org25)
        #sulphate
        sulf_01=d.so4_a01.values
        sulf_02=d.so4_a02.values
        sulf_03=d.so4_a03.values
        sulf_04=d.so4_a04.values
        sulf_05=d.so4_a05.values
        sulf_06=d.so4_a06.values
        sulf_07=d.so4_a07.values
        sulf_08=d.so4_a08.values
        sulf25=sulf_01+sulf_02+sulf_03+sulf_04+sulf_05+sulf_06+sulf_07+sulf_08
        out['sulf25'] = (('h_agl', 'y', 'x'), sulf25)
        #nitrate
        nitr_01=d.no3_a01.values
        nitr_02=d.no3_a02.values
        nitr_03=d.no3_a03.values
        nitr_04=d.no3_a04.values
        nitr_05=d.no3_a05.values
        nitr_06=d.no3_a06.values
        nitr_07=d.no3_a07.values
        nitr_08=d.no3_a08.values
        nitr25=nitr_01+nitr_02+nitr_03+nitr_04+nitr_05+nitr_06+nitr_07+nitr_08
        out['nitr25'] = (('h_agl', 'y', 'x'), nitr25)
        #Ammonium
        nh4_01=d.nh4_a01.values
        nh4_02=d.nh4_a02.values
        nh4_03=d.nh4_a03.values
        nh4_04=d.nh4_a04.values
        nh4_05=d.nh4_a05.values
        nh4_06=d.nh4_a06.values
        nh4_07=d.nh4_a07.values
        nh4_08=d.nh4_a08.values
        nh425=nh4_01+nh4_02+nh4_03+nh4_04+nh4_05+nh4_06+nh4_07+nh4_08
        out['nh425'] = (('h_agl', 'y', 'x'), nh425)
        #Chloride
        cl_01=d.cl_a01.values
        cl_02=d.cl_a02.values
        cl_03=d.cl_a03.values
        cl_04=d.cl_a04.values
        cl_05=d.cl_a05.values
        cl_06=d.cl_a06.values
        cl_07=d.cl_a07.values
        cl_08=d.cl_a08.values
        #24/2/20:LC: have just pm2.5, so up to 6th bin
        cl25=cl_01+cl_02+cl_03+cl_04+cl_05+cl_06
        out['cl25'] = (('h_agl', 'y', 'x'), cl25)
        #Black Carbon
        bc_01=d.bc_a01.values
        bc_02=d.bc_a02.values
        bc_03=d.bc_a03.values
        bc_04=d.bc_a04.values
        bc_05=d.bc_a05.values
        bc_06=d.bc_a06.values
        bc_07=d.bc_a07.values
        bc_08=d.bc_a08.values
        bc=bc_01+bc_02+bc_03+bc_04+bc_05+bc_06+bc_07+bc_08
        out['bc'] = (('h_agl', 'y', 'x'), bc)
        #print(d.PM2_5_DRY[1:10], d.PM10[1:10])
        #aaa=np.subtract(d.PM10.values, d.PM2_5_DRY.values)
        #out['aaa']=(('h_agl', 'y', 'x'), d.PM10[:,:,:]-d.PM2_5_DRY[:,:,:])
        #pm25=d.PM2_5_DRY.values
        #newpm25=pm25/4.
        #pm10=d.PM10.values
        #print(pm25[0:10,0:10,0:10])
        #print(pm10[0:10,0:10,0:10])
        #np.warnings.filterwarnings('ignore')
        #newpm10=(pm10-pm25+pm25/4.)
        #print(pm25.shape)
        #out['pm25'] = (
            #('h_agl', 'y', 'x'),d.PM2_5_DRY.values) 
        #out['pm10'] = (('h_agl', 'y', 'x'), (d.PM10-d.PM2_5_DRY).values+d.PM2_5_DRY.values/4)
        #out['pm10'] = (('h_agl', 'y', 'x'), d.PM10.values)
       
        #('h_agl', 'y', 'x'), (d.PM10.values-d.PM2_5_DRY.values/4.))
            #-d.PM2_5_DRY).values) 
##            wrf.x_to_yOm3(d.so2.values, (d.PB + d.P).values,
##                          d['T'].values, mm=64) 
##        )
        out['pb'] = (
            ('h_agl', 'y', 'x'), d.PB.values) 
        out['p_sl'] = (('y', 'x'), wrf.slp(h, p, t, q))
        out['rh'] = (('y', 'x'), wrf.rh(p, t, q)[0])
        out['swdown'] = (('y', 'x'), d.SWDOWN.values)
        out['cldfra'] = (('h_agl', 'y', 'x'), d.CLDFRA.values)
        out['qsnow'] = (('h_agl', 'y', 'x'), d.QSNOW.values)
        out['qgraup'] = (('h_agl', 'y', 'x'), d.QGRAUP.values)
        out.to_netcdf(
            os.path.join(config['output-wrf'],
                         '{domain}_{date}.nc'
                         .format(domain=domain,
                                 date=(time.strftime('%Y%m%d%H%M'))))
        )


if __name__ == '__main__':
    main('./config.json')
