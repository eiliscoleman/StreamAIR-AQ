#!/bin/bash

raw_dir=/mnt/raid/rong-ming/wrfchem/data/output-wrf-raw
opwrf_dir=/mnt/raid/rong-ming/wrfchem/data/output-wrf
run_dir=/mnt/raid/rong-ming/wrfchem/scripts/components
dir_out=/mnt/raid/rong-ming/web
dir_interim=/mnt/raid/rong-ming/web/temp
data_back_dir=/mnt/raid/rong-ming/wrfchem/data_back/output-wrf
#cp -rf $run_dir/earthjsonfromwrf*.* .
#cp -rf $run_dir/mygrid .
#cp -rf $run_dir/gen_plots.gs .
hourinc=0
for filename in $opwrf_dir/d01_*; do
cd $run_dir
    timestr=${filename:48:12}
    echo $timestr
    year=${timestr:0:4}
    month=${timestr:4:2}
    day=${timestr:6:2}
    hour=${timestr:8:2}
    datestr=${year}-${month}-${day}-${hour}:00:00
    echo $datestr
    if [ $hourinc -lt 19 ]; then
#      This section discards first 19hours of data as spin-up
       echo scrap this data $hourinc and date $datestr
    else
       echo 'this data is being processed' $filename
       cp $filename $run_dir/d01_org.nc
       datet="$year-$month-$day"
       timet="$hour:00"
       timestamp=$datet' '$timet
       cdo settime,$timet $run_dir/d01_org.nc d01_timeset.nc
       cdo setdate,$datet d01_timeset.nc d01_timestamped
       ncatted -O -a standard_name,lon,o,c,longitude d01_timestamped
       ncatted -O -a standard_name,lat,o,c,latitude d01_timestamped
       cdo remapbil,mygrid d01_timestamped d01_reg
       # from here, plot all the png files required directly from d file using python script produce_plots_args.py
       python produce_plots_args.py d01_reg
       mv -f test_t2.png $dir_out/$datestr-t2-wrf-chem-3.5.png
       mv -f test_wind.png $dir_out/$datestr-wind-10m-wrf-chem-3.5.png
       mv -f test_rain_psl.png $dir_out/$datestr-rain_psl-wrf-chem-3.5.png
       mv -f test_rh.png $dir_out/$datestr-rh-wrf-chem-3.5.png
       mv -f test_o3.png $dir_out/$datestr-o3-wrf-chem-3.5.png
       mv -f test_nox.png $dir_out/$datestr-nox-wrf-chem-3.5.png
       mv -f test_so2.png $dir_out/$datestr-so2-wrf-chem-3.5.png
       mv -f test_pm10.png $dir_out/$datestr-pm10-wrf-chem-3.5.png
       mv -f test_pm25.png $dir_out/$datestr-pm25-wrf-chem-3.5.png
# Create json files for all parameters
       echo  "****************************************************JSONS"
# uv 
       cdo selvar,u10,v10 d01_reg d01_uv
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_uv d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_uv.py d01_reg2 .
       mv 9999-wind-surface-level-gfs-1.0.json $dir_out/$datestr-wind-10m-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#t2 json
       cdo selvar,t2 d01_reg d01_t2
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_t2 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_t2.py d01_reg2 .       
       mv 9999-t2-surface-level-gfs-1.0.json $dir_out/$datestr-t2-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#rain json
       cdo selvar,rain d01_reg d01_rain
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_rain d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_rain.py d01_reg2 . 
       mv 9999-rain-surface-level-gfs-1.0.json $dir_out/$datestr-rain-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#psl json
       cdo selvar,p_sl d01_reg d01_psl
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_psl d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_psl.py d01_reg2 . 
       mv 9999-psl-surface-level-gfs-1.0.json $dir_out/$datestr-psl-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#rh json
       cdo selvar,rh d01_reg d01_rh
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_rh d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_rh.py d01_reg2 . 
       mv 9999-rh-surface-level-gfs-1.0.json $dir_out/$datestr-rh-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
# CHEMICAL CONSTITUENTS
cdo sellevel,29 d01_reg d01_surf
#o3 json
       cdo selvar,o3_concentration d01_surf d01_surf_o3
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_o3 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_o3.py d01_reg2 . 
       mv 9999-o3-surface-level-gfs-1.0.json $dir_out/$datestr-o3-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2

#nox json
       cdo selvar,nox_concentration d01_surf d01_surf_nox
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_nox d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_nox.py d01_reg2 . 
       mv 9999-nox-surface-level-gfs-1.0.json $dir_out/$datestr-nox-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2


#so2 json
       cdo selvar,so2_concentration d01_surf d01_surf_so2
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_so2 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_so2.py d01_reg2 . 
       mv 9999-so2-surface-level-gfs-1.0.json $dir_out/$datestr-so2-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#pm10 json
       cdo selvar,pm10 d01_surf d01_surf_pm10
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_pm10 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_pm10.py d01_reg2 . 
       mv 9999-pm10-surface-level-gfs-1.0.json $dir_out/$datestr-pm10-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
#pm25 json
       cdo selvar,pm25 d01_surf d01_surf_pm25
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_pm25 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_pm25.py d01_reg2 . 
       mv 9999-pm25-surface-level-gfs-1.0.json $dir_out/$datestr-pm25-wrf-chem-3.5.json
       rm d01_reg1 d01_reg2
   fi
#   mv -f ${filename} ${data_back_dir}/
# Clean up
   rm -f d01*  
   hourinc=$((hourinc+1))

done
rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/

