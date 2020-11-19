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
#    if [ $hourinc -gt 19 ]; then
#       mv ${filename} ${data_back_dir}/
#       echo scrap this data $hourinc and date $datestr
#    else
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
       # from here, plot all the png files required directly from d file: 
#       python produce_plots_args.py d01_reg
#       mv -f test_t2.png $dir_out/$datestr-t2.png
#       mv -f test_rain_psl.png $dir_out/$datestr-rain_psl.png
#       mv -f test_rh.png $dir_out/$datestr-rh.png
#       mv -f test_o3.png $dir_out/$datestr-o3.png
#       mv -f test_nox.png $dir_out/$datestr-nox.png
#       mv -f test_so2.png $dir_out/$datestr-so2.png
#       mv -f test_pm10.png $dir_out/$datestr-pm10.png
#       mv -f test_pm25.png $dir_out/$datestr-pm25.png
       echo  "****************************************************UV"
# uv 
       cdo selvar,u10,v10 d01_reg d01_uv
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_uv d01_reg1
       python earthjsonfromwrf_uv.py d01_reg1 .
       mv 9999-wind-surface-level-gfs-1.0.json $dir_out/$datestr-wind-10m-wrf-chem-3.5.json
#    fi
#   mv -f ${filename} ${data_back_dir}/
   rm -f d01*  
   hourinc=$((hourinc+1))

done
rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/

