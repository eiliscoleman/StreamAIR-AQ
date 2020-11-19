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
for f in $raw_dir/wrfout_d01_*; do
cd $run_dir
#for f in $raw_dir/wrfout_d01_2017-07-28*; do
    time=${f:59:19}
    year=${time:0:4}
    month=${time:5:2}
    day=${time:8:2}
    hour=${time:11:2}
    datestr=${year}-${month}-${day}-${hour}:00:00
    filedatestr=${year}${month}${day}${hour}00
    fileroot=${f:55:4}$filedatestr
    filename=${opwrf_dir}/${fileroot}.nc
    echo $filename
    if [ $hourinc -lt 19 ]; then
       mv ${filename} ${data_back_dir}/
       echo scrap this data $hourinc and date $datestr
    else
       cp $f $run_dir/wrfout_d01_org.nc
       cp $filename $run_dir/d01_org.nc
# on d file: set time, lon and u10, v10, see if we can use earthjsonfromwrf python
       cdo showtime $run_dir/wrfout_d01_org.nc > $run_dir/filetime
       cdo showdate $run_dir/wrfout_d01_org.nc > $run_dir/filedate
       readdatefile=$run_dir/filedate
       readtimefile=$run_dir/filetime
       echo 'flag 1'
       while read -r line; do
       datet="$line"
       done < "$readdatefile"
       while read -r line;do
       timet="$line"
       done < "$readtimefile"
       timestamp=$datet' '$timet
       cdo settime,$timet $run_dir/d01_org.nc d01_timeset.nc
       cdo setdate,$datet d01_timeset.nc d01_timestamped
       rm $readdatefile $readtimefile
       echo  "****************************************************UV"
# uv 
       cdo selvar,u10,v10 d01_timestamped d01_vars
       ncatted -O -a standard_name,lon,o,c,longitude d01_vars
       ncatted -O -a standard_name,lat,o,c,latitude d01_vars
       cdo remapbil,mygrid d01_vars d01_reg
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python earthjsonfromwrf_uv.py d01_reg2 .
       mv 9999-wind-surface-level-gfs-1.0.json $dir_out/$time-wind-10m-wrf-chem-3.5.json

echo  "****************************************************T2"

# t2 
    cdo selvar,t2 d01_timestamped d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_t2.py d01_reg1 .
    mv 9999-t2-surface-level-gfs-1.0.json $dir_out/$time-t2-wrf-chem-3.5.json
    echo  "****************************************************RAIN"


# rain 
    cdo selvar,rain d01_timestamped d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_rain.py d01_reg1 .
    mv 9999-rain-surface-level-gfs-1.0.json $dir_out/$time-rain-wrf-chem-3.5.json
    echo  "****************************************************PSL"

# psl
    cdo selvar,p_sl d01_timestamped d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_psl.py d01_reg1 .
    mv 9999-psl-surface-level-gfs-1.0.json $dir_out/$time-psl-wrf-chem-3.5.json
    echo  "****************************************************RH"

# rh
    cdo selvar,rh d01_timestamped d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_rh.py d01_reg1 .
    mv 9999-rh-surface-level-gfs-1.0.json $dir_out/$time-rh-wrf-chem-3.5.json
    echo  "****************************************************SO2"
# so2
    cdo sellevel,29 d01_timestamped d01_surf
    cdo selvar,so2_concentration d01_surf d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_so2.py d01_reg1 .
    mv 9999-so2-surface-level-gfs-1.0.json $dir_out/$time-so2-wrf-chem-3.5.json
    echo  "****************************************************PM25"

# pm25
    cdo selvar,pm25 d01_surf d01_vars
    ncatted -O -a standard_name,lon,o,c,longitude d01_vars
    ncatted -O -a standard_name,lat,o,c,latitude d01_vars
    cdo remapbil,mygrid d01_vars d01_reg
    ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_reg d01_reg1
    python earthjsonfromwrf_pm25.py d01_reg1 .
    mv 9999-pm25-surface-level-gfs-1.0.json $dir_out/$time-pm25-wrf-chem-3.5.json

    echo  "****************************************************JSONS done"
#  i=$((i+1))
   cdo selvar,u10,v10,rain,t2,p_sl,rh,pm25,so2_concentration d01_timestamped d01_allvars
   ncatted -O -a standard_name,lon,o,c,longitude d01_allvars
   ncatted -O -a standard_name,lat,o,c,latitude d01_allvars
   cdo remapbil,mygrid d01_allvars d01_allvars_remapped
   ieg_file=${dir_interim}/${fileroot}_remapped.ieg
   cdo -f ieg copy d01_allvars_remapped $ieg_file
   cdo gradsdes $ieg_file
   ctl_file=${fileroot}_remapped.ctl
   newline="ctlfile='${fileroot}'"
   echo $newline
   sed -i '14s/.*/'"${newline}"'/'  gen_plots.gs
   grads -a 1.1873166 -pbcx gen_plots.gs
   new_windspeed_png=w10_${timet}.png
   mv -f ${dir_out}/windspeed_today.png  ${dir_out}/${datestr}-wind-10m.png
   mv -f ${dir_out}/rain_today.png  ${dir_out}/${datestr}-rain.png
   mv -f ${dir_out}/t2_today.png  ${dir_out}/${datestr}-t2.png
   mv -f ${dir_out}/psl_today.png  ${dir_out}/${datestr}-psl.png
   mv -f ${dir_out}/rh_today.png  ${dir_out}/${datestr}-rh.png
   mv -f ${dir_out}/so2_today.png  ${dir_out}/${datestr}-so2.png
   mv -f ${dir_out}/pm_today.png  ${dir_out}/${datestr}-pm25.png
   #clean up 
   rm -f $ieg_file
   rm -f ${dir_interim}'/'$ctl_file
   # move the original d01* file to 'processed' folder
   mv -f ${filename} ${data_back_dir}/
   rm -f d01* 
   rm -f wrfout* 
   fi
   hourinc=$((hourinc+1))

done
rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/

#done
