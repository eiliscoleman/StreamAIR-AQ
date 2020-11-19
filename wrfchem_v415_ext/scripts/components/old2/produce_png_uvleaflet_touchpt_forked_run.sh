#!/bin/bash

raw_dir=/mnt/raid/rong-ming/wrfchem/data/output-wrf-raw
opwrf_dir=/mnt/raid/rong-ming/wrfchem/data/output-wrf
run_dir=/mnt/raid/rong-ming/wrfchem/scripts/components
dir_out=/mnt/raid/rong-ming/web/temp
dir_finaldestination=/mnt/raid/rong-ming/web
data_back_dir=/mnt/raid/rong-ming/wrfchem/data_back/output-wrf
#process all files in a tempoary folder
cd $run_dir

hourlist='00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23' 
#loop to be parallelised
process_pgm () {
local runhr=$1
for filename in $opwrf_dir/d01_*${runhr}"00.nc"; do
echo '1'



Pathname="${filename}"
timestr="${Pathname:48:12}"
if [ -e "$filename" ]; then
    echo "File exists" $filename
runhrdir=$run_dir/run_$runhr
mkdir $runhrdir
cd $runhrdir
cp -rf $run_dir/mygrid .
    timestr=${filename:48:12}
    echo $timestr
    year=${timestr:0:4}
    month=${timestr:4:2}
    day=${timestr:6:2}
    hour=${timestr:8:2}
    datestr=${year}-${month}-${day}-${hour}:00:00
    echo $datestr
#     if [ $hourinc -lt 19 ]; then
#      This section discards first 19hours of data as spin-up
#        echo scrap this data $hourinc and date $datestr
#     else
       echo 'this data is being processed' $filename
       cp $filename d01_org.nc
       datet="$year-$month-$day"
       timet="$hour:00"
       timestamp=$datet' '$timet
       cdo settime,$timet d01_org.nc d01_timeset.nc
       cdo setdate,$datet d01_timeset.nc d01_timestamped
       ncatted -O -a standard_name,lon,o,c,longitude d01_timestamped
       ncatted -O -a standard_name,lat,o,c,latitude d01_timestamped
       cdo remapbil,mygrid d01_timestamped d01_reg
       find d* -type f -not -name d01_reg -delete


       # from here, plot all the png files required directly from d file using python script produce_plots_args.py
#        cp d01_reg $run_dir
       python $run_dir/produce_plots_args.py d01_reg $runhrdir
       mv -f test_t2.png $dir_out/$datestr-t2.png
       mv -f test_wind.png $dir_out/$datestr-wind-10m.png
       mv -f test_rain_psl.png $dir_out/$datestr-rain_psl.png
       mv -f test_rh.png $dir_out/$datestr-rh.png
       mv -f test_o3.png $dir_out/$datestr-o3.png
       mv -f test_nox.png $dir_out/$datestr-nox.png
       mv -f test_so2.png $dir_out/$datestr-so2.png
       mv -f test_pm10.png $dir_out/$datestr-pm10.png
       mv -f test_pm25.png $dir_out/$datestr-pm25.png
# Create json files for all parameters

# Create json files for all parameters
       echo  "****************************************************JSONS"
       echo  "uv"
# uv for leaflet javascript framework
       cdo selvar,u10,v10 d01_reg d01_uv
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_uv d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/earthjsonfromwrf_uv.py d01_reg2 .
       mv 9999-wind-surface-level-gfs-1.0.json $dir_out/$datestr-wind-10m-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
# 10mwind json
       cdo selvar,u10,v10 d01_reg d01_uv
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_uv d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_w10.py d01_reg2 .
       mv 9999-w10-surface-level-gfs-1.0.json $dir_out/$datestr-w10-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
#t2 json
       echo  "t2"
       cdo selvar,t2 d01_reg d01_t2
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_t2 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_t2.py d01_reg2 .        
       mv 9999-t2-surface-level-gfs-1.0.json $dir_out/$datestr-t2-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
#rain json
       echo  "rain"
       cdo selvar,rain d01_reg d01_rain
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_rain d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_rain.py d01_reg2 . 
       mv 9999-rain-surface-level-gfs-1.0.json $dir_out/$datestr-rain-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
#psl json
       echo  "psl"
       cdo selvar,p_sl d01_reg d01_psl
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_psl d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_psl.py d01_reg2 . 
       mv 9999-psl-surface-level-gfs-1.0.json $dir_out/$datestr-psl-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
#rh json
       echo  "rh"
       cdo selvar,rh d01_reg d01_rh
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_rh d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_rh.py d01_reg2 . 
       mv 9999-rh-surface-level-gfs-1.0.json $dir_out/$datestr-rh-wrf-chem-3.5.json
       find d* -type f -not -name d01_reg -delete
# CHEMICAL CONSTITUENTS
       cdo sellevel,1 d01_reg d01_surf
#o3 json
       echo  "o3"
       cdo selvar,o3_concentration d01_surf d01_surf_o3
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_o3 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_o3.py d01_reg2 . 
       mv 9999-o3-surface-level-gfs-1.0.json $dir_out/$datestr-o3-wrf-chem-3.5.json
       find d* -type f -not -name d01_surf -delete
#nox json
       echo  "nox"
       cdo selvar,nox_concentration d01_surf d01_surf_nox
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_nox d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_nox.py d01_reg2 . 
       mv 9999-nox-surface-level-gfs-1.0.json $dir_out/$datestr-nox-wrf-chem-3.5.json
       find d* -type f -not -name d01_surf -delete
#so2 json
       echo  "so2"
       cdo selvar,so2_concentration d01_surf d01_surf_so2
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_so2 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_so2.py d01_reg2 . 
       mv 9999-so2-surface-level-gfs-1.0.json $dir_out/$datestr-so2-wrf-chem-3.5.json
       find d* -type f -not -name d01_surf -delete
#pm10 json
       echo  "pm10"
       cdo selvar,pm10 d01_surf d01_surf_pm10
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_pm10 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_pm10.py d01_reg2 . 
       mv 9999-pm10-surface-level-gfs-1.0.json $dir_out/$datestr-pm10-wrf-chem-3.5.json
       find d* -type f -not -name d01_surf -delete
#pm25 json
       echo  "pm25"
       cdo selvar,pm25 d01_surf d01_surf_pm25
       ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_surf_pm25 d01_reg1
       ncpdq -O -h -a -lat d01_reg1 d01_reg2
       python $run_dir/touchpoint_pm25.py d01_reg2 . 
       mv 9999-pm25-surface-level-gfs-1.0.json $dir_out/$datestr-pm25-wrf-chem-3.5.json
       find d* -type f -not -name d01_surf -delete
cd ../
rm -r run_$runhr
mv -f ${filename} ${data_back_dir}/
echo 'process complete'





else 
    echo "File does not exist" $filename
fi 

done
}

for runhr in $hourlist; do process_pgm "$runhr" &  done
wait
# 
# #now, move only relevant files to the output web folder
hourinc=0

 for fin_fil in ${dir_out}"/2017"*"t2.png"; do
# # echo "${fin_fil}"
 yyyy=${fin_fil#*temp/}
 echo ${yyyy:0:13}
 if [ $hourinc -lt 19 ]; then
# #        mv ${filename} ${data_back_dir}/
        echo scrap this data ${yyyy:0:13} $hourinc
        mv -f ${dir_out}/${yyyy:0:13}* ${dir_out}/scrap
 else
        echo keep this data ${yyyy:0:13} $hourinc
        mv -f ${dir_out}/${yyyy:0:13}* $dir_finaldestination
 fi
 hourinc=$((hourinc+1))
 done
 
# echo '5'
 rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/
