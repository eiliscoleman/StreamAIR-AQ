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
       python $run_dir/produce_plots_args_rain_higherres.py d01_reg $runhrdir
       
       mv -f test_rain_psl.png $dir_out/$datestr-rain_psl.png
           mv -f test_nox.png $dir_out/$datestr-nox.png 


cd ../
rm -r run_$runhr
mv -f ${filename} ${data_back_dir}/
echo 'process complete'





else 
    echo '2'
    echo "File does not exist" $filename
fi 

done
}

for runhr in $hourlist; do process_pgm "$runhr" &  done
# 
# echo '3'
# #now, move only relevant files to the output web folder
# hourinc=0
# echo '4'
# for fin_fil in ${dir_out}"/2017"*"t2.png"; do
# # echo "${fin_fil}"
# yyyy=${fin_fil#*temp/}
# echo ${yyyy:0:13}
# if [ $hourinc -lt 19 ]; then
# #        mv ${filename} ${data_back_dir}/
#        echo scrap this data ${yyyy:0:13} $hourinc
# #        mv ${dir_out}/${yyyy:0:13}* ${dir_out}/scrap
# else
#        echo keep this data ${yyyy:0:13} $hourinc
# #        mv ${dir_out}/${yyyy:0:13}* $dir_finaldestination
# fi
# hourinc=$((hourinc+1))
# done
# 
# echo '5'
# rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/
