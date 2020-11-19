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
    timestr=${filename:48:12}
    echo $timestr
    year=${timestr:0:4}
    month=${timestr:4:2}
    day=${timestr:6:2}
    hour=${timestr:8:2}
    datestr=${year}-${month}-${day}-${hour}:00:00
    echo $datestr
    mv ${filename} $data_back_dir
echo $filename ' moved to back directory'





else 
    echo '2'
    echo "File does not exist" $filename
fi 

done
#Check to see if the data is finished processing (data finished if opwrf_dir is empty): if so, kill this job
jobs -p
#aaa=$(printf "Empty:now we can move on\r\n")
#[ "$(ls -A $opwrf_dir)" ] && echo "Still Processing" || echo "$aaa"

}

for runhr in $hourlist; do process_pgm "$runhr"  & done
wait
echo 'this should be the end'
#

 
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
