#!/bin/bash
basdir=/mnt/raid/wrf-chem/wrfchem_v39
raw_dir=$basdir/data/output-wrf-raw
opwrf_dir=$basdir/data/output-wrf
remapped_dir=$basdir/data/output-wrf-remapped
data_back_remapped_dir=$basdir/data_back/output-wrf-remapped
# opwrf_dir=/mnt/raid/rong-ming/wrfchem/data_in_dfiles
run_dir=$basdir/scripts/components
dir_out=$basdir/web7
# dir_finaldestination=$basdir/web
data_back_dir=$basdir/data_back/output-wrf
# spinup_data_back_dir=$basdir/data_back/output-wrf/spinup
rainfiledir=$run_dir/rainfilesmm
snowfiledir=$run_dir/snowfilesmm
pslfiledir=$run_dir/psl
#process all files in a tempoary folder
cd $run_dir
# ==============================================
# if any files are stuck in remapped dir, copy them to the data_back remapped dir
cp $remapped_dir/* $data_back_remapped_dir

hourlist='01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 00' 
#loop to be parallelised
process_pgm () {
local runhr=$1
# for filename in $raw_dir/wrfout_d01_????-??-??_*${runhr}":00:00"; do



for filename in $raw_dir/wrfout_d01_????-??-??_*${runhr}":00:00"; do
Pathname="${filename}"
if [ -e "$filename" ]; then
    echo "File exists" $filename
runhrdir=$run_dir/run_$runhr
mkdir $runhrdir
cd $runhrdir
cp -rf $run_dir/mygrid .
cp $filename wrfout_org
# Extract a variable from wrffile to define grid
cdo selvar,T2 wrfout_org t2
cdo griddes t2 > griddes.txt
# Apply this grid to the dfile
datestr=${filename#*d01_}
echo $datestr
year=${datestr:0:4}
month=${datestr:5:2}
day=${datestr:8:2}
hr=${datestr:11:2}

datet="$year-$month-$day"
timet="$hr:00"
dateid="$datet-$timet"
dfilepath="$opwrf_dir/d01_$year$month$day${hr}00.nc"
echo $dfilepath
cp $dfilepath d01_org
cdo setgrid,griddes.txt d01_org d01_gridset
cdo settime,$timet d01_gridset d01_timeset.nc
cdo setdate,$datet d01_timeset.nc d01_timestamped
ncatted -O -a standard_name,lon,o,c,longitude d01_timestamped
ncatted -O -a standard_name,lat,o,c,latitude d01_timestamped
cdo remapbil,mygrid d01_timestamped d01_reg
remapped_f=$remapped_dir/d01reg_$year$month$day${hr}
echo hupdedoo $remapped_f
cp d01_reg $remapped_f
find d* -type f -not -name d01_reg -delete


# from here, plot all the png files required directly from d file using 
# =====================================
# python $run_dir/produce_merc_plots.py d01_reg $runhrdir $dir_out $dateid
# python $run_dir/produce_pm25_plots.py d01_reg $runhrdir $dir_out $dateid


# 
# #deal with rain files
#         echo "extracting raindata and moving to rainfile directory"
#         cdo selvar,rain d01_reg d01_rain
#         ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_rain d01_reg1
# #        This following step is  unnecceary
# # ncpdq -O -h -a -lat d01_reg1 d01_reg2 
#         mv d01_reg1 $rainfiledir/d01_$year$month$day${hr}00_rainmm
#         #deal with snow files
#         echo "extracting snowdata and moving to rainfile directory"
#         cdo selvar,snow d01_reg d01_snow
#         ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_snow d01_reg1
# #        This following step (now commented out) is  unnecceary:
# # ncpdq -O -h -a -lat d01_reg1 d01_reg2 
#         mv d01_reg1 $snowfiledir/d01_$year$month$day${hr}00_snowmm
# # #psl json
#        echo  "psl"
#        cdo selvar,p_sl d01_reg d01_psl
#        ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_psl d01_reg1
# #        ncpdq -O -h -a -lat d01_reg1 d01_reg2
#        cp d01_reg1 $pslfiledir/d01_$year$month$day${hr}00_psl
#        python $run_dir/touchpoint_psl_reduced.py d01_reg1 . 
#        mv 9999-psl-surface-level-gfs-1.0.json $dir_out/${dateid}:00-psl-wrf-chem-3.5.json
#        find d* -type f -not -name d01_reg -delete
cd ../
rm -r run_$runhr

echo 'process complete'





else 
    echo "File does not exist" $filename
fi 

done
}


for runhr in $hourlist; do process_pgm "$runhr" &  done
wait
# # 
#-----------------------------------------
# #run python file for making rain pngs and jsons
# python3 produce_rain_snow_rate.py  ${rainfiledir} ${snowfiledir} $dir_out $pslfiledir
# #================================================
# rm ${rainfiledir}/* ${pslfiledir}/* ${snowfiledir}/*
# 
# # ===================================================
# # overwrite the generated aqi files
python3 produce_aqi.py $basdir
# # #now, move only relevant files to the output web folder
# 
# 
# 
# # hourinc=0
# # echo 'sorting out the stuff'
#  for fin_fil in ${dir_out}"/"*"t2.png"; do
# # # echo "${fin_fil}"
# yyyy=${fin_fil#*web/}
# datestr=${yyyy:0:13}
# yrstr=${datestr:0:4}
# mnstr=${datestr:5:2}
# daystr=${datestr:8:2}
# hrstr=${datestr:11:2}
# dfilename=$opwrf_dir"/d01_"$yrstr$mnstr$daystr$hrstr"00.nc"
# remapped_f=$remapped_dir"/d01reg_"$yrstr$mnstr$daystr${hrstr}
# echo 'flat1'
# echo $dfilename
# mv -f ${dfilename} ${data_back_dir}/
# echo 'flat2'
# # mv -f ${remapped_f} $data_back_remapped_dir/
# 
#  hourinc=$((hourinc+1))
# 
#  done
# ========================================== 
# echo '5'
# rsync -zavu -e "ssh -p 8222"  --include ${year}"*.json"  --include ${year}"*.png" $dir_out/ update@140.203.204.132:/home/www/html/rt/weather/
# # rsync -zavu -e "ssh -p 8222"  --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" $dir_out/* update@140.203.204.132:/home/www/html/rt/weather/
# # Housekeeping: delete old files
# # find $dir_finaldestination/*.png -mtime +30 -exec rm -f {} \;
# # find $dir_finaldestination/*.json -mtime +30 -exec rm -f {} \;
