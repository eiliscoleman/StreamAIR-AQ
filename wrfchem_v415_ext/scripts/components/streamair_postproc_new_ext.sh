#!/bin/bash
basdir=/mnt/raid/wrf-chem/wrfchem_v415_ext
raw_dir=$basdir/data/output-wrf-raw
opwrf_dir=$basdir/data/output-wrf
remapped_dir=$basdir/data/output-wrf-remapped
data_back_remapped_dir=$basdir/data_back/output-wrf-remapped
# opwrf_dir=/mnt/raid/rong-ming/wrfchem/data_in_dfiles
run_dir=$basdir/scripts/components
dir_out=$basdir/web
# dir_finaldestination=$basdir/web
gen_data_back=$basdir/data_back
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

hourlist='00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23' 
# hourlist='12' 
#loop to be parallelised
process_pgm () {
local runhr=$1
# for filename in $raw_dir/wrfout_d01_????-??-??_*${runhr}":00:00"; 
# do



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
cdo selvar,COSALPHA wrfchem_org cosalpha
cdo selvar,SINALPHA wrfchem_org cosalpha
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
# python $run_dir/produce_newplots.py d01_reg $runhrdir $dir_out $dateid
 python $run_dir/produce_merc_plots_jsons.py d01_reg $runhrdir $dir_out $dateid 

# ==============================
# # Create json files for all parameters
#        echo  "****************************************************JSONS"
       echo  "uv"
# uv for leaflet javascript framework
cdo selvar,u10,v10 d01_reg d01_uv
ncrename -O -v lon,XLONG -v lat,XLAT -v time,Times d01_uv d01_reg1
ncpdq -O -h -a -lat d01_reg1 d01_reg2
python $run_dir/earthjsonfromwrf_uv.py d01_reg2 .
mv 9999-wind-surface-level-gfs-1.0.json $dir_out/${dateid}:00-wind-10m-wrf-chem-3.5.json
find d* -type f -not -name d01_reg -delete


cd ../
rm -r run_$runhr

echo 'process complete'

else 
echo "File does not exist" $filename

fi 

done
}
# 
# 
for runhr in $hourlist; do process_pgm "$runhr" &  done
wait
# 
#-----------------------------------------
#run python file for making rain pngs and jsons
python3 produce_rain.py 
#================================================
rm ${rainfiledir}/* ${pslfiledir}/* ${snowfiledir}/*

# # ===================================================
# # overwrite the generated aqi files
python3 produce_aq_who_sp_plots.py $basdir
python3 produce_aqi_eu.py $basdir
# #now, move only relevant files to the output web folder



# hourinc=0
# echo 'sorting out the stuff'
# make directory in data_back_dir, named today's date
bkupdirbranch=$(date +'%Y%m%d')
bkupdir=$gen_data_back"/"$bkupdirbranch
mkdir $bkupdir
mkdir $bkupdir"/output-wrf"
file_dir=$opwrf_dir"/*"

for fin_fil in $file_dir; do
yyyy=${fin_fil#*"d01_"}
datestr=${yyyy:0:12}
dfilename=$opwrf_dir"/d01_"$datestr".nc"
mv -f ${dfilename} $bkupdir"/output-wrf"
done
 
# ========================================== 
# echo '5'
rsync -zavu -e "ssh -p 8222"  --include ${year}"*.json"  --include ${year}"*.png" $dir_out/ update@140.203.204.132:/home/www/html/rt/weather/


