#!/usr/bin/env sh
#
# go to relevant directory 
emit_dir="/mnt/raid/wrf-chem/emit"
run_dir="/mnt/raid/wrf-chem/WRFV3/run"
ref_emit_dir="/mnt/raid/wrf-chem/emit/ref_emit_2018"
cd $emit_dir


#copy
cp $run_dir/wrfinput_d01 .

echo this file to be processed: $filename

#Define times(=simulation start time): same for all days of simulation
ncks -A wrfinput_d01 -v Times $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00
ncks -A wrfinput_d01 -v Times $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00
ncks -A wrfinput_d01 -v Times $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00
d0=$(date +%Y-%m-%d)
d1=$(date +%Y-%m-%d --date=+1day)
d2=$(date +%Y-%m-%d --date=+2day)
d3=$(date +%Y-%m-%d --date=+3day)
d4=$(date +%Y-%m-%d --date=+4day)
echo $d $d0 $d1 $d2 $d3
#d0
filename_wrfchem0=wrfchemi_d01_${d0}_00:00:00
filename_wrffire0=wrffirechemi_d01_${d0}_00:00:00
filename_wrfgocart0=wrfchemi_gocart_bg_d01_${d0}_00:00:00
cp $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrffire0
cp $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfchem0
cp $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfgocart0
#d1
filename_wrfchem1=wrfchemi_d01_${d1}_00:00:00
filename_wrffire1=wrffirechemi_d01_${d1}_00:00:00
filename_wrfgocart1=wrfchemi_gocart_bg_d01_${d1}_00:00:00
cp $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrffire1
cp $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfchem1
cp $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfgocart1
#d2
filename_wrfchem2=wrfchemi_d01_${d2}_00:00:00
filename_wrffire2=wrffirechemi_d01_${d2}_00:00:00
filename_wrfgocart2=wrfchemi_gocart_bg_d01_${d2}_00:00:00
cp $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrffire2
cp $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfchem2
cp $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfgocart2
#d3
filename_wrfchem3=wrfchemi_d01_${d3}_00:00:00
filename_wrffire3=wrffirechemi_d01_${d3}_00:00:00
filename_wrfgocart3=wrfchemi_gocart_bg_d01_${d3}_00:00:00
cp $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrffire3
cp $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfchem3
cp $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfgocart3

#d4
filename_wrfchem4=wrfchemi_d01_${d4}_00:00:00
filename_wrffire4=wrffirechemi_d01_${d4}_00:00:00
filename_wrfgocart4=wrfchemi_gocart_bg_d01_${d4}_00:00:00
cp $ref_emit_dir/wrffirechemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrffire4
cp $ref_emit_dir/wrfchemi_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfchem4
cp $ref_emit_dir/wrfchemi_gocart_bg_d01_2018-03-05_00:00:00 $emit_dir/$filename_wrfgocart4

cd $run_dir

ln -sf $emit_dir/wrf* .


