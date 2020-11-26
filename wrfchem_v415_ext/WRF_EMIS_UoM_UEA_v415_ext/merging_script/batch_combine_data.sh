#!/bin/bash --login
#
#PBS -l select=serial=true:ncpus=1
#PBS -l walltime=12:00:00
#PBS -A n02-weat

cd $PBS_O_WORKDIR
module load ncl


# first read in the list of files to modify
k=0
while read temp
do
	#temp=$line
	echo "Line # $k: $temp"
	giles[$k]=$temp
	((k++))
done < file_list.txt

# then go through the list of files to modify
#   IMPORTANT: Don't combine these loops - over ncl will fuck it all up by
#              trying to read file_list.txt after it's read the proper script.
#              And I don't know how to stop this happening, other than by keeping the loops separate.
for (( g=0; g<$k; g++ )); do
	echo "Line # $g: ${giles[$g]}"
	sed -e '33 c\ file_name = "'${giles[$g]}'"'  combine_emissions.ncl > combine_emissions.temp.ncl ;
	ncl combine_emissions.temp.ncl ;
done