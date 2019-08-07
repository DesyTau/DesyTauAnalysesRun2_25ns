#!/bin/sh

# This script loops over all folders in directory and checks if there is a *.root file in these folders that is corrupted. If it is it is resubmitted.

# 1.) Find all jobs that should have run
find *_files/. -type f -name "*.sh" -not -name "HTC_submit.sh" > filelist_to_check.txt

# 2.) Split this list into four
n=$(cat filelist_to_check.txt | wc -l)
part=$((n/4+1))
split -l $part -d filelist_to_check.txt filelist_to_check_

path=$PWD
echo $path

# 3.) 
for files in ./filelist_to_check_*
do
    echo $files   
    while read line
    do

	searchstring="/./"
	dir=${line%$searchstring*}
	#echo $dir
	file=${line#*$searchstring} 
	#echo $file
	cd $dir
	list=$(echo $file | sed 's/.sh//g')
	checkRoot ${list}.root
	exit=$?
	if [[ $exit -ne 0 ]]; then
	    echo reason $exit
	    echo "resubmit" $list
	    condor_submit $list.submit
	    echo ""
	fi
	cd $path
	done< $files&
done

wait

echo ""
echo "Done"
echo ""