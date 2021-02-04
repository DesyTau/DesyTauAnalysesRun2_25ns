#!/bin/sh

# This script loops over all folders in directory and checks if there is a *.root file in these folders that is corrupted. If it is it is resubmitted.

for dir in $(find -maxdepth 1 -type d -name "*_files")
do
    if [[ $dir == *"xxxx"* ]]; then
	echo ".......................... don't process this: " $dir
	echo ""
	continue
    fi
    echo $dir
    cd $dir

    for file in $(ls ${dir%??????}*.sh)
    do

	list=$(echo $file | sed 's/.sh//g')
	checkRoot ${list}.root
	exit=$?
	if [[ $exit -ne 0 ]]; then
	        echo "resubmit" $list
		    echo reason $exit
                        condor_submit ${list}.submit

			fi
    done
    echo ""
    cd ../
done
