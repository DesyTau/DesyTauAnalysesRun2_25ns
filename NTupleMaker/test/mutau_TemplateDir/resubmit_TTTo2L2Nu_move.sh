#!/bin/sh

# This script loops over all folders in directory and checks if there is a *.root file in these folders that is corrupted. If it is it is resubmitted.

for dir in $(find -maxdepth 1 -type d -name "TTTo2L2Nu*_files")
do
    echo ''
    echo ''
    echo 'dir = ' $dir
    cd $dir
    samplename=${dir:2:${#dir}-8}
    echo 'sample name = ' $samplename
    for file in $(ls ${samplename}*.sh)
    do
	list=$(echo $file | sed 's/.sh//g')
	piet="_0_mt_Sync"
	aap=${list}${piet}
	#echo "checking for ${aap}.root"
	checkRoot ${aap}.root
	exit=$?
	if [[ $exit -ne 0 ]]; then
	    echo reason $exit
#	    echo "resubmit" $list
#	    condor_submit $list.submit
	    mv ${aap}.root ${aap}___
	    echo ""
	else
	   echo "file 	${aap}.root looks in order"
	fi
    done
    cd ../
done

