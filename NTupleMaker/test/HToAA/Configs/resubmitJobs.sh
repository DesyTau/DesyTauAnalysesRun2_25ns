#!/bin/zsh

for dir in `ls | grep $1_files*` 
do
    cd $dir
    echo Inspecting directory $dir -------------
    for file in `ls *.submit`
    do
#	echo $file
	list=`echo $file | sed 's/.submit//g'`
#	echo checking root file ${list}.root
	checkRoot ${list}.root
	exit=$?
	if [[ $exit -ne 0 ]]; then
	   echo "resubmit" $list
	   echo reason $exit
	   condor_submit ${list}.submit
	fi
    done
    cd ../
done
   
