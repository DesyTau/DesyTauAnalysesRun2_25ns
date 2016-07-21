#!/bin/zsh

for dir in `ls | grep _files` 
do
    cd $dir
    if [ ! -z $1 ]
       then
       cp $1 ./qsub.sh
    fi
    echo $dir
    for file in `ls *.zsh`
    do
	list=`echo $file | sed 's/.zsh//g'`
	checkRoot ${list}_*.root
	exit=$?
	if [[ $exit -ne 0 ]]; then
	   echo "resubmit" $list
	   echo reason $exit

	   params=`while read line; do; if [[ -n $line ]] then; if [[ $line[1] != "#" ]] then; echo $line; break; fi; fi; done < $file`
	   ./qsub.sh $params
	fi
    done
    cd ../
done
   
