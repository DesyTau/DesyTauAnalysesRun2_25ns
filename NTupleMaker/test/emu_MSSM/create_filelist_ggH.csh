#!/bin/bash
masses=(80
100
120
125
130
140
160
180
200
250
300
350
400
450
500
600
700
800
900
1000
1200
1400
1600
1800
2000
2300
2600
2900
3200 
3500)

i=0
for dataset in $(less ggh_2018);
do
    mass=${masses[$i]}
    python read_filelist_from_das.py --nick SUSYGluGluToHToTauTau_M${mass} --query "${dataset}" --outputfile 2018/SUSYGluGluToHToTauTau_powheg_M${mass}.list
    i=`expr $i + 1` 
done

j=0
for dataset in $(less ggh_2017);
do
    mass=${masses[$j]}
    python read_filelist_from_das.py --nick SUSYGluGluToHToTauTau_M${mass} --query "${dataset}" --outputfile 2017/SUSYGluGluToHToTauTau_powheg_M${mass}.list
    j=`expr $j + 1` 
done

k=0
for dataset in $(less ggh_2016);
do
    mass=${masses[$k]}
    python read_filelist_from_das.py --nick SUSYGluGluToHToTauTau_M${mass} --query "${dataset}" --outputfile 2016/SUSYGluGluToHToTauTau_powheg_M${mass}.list
    k=`expr $k + 1` 
done
