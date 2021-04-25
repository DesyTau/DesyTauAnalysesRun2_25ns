#!/bin/bash
masses=(60
80
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
for dataset in $(less bbH_2018);
do
    mass=${masses[$i]}
    python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_M${mass} --query "${dataset}" --outputfile 2018/SUSYGluGluToBBHToTauTau_powheg_M${mass}.list
    i=`expr $i + 1` 
done
