#!/bin/bash
# trigger option
# $1 - trigger
trigger=$1

for era in 2016 2017 2018
do
    for dataset in Data EMB TTbar DYJets EWK SMHiggs MSSMHiggs
    do
	for category in em_inclusive em_ttbar_control em_btag_highdzeta em_btag_mediumdzeta em_btag_lowdzeta em_nobtag_highdzeta em_nobtag_mediumdzeta em_nobtag_lowdzeta em_nobtag_highmsv_highdzeta em_nobtag_highmsv_mediumdzeta em_nobtag_highmsv_lowdzeta em_btag em_nobtag em_nobtag_highmsv em_ttbar_btag em_ttbar_nobtag
	do 
	    ./Datacards_submit.sh ${era} ${dataset} ${category} ${trigger}
	done
    done 
done
