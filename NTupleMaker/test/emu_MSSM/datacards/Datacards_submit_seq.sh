
#!/bin/bash
for era in 2016 2017 2018
do
    for dataset in Data EMB TTbar DYJets EWK SMggH SMothers bbH ggh_t ggh_b ggh_i ggH_t ggH_b ggH_i ggA_t ggA_b ggA_i
    do
	for category in em_inclusive em_NbtagGt1 em_Nbtag0 em_DZetaLtm35 em_NbtagGt1_DZetaLtm35 em_Nbtag0_DZetaLtm35 em_Nbtag0_DZetaGt30 em_Nbtag0_DZetam10To30 em_Nbtag0_DZetam35Tom10 em_Nbtag0_DZetaGt30_MHGt250 em_Nbtag0_DZetam10To30_MHGt250 em_Nbtag0_DZetam35Tom10_MHGt250 em_NbtagGt1_DZetaGt30 em_NbtagGt1_DZetam10To30 em_NbtagGt1_DZetam35Tom10 em_Nbtag0_NjetLt1_DZetaGt30 em_Nbtag0_NjetLt1_DZetam10To30 em_Nbtag0_NjetLt1_DZetam35Tom10 em_NbtagGt1_NjetLt1_DZetaGt30 em_NbtagGt1_NjetLt1_DZetam10To30 em_NbtagGt1_NjetLt1_DZetam35Tom10
	do 
	    ./Datacards_submit.sh ${era} ${dataset} ${category} 
	done
    done 
done
