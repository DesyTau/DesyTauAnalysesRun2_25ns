#!/bin/bash


cp analysisMacroSynch_lept_tt_MC17.conf analysisMacroSynch_lept_tt_backup.conf

sed 's/pileUpforMC =/pileUpforMC = DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY1JetsToLL.conf

sed 's/pileUpforMC =/pileUpforMC = DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY1JetsToLL_ext1.conf

sed 's/pileUpforMC =/pileUpforMC = DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY2JetsToLL.conf

sed 's/pileUpforMC =/pileUpforMC = DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY2JetsToLL_ext1.conf

sed 's/pileUpforMC =/pileUpforMC = DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY3JetsToLL.conf

sed 's/pileUpforMC =/pileUpforMC = DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY3JetsToLL_ext1.conf

sed 's/pileUpforMC =/pileUpforMC = DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DY4JetsToLL.conf

sed 's/pileUpforMC =/pileUpforMC = DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DYJetsToLL.conf

sed 's/pileUpforMC =/pileUpforMC = DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_DYJetsToLL_ext1.conf

sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ST_t-channel_antitop.conf

sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ST_t-channel_top.conf

sed 's/pileUpforMC =/pileUpforMC = ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ST_tW_antitop.conf

sed 's/pileUpforMC =/pileUpforMC = ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ST_tW_top.conf

sed 's/pileUpforMC =/pileUpforMC = TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_TTTo2L2Nu.conf

sed 's/pileUpforMC =/pileUpforMC = TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_TTToHadronic.conf

sed 's/pileUpforMC =/pileUpforMC = TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_TTToSemiLeptonic.conf

sed 's/pileUpforMC =/pileUpforMC = W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_W1JetsToLNu.conf

sed 's/pileUpforMC =/pileUpforMC = W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_W2JetsToLNu.conf

sed 's/pileUpforMC =/pileUpforMC = W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_W3JetsToLNu.conf

sed 's/pileUpforMC =/pileUpforMC = W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_W4JetsToLNu.conf

sed 's/pileUpforMC =/pileUpforMC = WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_WJetsToLNu.conf

sed 's/pileUpforMC =/pileUpforMC = WW_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_WW.conf

sed 's/pileUpforMC =/pileUpforMC = WZ_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_WZ.conf

sed 's/pileUpforMC =/pileUpforMC = ZZ_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ZZ.conf

sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_ggH125.conf

sed 's/pileUpforMC =/pileUpforMC = VBFHToTauTau_M125_13TeV_powheg_pythia8_pileup/' analysisMacroSynch_lept_tt_MC17.conf > analysisMacroSynch_lept_tt_VBF125.conf


for file in $(ls *.conf)
do
    echo "**********************************"
    echo $file
    echo "**********************************"

    diff analysisMacroSynch_lept_tt_MC17.conf $file
done