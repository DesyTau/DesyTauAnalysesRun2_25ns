#!/bin/bash


cp analysisMacroSynch_lept_mt_MC17.conf analysisMacroSynch_lept_mt_backup.conf

sed 's/pileUpforMC =/pileUpforMC = DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = MC_PU2017_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_DYJetsToLL_M-5to50_13TeV-12Apr2018.conf

sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf

sed 's/pileUpforMC =/pileUpforMC = TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf

sed 's/pileUpforMC =/pileUpforMC = W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf

sed 's/pileUpforMC =/pileUpforMC = WW_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_WW_TuneCP5_13TeV-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = WZ_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_WZ_TuneCP5_13TeV-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = ZZ_TuneCP5_13TeV-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_ZZ_TuneCP5_13TeV-pythia8.conf

sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = VBFHToTauTau_M125_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_VBF125.conf

#Tauspinner samples
sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8_TauSpinner.conf
sed -i 's/applyTauSpinnerWeights=false/applyTauSpinnerWeights=true/g' analysisMacroSynch_lept_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8_TauSpinner.conf

sed 's/pileUpforMC =/pileUpforMC = VBFHToTauTau_M125_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_VBF125_TauSpinner.conf
sed -i 's/applyTauSpinnerWeights=false/applyTauSpinnerWeights=true/g' analysisMacroSynch_lept_mt_VBF125_TauSpinner.conf


sed 's/pileUpforMC =/pileUpforMC = MC_PU2017_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_WGToLNuG.conf

sed 's/pileUpforMC =/pileUpforMC = EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_EWKWMinus.conf
sed 's/pileUpforMC =/pileUpforMC = EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_EWKWPlus.conf
sed 's/pileUpforMC =/pileUpforMC =  EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_EWKZ2Jets_ZToLL.conf
sed 's/pileUpforMC =/pileUpforMC =  EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_pileup/' analysisMacroSynch_lept_mt_MC17.conf > analysisMacroSynch_lept_mt_EWKZ2Jets_ZToNuNu.conf
