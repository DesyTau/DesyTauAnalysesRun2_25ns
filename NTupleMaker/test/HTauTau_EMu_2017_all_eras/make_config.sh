#!/bin/bash

sed 's/IsEmbedded = false/IsEmbedded = true/' analysisMacroSynch_em_DATA.conf > analysisMacroSynch_em_Embedded.conf
sed -i 's/IsData = true/IsData = false/g' analysisMacroSynch_em_Embedded.conf
sed 's/IsW = false/IsW = true/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_W.conf
sed 's/IsDY = false/IsDY = true/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_DY.conf
sed 's/IsSignal = false/IsSignal = true/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_Signal.conf

cp analysisMacroSynch_em_Signal.conf analysisMacroSynch_em_Signal_VBF.conf
cp analysisMacroSynch_em_Signal.conf analysisMacroSynch_em_Signal_ggZH.conf
sed -i 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal_VBF.conf
sed -i 's/ApplyVBFUncertainties = false/ApplyVBFUncertainties = true/g' analysisMacroSynch_em_Signal_VBF.conf
sed -i 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal_ggZH.conf


#

sed 's/SampleNameForPUHist =/SampleNameForPUHist = DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = DYJetsToLL_M-10to50_13TeV-12Apr2018_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_DYJetsToLL_M-10to50_13TeV-12Apr2018.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
#sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_DYJetsToLL_M-5to50_13TeV-12Apr2018.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = WW_TuneCP5_13TeV-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_WW_TuneCP5_13TeV-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = WZ_TuneCP5_13TeV-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_WZ_TuneCP5_13TeV-pythia8.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = ZZ_TuneCP5_13TeV-pythia8_/' analysisMacroSynch_em_MC.conf > analysisMacroSynch_em_ZZ_TuneCP5_13TeV-pythia8.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_GGH_Htautau_M125.conf
sed -i 's/ApplygghReweighting = false/ApplygghReweighting = true/g' analysisMacroSynch_em_Signal_GGH_Htautau_M125.conf 
sed -i 's/ApplygghUncertainties = false/ApplygghUncertainties = true/g' analysisMacroSynch_em_Signal_GGH_Htautau_M125.conf 
sed 's/SampleNameForPUHist = GluGluHToTauTau_M125_13TeV_powheg_pythia8_/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal_GGH_Htautau_M125.conf > analysisMacroSynch_em_Signal_ggh.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_VBF_Htautau_M125.conf
sed -i 's/ApplyVBFUncertainties = false/ApplyVBFUncertainties = true/g' analysisMacroSynch_em_Signal_VBF_Htautau_M125.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_GGH_HWW_M125.conf
sed -i 's/ApplygghReweighting = false/ApplygghReweighting = true/g' analysisMacroSynch_em_Signal_GGH_HWW_M125.conf
sed -i 's/ApplygghUncertainties = false/ApplygghUncertainties = true/g' analysisMacroSynch_em_Signal_GGH_HWW_M125.conf  
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_VBF_HWW_M125.conf
sed -i 's/ApplyVBFUncertainties = false/ApplyVBFUncertainties = true/g' analysisMacroSynch_em_Signal_VBF_HWW_M125.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_ZH.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = WminusHToTauTau_M125_13TeV_powheg_pythia8_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_WMinusH.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = WplusHToTauTau_M125_13TeV_powheg_pythia8_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_WPlusH.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_WGToLNuG.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_Signal.conf > analysisMacroSynch_em_Signal_ttH.conf

sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_EWKWMinus.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist = MC_PU2017_/' analysisMacroSynch_em_W.conf > analysisMacroSynch_em_EWKWPlus.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist =  MC_PU2017_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_EWKZ2Jets_ZToLL.conf
sed 's/SampleNameForPUHist =/SampleNameForPUHist =  MC_PU2017_/' analysisMacroSynch_em_DY.conf > analysisMacroSynch_em_EWKZ2Jets_ZToNuNu.conf



