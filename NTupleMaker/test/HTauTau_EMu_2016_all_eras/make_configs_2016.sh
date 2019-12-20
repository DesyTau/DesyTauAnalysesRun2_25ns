#!bin/bash

cp analysisMacroSynch_em_mc.conf analysisMacroSynch_em_W.conf
sed -i 's/IsW = false/IsW = true/g' analysisMacroSynch_em_W.conf

cp analysisMacroSynch_em_mc.conf analysisMacroSynch_em_DY.conf
sed -i 's/IsDY = false/IsDY = true/g' analysisMacroSynch_em_DY.conf

cp analysisMacroSynch_em_mc.conf analysisMacroSynch_em_Signal_VBF.conf
cp analysisMacroSynch_em_mc.conf analysisMacroSynch_em_Signal.conf
sed -i 's/IsSignal = false/IsSignal = true/g' analysisMacroSynch_em_Signal_VBF.conf
sed -i 's/IsSignal = false/IsSignal = true/g' analysisMacroSynch_em_Signal.conf

cp analysisMacroSynch_em_Signal_VBF.conf analysisMacroSynch_em_Signal_ggh.conf
sed -i 's/ApplygghReweighting = false/ApplygghReweighting = true/g' analysisMacroSynch_em_Signal_ggh.conf
sed -i 's/ApplygghUncertainties = false/ApplygghUncertainties = true/g' analysisMacroSynch_em_Signal_ggh.conf
sed -i 's/ApplyVBFUncertainties = false/ApplyVBFUncertainties = true/g' analysisMacroSynch_em_Signal_VBF.conf

cp analysisMacroSynch_em_DATA.conf analysisMacroSynch_em_Embedded.conf
sed -i 's/IsEmbedded = false/IsEmbedded = true/g' analysisMacroSynch_em_Embedded.conf
sed -i 's/IsData = true/IsData = false/g' analysisMacroSynch_em_Embedded.conf

cp analysisMacroSynch_em_W.conf analysisMacroSynch_em_WG.conf
sed -i 's/RemoveGammaStar = false/RemoveGammaStar = true/g' analysisMacroSynch_em_WG.conf

cp analysisMacroSynch_em_DATA.conf analysisMacroSynch_em_GH.conf
sed -i 's/ApplyICHEPMuonId = true/ApplyICHEPMuonId = false/g' analysisMacroSynch_em_GH.conf
sed -i 's/ApplyDzFilterMatch = false/ApplyDzFilterMatch = true/g' analysisMacroSynch_em_GH.conf

