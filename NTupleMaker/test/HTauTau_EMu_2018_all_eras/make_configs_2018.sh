#!bin/bash

cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_W.conf
sed -i 's/IsW = false/IsW = true/g' analysisMacroSynch_em_W.conf
cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_DY.conf
sed -i 's/IsDY = false/IsDY = true/g' analysisMacroSynch_em_DY.conf
cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_Signal.conf
sed -i 's/IsSignal = false/IsSignal = true/g' analysisMacroSynch_em_Signal.conf
cp analysisMacroSynch_em_Signal.conf analysisMacroSynch_em_Signal_ggh.conf
sed -i 's/ApplygghReweighting = false/ApplygghReweighting = true/g' analysisMacroSynch_em_Signal_ggh.conf
sed -i 's/ApplygghUncertainties = false/ApplygghUncertainties = true/g' analysisMacroSynch_em_Signal_ggh.conf
cp analysisMacroSynch_em_Signal.conf analysisMacroSynch_em_Signal_vbf.conf
sed -i 's/ApplyVBFUncertainties = false/ApplyVBFUncertainties = true/g' analysisMacroSynch_em_Signal_vbf.conf
cp analysisMacroSynch_em_DATA.conf analysisMacroSynch_em_Embedded.conf
sed -i 's/IsEmbedded = false/IsEmbedded = true/g' analysisMacroSynch_em_Embedded.conf
sed -i 's/IsData = true/IsData = false/g' analysisMacroSynch_em_Embedded.conf