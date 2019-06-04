#!bin/bash

cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_W.conf
sed -i 's/IsW = false/IsW = true/g' analysisMacroSynch_em_W.conf
cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_DY.conf
sed -i 's/IsDY = false/IsDY = true/g' analysisMacroSynch_em_DY.conf
cp analysisMacroSynch_em_MC.conf analysisMacroSynch_em_Signal.conf
sed -i 's/IsSignal = false/IsSignal = true/g' analysisMacroSynch_em_Signal.conf