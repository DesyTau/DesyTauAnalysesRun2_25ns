#!/bin/bash

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DYJetsToLL_M-50 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY1JetsToLL_M-50 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY2JetsToLL_M-50 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY3JetsToLL_M-50 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY4JetsToLL_M-50 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTTo2L2Nu 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTToSemiLeptonic 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTToHadronic 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WJetsToLNu 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W1JetsToLNu 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W2JetsToLNu 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W3JetsToLNu 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W4JetsToLNu 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_t-channel_antitop_4f 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_t-channel_top_4f 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_tW_antitop_5f 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_tW_top_5f 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WW 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WZ 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ZZ 20


./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_Data.conf SingleEle_Run2016 10
