#!/bin/bash


./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf TTToHadronic_TuneCP5_13TeV-powheg-pythia8 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 20

./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WW_TuneCP5_13TeV-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf WZ_TuneCP5_13TeV-pythia8 20
./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_MC.conf ZZ_TuneCP5_13TeV-pythia8 20


./HTC_submit_seq.sh EleTauFR analysisMacro_etauFR_Data.conf EGamma_Run2018 10
