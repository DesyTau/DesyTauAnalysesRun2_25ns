#!/bin/bash

#DY Samples
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_DY.conf  DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_DY.conf  DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_DY.conf  DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_DY.conf  DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_DY.conf  DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 10

#WJets
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_W.conf  WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_W.conf  W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_W.conf  W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_W.conf  W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_W.conf  W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 10


#SingleTop 
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 10

#TTbar
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  TTToHadronic_TuneCP5_13TeV_powheg_pythia8 10

#Diboson 
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  WW_TuneCP5_13TeV-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  WZ_TuneCP5_13TeV-pythia8 10
./HTC_submit_seq.sh SynchNTupleProducer_em_allEras  analysisMacroSynch_em_MC.conf  ZZ_TuneCP5_13TeV-pythia8 10

