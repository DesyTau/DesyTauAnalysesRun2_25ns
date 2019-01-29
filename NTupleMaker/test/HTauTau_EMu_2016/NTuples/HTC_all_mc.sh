#!/bin/bash

#DY Samples
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DY1JetsToLL_M-50    5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DY2JetsToLL_M-50    5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DY3JetsToLL_M-50    5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DY4JetsToLL_M-50   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DYJetsToLL_M-50   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf  DYJetsToLL_M-10to50  5

#WJets
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf  W1JetsToLNu   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf  W2JetsToLNu   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf  W3JetsToLNu   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf  W4JetsToLNu   5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf   WJetsToLNu   5


#SingleTop 
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ST_t-channel_antitop  4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ST_t-channel_top  4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ST_tW_top   4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ST_tW_antitop   4

#TTbar
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Top.conf  TTbar   4

#Diboson 
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  VVTo2L2Nu  1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  WWToLNuQQ  1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  WZJToLLLNu   1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  WZTo1L1Nu2Q  1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  WZTo1L3Nu  1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  WZTo2L2Q   1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ZZTo2L2Q   1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_mc.conf  ZZTo4L    1


#Signal
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_ggh.conf  GluGluHToTauTau_M125  4
./HTC_submit_seq.sh SynchNTupleProducer_em analysisMacroSynch_em_Signal_VBF.conf VBFHToTauTau_M125 4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  WminusHToTauTau_M125 1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  WplusHToTauTau_M125 1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  ZHToTauTau_M125 1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  GluGluHToWWTo2L2Nu_M125 1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  VBFHToWWTo2L2Nu_M125 1
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Signal_VBF.conf  ttHJetToTT_M125 1

#WG and WGstar 
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_WG.conf  WGToLNuG   4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf  WGstarToLNuEE   4
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf   WGstarToLNuMuMu  4

#EWK
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf EWKWMinus2Jet 5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_W.conf EWKWPlus2Jets 5
./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_DY.conf EWKZ2Jets  5

./HTC_submit_seq.sh SynchNTupleProducer_em  analysisMacroSynch_em_Embedded.conf Embedding_Run2016  1
