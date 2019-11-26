#!/bin/bash

# Data
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016B mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016C mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016D mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016E mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016F mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016G mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_data.conf SingleMuon_Run2016H mt 10

# Signals
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf GluGluHToTauTau_M125 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf VBFHToTauTau_M125 mt 10

# DY
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DYJetsToLL_M-50 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DY1JetsToLL_M-50 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DY2JetsToLL_M-50 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DY3JetsToLL_M-50 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DY4JetsToLL_M-50 mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf DYJetsToLL_M-10to50 mt 10

# W+jets
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf WJetsToLNu mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf W1JetsToLNu mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf W2JetsToLNu mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf W3JetsToLNu mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf W4JetsToLNu mt 10

# VV
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf WW mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf WZ mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf ZZ mt 10

# TT
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf TT mt 10

# Single Top
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf ST_t-channel_antitop_4f mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf ST_t-channel_top_4f mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf ST_tW_antitop_5f mt 10
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_16_MC.conf ST_tW_top_5f mt 10
