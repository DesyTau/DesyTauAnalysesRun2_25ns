#!/bin/bash

# Data
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018A mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018B mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018C mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018D mt 3

# Signals
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf GluGluHToTauTau_M125 mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf VBFHToTauTau_M125 mt 5

# DY
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DYJetsToLL_M-50 mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DY1JetsToLL_M-50 mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DY2JetsToLL_M-50 mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DY3JetsToLL_M-50 mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DY4JetsToLL_M-50 mt 3
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf DYJetsToLL_M-10to50 mt 3

# W+jets
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf WJetsToLNu mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf W1JetsToLNu mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf W2JetsToLNu mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf W3JetsToLNu mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf W4JetsToLNu mt 5

# VV
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf WW mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf WZ mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf ZZ mt 5

# TT
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf TTTo2L2Nu mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf TTToHadronic mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf TTToSemiLeptonic mt 5

# Single Top
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf ST_t-channel_antitop_4f mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf ST_t-channel_top_4f mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf ST_tW_antitop_5f mt 5
./condorsub_seq_leptau.sh SynchNTupleProducer_Run2 analysisMacroSynch_mt_18_MC.conf ST_tW_top_5f mt 5
