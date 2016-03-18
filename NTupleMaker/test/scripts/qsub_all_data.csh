#!/bin/csh
./qsub_seq.csh AnalysisMacro_dimuons analysisMacro_dimuons.conf /nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v2 SingleMuon_2015D_PRv4 10 99
./qsub_seq.csh AnalysisMacro_dimuons analysisMacro_dimuons.conf /nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v2 SingleMuon_2015D_05Oct 10 99
./qsub_seq_x.csh AnalysisMacro_dimuons analysisMacro_dimuons.conf /nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v2 SingleMuon_2015D_PRv4 0 9
./qsub_seq_x.csh AnalysisMacro_dimuons analysisMacro_dimuons.conf /nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v2 SingleMuon_2015D_05Oct 0 9
