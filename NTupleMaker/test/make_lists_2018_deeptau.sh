#!/bin/sh
dirData=/nfs/dust/cms/group/higgs-kit/2018_v2/data
dirMC=/nfs/dust/cms/group/higgs-kit/2018_v2/mc/RunIIAutumn18MiniAOD


ls $dirData/SingleMuon/SingleMuon_Run2018A-17Sep2018-v2/*root > SingleMuon_Run2018A
ls $dirData/SingleMuon/SingleMuon_Run2018B-17Sep2018-v1/*root > SingleMuon_Run2018B
ls $dirData/SingleMuon/SingleMuon_Run2018C-17Sep2018-v1/*root > SingleMuon_Run2018C
ls $dirData/SingleMuon/SingleMuon_Run2018D-22Jan2019-v2/*root > SingleMuon_Run2018D

cat SingleMuon_Run2018A SingleMuon_Run2018B SingleMuon_Run2018C SingleMuon_Run2018D > SingleMuon_Run2018

ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8


ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8


ls $dirMC/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8
ls $dirMC/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/*root > WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8
ls $dirMC/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/*root > WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8
ls $dirMC/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/*root > ZZTo4L_TuneCP5_13TeV_powheg_pythia8
ls $dirMC/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8
ls $dirMC/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8


ls $dirMC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
ls $dirMC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic_TuneCP5_13TeV-powheg-pythia8



ls $dirMC/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ls $dirMC/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ls $dirMC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8