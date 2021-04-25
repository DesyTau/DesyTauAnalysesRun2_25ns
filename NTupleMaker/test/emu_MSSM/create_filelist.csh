#!/bin/csh 

python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_powheg_M${1} --query "/SUSYGluGluToBBHToTauTau_M${1}_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" --outputfile 2017/SUSYGluGluToBBHToTauTau_powheg_M${1}.list
python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_powheg_M${1} --query "/SUSYGluGluToBBHToTauTau_M${1}_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" --outputfile 2016/SUSYGluGluToBBHToTauTau_powheg_M${1}.list
#python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_M${1} --query "/SUSYGluGluToBBHToTauTau_M${1}_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v*/MINIAODSIM" --outputfile 2018/SUSYGluGluToBBHToTauTau_${1}.list
