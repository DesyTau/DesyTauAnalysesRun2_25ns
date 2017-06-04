#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=0:35:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=1500M
#
#(use hh site)
#$ -l site=hh
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL5
#$ -l os=sld6
#
# use current dir and current environment
#$ -cwd
#$ -V
#

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/plots_Wtemplate;eval `scramv1 runtime -sh` ;

root -l -q -b OverlapQCDNominal.C
mv QCD_Nominal_B.root QCD_DataDriven_Nominal_B.root
