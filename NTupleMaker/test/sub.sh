#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:39:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=10000M
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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_7_4_14/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;

SUSYmutau analysisMacroSUSY.conf 25ns
