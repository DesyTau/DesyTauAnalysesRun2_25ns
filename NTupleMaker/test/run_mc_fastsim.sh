#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=10:29:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=5000M
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

wdir=`pwd`

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test ;  eval `scramv1 runtime -sh` ;

#find /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_4/src/FastSim/Output/*/*/ -type f -name "*.root" > filesMC

#ls FastSim_Stau > filesMC
#ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_8_0_12/src/miniAOD/Output/StauScanv2/*stau1* > filesMC
#ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_8_0_12/src/miniAOD/Output/StauScanv2/*stau2* >> filesMC
ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_6_3/src/miniAOD/Output/Scan/$1*.root > filesMC
#ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_4/src/FastSim/Output/f > filesMC
#cat filesMC | head -1  > t 
#mv t filesMC
while read line 
do

	fname=$(basename $line)
	echo the name is $fname
	if [[ ! -f $PWD/FastSim_Stau_Out/${fname} ]]  && [[ ! -f $PWD/FastSim_Stau_Out/RootFiles/${fname} ]]; then 
	cat bssFast > ${fname}_exec.sh

	echo cmsRun TreeProducerFromMiniAOD_76x_MC25ns_signal_${fname}.py >> ${fname}_exec.sh


	cp TreeProducerFromMiniAOD_76x_MC25ns_signal_template.py TreeProducerFromMiniAOD_76x_MC25ns_signal_${fname}.py


	sed -i 's/FILEIN/'${fname}'/g' TreeProducerFromMiniAOD_76x_MC25ns_signal_${fname}.py
	sed -i 's/FILEOUT/'${fname}'/g' TreeProducerFromMiniAOD_76x_MC25ns_signal_${fname}.py
	mv TreeProducerFromMiniAOD_76x_MC25ns_signal_${fname}.py FastSim_Stau
	echo $fname
	qsub ${fname}_exec.sh

fi

done <filesMC
