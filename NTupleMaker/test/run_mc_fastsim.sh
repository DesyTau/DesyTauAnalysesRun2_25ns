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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_7_4_14/src/DesyTauAnalyses/NTupleMaker/test ;  eval `scramv1 runtime -sh` ;

#find /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_4/src/FastSim/Output/*/*/ -type f -name "*.root" > filesMC

#ls FastSim_Stau > filesMC
ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_14/src/miniAOD/Output/StauScan > filesMC
#ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_4/src/FastSim/Output/f > filesMC
#cat filesMC | head -1  > t 
#mv t filesMC
while read line 
do

	fname=$(basename $line)
	echo the name is $fname
	if [[ ! -f $PWD/FastSim_Stau_Out/${fname} ]] ; then 
	cat bssFast > ${fname}_exec.sh

	echo cmsRun TreeProducerFromMiniAOD_74x_MC25ns_V6_${fname}.py >> ${fname}_exec.sh


	cp TreeProducerFromMiniAOD_74x_MC25ns_V6_template2.py TreeProducerFromMiniAOD_74x_MC25ns_V6_${fname}.py


	sed -i 's/FILEIN/'${fname}'/g' TreeProducerFromMiniAOD_74x_MC25ns_V6_${fname}.py
	sed -i 's/FILEOUT/'${fname}'/g' TreeProducerFromMiniAOD_74x_MC25ns_V6_${fname}.py
	mv TreeProducerFromMiniAOD_74x_MC25ns_V6_${fname}.py FastSim_Stau
	echo $fname
	qsub ${fname}_exec.sh

fi

done <filesMC
