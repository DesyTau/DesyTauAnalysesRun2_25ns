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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test ;  eval `scramv1 runtime -sh` ;

#ls FastSim_Stau > filesMC
\ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_8_0_21/src/FastSim/Output_mAOD/files/$1*.root > filesMC

#cat filesMC | head -1  > t 
#mv t filesMC
while read line 
do

	fname=$(basename $line)
	echo the name is $fname
	if [[ ! -f $PWD/FastSim_Stau_Out/${fname} ]]  && [[ ! -f $PWD/Skims/${fname} ]]; then 
	cat bssFast > ${fname}_exec.sh

	echo cmsRun FastSim_Stau/TreeProducerFromMiniAOD_80x_MC25ns_signal_${fname}.py >> ${fname}_exec.sh


	cp TreeProducerFromMiniAOD_80x_MC25ns_FastSim_noHLT_template.py TreeProducerFromMiniAOD_80x_MC25ns_signal_${fname}.py


	sed -i 's/FILEIN/'${fname}'/g' TreeProducerFromMiniAOD_80x_MC25ns_signal_${fname}.py
	sed -i 's/FILEOUT/'${fname}'/g' TreeProducerFromMiniAOD_80x_MC25ns_signal_${fname}.py
	mv TreeProducerFromMiniAOD_80x_MC25ns_signal_${fname}.py FastSim_Stau/.
	echo $fname
	qsub ${fname}_exec.sh

fi

done <filesMC
