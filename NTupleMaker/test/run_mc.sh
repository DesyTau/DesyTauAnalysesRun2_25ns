#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_rt=2:45:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=2000M
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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test;  eval `scramv1 runtime -sh` ;

wdir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test"
channel=$2 

##type MC or Data
type=MC




dir=${2}


cp *.conf Jobs/.

while read line
do

unset xsec
xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
#xsec=1
#echo FOUND XSEC for ${line} to be $xsec
unset f
while read f
	do

#echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f $dir/$bas.root ] 
then
echo $f > $dir/$bas

#	echo " "$bas $xsec >> xsecs


	if [ -f Jobs/job${line}$channel$dir${bas}_B.sh ] ; then
rm Jobs/job${line}$channel$dir${bas}_B.sh
fi


cat bss > Jobs/job${line}$channel$dir${bas}_B.sh




if [ ! -f ${dir}/${bas}_B_OS.root ] ;then

echo $bas $xsec $dir
echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} ${channel} 1 Nominal>> Jobs/job${line}$channel$dir${bas}_B.sh
chmod u+x $wdir/Jobs/job${line}$channel$dir${bas}_B.sh
 qsub -l h_rt=0:45:00 -l h_cpu=1500M $wdir/Jobs/job${line}$channel$dir${bas}_B.sh 
fi


fi

done<$dir/${line}
done<$1

