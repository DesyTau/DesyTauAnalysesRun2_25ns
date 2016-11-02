#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_rt=0:45:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=3000M
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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test ;  eval `scramv1 runtime -sh` ;

wdir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test"
channel=$2 
dir=$2




##type MC or Data
type=MC



cp *.conf Jobs/.

while read line
do

unset xsec
xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
#xsec=1
echo FOUND XSEC for ${line} to be $xsec
unset f
while read f
	do

echo $f
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

echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} $dir 1 1 1 1>> Jobs/job${line}$channel$dir${bas}_B.sh
#echo SUSY${channel} analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh



#echo $bas $xsec >> xsecs
echo $bas $xsec 

if [ ! -f $dir/${bas}_B_OS.root ] ;then

chmod u+x $wdir/Jobs/job${line}$channel$dir${bas}_B.sh
qsub $wdir/Jobs/job${line}$channel$dir${bas}_B.sh 
fi


fi

done<$dir/${line}
done<$1



