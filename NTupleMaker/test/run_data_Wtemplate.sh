#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=2:29:00
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
##$ -o NTuple_$1.out
#
##$ -e NTuple_$1.err
#mkdir $1

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test;cmsenv



channel=mutau
dir=Wtemplate

region=InvMET
#region=InvMuIso

if [ ! -d Jobs ]
then
	mkdir Jobs
fi

if [ ! -d $3 ]
then
	mkdir $3
fi


cp *.conf 25ns/.
cp *.conf Htt/.
cp *.conf Jobs/.

while read line
do

unset xsec
#xsec=`grep " $line " xsecs | cut -d " " -f2-3`	
#	cp $dir/$line input$line
xsec=1
echo FOUND XSEC for $line to be $xsec
unset f
while read f
	do

echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`


echo $bas $xsec 

if [ ! -f $dir/${bas}_B_OS_DataDriven.root ]
then
	echo $dir/${bas}_B_OS_DataDriven.root
echo $f > $dir/$bas
cat bss > Jobs/job$line$channel$dir${bas}_B.sh
echo SUSYWtemplate analysisMacroSUSY_Data_B.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_B.sh

chmod u+x Jobs/job$line$channel$dir${bas}_B.sh
qsub Jobs/job$line$channel$dir${bas}_B.sh 
fi






done<$dir/$line
rm input$line
done<$1



