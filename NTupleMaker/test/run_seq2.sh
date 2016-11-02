#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=0:29:00
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
##$ -o NTuple_$1.out
#
##$ -e NTuple_$1.err
#mkdir $1


cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test;cmsenv

channel=$2 
era=$3

era=Htt
era=25ns

channel=mutau

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
xsec=`grep "$line " xsecs | cut -d " " -f2-3`	
	cp $era/$line input$line
echo FOUND XSEC for $line to be $xsec
unset f
while read f
	do

echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f $era/$bas.root ] 
then
echo $f > $era/$bas
rm Jobs/job$line$channel$era$bas.sh
cat bss > Jobs/job$line$channel$era$bas.sh
#echo 	AnalysisMacro_dimuons  analysisMacro_dimuons.conf ${bas} Htt>> Jobs/job$line$channel$era$bas.sh
#echo 	AnalysisMacro_dimuonsMET  analysisMacro_dimuons.conf ${bas} $era>> Jobs/job$line$channel$era$bas.sh
#echo 	AnalysisMacro_dijet  analysisMacro_dimuons.conf ${bas} $era>> Jobs/job$line$channel$era$bas.sh

echo 	SUSY$channel analysisMacroSUSY_Data.conf ${bas} $era>> Jobs/job$line$channel$era$bas.sh

#echo 	SUSY$channel DataPRv1.conf ${bas} $era>> Jobs/job$line$channel$era$bas.sh
#echo 	SUSY$channel Data.conf ${bas} $era>> Jobs/job$line$channel$era$bas.sh

#echo 	SUSYQCD analysisMacroSUSY_Data.conf ${bas} 25ns>> Jobs/job$line$channel$era$bas.sh



echo " "$bas $xsec >> xsecs
echo $bas $xsec 
chmod u+x Jobs/job$line$channel$era$bas.sh
qsub Jobs/job$line$channel$era$bas.sh
fi

done<input$line
rm input$line
done<$1



