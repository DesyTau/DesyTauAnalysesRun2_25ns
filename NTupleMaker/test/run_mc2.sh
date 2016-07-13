#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:29:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=6000M
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

channel=$2 
dir=$2



#dir=InvMuIso
region=InvElIso

##type MC or Data
type=MC

if [ ! -d Jobs ]
then
	mkdir Jobs
fi

if [ ! -d $3 ]
then
	mkdir $3
fi


cp *.conf Jobs/.

while read line
do

unset xsec
xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
#	cp $dir/${line} input$line
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

	echo " "$bas $xsec >> xsecs

	if [ -f Jobs/job${line}$channel$dir${bas}_B.sh ] ; then
rm Jobs/job${line}$channel$dir${bas}_B.sh
fi
#rm Jobs/job${line}$channel$dir${bas}_A.sh
#rm Jobs/job${line}$channel$dir${bas}_C_${region}.sh
#rm Jobs/job${line}$channel$dir${bas}_D_${region}.sh

#analysisMacroSUSY_${type}_B.conf  analysisMacroSUSY_${type}_C_InvMu.conf  analysisMacroSUSY_${type}_D_InvMu.conf

cat bss > Jobs/job${line}$channel$dir${bas}_B.sh
#echo SUSYTtemplate analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh

echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} $dir 1 $3 $4>> Jobs/job${line}$channel$dir${bas}_B.sh
#echo SUSY${channel} analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh
#echo taufakerateMu analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh

#cat bss > Jobs/job${line}$channel$dir${bas}_A.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_A.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_A.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_A.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_A.sh
#echo 	SUSYTtemplate analysisMacroSUSY_${type}_A.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_A.sh

#cat bss > Jobs/job${line}$channel$dir${bas}_C_${region}.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_C_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_C_${region}.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_C_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_C_${region}.sh
#echo 	SUSYTtemplate analysisMacroSUSY_${type}_C_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_C_${region}.sh

#cat bss > Jobs/job${line}$channel$dir${bas}_D_${region}.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_D_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_D_${region}.sh
#echo 	SUSY${channel} analysisMacroSUSY_${type}_D_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_D_${region}.sh
#echo 	SUSYTtemplate analysisMacroSUSY_${type}_D_${region}.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_D_${region}.sh


#echo $bas $xsec >> xsecs
echo $bas $xsec 

if [ ! -f $dir/${bas}_B_OS.root ] ;then

chmod u+x Jobs/job${line}$channel$dir${bas}_B.sh
qsub Jobs/job${line}$channel$dir${bas}_B.sh 
fi


#if [ ! -f $dir/${bas}_${region}__C_OS.root ] ;then
#chmod u+x Jobs/job${line}$channel$dir${bas}_C_${region}.sh
#qsub Jobs/job${line}$channel$dir${bas}_C_${region}.sh 
#fi

#if [ ! -f $dir/${bas}_A_SS.root ] ;then
#chmod u+x Jobs/job${line}$channel$dir${bas}_A.sh
#qsub Jobs/job${line}$channel$dir${bas}_A.sh 
#fi

#if [ ! -f $dir/${bas}_${region}__D_SS.root ] ;then
#chmod u+x Jobs/job${line}$channel$dir${bas}_D_${region}.sh
#qsub Jobs/job${line}$channel$dir${bas}_D_${region}.sh 
#fi

fi

done<$dir/${line}
rm input${line}
done<$1



