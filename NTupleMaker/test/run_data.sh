#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=2:29:00
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
##$ -o NTuple_$1.out
#
##$ -e NTuple_$1.err
#mkdir $1

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test;cmsenv


channel=$2 
dir=$2


#dir=InvMuIso
region=InvElIso

#dir=Ttemplate
#dir=eltau
#channel=eltau
if [ ! -d Jobs ]
then
	mkdir Jobs
fi

if [ ! -d $3 ]
then
	mkdir $3
fi


cp *.conf 25ns/.
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


#echo $bas $xsec >> xsecs
echo $bas $xsec 
#25ns/SingleMuon_2015D_05Oct_0_A_SS_DataDriven.root  25ns/SingleMuon_2015D_05Oct_0_B_OS_DataDriven.root  25ns/SingleMuon_2015D_05Oct_0_InvMuIso__D_SS_DataDriven.root

if [ ! -f $dir/${bas}_B_OS_DataDriven.root ]
then
	echo $dir/${bas}_B_OS_DataDriven.root
echo $f > $dir/$bas
cat bss > Jobs/job$line$channel$dir${bas}_B.sh
#echo SUSYTtemplate analysisMacroSUSY_Data_B.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_B.sh
echo SUSY$channel analysisMacroSUSY_Data_B.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_B.sh
#echo SUSYeltau analysisMacroSUSY_Data_B.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_B.sh
#echo taufakerateMu analysisMacroSUSY_Data_B.conf $bas $dir>> Jobs/job$line$channel$dir${bas}_B.sh
#echo taufakerateJets analysisMacroSUSY_DataJets_TauFakeRate.conf $bas $dir>> Jobs/job$line$channel$dir${bas}_B.sh

chmod u+x Jobs/job$line$channel$dir${bas}_B.sh
qsub Jobs/job$line$channel$dir${bas}_B.sh 
fi



if [ ! -f $dir/${bas}_A_SS_DataDriven.root ]
then
echo $f > $dir/$bas
cat bss > Jobs/job$line$channel$dir${bas}_A.sh
echo 	SUSY${channel} analysisMacroSUSY_Data_A.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_A.sh
#echo 	SUSYeltau analysisMacroSUSY_Data_A.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_A.sh
#echo 	SUSYTtemplate analysisMacroSUSY_Data_A.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_A.sh
chmod u+x Jobs/job$line$channel$dir${bas}_A.sh
#qsub Jobs/job$line$channel$dir${bas}_A.sh 
fi

if [ ! -f $dir/${bas}_${region}__C_OS_DataDriven.root ]
then
echo $f > $dir/$bas
cat bss > Jobs/job$line$channel$dir${bas}_C_${region}.sh
echo 	SUSY${channel} analysisMacroSUSY_Data_C_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_C_${region}.sh
#echo 	SUSYeltau analysisMacroSUSY_Data_C_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_C_${region}.sh
#echo 	SUSYTtemplate analysisMacroSUSY_Data_C_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_C_${region}.sh
chmod u+x Jobs/job$line$channel$dir${bas}_C_${region}.sh
#qsub Jobs/job$line$channel$dir${bas}_C_${region}.sh 
fi


if [ ! -f $dir/${bas}_${region}__D_SS_DataDriven.root ]
then
echo $f > $dir/$bas
cat bss > Jobs/job$line$channel$dir${bas}_D_${region}.sh
echo 	SUSY${channel} analysisMacroSUSY_Data_D_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_D_${region}.sh
#echo 	SUSYeltau analysisMacroSUSY_Data_D_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_D_${region}.sh
#echo 	SUSYTtemplate analysisMacroSUSY_Data_D_${region}.conf ${bas} $dir>> Jobs/job$line$channel$dir${bas}_D_${region}.sh
chmod u+x Jobs/job$line$channel$dir${bas}_D_${region}.sh
#qsub Jobs/job$line$channel$dir${bas}_D_${region}.sh 
fi


done<$dir/$line
rm input$line
done<$1



