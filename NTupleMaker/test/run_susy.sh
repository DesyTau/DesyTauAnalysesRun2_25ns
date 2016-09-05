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


cp *.conf Jobs/.

while read line
do


while read line2
do

	unset mass
	unset lsp
mass=`echo $line2 | awk -F " " '{print $1}'`
lsp=`echo $line2 | awk -F " " '{print $2}'`

unset f
while read f
	do

echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f $dir/$bas.root ] 
then
echo $f > $dir/$bas


	if [ -f Jobs/job${line}$channel$dir${bas}_B.sh ] ; then
rm Jobs/job${line}$channel$dir${bas}_B.sh
fi



cat bss > Jobs/job${line}$channel$dir${bas}_B.sh
#echo SUSYTtemplate analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh


#echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} $dir 1 $mass $lsp

echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} $dir 1 $mass $lsp>> Jobs/job${line}$channel$dir${bas}_B.sh


if [ ! -f $dir/${bas}_B_OS_${mass}_${lsp}.root ] ;then

chmod u+x Jobs/job${line}$channel$dir${bas}_B.sh
qsub Jobs/job${line}$channel$dir${bas}_B.sh 
fi


fi


done<$dir/${line}
done<$3
done<$1



