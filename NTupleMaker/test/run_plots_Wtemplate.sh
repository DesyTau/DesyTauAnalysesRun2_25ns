#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=0:15:00
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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_7_4_14/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;


dir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_7_4_14/src/DesyTauAnalyses/NTupleMaker/test"

while read line
do

unset file
file=`echo $line | cut -d '/' -f2`
	
mkdir dirW_$file
cd dirW_$file
echo ==============================================
pwd
echo ==============================================

cp $dir/analyzer_h .
cp $dir/analyzer_Wtemplate_C .


cp $dir/runme.C .
cp $dir/plots.h .


unset fileA
fileA=`echo $file | awk -F "_A_SS" '{print $1}'`
unset fileB
fileB=`echo $file | awk -F "_B_OS" '{print $1}'`

echo $line , $fileA , $fileB
if [[  $file == *"B_OS"* &&  ! -f $dir/plotsW/${fileB}.root ]]  ; then
cp analyzer_h analyzer.h
cp analyzer_Wtemplate_C analyzer.C

echo the filein for main region : $file , the fileout : $file
sed -i 's/FILEIN/'$file'/g' analyzer.h
sed -i 's/FILEIN/'$file'/g' analyzer.C

rm plots.root
root -l -q -b runme.C 

mv plots.root $dir/plotsW/${fileB}.root
fi


cd $dir
rm -fr dirW_$line
done<$1
