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
	
mkdir dir_$file
cd dir_$file
echo ==============================================
pwd
echo ==============================================

if [[ $line == *"stau"* ]] ; then
lsp=`echo $line | awk -F "_LSP" '{print $2}' | cut -d '_' -f1`
else
lsp=0

fi


cp $dir/analyzer_top_h .

cp $dir/analyzer_top_C .


sed -i 's/CHIMASSS/'$lsp'/g' analyzer*C

cp $dir/runme.C .
cp $dir/plots.h .


unset fileA
fileA=`echo $file | awk -F "_A_SS" '{print $1}'`
unset fileB
fileB=`echo $file | awk -F "_B" '{print $1}'`

echo $line , $fileA , $fileB
if [[  $file == *"_B"* &&  ! -f $dir/plotsT/${fileB}_B.root ]]  ; then
cp analyzer_top_h analyzer.h
cp analyzer_top_C analyzer.C

echo the filein for main region : $file , the fileout : $file
sed -i 's/FILEIN/'$file'/g' analyzer.h
sed -i 's/FILEIN/'$file'/g' analyzer.C

rm plots.root
root -l -q -b runme.C 

mv plots.root $dir/plotsT/${fileB}_B.root

fi



cd $dir
#rm -fr dir_$line
done<$1
