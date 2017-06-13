#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:45:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=3500M
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

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;


wdir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test"
channel=$2 

##type MC or Data
type=MC
#systematics="Nominal"
#systematics="Nominal"

systematics="$3"

if [[ $3 == "all" || $3 == "All"  || $3 == "list" ]];then
#systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
#systematics="Nominal JetEnUp JetEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
#systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown"

fi

if [[ -z $3  ]];then
systematics="Nominal"
syst="1"
fi

cp *.conf Jobs/.

for syst in $systematics
do


export dir=${2}_${syst}


if [[ $syst == "Nominal" ]] || [[ -z $3  ]];then
export	dir=${channel}
fi


while read line
do

if [[ $line != *"C1"* ]] ; then
ct=`ls ${dir}/${line}_[0-9]*_B_OS.root | wc -l`

else
ct=`ls ${dir}/${line}_B_OS.root | wc -l`

fi


ctt=`cat ${dir}/${line} | wc -l`


echo There are  $ct out of $ctt for $line in $dir dir for systematic $syst

if [[ $ct -ge $ctt ]] ;then
	continue;
fi

#unset xsec
#xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
xsec=1
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


	if [ -f Jobs/job${channel}${line}$dir${bas}_B${syst}.sh ] ; then
rm Jobs/job${channel}${line}$dir${bas}_B${syst}.sh
fi


cat bss > Jobs/job${channel}${line}$dir${bas}_B${syst}.sh




if [ ! -f ${dir}/${bas}_B_OS.root ] ;then

echo $bas $xsec $dir
echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} ${channel} 1 $syst>> Jobs/job${channel}${line}$dir${bas}_B${syst}.sh
chmod u+x $wdir/Jobs/job${channel}${line}$dir${bas}_B${syst}.sh
 qsub -l h_rt=1:45:00 -l h_cpu=3500M -e /dev/null -o /dev/null $wdir/Jobs/job${channel}${line}$dir${bas}_B${syst}.sh 
fi


fi

done<$dir/${line}
done<$1
done
