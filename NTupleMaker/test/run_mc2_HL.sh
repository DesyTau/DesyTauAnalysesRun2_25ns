#!/bin/sh
#


################3 CHANGE THIS TO YOUR WORKING DIR
wdir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test"

cd $wdir/test;eval `scramv1 runtime -sh` ;


channel=$2 

##type MC or Data
type=MC

systematics="$3"

if [[ $3 == "all" || $3 == "All"  || $3 == "list" ]];then
#systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
#systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
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

if [[ $line != *"C1"*  && $line != *"stau"*  && $line != *"Chi"* ]] ; then
ct=`ls ${dir}/${line}_[0-9]*_B_OS.root | wc -l`

else
ct=`ls ${dir}/${line}_B_OS.root | wc -l`

fi


ctt=`cat ${dir}/${line} | wc -l`


echo There are  $ct out of $ctt for $line in ${dir} dir for systematic $syst

if [[ $ct -ge $ctt ]] ;then
	continue;
fi

#unset xsec
xsec=1
unset f
while read f
	do

#echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f ${dir}/$bas.root ] 
then
echo $f > ${dir}/$bas

#	echo " "$bas $xsec >> xsecs


if [ -f Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh ] ; then
rm Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh
fi


cat bss > Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh




if [ ! -f ${dir}/${bas}_B_OS.root ] ;then

echo $bas $xsec ${dir}


	echo SUSY$channel analysisMacroSUSY_MC_B_HL.conf ${bas} ${channel} 1 $syst>> Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh

 qsub -l h_rt=0:30:00 -l h_cpu=1500M -e /dev/null -o /dev/null $wdir/Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh 
 #. $wdir/Jobs/job${channel}_${line}_${dir}_${bas}_B${syst}.sh 
fi


fi

done<${dir}/${line}
done<$1
done
