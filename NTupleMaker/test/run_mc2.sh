#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
# use current dir and current environment
#$ -cwd
#$ -V
#


wdir="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test"
channel=$2 

##type MC or Data
type=MC

systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

systematics="Nominal"

syst=$3


export dir=${2}_${3}

if [[ -z $3  ]];then
systematics="Nominal"
fi

if [[ $syst == "Nominal" ]] || [[ -z $3  ]];then
export	dir=${channel}
fi



cp *.conf Jobs/.

while read line
do

ct=`ls ${dir}/${line}*.root | wc -l`
ctt=`cat ${dir}/${line} | wc -l`


echo There are  $ct out of $ctt for $line in $dir dir for systematic $syst

if [[ $ct -ge $ctt ]] ;then
	continue;
fi

unset xsec
xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
#xsec=1
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


	if [ -f Jobs/job${line}$channel$dir${bas}_B${syst}.sh ] ; then
rm Jobs/job${line}$channel$dir${bas}_B${syst}.sh
fi


cat bss > Jobs/job${line}$channel$dir${bas}_B${syst}.sh




if [ ! -f ${dir}/${bas}_B_OS.root ] ;then

echo $bas $xsec $dir
echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} ${channel} 1 $syst>> Jobs/job${line}$channel$dir${bas}_B${syst}.sh
chmod u+x $wdir/Jobs/job${line}$channel$dir${bas}_B${syst}.sh
 qsub -l h_rt=1:45:00 -l h_cpu=1500M $wdir/Jobs/job${line}$channel$dir${bas}_B${syst}.sh 
fi


fi

done<$dir/${line}
done<$1

