#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:59:00
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


cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;

dir=MET_filters_new

filters="csc2015 hbheiso hbher1nozeros hbher1 hbher2l hbher2t"

filters="nofilters csc2015 hbheiso hbher1nozeros hbher1 hbher2l hbher2t custom"

filters="nofilters SingleMuon_badMuonTrack  SingleMuon_badTrack SingleMuon_ecalscn1043093 SingleMuon_HBHE_CSC_BADTRK_BADMUON"

filters="SingleMuon_nofilters SingleMuon_HBHE_CSC_BADTRK_BADMUON"

filters="MET_nofilters MET_badMuonTrack  MET_badTrack MET_ecalscn1043093 MET_HBHE_CSC_BADTRK_BADMUON"

filters="csc2015_Dec01"

#filters="MET_badMuonTrack  MET_badTrack"


#filters="allfilters nofilters"

#filters="csc2015"


#filters="nofilters csc2015 hbheiso hbher1nozeros hbher1 hbher2l hbher2t"
##filters="all badResolutionTrack_Jan13 csc2015_Dec01  ecalscn1043093_Dec01   muonBadTrack_Jan13  nofilters"

#filters to be run after applied the noise
#filters="all_and_noise badResolutionTrack_Jan13 csc2015_Dec01  ecalscn1043093_Dec01   muonBadTrack_Jan13 noise"


filters="all_and_noise noise"

### filters to be run on the no-noise applied ie 1-by-1
#filters="badResolutionTrack_Jan13 csc2015_Dec01  ecalscn1043093_Dec01   muonBadTrack_Jan13 nofilters"


#filters="noise"

if [ ! -d Jobs ]
then
	mkdir Jobs
fi

if [ ! -d $3 ]
then
	mkdir $3
fi


while read line
do

unset xsec
xsec=`grep "$line " xsecs | cut -d " " -f2-3`	
	cp $dir/$line input$line
xsec=1;
	
echo FOUND XSEC for $line to be $xsec
unset f
while read f
	do




for filt in $filters
do

if [ ! -d ${dir}/$filt ] ; then
mkdir ${dir}/$filt
fi

echo $f and filter $filt
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f ${dir}/$filt/$bas.root ] 
then
echo $f > ${dir}/$filt/$bas
rm Jobs/job$line$channel$filt$bas.sh

cat bss > Jobs/job$line$channel$filt${bas}_dijet.sh
echo 	AnalysisMacro_dijet  analysisMacro_dimuons.conf $bas ${dir}/$filt $filt >> Jobs/job$line$channel$filt${bas}_dijet.sh

#cat bss > Jobs/job$line$channel$filt${bas}_dimuon.sh
#echo 	AnalysisMacro_dimuonsMET analysisMacro_dimuons.conf $bas ${dir}/$filt $filt >> Jobs/job$line$channel$filt${bas}_dimuon.sh



echo " "$bas $xsec >> xsecs
echo $bas $xsec 
if [ -f Jobs/job$line$channel$filt${bas}_dimuon.sh ] ;then

chmod u+x Jobs/job$line$channel$filt${bas}_dimuon.sh
#qsub Jobs/job$line$channel$filt${bas}_dimuon.sh
fi

if [ -f Jobs/job$line$channel$filt${bas}_dijet.sh ] ;then

chmod u+x Jobs/job$line$channel$filt${bas}_dijet.sh
qsub Jobs/job$line$channel$filt${bas}_dijet.sh

fi

fi

done

done<input$line
rm input$line
done<$1



