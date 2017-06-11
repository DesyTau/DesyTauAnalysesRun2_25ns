#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:29:00
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

systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown TopPtUp TopPtDown ZPtUp ZPtDown BTagUp BTagDown"
systematics="Nominal JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown genMET ScalesDown ScalesUp PDFUp PDFDown BTagUp BTagDown METRecoilUp METRecoilDown"
#systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown ZPtUp ZPtDown"

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;

channel=$2

	if [[  -z "$2" ]] ;then

		echo you must provide a channel....
		return 1
	fi

while read line
do

if [[  -z "$3" || $3 == "Nominal" ]] ;then
systematics="Nominal"
fi

if [[  $3 == "list" ]] ;then
systematics="list"
fi

if [[  $3 == "listTT" ]] ;then
#systematics="Nominal JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown BTagUp BTagDown METRecoilUp METRecoilDown"
systematics="TopPtUp TopPtDown ZPtUp ZPtDown ScalesDown ScalesUp PDFUp PDFDown  METRecoilUp METRecoilDown"
#systematics="JetEnUp JetEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
fi

if [[  $3 == "listDY" ]] ;then
systematics="Nominal JetEnUp JetEnDown ZPtUp ZPtDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown METRecoilUp METRecoilDown"
#systematics="ZPtUp ZPtDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown METRecoilUp METRecoilDown"
systematics="ZPtUp ZPtDown MuEnUp MuEnDown UnclEnUp UnclEnDown JetEnUp JetEnDown ElEnUp ElEnDown"
fi



if [[  $3 == "listWJ" ]] ;then
systematics="Nominal JetEnUp JetEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown METRecoilUp METRecoilDown BTagUp BTagDown ZPtUp ZPtDown TopPtUp TopPtDown"
#systematics="TopPtUp TopPtDown ZPtUp ZPtDown"
#systematics="JetEnUp JetEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown BTagUp BTagDown"
systematics="BTagDown BTagUp"
fi

if [[  $3 == "top" ]] ;then
systematics="TopPtUp TopPtDown"
fi

if [[  $3 == "Tau" ]] ;then
systematics="TauEnUp TauEnDown"
fi
if [[  $3 == "El" ]] ;then
systematics="ElEnUp ElEnDown"
fi

if [[  $3 == "Jet" ]] ;then
systematics="JetEnUp JetEnDown"
fi

if [[  $3 == "new" ]] ;then
systematics="PDFUp PDFDown ScalesUp ScalesDown BTagUp BTagDown"
systematics="PDFUp PDFDown ScalesUp ScalesDown"
systematics="PDFUp PDFDown ScalesUp ScalesDown"
fi

if [[  $3 == "MET" ]] ;then
systematics="METRecoilUp METRecoilDown"
fi

lt=`echo $line | cut -d '/' -f2`


	echo $lt > list_$lt
	
	for syst in $systematics
	do
	

		echo  plots for channel $2 and syst $syst and $lt 
	 	qsub -N p$2 -l h_rt=1:30:00 -l h_cpu=4500M run_plots_new.sh list_$lt $2 $syst


	done
done<$1

