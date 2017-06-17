
channel=$1

rg=$2
unset systematics
systematics="JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown BTagUp BTagDown METRecoilUp METRecoilDown TFRJetEnUp TFRJetEnDown TFRMuEnUp TFRMuEnDown TFRTauEnUp TFRTauEnDown"
#systematics="JetEnUp JetEnDown ZPtUp ZPtDown TopPtUp TopPtDown TauEnUp TauEnDown MuEnUp MuEnDown ElEnUp ElEnDown UnclEnUp UnclEnDown"
systematics="MuEnUp MuEnDown"
#systematics="Nominal"
unset syst

for syst in $systematics
do
	echo will set up now $syst

	cp QCDsamples QCDsamples${syst}
	cp OverlapQCD_C OverlapQCD${syst}.C

	sed -i 's/_A /_'${syst}'_A /g' QCDsamples${syst}
	sed -i 's/_C /_'${syst}'_C /g' QCDsamples${syst}
	sed -i 's/_D /_'${syst}'_D /g' QCDsamples${syst}
	sed -i 's/35invfb_'${syst}'_A/35invfb_A/g' QCDsamples${syst}
	sed -i 's/35invfb_'${syst}'_C/35invfb_C/g' QCDsamples${syst}
	sed -i 's/35invfb_'${syst}'_D/35invfb_D/g' QCDsamples${syst}
	sed -i 's/SYST/'${syst}'/g' OverlapQCD${syst}.C
	sed -i 's/CHANNELHERE/'${channel}'/g' OverlapQCD${syst}.C

	cat bss > run_QCD${syst}.sh
	echo root -l -q -b OverlapQCD${syst}.C >> run_QCD${syst}.sh
	#echo mv QCD_${syst}_B.root QCD_DataDriven_${rg}_${syst}_B.root >> run_QCD${syst}.sh
	echo mv QCD_${syst}_B.root QCD_DataDriven_${syst}_B.root >> run_QCD${syst}.sh
	qsub run_QCD${syst}.sh
done
