

region="SR_CR1"

systematics="METEn JetEn"
systematics="JetEn UnclEn MuEn ElEn BTag TFRJetEn TFRMuEn TFRTauEn"
systematics="TauEn"

processes="tt wj dyj ztt dib ww qcd sT ttx "
#processes="sT"

unset fullproc
unset systt

for syst in $systematics
do

systt=$syst

if [[ $syst == "JetEn" ]] ; then
	systt="JES"
fi

if [[ $syst == "METen" ]] ; then
	systt="MET"
fi

for proc in $processes 
do


if [[ $proc == "tt" ]] ; then
fullproc="TTJets"
	echo will use the $proc, $fullproc
fi

if [[ $proc == "wj" ]] ; then
fullproc="WJets"

if [[ $syst == "BTag" ]] ; then
	syst="misBTag"
fi


fi

if [[ $proc == "ztt" ]] ; then
fullproc="DY#rightarrow#tau#tau"

if [[ $syst == "BTag" ]] ; then
	syst="misBTag"
fi

fi

if [[ $proc == "dyj" ]] ; then
fullproc="DY#rightarrowll"

if [[ $syst == "BTag" ]] ; then
	syst="misBTag"
fi

fi


if [[ $proc == "ww" ]] ; then
fullproc="WW"

if [[ $syst == "BTag" ]] ; then
	syst="misBTag"
fi

fi

if [[ $proc == "dib" ]] ; then
fullproc="VV(X)"

if [[ $syst == "BTag" ]] ; then
	syst="misBTag"
fi

fi


if [[ $proc == "ttx" ]] ; then
fullproc="TTXJets"


fi

if [[ $proc == "sT" ]] ; then
fullproc="singleTop"
fi





	cp PlotSys_C PlotSys.C

#fileNom=Templ_met_MT2lester_DZeta01J1D_17_35invfb_mt_C1N2_SR_CR1.root
#fileSystUp=Templ_met_MT2lester_DZeta01J1D_17_35invfb_mt_C1N2_SR_CR1_${syst}Up.root
#fileSystDown=Templ_met_MT2lester_DZeta01J1D_17_35invfb_mt_C1N2_SR_CR1_${syst}Up.root


	sed  -i 's@SYSTEMATICS@'${systt}'@g'  PlotSys.C
	sed  -i 's@PROCESS@'${proc}'@g'  PlotSys.C
	sed  -i 's@SYSTEMATIC@'$syst'@g'  PlotSys.C
	sed  -i 's@FILEHERE@'$1'@g'  PlotSys.C
	sed  -i 's@FILEHERE@'$1'@g'  PlotSys.C
	sed  -i 's@REGION@'$region'@g'  PlotSys.C
	sed  -i 's@NAMEHERE@'${fullproc}'@g' PlotSys.C
	sed  -i 's@FILEHERE@'$1'@g' PlotSys.C

	root -l -q -b PlotSys.C
done
done


