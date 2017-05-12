systematics="JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"

regions="SR CRA CRB SR_CR"


for syst in $systematics
do
	for rg in $regions
	do



 mv Templ_met_MT2lester_DZeta01J1D_17_35invfb_mt_${1}_${syst}_${rg}.root	Templ_met_MT2lester_DZeta01J1D_17_35invfb_mt_${1}_${rg}_${syst}.root	
 mv Templ_met_MT2lester_DZeta01J1D_18_35invfb_mt_${1}_${syst}_${rg}.root	Templ_met_MT2lester_DZeta01J1D_18_35invfb_mt_${1}_${rg}_${syst}.root	

done
done
