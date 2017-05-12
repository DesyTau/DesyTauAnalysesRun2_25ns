
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
systematics="Nominal JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
#systematics="ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
#systematics="Nominal JetEnUp"
#systematics="JetEnDown JetEnUp UnclEnUp UnclEnDown"
systematics="Nominal"
#systematics="TopPtUp TopPtDown ZPtUp ZPtDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="JetEnDown JetEnUp"


#region="SR CRA CRB SR_CR"
region="SR_CR1"

syst="Nominal"
for rg in $region
do


		root -l -q -b 'Extract1DFakes.C("'$syst'", "'$rg'")' 
done
