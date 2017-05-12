
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
systematics="Nominal JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
#systematics="ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
#systematics="Nominal JetEnUp"
#systematics="JetEnDown JetEnUp UnclEnUp UnclEnDown"
#systematics="Nominal"
#systematics="TopPtUp TopPtDown ZPtUp ZPtDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="TauEnDown"


#region="SR CRA CRB SR_CR"
region="SR_CR CRB CRC SR_CR1"
region="CRB"

for rg in $region
do

for syst in $systematics
do

		root -l -q -b 'Overlap1DMod.C("'$syst'", "'$rg'")' 
done
done
