
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
systematics="Nominal JetEnUp JetEnDown"
#systematics="Nominal"
for syst in $systematics
do

		root -l -q -b 'Overlap1DMod.C("'$syst'")' 
done
