
systematics="Nominal JetEnUp JetEnDown ZPtUp ZPtDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown METRecoilUp METRecoilDown"
systematics="Nominal"

for syst in $systematics
do

	                root -l -q -b 'Overlap.C("'$syst'")'
			
		done

. ren.sh
for syst in $systematics
do
	cp cutflow_py cutflow.py
	if [[ $syst == "Nominal" ]] ; then
		syst=""
	sed -i 's/SYSTEMATICS//g' cutflow.py
	fi

	if [[ $syst != "Nominal" ]] ; then
	sed -i 's/SYSTEMATICS/_'$syst'/g' cutflow.py
	fi

	python cutflow.py > LabelList$syst
#	sed -i '4d' LabelList$syst
	sed -i '16d' LabelList$syst
	sed -i '4d' LabelList$syst
	sed -i '1,2d' LabelList$syst
	sed -i 's/JetsCuts/SRCuts/g'  LabelList$syst
	sed -i 's/MET>/\#it{p}_{T}^{miss}>/g'  LabelList$syst

done
