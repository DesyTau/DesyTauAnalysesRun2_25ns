
systematics="JetEnUp JetEnDown TopPtUp TopPtDown ZPtUp ZPtDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown ScalesDown ScalesUp PDFUp PDFDown BTagUp BTagDown METRecoilUp METRecoilDown"
while read line
do
	for syst in $systematics
	do
		file=${line}_B.root
		cp $file ${line}_${syst}_B.root
	done
done<qcd
