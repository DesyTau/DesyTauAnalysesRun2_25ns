dir="/nfs/dust/cms/group/higgs-kit/80x_MC_SUSY_v1/"

sources="TChipmStauSnu   TChiStauStau   TStauStau"

sources="TStauStau"
#sources="TChiStauStau"
#sources="TChipmStauSnu"



for model in $sources
do

unset out

if [ $model == "TStauStau" ] ;
then


out="stau-stau"
fi

if [ $model == "TChiStauStau" ] ;
then


out="C1N2"
fi

if [ $model == "TChipmStauSnu" ] ;
then


out="C1C1"
fi

COUNTER=0

while read line
do


	susy=`echo $line  | awk -F " " '{print $1}'`
	lsp=`echo $line  | awk -F " " '{print $2}'`


if [ $lsp == "0" ];then

	lsp=1
fi
outfile=${model}_${susy}_LSP${lsp}_${COUNTER}.root
if [ ! -f ${out}_${susy}_LSP${lsp}.root ] ; then



echo "hadd -f -k ${out}_${susy}_LSP${lsp}.root ${model}_${susy}_LSP${lsp}_*.root"


fi
let COUNTER=COUNTER+1 


done<points_$model

done
