rm m
dir="/nfs/dust/cms/group/higgs-kit/80x_MC_SUSY_v1/"

sources="TChipmStauSnu   TChiStauStau   TStauStau"

sources="TChiStauStau"
sources="TChipmStauSnu"
#sources="TStauStau"
sources="TChipmWW"
sources="TStauStau_left"
sources="TChiStauStaux0p05"


for model in $sources
do

unset out

if [[ $model == "TStauStau"  ]] ;
then

out="stau-stau"
fi

if [[ $model == "TStauStau_left" ]] ;then
out="stau-stau_left"
fi

if [[ $model == "TStauStau_right" ]] ;then
out="stau-stau_right"
fi

if [[ $model == "TChiStauStau" || $model == "TChiStauStau_ext" || $model == "TChiStauStau_ext2" ]] ;
then


out="C1N2"
fi

if [[ $model == "TChiStauStaux0p95" ]] ;
then


out="C1N2x0p95"
fi

if [[ $model == "TChiStauStaux0p05" ]] ;
then


out="C1N2x0p05"
fi


if [ $model == "TChipmStauSnu" ] ;
then


out="C1C1"
fi

if [ $model == "TChiWZ" ] ;
then


out="ChiWZ"
fi

if [ $model == "TChiWH" ] ;
then


out="ChiWH"
fi

if [ $model == "TChipmWW" ] ;
then


out="C1C1WW"
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


#cat bss > jobmerge${model}_${susy}_LSP${lsp}.sh

#echo "hadd -f -k ${out}_${susy}_LSP${lsp}.root ${model}_${susy}_LSP${lsp}_*.root  ${model}_ext_${susy}_LSP${lsp}_*.root " >> jobmerge${model}_${susy}_LSP${lsp}.sh
echo "hadd -f -k ${out}_${susy}_LSP${lsp}.root ${model}_${susy}_LSP${lsp}_*.root" >> m
#echo "hadd -f -k ${out}_${susy}_LSP${lsp}.root ${model}_${susy}_LSP${lsp}_*.root  ${model}_ext_${susy}_LSP${lsp}_*.root " >> m


#qsub jobmerge${model}_${susy}_LSP${lsp}.sh


fi
let COUNTER=COUNTER+1 


done<points_$model

done
