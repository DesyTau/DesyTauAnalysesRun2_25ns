dir="/nfs/dust/cms/group/higgs-kit/80x_MC_SUSY_v1/"
dir="/nfs/dust/cms/group/susy-desy/Run2/Stau/8025_MC_SUSY_v1/"

sources="TChipmStauSnu   TChiStauStau   TStauStau"

#sources="TStauStau"
sources="TChiStauStau"
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
ls $dir/SMS-$model/*.root > listfiles_$model

### first loop for all .root files
while read file
do

## next for each point
while read line
do


	susy=`echo $line  | awk -F " " '{print $1}'`
	lsp=`echo $line  | awk -F " " '{print $2}'`


if [ $lsp == "0" ];then

	lsp=1
fi

outfile=${model}_${susy}_LSP${lsp}_${COUNTER}.root
if [ ! -f $outfile ] ; then


cp copytrees_C copytrees.C

chmod 777 copytrees.C
f2=$(basename $file)

sed -i 's/MODEL/'$model'/g' copytrees.C

sed -i 's/FILEIN/'$f2'/g' copytrees.C

sed -i 's/SUSYMASS/'$susy'/g' copytrees.C

sed -i 's/N1MASS/'$lsp'/g' copytrees.C

sed -i 's/NUMBER/'$COUNTER'/g' copytrees.C

cp bss job_skim${susy}_LSP${lsp}_${COUNTER}.sh
mv copytrees.C copytrees${model}_${susy}_LSP${lsp}_${COUNTER}.C

echo "root -l -q -b copytrees${model}_${susy}_LSP${lsp}_${COUNTER}.C " >>job_skim${susy}_LSP${lsp}_${COUNTER}.sh

 qsub job_skim${susy}_LSP${lsp}_${COUNTER}.sh

#mv skim.root ${out}${susy}_LSP${lsp}_${COUNTER}.root

fi
let COUNTER=COUNTER+1 


done<points_$model


done<listfiles_$model

done
