dir="/nfs/dust/cms/group/susy-desy/Run2/Stau/8025_MC_SUSY_v2/"

sources="TChipmStauSnu   TChiStauStau   TStauStau"

sources="TStauStau"
sources="TChiStauStau_ext2"
#sources="TStauStau"
#sources="TChipmStauSnu"
#sources="TStauStau"
#sources="TChiWZ"
#sources="TChiStauStau TChiStauStau_ext TChiStauStau_ext2"
#sources="TChipmWW"
#sources="SMS-TChipmWW"
sources="TChiStauStaux0p05 TChiStauStaux0p95 TChiStauStau TChiStauStau_ext TChiStauStau_ext2"
sources="TChiStauStaux0p05 TChiStauStaux0p95"
#sources="TChiStauStau TChiStauStau_ext TChiStauStau_ext2"
#sources="TChipmWW TChiWH TChipmStauNu TStauStau_left TStauStau_right"



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

if [ $model == "TChiWZ" ] ;
then


out="ChiWZ"
fi


COUNTER=0
ls $dir/SMS-$model/*.root > listfiles_$model

### first loop for all .root files
#while read file
#do

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


#cp copytreesext_C copytrees${model}.C
cp copytrees2_C copytrees${model}.C

chmod 777 copytrees${model}.C
#f2=$(basename $file)

sed -i 's/MODEL/'$model'/g' copytrees${model}.C

#sed -i 's/FILEIN/'$f2'/g' copytrees${model}.C

sed -i 's/SUSYMASS/'$susy'/g' copytrees${model}.C

sed -i 's/N1MASS/'$lsp'/g' copytrees${model}.C

sed -i 's/NUMBER/'$COUNTER'/g' copytrees${model}.C

cp bss job_skim${susy}_LSP${lsp}_${COUNTER}_${model}.sh
mv copytrees${model}.C copytrees${model}_${susy}_LSP${lsp}_${COUNTER}.C

echo "root -l -q -b copytrees${model}_${susy}_LSP${lsp}_${COUNTER}.C " >>job_skim${susy}_LSP${lsp}_${COUNTER}_${model}.sh

echo "rm copytrees${model}_${susy}_LSP${lsp}_${COUNTER}.C" >> job_skim${susy}_LSP${lsp}_${COUNTER}_${model}.sh
 qsub job_skim${susy}_LSP${lsp}_${COUNTER}_${model}.sh


fi
let COUNTER=COUNTER+1 


done<points_$model


#done<listfiles_$model

done
