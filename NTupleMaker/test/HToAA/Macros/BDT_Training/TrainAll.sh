#!/bin/csh

for i in 4 5 6 7 8 9 10 11 12 13 14 15 17 19 21
do

echo "  "

cd ma$i/
echo "Starting Training on mass point Ma = $i GeV"

#rm -r -f *

cp ../../SCDoMuFilter_Samples/SUSYGluGluToHToAA_AToTauTau_SCDoMuFilter_M-$i.root ../../DoubleMuon_Run2016.root ./
cp ../../SUSYVBFToHToAA_AToTauTau_M-$i.root ../../SUSYVH_HToAA_AToTauTau_M-$i.root ./
cp ../../SUSYttH_HToAA_AToTauTau_M-$i.root ../../SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i.root ./

cp ../trainBDT.py ./
sed -i "s/herereplacemasspointstring/$i/g" trainBDT.py
python trainBDT.py



cd ..

echo "  "

done
