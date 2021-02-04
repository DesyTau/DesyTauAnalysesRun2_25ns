#!/bin/sh

echo submitting job 3.6
cat > HTC_TheScript_3.6.sh <<EOF
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cd ${CMSSW_BASE}/src
cmsenv
cd -

echo "  "
echo " Runing over the samples between 3.6 and 4 the classification and creating corresponding datacards "
echo "  "
root -l -b -q 'CreateInputs.C(3.6,4.0)'
EOF
chmod u+x HTC_TheScript_3.6.sh
./HTC_qsub.sh HTC_TheScript_3.6.sh

for i in 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0
do
echo submitting job $i
cat > HTC_TheScript_$i.sh <<EOF
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cd ${CMSSW_BASE}/src
cmsenv
cd -

echo "  "
echo " Runing over the samples between $i and $(( $i+1.0 )) the classification and creating corresponding datacards "
echo "  "
root -l -b -q 'CreateInputs.C($i,$(( $i+1.0 )))'
EOF
chmod u+x HTC_TheScript_$i.sh
./HTC_qsub.sh HTC_TheScript_$i.sh
done


echo submitting job 21.0
cat > HTC_TheScript_21.0.sh <<EOF
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cd ${CMSSW_BASE}/src
cmsenv
cd -

echo "  "
echo " Runing over the samples between 21.0 and 21.0 the classification and creating corresponding datacards "
echo "  "
root -l -b -q 'CreateInputs.C(21.0,21.2)'
EOF
chmod u+x HTC_TheScript_21.0.sh
./HTC_qsub.sh HTC_TheScript_21.0.sh
