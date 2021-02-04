#!/bin/sh

mass=mass
massf=massf

for i in {3..21}

do

  for j in 0 2 4 6 8

  do

  if [[ "$i" -eq 3 && "$j" -le 4 ]]
  then
   continue
  fi

  if [[ "$i" -eq 21 && "$j" -ge 2 ]]
  then
   continue
  fi

  echo "  "
  echo "  "
  
  eval "mass=${i}p${j}"
  eval "massf=${i}.${j}"

   cat > HTC_TheScript_$mass.sh <<EOF
   source /cvmfs/cms.cern.ch/cmsset_default.sh
   export SCRAM_ARCH=slc6_amd64_gcc630
   cd ${CMSSW_BASE}/src
   cmsenv
   cd -

   echo "Computing limits for Ma = $mass GeV"
   combineCards.py _lep_lep=haa-13TeV_2018_ma$mass"_"lep_lep.txt _lep_had=haa-13TeV_2018_ma$mass"_"lep_had.txt _had_had=haa-13TeV_2018_ma$mass"_"had_had.txt > haa-13TeV_2018_ma$mass.txt

  combine -M AsymptoticLimits haa-13TeV_2018_ma$mass.txt -m $massf

EOF
   chmod u+x HTC_TheScript_$mass.sh
   ./HTC_qsub.sh HTC_TheScript_$mass.sh HTC_TheScript_$mass

  done

done
