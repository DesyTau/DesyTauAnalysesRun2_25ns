#!/bin/sh


mass=mass
massD=massD
year=(2016 2017 2018)

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
  eval "massD=${i}.${j}"
  

  echo "Copying files needed for limit calculation for Ma = $mass GeV"

    for k in {1..3}

    do

    cp /nfs/dust/cms/user/consuegs/Analyses/H2aa_2mu2tau/Run$year[$k]/Inputs/DataCards/haa-13TeV_$year[$k]_ma$mass.txt ./
    cp /nfs/dust/cms/user/consuegs/Analyses/H2aa_2mu2tau/Run$year[$k]/Inputs/DataCards/haa-13TeV_$year[$k]_ma$mass*.root ./

    done

  done


done

