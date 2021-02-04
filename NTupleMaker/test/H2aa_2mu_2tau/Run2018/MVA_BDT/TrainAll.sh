#!/bin/csh

mass=mass

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

  mkdir ma$mass
  
  cd ma$mass

  echo "Preparing files for training on mass point Ma = $mass GeV"

  cp ../trainBDT_lep_lep.py ./
  sed -i "s/herereplacemasspointstring/$mass/g" trainBDT_lep_lep.py
  echo "  "

  cp ../trainBDT_lep_had.py ./
  sed -i "s/herereplacemasspointstring/$mass/g" trainBDT_lep_had.py
  echo "  "

  cp ../trainBDT_had_had.py ./
  sed -i "s/herereplacemasspointstring/$mass/g" trainBDT_had_had.py
  echo "  "

  cp ../HTC_TheScript.sh ../TrainingBDT.submit ./

  echo "Submitintg Training to HTCondor"

  chmod u+x TrainingBDT.submit
  condor_submit TrainingBDT.submit


  cd ..

  echo "  "

  done

echo "  "

done
