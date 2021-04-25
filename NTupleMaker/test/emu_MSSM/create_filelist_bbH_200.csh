#!/bin/bash
python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_M200 --query "/bbh_m200_LHE/dwinterb-bbH_m200_2016_MINIAOD-53f8667ba4b240d5eafd36e71bf34742/USER" --phys03 --outputfile 2016/SUSYGluGluToBBHToTauTau_powheg_M200.list
python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_M200 --query "/bbh_m200_GENSIM-v3/dwinterb-bbH_m200_2017_MINIAOD-b63beb1ae05c0e254c43785544367ee5/USER" --phys03 --outputfile 2017/SUSYGluGluToBBHToTauTau_powheg_M200.list
python read_filelist_from_das.py --nick SUSYGluGluToBBHToTauTau_M200 --query "/bbh_m200_GENSIM-v2/dwinterb-bbH_m200_2018_MINIAOD-0bd58594e6ade05f64e0c3a8301c3139/USER" --phys03 --outputfile 2018/SUSYGluGluToBBHToTauTau_powheg_M200.list
