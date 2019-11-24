# PATH_TO_TUPLES="/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HiggsCP/Inputs/test/NTuples_mt_2017/"
PATH_TO_TUPLES="/nfs/dust/cms/user/filatovo/HTT/CMSSW_10_2_16/src/HiggsCP/Inputs/test/NTuples_mt_2018/"
PATH_FOR_OUTPUT="./Plots_v3_2018/"

WEIGHT="xsec_lumi_weight*puweight*effweight*mcweight*"
USE_DEEPTAU="true"
APPLY_FF="true"
SHOW_SIGNAL="false"
COMPARE_CP="false"

# eval "root -l -b -q 'Plot_lept_mutau_NNNTuples.C(${"pt_1","p_{T,#mu}[GeV]"},30,20.,80.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)'" 

### beware of spaces inside of parameter brackets! they might crush things 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("pt_1","p_{T,#mu}[GeV]",30,20.,80.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("eta_1","#eta_{#mu}",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("phi_1","#phi_{#mu}[rad]",30,-TMath::Pi(),TMath::Pi(),"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("pt_2","p_{T,#tau}[GeV]",30,20.,80.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("eta_2","#eta_{#tau}",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("phi_2","#phi_{#tau}[rad]",30,-TMath::Pi(),TMath::Pi(),"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("met","{E_{T}}^{mis}[GeV]",30,0.,150.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("njets","njets",6,0,6,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("jdeta","#Delta#eta_{jj}[GeV]",12,0,6,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("mjj","m_{jj}[GeV]",30,0,150.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("dijetpt","p_{T,jj}[GeV]",30,26.,146.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("jpt_1","p_{T,lead.j}[GeV]",30,26.,146.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("jpt_2","p_{T,trail.j}[GeV]",30,26.,146.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("pt_tt","p_{T,#tau#tau}[GeV]",30,0.,150.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("m_sv","m_{#tau#tau}[GeV]",50,0.,250.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("m_vis","m_{vis}[GeV]",30,0.,150.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("mt_1","m_{T}[GeV]",30,0.,150.,"%s","","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("acotautau_01","#phi_{CP}",11,0.,2*TMath::Pi(),"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 

### b-jets ###
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("nbtag","nbtags",6,0,6,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bpt_1","p_{T,lead.bj}[GeV]",30,16.,136.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bpt_2","p_{T,trail.bj}[GeV]",30,16.,136.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("beta_1","#eta_{lead.bj}",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bphi_1","#phi_{lead.bj}[rad]",30,-TMath::Pi(),TMath::Pi(),"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("beta_2","#eta_{trail.bj}",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bphi_2","#phi_{trail.bj}[rad]",30,-TMath::Pi(),TMath::Pi(),"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bcsv_1","DeepFlavour_probb_1",40,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("bcsv_2","DeepFlavour_probb_2",40,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 

### NEW ###
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("ipx_1","IP_helix_x_1",30,-0.025,0.035,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("ipy_1","IP_helix_y_1",30,-0.025,0.035,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("ipz_1","IP_helix_z_1",30,-0.025,0.035,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 

root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("mt_fast","mt_fast",30,0.,300.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("pt_fast","pt_fast",30,0.,300.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("eta_fast","eta_fast",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("phi_fast","phi_fast",40,-4.,4.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 

root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("deepTauVsEleRaw_2","deepTauVsEleRaw_2",20,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("deepTauVsJetRaw_2","deepTauVsJetRaw_2",20,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("deepTauVsMuRaw_2","deepTauVsMuRaw_2",20,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("MVADM2017v1DM10raw_2","MVADM2017v1DM10raw_2",20,0.,1.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 

root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("puppimt_1","puppimt_1",40,0.,400.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("puppimt_2","puppimt_2",40,0.,400.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("puppimet","puppimet",40,0.,400.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
root -l -b -q `printf 'Plot_lept_mutau_NNNTuples.C("puppimetphi","puppimetphi",30,-3.,3.,"%s","(mt_1<50)*","Events",-1,"%s","%s",2017,%s,%s,false,false,%s,%s)' $WEIGHT $PATH_TO_TUPLES $PATH_FOR_OUTPUT $USE_DEEPTAU $APPLY_FF $SHOW_SIGNAL $COMPARE_CP` 
