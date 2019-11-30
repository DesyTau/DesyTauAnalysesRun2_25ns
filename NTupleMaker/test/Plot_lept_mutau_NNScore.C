#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "DesyTauAnalyses/NTupleMaker/test/Plot_lept_mutau_NNNTuples.C"


void   Plot_lept_mutau_NNScore(TString directory = "/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HiggsCP/Outputs/nobveto/NTuples_mt_2017/predictions_2017/",
			       TString outputDir = "./Plots/",
			       int year=2017,
			       bool DeepTau = true, 
			       bool FFmethod = true,  
			       bool LargeScale = false,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = true,
			       int scaleSignal = 1.,
			       bool blindData = false,
			       bool FORCE = false
			       )
{
  
  const int nCategories = 8;
  const int nSigCategories = 2;
  
  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_00",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "xsec_lumi_weight*weight*",
			      "(puppimt_1<50&&pt_1>20)*(tau_decay_mode_2==0)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      DeepTau,
			      FFmethod,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }

  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_01",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "xsec_lumi_weight*weight*",
			      "(puppimt_1<50&&pt_1>20)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      DeepTau,
			      FFmethod,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
  
  bool _logY = false;
  bool _largeScale = false;
  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    if(categoryIndex<nSigCategories){
      blindData=true;
      _logY=true;
      _largeScale=true;
    }
    else {
      _logY=logY;
      _largeScale=LargeScale;
    }
    Plot_lept_mutau_NNNTuples("predicted_prob",
			      "NN Score",
			      10,0.,1.,
			      "xsec_lumi_weight*weight*",
			      "(puppimt_1<50&&pt_1>20)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      DeepTau,
			      FFmethod,  
			      _largeScale,  
			      _logY,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
}
