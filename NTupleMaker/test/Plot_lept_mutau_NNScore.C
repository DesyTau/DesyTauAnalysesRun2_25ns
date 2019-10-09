#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "DesyTauAnalyses/NTupleMaker/test/Plot_lept_mutau_NNNTuples.C"


void   Plot_lept_mutau_NNScore(TString directory = "HtautauCP_mutau/Outputs/NTuples_mt_2017/predictions_2017/",
			       TString outputDir = "./Plots/",
			       int year=2017,
			       bool FFmethod = true,  
			       bool LargeScale = true,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = false,
			       int scaleSignal = 1.,
			       bool blindData = false,
			       bool FORCE = false
			       )
{
  
  for(int categoryIndex=0;categoryIndex<2;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_01",
			      "#phi_{CP}",
			      5,0.,2*TMath::Pi(),
			      "xsec_lumi_weight*puweight*effweight*mcweight*",
			      "(mt_1<50)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      FFmethod,  
			      LargeScale,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }

  bool _logY = false;
  for(int categoryIndex=0;categoryIndex<8;categoryIndex++){
    if(categoryIndex<2) _logY=true;
    else _logY=logY;
    Plot_lept_mutau_NNNTuples("predicted_prob",
			      "NN Score",
			      10,0.,1.,
			      "xsec_lumi_weight*puweight*effweight*mcweight*",
			      "(mt_1<50)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      FFmethod,  
			      LargeScale,  
			      _logY,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
}
