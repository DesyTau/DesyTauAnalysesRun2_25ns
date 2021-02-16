
#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/CardsEMu.h"

using namespace std;

int main(int argc, char * argv[]) {
  
  // argument - config file
  if (argc!=5) {
    std::cout << "usage of the script : CreateCards_emu [era=2016,2017,2018] [Sample=Data,EWK,DYJets,EMB,TTbar,SMggH,SMothers,bbH,ggH_{t,b,i}] [category] [Systematics=0,1]" << std::endl;
    exit(-1);
  }

  int nBinsNoBTag = 27;
  double binsNoBTag[28] = 
    {40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,
     225,250,275,300,325,350,400,500,700,900,5000};
  
  int nBinsBTag = 15;
  double binsBTag[16] =
    {40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,5000};

  int nBinsTTbar = 12;
  double binsTTbar[13] = {40, 60, 80, 100,
			  120, 150, 200, 250, 300,
			  500, 700, 1000, 5000};

  TString era(argv[1]);
  TString sample(argv[2]);
  TString category(argv[3]);
  int systematics = atoi(argv[4]);
  int triggerOption = 0;

  TString inputDir("/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/");
  TString outputDir("/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/datacards/");
  TString variable("mt_tot");
  int nbins = 20;
  double xmin = 0;
  double xmax = 1000;
  bool runWithSystematics = true;
  if (systematics==0) runWithSystematics = false;
  bool runOnEmbedded = true;

  CardsEMu * cardsEMu = new CardsEMu(sample,
				     era,
				     category,
				     inputDir,
				     outputDir,
				     variable,
				     nbins,
				     xmin,
				     xmax,
				     triggerOption,
				     runWithSystematics,
				     runOnEmbedded
				     );
  
  cardsEMu->PrintSamples();
  cardsEMu->PrintShapeSystematics();
  cardsEMu->PrintWeightSystematics();
  if (category.Contains("_DZetaLtm35"))
  cardsEMu->SetVariableToPlot(variable,nBinsTTbar,binsTTbar);
  else {
    if (category.Contains("_NbtagGt1"))
      cardsEMu->SetVariableToPlot(variable,nBinsBTag,binsBTag);
    else 
      cardsEMu->SetVariableToPlot(variable,nBinsNoBTag,binsNoBTag);
  }
  cardsEMu->Run();
  cardsEMu->CloseFile();
  delete cardsEMu;


}
