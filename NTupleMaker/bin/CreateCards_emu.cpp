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
  if (argc!=6) {
    std::cout << "usage of the script : CreateCards_emu [era=2016,2017,2018] [Sample=Data,EWK,DYJets,EMB,TTbar,SMHiggs,MSSMHiggs] [category] [TriggerOption=0,1,2] [Systematics=0,1]" << std::endl;
    exit(-1);
  }
  int nBinsNoBTag = 29;
  double binsNoBTag[30] = 
    {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,
     225,250,275,300,325,350,400,500,700,900,4000};
  
  int nBinsBTag = 16;
  double binsBTag[17] =
    {20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};

  int nBinsTTbar = 13;
  double binsTTbar[14] = {0, 50, 60, 80, 100,
			  120, 150, 200, 250, 300,
			  500, 700, 1000, 5000};

  TString era(argv[1]);
  TString sample(argv[2]);
  TString category(argv[3]);
  TString trigger(argv[4]);
  int triggerOption = atoi(argv[4]);
  int systematics = atoi(argv[5]);

  TString inputDir("/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan1/");
  //  TString outputDir("/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan1/datacards_"+trigger+"/");
  TString outputDir("/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan10/datacards_"+trigger+"/");
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
  if (category.Contains("_ttbar_"))
  cardsEMu->SetVariableToPlot(variable,nBinsTTbar,binsTTbar);
  else {
    if (category.Contains("_btag_"))
      cardsEMu->SetVariableToPlot(variable,nBinsBTag,binsBTag);
    else 
      cardsEMu->SetVariableToPlot(variable,nBinsNoBTag,binsNoBTag);
  }
  cardsEMu->Run();
  cardsEMu->CloseFile();
  delete cardsEMu;


}
