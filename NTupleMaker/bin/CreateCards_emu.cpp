//#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
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

  int nBinsBTag2 = 7;
  double binsBTag2[8] =
    {40,80,120,160,200,300,500,5000};

  int nBinsTTbar = 12;
  double binsTTbar[13] = {40, 60, 80, 100,
			  120, 150, 200, 250, 
			  300, 400, 500, 1000,5000};

  //  std::vector<double> binsVector = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000};

  std::vector<double> binsNoBTagVector = {0,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,450,500,600,700,800,900,1100,1300,1500,1700,1900,2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100,4300,4500,4700,5000};

  std::vector<double> binsBTagVector = {0,60,80,100,120,140,160,180,200,250,300,350,400,500,600,700,800,900,1100,1300,1500,1700,1900,2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100,4300,4500,4700,5000};

  //  int nBinsFine = binsVector.size()-1;
  nBinsBTag = binsBTagVector.size()-1;
  nBinsNoBTag = binsNoBTagVector.size()-1;

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
    cardsEMu->SetVariableToPlot(variable,nBinsNoBTag,binsNoBTagVector);
  else {
    if (category.Contains("_NbtagGt1")||
	category.Contains("_Nbtag1"))
      cardsEMu->SetVariableToPlot(variable,nBinsBTag,binsBTagVector);
    else 
      cardsEMu->SetVariableToPlot(variable,nBinsNoBTag,binsNoBTagVector);
  }
  cardsEMu->Run();
  cardsEMu->CloseFile();
  delete cardsEMu;

}
