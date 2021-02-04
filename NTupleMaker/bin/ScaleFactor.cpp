#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"

using namespace std;

int main (int argc, char * argv[]) {

  std::cout<< " test scale factors " << std::endl;

  TString era(argv[1]);
  string cmsswBase = (getenv("CMSSW_BASE"));

  TString workspace_filename = TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/CorrectionWS_IC/htt_scalefactors_legacy_"+era+".root";
  TFile *f_workspace = new TFile(workspace_filename, "read");
  if (f_workspace->IsZombie()) {
    std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
     exit(-1);
   }

  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");

  int nbins = 10;
  float bins[11] = {25,30,35,40,45,50,55,60,70,80,100};
  
  int nbinsEta = 4;
  float etaBins[5] = {0.0,0.9,1.2,2.1,2.4};

  TString EtaBins[4] = {"EtaLt0p9","Eta0p9to1p2","Eta1p2to2p1","EtaGt2p1"};
  
  TFile * file = new TFile("MuonIdIso_"+era+".root","recreate");
  file->cd("");
  for (int iEta=0; iEta<nbinsEta; ++iEta) {
    file->cd("");
    TH1D * histData = new TH1D("data_"+EtaBins[iEta],"",nbins,bins);
    TH1D * histMC   = new TH1D("mc_"+EtaBins[iEta],"",nbins,bins);
    for (int iPt=0; iPt<nbins; ++iPt) {
      double eta = 0.5*(etaBins[iEta]+etaBins[iEta+1]);
      double pt = 0.5*(bins[iPt]+bins[iPt+1]);
      w->var("m_pt")->setVal(pt);
      w->var("m_eta")->setVal(eta);
      w->var("m_iso")->setVal(0.0);
      double eff_data = w->function("m_idiso_ic_data")->getVal();
      double eff_mc = w->function("m_idiso_ic_mc")->getVal();
      histData->SetBinContent(iPt+1,eff_data);
      histMC->SetBinContent(iPt+1,eff_mc);
    }
    file->cd("");
    histData->Write("data_"+EtaBins[iEta]);
    histMC->Write("mc_"+EtaBins[iEta]);
  }
  file->Write();
  file->Close();

  return 0;
}





