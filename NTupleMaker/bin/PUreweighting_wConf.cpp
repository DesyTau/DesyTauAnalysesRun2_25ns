#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"

#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUpSyst.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUpInit.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"

using namespace std;

int main () {


//read the MC histo from TTree
TFile * fileMC = new TFile ("../test/PileUp/MC_pileup_fake.root");
TTree * treeMC = (TTree*)fileMC->Get("MCpileup");

double pu;
treeMC->SetBranchAddress("pu", &pu);

// MC reweighted distributions 
int nbins = 50; double xmin = 0.; double xmax = 50.;
TH1D * h_systUp = new TH1D("h_MC_rw_up", "MC reweighted - Syst UP",  nbins, xmin, xmax);
TH1D * h_systDown = new TH1D("h_MC_rw_down", "MC reweighted - Syst DOWN", nbins, xmin, xmax);
TH1D * h_systCentral = new TH1D("h_MC_rw_central", "MC reweighted - CENTRAL", nbins, xmin, xmax);

// systematics 


Config configFile("../test/PileUp/PUconfig.conf");
PileUpSyst * PUsyst = 0;
PUsyst = PileUpInit::initPU(configFile);

double weightUp, weightDown, weightCentral;


for (int i=0; i<treeMC->GetEntries(); i++) 
	{
		treeMC->GetEntry(i);
		weightUp = PUsyst->get_PUweight(pu, "up");
		h_systUp->Fill(pu,weightUp);
		weightDown = PUsyst->get_PUweight(pu, "down");
		h_systDown->Fill(pu, weightDown);
		weightCentral = PUsyst->get_PUweight(pu, "central");
		h_systCentral->Fill(pu, weightCentral);
	}


TFile * fileOut = new TFile ("../test/PileUp/reweighting_syst_wCentral.root","recreate");
fileOut->cd();
h_systUp->Write("h_MC_rw_up");
h_systDown->Write("h_MC_rw_down");
h_systCentral ->Write("h_MC_rw_central");

fileMC->Close();
fileOut->Close();

std::cout << "----------------------------- done PUreweighting_wConf.cpp ---------------------" << std::endl;
std::cout << "fileOut : ../test/PileUp/reweighting_syst_wCentral.root " << std::endl;

return 0;

}


