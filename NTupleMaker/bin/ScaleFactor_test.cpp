#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TROOT.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"

using namespace std;

int main () {

std::cout<< " test fake rates " << std::endl;

float tauEta [3]= {0.1, 1.5, 2.0};
std::string WP [6] = {"VLoose","Loose","Medium","Tight","VTight"};

LepTauFakeRate *LepTauFakeRate_test = new LepTauFakeRate();
LepTauFakeRate_test->Init();


std::cout << " --------------------------------------- " << std::endl;
std::cout << "putting gen_match == 1" << std::endl;
std::cout << " --------------------------------------- " << std::endl;

for (int i=0; i<3; i++){
	for (int j=0; j<5; j++){ 
		float SF = LepTauFakeRate_test->get_fakerate("electron",WP[j],tauEta[i],1);
		std::cout<<"WP : "<< WP[j] << " | eta " << tauEta[i] << " | SF " <<  SF << std::endl;
	}
}

std::cout << " --------------------------------------- " << std::endl;
std::cout << "putting gen_match == 3" << std::endl;
std::cout << " --------------------------------------- " << std::endl;

for (int i=0; i<3; i++){
	for (int j=0; j<5; j++){ 
		float SF = LepTauFakeRate_test->get_fakerate("electron",WP[j],tauEta[i],3);
		std::cout<<"WP : "<< WP[j] << " | eta " << tauEta[i] << " | SF " <<  SF << std::endl;
	}
}


std::string WP_mu [2] = {"Loose","Tight"};

std::cout << " --------------------------------------- " << std::endl;
std::cout << "putting gen_match == 2" << std::endl;
std::cout << " --------------------------------------- " << std::endl;
for (int i=0; i<3; i++){
	for (int j=0; j<2; j++){ 
		float SF = LepTauFakeRate_test->get_fakerate("muon",WP_mu[j],tauEta[i],2);
		std::cout<<"WP : "<< WP_mu[j] << " | eta " << tauEta[i] << " | SF " <<  SF << std::endl;
	}
}


std::cout << " --------------------------------------- " << std::endl;
std::cout << "putting gen_match == 4" << std::endl;
std::cout << " --------------------------------------- " << std::endl;
for (int i=0; i<3; i++){
	for (int j=0; j<2; j++){ 
		float SF = LepTauFakeRate_test->get_fakerate("muon",WP_mu[j],tauEta[i],4);
		std::cout<<"WP : "<< WP_mu[j] << " | eta " << tauEta[i] << " | SF " <<  SF << std::endl;
	}
}

return 0;

}

/*
ScaleFactor * SF_test = new ScaleFactor();
SF_test->init_ScaleFactor("/nfs/dust/cms/user/bottav/Area76X/CMSSW_7_6_3_patch2/src/HTT-utilities/LepEffInterface/data/Electron/Electron_IdIso0p15_fall15.root");

double pt [10] = {5, 10, 15, 30, 57, 88};
double eta [3] = {0.5, -2.6};


for (int i=0; i<6; i++) {
	for (int j=0; j<2; j++) {
		cout << "pt : "<< pt[i] << " | eta : " << eta[j] << std::endl;
		cout << "Data : " << SF_test->get_EfficiencyData(pt[i], eta[j]) << " +- " << SF_test->get_EfficiencyDataError(pt[i],eta[j]) <<std::endl;
		cout <<"MC : " << SF_test->get_EfficiencyMC(pt[i], eta[j]) << " +- " << SF_test->get_EfficiencyMCError(pt[i],eta[j]) << std::endl;
		cout << "SF = " << SF_test->get_ScaleFactor(pt[i], eta[j])<< " +- " << SF_test->get_ScaleFactorError(pt[i],eta[j]) << std::endl;
		cout << " --------------------------------------------------------------- " << std::endl;
	}
}

return 0;
}
*/




