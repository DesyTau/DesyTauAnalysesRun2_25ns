#ifndef LepTauFakeRate_h
#define LepTauFakeRate_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TString.h"

using namespace std;

class LepTauFakeRate {

std::map<std::string, std::vector<float>> ETAUfakerates;
std::map<std::string, std::vector<float>> MUTAUfakerates;

public: 

	LepTauFakeRate(){};
	~LepTauFakeRate(){};

    void Init(){
		// e->tau fake rates
		// Scale Factor values corresponding to eta<1.46, 1.46<eta<1.558, eta>1.558

		// ICHEP 2016 dataset
		/*std::vector<float> FakeRates_VLoose_ele = {1.292, 1., 1.536}; 
		std::vector<float> FakeRates_Loose_ele =  {1.431, 1., 1.722}; 
		std::vector<float> FakeRates_Medium_ele = {1.549, 1., 1.678}; 
		std::vector<float> FakeRates_Tight_ele =  {1.505, 1., 1.994}; 
		std::vector<float> FakeRates_VTight_ele = {1.420, 1., 2.017}; */		

		// full 2016 dataset
		// no number available for loose and very loose. Set to 1. 
		std::vector<float> FakeRates_VLoose_ele = {1., 1., 1.}; 
		std::vector<float> FakeRates_Loose_ele =  {1., 1., 1.}; 
		std::vector<float> FakeRates_Medium_ele = {1.648, 1., 1.241}; 
		std::vector<float> FakeRates_Tight_ele =  {1.867, 1., 1.456}; 
		std::vector<float> FakeRates_VTight_ele = {1.967, 1., 1.401}; 
		ETAUfakerates.insert(std::make_pair("VLoose", FakeRates_VLoose_ele));
		ETAUfakerates.insert(std::make_pair("Loose",  FakeRates_Loose_ele));
		ETAUfakerates.insert(std::make_pair("Medium", FakeRates_Medium_ele));
		ETAUfakerates.insert(std::make_pair("Tight",  FakeRates_Tight_ele));
		ETAUfakerates.insert(std::make_pair("VTight", FakeRates_VTight_ele));

	
		// mu->tau fake rates
	    // Scale Factor values corresponding to eta={0, 0.4, 0.8, 1.2, 1.7, >1.7} 
			
		// full 2016 from Artur
		/*std::vector<float> FakeRates_Loose_mu = {1.09, 1.05, 1.1, 1.03, 1.2};
		std::vector<float> FakeRates_Tight_mu = {1.37, 1.2, 1.14, 1.0, 1.8};*/

		// full 2016 from Yiwen on 28.02.17
		std::vector<float> FakeRates_Loose_mu = {1.012, 1.007, 0.870, 1.154, 2.281};
		std::vector<float> FakeRates_Tight_mu = {1.263, 1.364, 0.854, 1.712, 2.324};
		MUTAUfakerates.insert(std::make_pair("Loose", FakeRates_Loose_mu));
		MUTAUfakerates.insert(std::make_pair("Tight", FakeRates_Tight_mu));			
	}

    float get_fakerate(string lepton, string AntiLeptonDiscriminatorWP, float tauEta, int tau_gen_match ){
		float fakerate = 1.;

		if (lepton == "electron" && (tau_gen_match == 1 || tau_gen_match == 3) ){
			if (fabs(tauEta)<1.46 )                       fakerate = ETAUfakerates.find( AntiLeptonDiscriminatorWP )->second[0];
			if (fabs(tauEta)>=1.46 && fabs(tauEta)<1.558) fakerate = ETAUfakerates.find( AntiLeptonDiscriminatorWP )->second[1];
			if (fabs(tauEta)>=1.558 )                     fakerate = ETAUfakerates.find( AntiLeptonDiscriminatorWP )->second[2];
        	}
		else if (lepton == "muon" && (tau_gen_match == 2 || tau_gen_match == 4) ) {
				if (fabs(tauEta)<0.4 )                      fakerate = MUTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[0];
				if (fabs(tauEta)>=0.4 && fabs(tauEta)<0.8)  fakerate = MUTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[1];
				if (fabs(tauEta)>=0.8 && fabs(tauEta)<1.2 ) fakerate = MUTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[2];
				if (fabs(tauEta)>=1.2 && fabs(tauEta)<1.7 ) fakerate = MUTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[3];
				if (fabs(tauEta)>=1.7 )                     fakerate = MUTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[4];
			}
		else fakerate = 1.;

		return fakerate;

	}

}; // end of class 
#endif
