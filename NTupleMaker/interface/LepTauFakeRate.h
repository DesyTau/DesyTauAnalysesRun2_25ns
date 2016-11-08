#ifndef LepTauFakeRate_h
#define LepTauFakeRate_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TString.h"

using namespace std;

class LepTauFakeRate {

std::map<std::string, std::vector<float>> LTAUfakerates;
//std::map<std::string, std::vector<float>> MUTAUfakerates;

public: 

	LepTauFakeRate(){};
	~LepTauFakeRate(){};

    void Init(TString channel){
		if (channel == "et"){
			// Scale Factor values corresponding to eta<1.46, 1.46<eta<1.558, eta>1.558
			std::vector<float> FakeRates_VLoose = {1.292, 1., 1.536}; 
			std::vector<float> FakeRates_Loose =  {1.431, 1., 1.722}; 
			std::vector<float> FakeRates_Medium = {1.549, 1., 1.678}; 
			std::vector<float> FakeRates_Tight =  {1.505, 1., 1.994}; 
			std::vector<float> FakeRates_VTight = {1.420, 1., 2.017}; 
			LTAUfakerates.insert(std::make_pair("VLoose", FakeRates_VLoose));
			LTAUfakerates.insert(std::make_pair("Loose",  FakeRates_Loose));
			LTAUfakerates.insert(std::make_pair("Medium", FakeRates_Medium));
			LTAUfakerates.insert(std::make_pair("Tight",  FakeRates_Tight));
			LTAUfakerates.insert(std::make_pair("VTight", FakeRates_VTight));			 
	    }

        if (channel == "mt") {
	    	// Scale Factor values corresponding to eta={0, 0.4, 0.8, 1.2, 1.7, >1.7} 
			std::vector<float> FakeRates_Loose = {1.137, 1.067, 1.221, 1.460, 1.609};
			std::vector<float> FakeRates_Tight = {1.418, 1.134, 1.260, 1.660, 1.205};
			LTAUfakerates.insert(std::make_pair("Loose", FakeRates_Tight));
			LTAUfakerates.insert(std::make_pair("Tight", FakeRates_Loose));			
		}
		else {
			std::cout << "Error initialising LepTauFakeRate , channel "<< channel << " not known. Exiting. " << std::endl;
			exit(1); 	
		}

	}

    float get_fakerate(TString channel, string AntiLeptonDiscriminatorWP, float tauEta, int tau_gen_match ){
		float fakerate;
		if (tau_gen_match != 1 )
			fakerate = 1.;
	    else {
			if (channel == "et"){
				if (fabs(tauEta)<1.46 )                       fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[0];
				if (fabs(tauEta)>=1.46 && fabs(tauEta)<1.558) fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[1];
				if (fabs(tauEta)>=1.558 )                     fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[2];
        	}
			if (channel == "mt") {
				if (fabs(tauEta)<0.4 ) fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[0];
				if (fabs(tauEta)>=0.4 && fabs(tauEta)<0.8)  fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[1];
				if (fabs(tauEta)>=0.8 && fabs(tauEta)<1.2 ) fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[2];
				if (fabs(tauEta)>=1.2 && fabs(tauEta)<1.7 ) fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[3];
				if (fabs(tauEta)>=1.7 )                     fakerate = LTAUfakerates.find( AntiLeptonDiscriminatorWP )->second[4];
			}
        }  

		return fakerate;
		//else if (channel == "mt")
	}

}; // end of class 
#endif
