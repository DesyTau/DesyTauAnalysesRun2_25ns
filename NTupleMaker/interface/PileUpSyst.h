#ifndef PileUpSyst_h
#define PileUpSyst_h

#include "TROOT.h"
#include "TH1D.h" 
#include "PileUp.h"
#include <iostream>
#include <map>


//calculate sys shifted pu weights (shifting data) wrt. given central data and MC PU distributions (passed to the constructor)
class PileUpSyst {  


	private: 
		std::map<std::string, PileUp *> PUsyst;
		TH1D * h_MC_central;
		TH1D * h_data_central;

	public: 
		PileUpSyst(PileUp*); //constructor

		TH1D * get_h_MC_central();
		TH1D * get_h_data_central();
		void set_shifted_data(TH1D *, std::string);
		double get_PUweight(double, std::string);

		~ PileUpSyst() {}; 


};



#endif
