#ifndef PileUp_h
#define PileUp_h

#include "TROOT.h"
#include "TH1D.h" 
#include <iostream>


class PileUp {

	private:
		TH1D * h_data;
		TH1D * h_MC;
		TH1D * h_ratio;

	protected:

		void make_ratio(); 

	public:
		PileUp() {h_data = 0; h_MC=0; h_ratio=0; }; //constructor

		void set_h_data (TH1D*);
		void set_h_MC (TH1D*);

		TH1D* get_h_data_norm();
		TH1D* get_h_MC_norm();

		double get_PUweight (double) ;
	
		~ PileUp() {}; //destructor
};

#endif


