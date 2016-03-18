#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include <cmath>

void PileUp::set_h_data(TH1D * user_data_h) {
 
	h_data = new TH1D (*user_data_h);
	h_data->SetDirectory(0);
	if (h_data ==0 ) 
		std::cout << "Error in  PileUp::set_h_data : failed reading the data histogram" << std::endl;
	h_data->Scale( (1/(h_data->Integral()) ) );
	//std::cout << "----------in PileUp::set_h_data --- " << std::endl;
	//std::cout << "Normalized histogram, area = " << h_data->Integral() << std::endl;

	if (h_data != 0 and h_MC != 0) 
		this->make_ratio();
	
}

void PileUp::set_h_MC(TH1D * user_MC_h) {

	h_MC = new TH1D (*user_MC_h);
	h_MC->SetDirectory(0);
	if (h_MC ==0 ) 
		std::cout << "Error in  PileUp::set_h_MC : failed reading the Monte Carlo histogram" << std::endl;
	h_MC->Scale( (1/(h_MC->Integral()) ) );
	if (h_data != 0 and h_MC != 0) 
		this->make_ratio();
}


TH1D* PileUp::get_h_data_norm () {
	
	return h_data;
}


TH1D* PileUp::get_h_MC_norm () {

	return h_MC;
}


void 	PileUp::make_ratio () {

	
	if ( fabs( (h_data->Integral()) -1) > 0.01 or fabs( (h_MC->Integral()) -1 ) > 0.01) {
		std::cout << "Error in PileUp::get_ratio. Distributions not normalized to 1 " << std::endl;
		std::cout << " Integral of h_data is " << h_data->Integral() << std::endl;
		std::cout << " Integral of h_MC is " << h_MC->Integral() << std::endl;		
		exit(1);
	}
	
	h_ratio = new TH1D(*h_data); h_ratio->SetDirectory(0);
	h_ratio->Divide(h_MC);

}



double PileUp::get_PUweight (double Nvertex_MC) {

	if (h_MC == 0 or h_data == 0) {
	std::cout << "Error in PileUp::get_PUweight : missing MC and/or data histogram! " << std::endl;
	exit(1);
	}

	if (h_ratio == 0) {
		std::cout << "Error : missing h_ratio (h_ratio = 0) " << std::endl;
		exit(1);
	}

	int binNumber = h_ratio->FindBin(Nvertex_MC);
	double PUweight = h_ratio->GetBinContent(binNumber); 
	
	return PUweight;

}


