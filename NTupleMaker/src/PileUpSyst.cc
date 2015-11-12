//#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUpSyst.h"


PileUpSyst::PileUpSyst(PileUp* PU_central) {	//constructor

		h_MC_central = 0;
		h_data_central = 0;
		h_MC_central   = new TH1D(*PU_central->get_h_MC_norm());
		h_data_central = new TH1D(*PU_central->get_h_data_norm());

		h_MC_central->SetDirectory(0);
		h_data_central->SetDirectory(0);
		
		PUsyst["central"] = new PileUp();
		std::cout<<h_MC_central->GetNbinsX()<<std::endl;
		PUsyst["central"]->set_h_MC(h_MC_central);   // the MC distribution is always the same 
		PUsyst["central"]->set_h_data(h_data_central); // pass the shifted data distribution

		// controlli su esistenza, normalizzazione
}

TH1D * PileUpSyst::get_h_MC_central() { return h_MC_central; }
TH1D * PileUpSyst::get_h_data_central() { return h_data_central; }


void PileUpSyst::set_shifted_data(TH1D * h_data_shifted, std::string name) {

	PileUp * PileUp_shifted = new PileUp();   // creating object of type PileUp, containg MC(central) and data(shifted) distrib	
	PileUp_shifted->set_h_MC(h_MC_central);   // the MC distribution is always the same 
	std::cout << "-----------in set_shifted_data --- " << std::endl;
	PileUp_shifted->set_h_data(h_data_shifted); // pass the shifted data distribution

	PUsyst[name] = PileUp_shifted;		    // create an element in the map

	std::cout << "Set shifted data, name =  " << name << " area = " << h_data_shifted->Integral() << std::endl;
	std::cout << "Check central MC : area = " << h_MC_central->Integral() << std::endl;
}

double PileUpSyst::get_PUweight(double Nvertex_MC, std::string name ) {

	std::map<std::string, PileUp *>::iterator it;

 	it = PUsyst.find(name);  

	if ( it == PUsyst.end() ) { 
	std::cout << "Error in PileUpSyst::get_PUweight(double Nvertex_MC, std::string name ): no PU objetc corresponding to name " << name << std::endl;
	exit(1);
	}

	double weight = PUsyst[name]->get_PUweight(Nvertex_MC);
	return weight;

}

