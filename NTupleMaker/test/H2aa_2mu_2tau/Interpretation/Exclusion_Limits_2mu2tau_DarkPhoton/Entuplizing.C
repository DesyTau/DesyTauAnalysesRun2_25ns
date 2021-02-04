#include "Riostream.h"
#include <string>

void Entuplizing()
{

	std::vector<Double_t> mass_values;

	auto f = TFile::Open("BR.root", "RECREATE");

	TNtuple ntuple("Ntuple_DarkPhoton", "DarkPhoton", "m_a:BR_a_ll");

	TString dir = "BR/";

	ifstream in ; in .open(Form("%sHiggsedDarkPhoton_BrTableData.txt", dir.Data()));

	Float_t x1, x2, x3, x4, x5;
	Int_t nlines = 0;

	while (1)
	{
		in >> x1 >> x2 >> x3 >> x4 >> x5;
		if (! in .good()) break;
		if (nlines < 15) printf("x1=%8f, x3=%8f\n", x1, x3);
		ntuple.Fill(x1, x3);

		mass_values.push_back(x1);

		nlines++;
	}

	printf(" found %d points\n", nlines);

	in .close();

	f->Write();

	std::cout << mass_values.size() << std::endl;

	f->WriteObject(&mass_values, "mass_values");

}
