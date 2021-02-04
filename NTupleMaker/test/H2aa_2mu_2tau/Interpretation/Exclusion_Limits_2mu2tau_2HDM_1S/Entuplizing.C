#include "Riostream.h"
#include <string>

void Entuplizing()
{

	int nType = 4;
	TString HDM_Type[] = { "_I", "_II", "_III", "_IV" };
	int nTan = 24;
	TString Tan_Beta[] = { "_0.5", "_0.6", "_0.7", "_0.8", "_0.9", "_1.0", "_1.5", "_2.0", "_2.5", "_3.0", "_3.5", "_4.0", "_4.5", "_5.0", "_5.5", "_6.0", "_6.5", "_7.0", "_7.5", "_8.0", "_8.5", "_9.0", "_9.5", "_10.0", "_" };

	std::vector<Double_t> mass_values;
	std::vector<Double_t> tanbeta_values;

	auto f = TFile::Open("BR.root", "RECREATE");

	for (int i = 0; i < nType; i++)
	{

		TNtuple ntuple("Ntuple" + HDM_Type[i], "Ntuple" + HDM_Type[i], "m_a:tan_beta:BR_a_tautau:BR_a_mumu");

		for (int j = 0; j < nTan; j++)
		{

			if (i == 0) j = nTan;
			// each file has 11 columns of float data

			TString dir = "BR/";

			ifstream in ; in .open(Form("%sBR" + HDM_Type[i] + Tan_Beta[j] + ".dat", dir.Data()));

			Float_t x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
			Int_t nlines = 0;

			while (1)
			{
				in >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11;
				if (! in .good()) break;
				if (nlines < 15) printf("x1=%8f, x2=%8f, x3=%8f\n", x1, x2, x3);
				ntuple.Fill(x2, x3, x6, x7);

				if (HDM_Type[i] == "_II")
				{
					if (x3 == 0.5) mass_values.push_back(x2);
					if (x2 == 1.) tanbeta_values.push_back(x3);
				}

				nlines++;
			}
			printf(" found %d points\n", nlines);

			in .close();
		}

		f->Write();
	}

	std::cout << mass_values.size() << std::endl;
	std::cout << tanbeta_values.size() << std::endl;

	f->WriteObject(&mass_values, "mass_values");
	f->WriteObject(&tanbeta_values, "tanbeta_values");

}
