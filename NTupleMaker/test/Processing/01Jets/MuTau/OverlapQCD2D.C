#include <cmath>
#include <sstream>
#include <iomanip>
#include "TChain.h"
#include "TH1.h"
#include "THStack.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "TList.h"
#include <string>
#include "TObject.h"
#include "TBranch.h"
#include <functional>
#include "TAxis.h"
#include "TChain.h"
#include "TMath.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include "tdrstyle.C"


using namespace std;

//void Impose( TDirectory *ttarget, TList *ssourcelist, string &np_legend , vector<string> titles_ );
void Impose( TDirectory *ttarget, TList *ssourcelist, string &np_legend , vector<string> titles_ ,vector<float> xsecs);

// void MergeRootfile( TDirectory *target, TList *sourcelist );


void OverlapQCD2D()
{

	gROOT->SetStyle ("Plain");
	gStyle->SetPalette (1);
	gStyle->SetTextFont(22) ;
	gStyle->SetTitleFont(22,"xyz") ;
	gStyle->SetLabelFont(22,"xyz") ;


	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
	lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
	lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	int iPeriod = 3;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

	// second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
	// iPos=11 : top-left, left-aligned
	// iPos=33 : top-right, right-aligned
	// iPos=22 : center, centered
	// mode generally : 
	//   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

	//example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
	//  example_plot( iPeriod, 11 );  // left-aligned
	//  example_plot( iPeriod, 33 );  // right-aligned

	//  writeExtraText = false;       // remove Preliminary

	//  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

	//  example_plot( iPeriod, 11 );  // default: left-aligned
	//  example_plot( iPeriod, 22 );  // centered
	//  example_plot( iPeriod, 33 );  // right-aligned  

	vector <string> titles;
	// void MergeRootfile( TDirectory *target, TList *sourcelist );
	//TTrees
	TList *FileList;
	TFile *Target;
	titles.clear();
	int np=1;

	Float_t value=0;
	vector<float> xsecs_;
	ifstream ifs("QCDsamples");
	string line;
	while(std::getline(ifs, line)) // read one line from ifs
	{
		istringstream iss(line); // access line as a stream
		string dataname;
		float XSec;	
		float xs,fact,fact2;
		xs=0;fact=1;fact2=1;
		iss >> dataname >> xs >> fact >> fact2;
		titles.push_back(dataname+".root");
		XSec= xs*fact*fact2;
		cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dataname<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
		xsecs_.push_back(XSec);
	}

	string fout = "QCD_DataDriven_B.root";


	FileList = new TList ();


	for (unsigned int i=0; i <titles.size();++i){
		//string ext=".root";
		cout<<" loading dataset "<<titles[i]<<endl;
		//string file=titles[i]+".root";
		string file=titles[i];
		FileList->Add (TFile::Open (file.c_str()));
	}


	//return;
	Target = TFile::Open (fout.c_str (), "update");

	string np_title = titles[0];
	Impose (Target, FileList, np_title,titles,xsecs_);


	delete FileList;
	delete Target;
}


	void
Impose (TDirectory * target, TList * sourcelist, string & np_title_, vector<string> titles,vector<float> xsecs)
{
	cout << "	" << "========================================================" << endl;
	cout << "	" << "This is a macro to superimpose plots of different root files." << endl;
	cout << "	" << "Only TH1Dobjects are superimposed." << endl;
	cout << "	" << "Target path: " << target->GetPath () << endl;
	TString path ((char *) strstr (target->GetPath (), ":"));
	path.Remove (0, 2);



	float Lumi=1;
	Lumi = 15712.;
	bool norm_=false;
	cout<<titles[0]<<"   "<<titles.size()<<endl;

	//not really useful if plots already weighted to lumi - usefull is plots are in a.u.
	vector <float > lumiweights;
	lumiweights.clear();

	TFile *first_source = (TFile *) sourcelist->First ();
	first_source->cd ("mutau");
	float norm=1.;
	lumiweights.push_back(float(norm));


	// cout<< " for first source file, there where "<<nGen<<" events with xsec "<<xsec<<" SumOfWeights "<<nGen<<"  weight "<<lumiweights[0]<<" weight "<<float(xsec*Lumi/nGen)<<"  norm "<<norm<<endl;

	TDirectory *current_sourcedir = gDirectory;
	//gain time, do not add the objects in the list in memory
	Bool_t status = TH1::AddDirectoryStatus ();
	TH1::AddDirectory (kFALSE);
	// loop over all keys in this directory
	TChain *globChain = 0;
	TIter nextkey (current_sourcedir->GetListOfKeys ());
	//TIter nextkey (((TDirectory *) current_sourcedir->Get ("ana"))->GetListOfKeys ());
	TKey *key, *oldkey = 0;
	int count=0;
	while ((key = (TKey *) nextkey ())) {
		count++;
		//if (count>200) break;
		//keep only the highest cycle number for each key
		//        if (oldkey && !strcmp (oldkey->GetName (), key->GetName ()))
		//            continue;

		// read object from first source file and create a canvas
		first_source->cd ("mutau");
		TObject *obj = key->ReadObj ();
		string nn = obj->GetName();
/*
TH2D *hmet_MT[CutN];
TH2D *hmet_MTsum[CutN];
TH2D *hmet_DZeta[CutN];
TH2D *hmet_MCTb[CutN];
TH2D *hDZeta_MCTsum[CutN];

TH2D *hDZeta_MCTb[CutN];

TH2D *hMTsum_MCTb[CutN];
*/
		if ( string::npos == nn.find("met_")   && string::npos == nn.find("DZeta_M")  && string::npos == nn.find("MTsum_M") && string::npos == nn.find("MTtot_M") && string::npos == nn.find("MT2lester_D") ) continue;



		//if ( string::npos == nn.find("CutFlowUnW")  ) continue;
                //if (string::npos == nn.find("_15") ) flcont=false;
		//if (obj->GetName () != "CutFlow") continue;
		//if (obj->GetName () != "CutFlowUnW" ) continue;
		string nn = obj->GetName();
		//if ( nn != "MET_" ) continue;
		//if (std::string::npos == nn.find("_5")) continue; 
		//cout<<obj->GetName()<<endl;
		TCanvas *c1 = new TCanvas ("c1",obj->GetName (),0,22,600,600);


//		if (obj->IsA ()->InheritsFrom ("TH1")  || obj->IsA ()->InheritsFrom ("TTree")) cout<<" TH1 is not valid "<<endl;

		if (obj->IsA ()->InheritsFrom ("TH2D") ) {


			TH2D* hh[1000];
			TH2D* allRegA;  
			TH2D* allbkg;  
			TH2D *allRegC; 
			TH2D *allRegD;
			TH2D* h1 = (TH2D*) obj;


			TLegend *legend_c1 = new TLegend (0.6, 0.97, 0.98, 0.6);
			legend_c1->SetFillColor(1);
			legend_c1->SetFillStyle(0);
			legend_c1->SetLineColor(0);
			legend_c1->SetTextFont (22);
			legend_c1->SetTextSize (0.04);

			//legend_c1-> SetNColumns(2);     
			h1->Scale(lumiweights[0]);


			// loop over all source files and modify the correspondant
			// histogram to the one pointed to by "h1"
			TFile *nextsource = (TFile *) sourcelist->After (first_source);

			int cl=0;
			cl=1;
			h1->SetStats(000000);
			h1->SetLineWidth(5);

			//	cout<<" You have a cl  "<<cl<<" xsec "<<xsec<<" norm  "<<norm<<endl;
			hh[cl]=h1;
			//hh[1]->GetXaxis()->SetRange(hh[1]->FindFirstBinAbove(0),hh[1]->FindLastBinAbove(5));	

			if (cl==1){	
				allRegA  = (TH2D) hh[1]->Clone();
				allRegA->Reset(); 
				allRegA->Add(hh[1]); 
			}

			string regA, regB,regC,regD;

			string sn="stau";string sdata="Single";
			string qcd="QCD_DataDriven";
			regA="_A.root";
			regC="_C.root";
			regD="_D.root";
			while (nextsource) {

				string fname= nextsource->GetName();

				bool flagg= false;

				if (std::string::npos != fname.find(qcd) || std::string::npos != fname.find(sdata)    ) 	flagg=true;

				// make sure we are at the correct directory level by cd'ing to path
				nextsource->cd("mutau");
				TH1D* eventCountt = (TH1D*)nextsource->Get("mutau/histWeightsH");
				float normm =1.;
				TH1D* hxsecc = (TH1D*)nextsource->Get("mutau/xsec");
				float xsecc = xsecs[cl];


				float nGenn = eventCountt->GetSumOfWeights();


				normm = float(xsecc*Lumi) / float(nGenn)  ;

				if (flagg) { xsecc=1;normm =1.;}	
				lumiweights.push_back(normm);

				TKey *key2 = (TKey *) gDirectory->GetListOfKeys ()->FindObject (h1->GetName ());

				if (key2) {
					cl++;

					TH2D *h2;

					h2 = (TH2D*) key2->ReadObj ();
					h2->SetLineWidth(4);
					//h2->GetXaxis()->SetRange(hh[1]->FindFirstBinAbove(0),hh[1]->FindLastBinAbove(5));	
					h2->Scale(lumiweights[cl-1]);
					h2->SetStats(0);
					hh[cl] = h2;



					if (cl==2){	
						allRegC  = (TH2D*) h2->Clone();
					}
					if (cl==3){	
						allRegD  = (TH2D*) h2->Clone();
					}
					if (cl==4){	
						allbkg  = (TH2D*) h2->Clone();
					}

					if (cl>3)  {

						if (std::string::npos != fname.find("_A") ){ 
							allRegA->Add(h2,-1);
						}

						if (std::string::npos != fname.find("_C") ) {
							allRegC->Add(h2,-1);
							}

						if (std::string::npos != fname.find("_D") ) {
							allRegD->Add(h2,-1);
							}

					}

				}
				nextsource = (TFile *) sourcelist->After (nextsource);
				}				// while ( nextsource )
			}


			if (obj) {



				TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);


				pad1->SetBottomMargin(0.07);
				pad1->SetRightMargin(0.05);
				pad1->SetGridy();

				c1->cd();


				c1->Clear();

				pad1->SetLogy();
				pad1->Draw();


				pad1->cd();
				pad1->Clear();
				char namee[100];

				sprintf(namee,"%s",key->GetName ());

				double nd = allRegD->Integral();
				double nc = allRegC->Integral();
				                                
				allRegA->Scale(double(nc/nd));
			
				string nn = obj->GetName();
			        if ( string::npos != nn.find("CutFlow")  ) {				char f[100];char ff[100];
				if (!norm_)sprintf(f,"Plots/%s.pdf",namee);
				else sprintf(f,"Plots/%s_Norm.",namee);

				if (!norm_)sprintf(ff,"Plots/%s_Log.pdf",namee);
				else sprintf(ff,"Plots/%s_Log_Norm.pdf",namee);

				//if (!norm_)sprintf(f,"Plots/%s.pdf",namee);

				c1->SaveAs (f);
			}

				target->cd ("mutau");
 				
				allRegA->Write();



			}

		} 			// while ( ( TKey *key = (TKey*)nextkey() ) )





		target->SaveSelf (kTRUE);
		target->Write();
		TH1::AddDirectory (status);
		cout << "	" << "========================================================" << endl;
		cout<< " Ending extracting 2D files.... " <<endl;



	}






