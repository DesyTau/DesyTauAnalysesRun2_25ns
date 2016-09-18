#include <cmath>
#include <sstream>
#include <iomanip>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include <string>
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
#include "TPie.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include "tdrstyle.C"
#include "THStack.h"
#include <sstream>

using namespace std;

void Impose( TDirectory *ttarget, TList *ssourcelist, string &np_legend , vector<string> titles_ ,vector<float> xsecs);
//void Impose( TList *ssourcelist, string &np_legend , vector<string> titles_ ,vector<float> xsecs);
void ModifyHist (TH2D* &h, int cl,vector <string> title);
void ModifyHist (TH2D* &h, int cl,vector <string> title, TLegend * & tleg);
void ModifyHist (TH2D* &h, int cl);
void OverFlow (TH1D* &h, int bin);
//void ModifyHist (TH2D* &h, int cl_ ,float & lumi,float & weight,string & title_, bool norm_=false);
//TH1D *Unroll(TH2D *& hist2D,TString histName);
void Unroll(TH1D *&hist,TH2D *& hist2D,char *& histName);

//TH1D * unrolled(TH2D * hist2D,TString name);

// void MergeRootfile( TDirectory *target, TList *sourcelist );

vector <string > signal_names;
TString signalnames[1500];
TString variable;

int mycolor=TColor::GetColor("#ffcc66");
int mycolorvv=TColor::GetColor("#6F2D35");
//int mycolorvv=TColor::GetColor("#FF6633");
int mycolorqcd=TColor::GetColor("#ffccff");
int mycolortt=TColor::GetColor("#9999cc");
//int mycolorttx=TColor::GetColor("#bbccdd");
int mycolorttx=TColor::GetColor("#33CCFF");
int mycolorwjet=TColor::GetColor("#de5a6a");
//int mycolorwjet=TColor::GetColor("#66CC66");
int mycolordyj=TColor::GetColor("#ffcc66");


void Overlap2D()
{

	gROOT->SetStyle ("Plain");
	gStyle->SetPalette (1);
	gStyle->SetTextFont(22) ;
	gStyle->SetTitleFont(22,"xyz") ;
	gStyle->SetLabelFont(22,"xyz") ;
	setTDRStyle();

	vector <string> titles;
	//TTrees
	TList *FileList;
	TFile *Target;
	titles.clear();
	int np=1;
	Float_t value=0;
	vector<float> xsecs_;
	signal_names.clear();
	ifstream ifs("datasets2D_stau-stau");
	string line;
	while(std::getline(ifs, line)) // read one line from ifs
	{
		istringstream iss(line); // access line as a stream
		string dataname;
		float XSec;	
		float xs,fact,fact2;
		xs=0;fact=1;fact2=1;
		iss >> dataname >> xs >> fact >> fact2;
		titles.push_back(dataname+"_B.root");
		XSec= xs*fact*fact2;
		cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dataname<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
		xsecs_.push_back(XSec);
	}

	//string fout = "mergedplotsTau.root";

	FileList = new TList ();


	for (unsigned int i=0; i <titles.size();++i){
		//string ext=".root";
		cout<<" loading dataset "<<titles[i]<<endl;
		//string file=titles[i]+".root";
		string file=titles[i];
		FileList->Add (TFile::Open (file.c_str()));

	}

	string sign_="C1";
	string sign2_="stau";
	string n_;
	//signal_names.clear();
	for (unsigned int k=0; k<titles.size();k++)
	{
		if (   std::string::npos != titles[k].find(sign_) ||  std::string::npos != titles[k].find(sign2_)){
			n_ = titles[k];
			n_.erase(n_.length()-5);
			//string nn_ = n_+"_out.root";
			signal_names.push_back(n_.c_str());
			signalnames[k] = n_.c_str();
		}
	}

	//return;
	//Target = TFile::Open (fout.c_str (), "RECREATE");

	string np_title = titles[0];
	vector <string>  variables;
	variables.push_back("met_MT_17");
	
	variables.push_back("met_MTsum_17");
	variables.push_back("met_MTtot_17");
	variables.push_back("met_DZeta_17");
	variables.push_back("met_MCTb_17");

	variables.push_back("DZeta_MTsum_17");
	variables.push_back("DZeta_MTtot_17");
	variables.push_back("DZeta_MCTb_17");

	variables.push_back("MTsum_MCTb_17");
	variables.push_back("MTtot_MCTb_17");
	variables.push_back("MTsum_MCTb_17");

	variables.push_back("MT2lester_DZeta_17");
	variables.push_back("MT2lester_met_17");


	for (int vr=0;vr<variables.size();vr++){

	variable=variables[vr].c_str();

	cout<< " the variable name  "<<variable<<endl;
	Impose (FileList, np_title,titles,xsecs_, variable);
	}
	delete FileList;
	delete Target;
}


void
//Impose (TDirectory * target, TList * sourcelist, string & np_title_, vector<string> titles,vector<float> xsecs)
Impose (TList * sourcelist, string & np_title_, vector<string> titles,vector<float> xsecs, TString &variable)
{
	cout << "	" << "========================================================" << endl;
	cout << "	" << "This is a macro to superimpose plots of different root files." << endl;
	cout << "	" << "Only TH2Dobjects are superimposed." << endl;
	float Lumi=1;


	Lumi = 15939;
	bool norm_=false;
        int MaxEventsBin = 10;
	cout<<titles[0]<<"   "<<titles.size()<<endl;

	//not really useful if plots already weighted to lumi - usefull is plots are in a.u.
	vector <float > lumiweights;
	lumiweights.clear();
//	for (unsigned int kk=0; kk<signal_names.size();kk++){
//		cout<<" HERE is some signal  ===============================================  "<<signal_names[kk]<<"  "<<signalnames[kk]<<endl;
//	}

	TH2D* allbkg, *htt,*hstop,*hwj,*hdyj,*hrare,*hdib,*hqcd,*httx, *hrest;

	TFile *first_source = (TFile *) sourcelist->First ();
	first_source->cd ("mutau");

	TH1D* eventCount = (TH1D*)first_source->Get("mutau/histWeightsH");
	//TH1D* eventCount = (TH1D*)first_source->Get("mutau/inputEventsH");
	//TH1D* hxsec = (TH1D*)first_source->Get("mutau/xsec");
	float nGen = eventCount->GetSumOfWeights();
	float xsec = 1;//hxsec->GetMean();
	float norm = xsec*Lumi/nGen;

	norm =1;
	lumiweights.push_back(float(norm));


	//cout<< " for first source file, there where "<<nGen<<" events with xsec "<<xsec<<" weight "<<lumiweights[0]<<endl;//" weight "<<float(xsec*Lumi/nGen)<<"  norm "<<norm<<endl;

	TDirectory *current_sourcedir = gDirectory;
	//gain time, do not add the objects in the list in memory
	Bool_t status = TH1::AddDirectoryStatus ();
	TH1::AddDirectory (kFALSE);
	// loop over all keys in this directory
	TChain *globChain = 0;
	TIter nextkey (current_sourcedir->GetListOfKeys ());
	//TIter nextkey (((TDirectory *) current_sourcedir->Get ("ana"))->GetListOfKeys ());
	TKey *key, *oldkey = 0;
	while ((key = (TKey *) nextkey ())) {
	//variable="met_MTsum_16";

	int count=0;
		count++;
		//if (count>20) break;
		//keep only the highest cycle number for each key
		//        if (oldkey && !strcmp (oldkey->GetName (), key->GetName ()))
		//            continue;

		// read object from first source file and create a canvas
		// first_source->cd (path);
		first_source->cd ("mutau");
		TObject *obj = key->ReadObj ();

		//string nn = obj->GetName();
		// if (std::string::npos == nn.find("Cut")) continue;
		//cout<<obj->GetName()<<endl;



		string nn = obj->GetName();
		bool flcont = true;
		bool NormTT = false;
		if (string::npos != nn.find("_9") || string::npos != nn.find("_10") || string::npos != nn.find("_11") || string::npos != nn.find("_12") || string::npos != nn.find("_13") || string::npos != nn.find("_14") || string::npos != nn.find("_15") || string::npos != nn.find("_16") || string::npos != nn.find("_17") || string::npos != nn.find("_18") || string::npos != nn.find("_19")) NormTT = true;
		//if ( string::npos != nn.find("_11") || string::npos != nn.find("_12") || string::npos != nn.find("_13") || string::npos != nn.find("_14")) NormTT = true;
		NormTT=false;
		//if ( string::npos == nn.find("CutFlowUnW")) flcont=false;
		//if (string::npos == nn.find(""+variable) ) flcont=false;
		if (string::npos == nn.find(variable) ) continue;
		//		if (!flcont) continue;
		if (obj->IsA ()->InheritsFrom ("TTree") ) continue;
		//	if (obj->IsA ()->InheritsFrom ("TH1") ) continue;
		if (obj->IsA ()->InheritsFrom ("TH2") ) {

			cout<<"=================================================== OK for variable "<<variable<<endl;


			TH2D* hh[1500];
			TH2D* hsignal[1500];
			TH2D* h1 = (TH2D*) obj;

			ModifyHist (h1,1,Lumi,lumiweights[0],titles[0],norm_);
			TFile *nextsource = (TFile *) sourcelist->After (first_source);

			int cl, countsignal;
			h1->SetStats(000000);
			cl=1;
			countsignal=1;

			hh[cl]=h1;

			THStack *hs = new THStack(h1->GetName(),h1->GetTitle());



			string sn="stau";
			string sdata="Single";
			string sdata2="MuonEG";
			string cc1="C1";



			while (nextsource) {

			string fname= nextsource->GetName();


				bool flagg= false;

				if (std::string::npos != fname.find(sn) || std::string::npos != fname.find(cc1) || std::string::npos != fname.find(sdata)    ) 	flagg=true;

				//if (flagg) cout<<"=============================================================== "<<fname<<endl;
				// make sure we are at the correct directory level by cd'ing to path
				nextsource->cd("mutau");
				TH1D* eventCountt ;
				
				  if  ( std::string::npos == fname.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8") || !NormTT) eventCountt = (TH1D*)nextsource->Get("mutau/histWeightsH");
				  if ( NormTT && std::string::npos != fname.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8") ) eventCountt = (TH1D*)nextsource->Get("mutau/histTopPt");

				
				TH1D* hxsecc = (TH1D*)nextsource->Get("mutau/xsec");
				float xsecc = xsecs[cl];

				float nGenn = eventCountt->GetSumOfWeights();
				float normm = float(xsecc*Lumi) / float(nGenn)  ;
				if (std::string::npos != fname.find("DataDriven")) normm=1.;

				lumiweights.push_back(normm);

				TKey *key2 = (TKey *) gDirectory->GetListOfKeys ()->FindObject (h1->GetName ());

				if (key2) {
					cl++;
					countsignal++;
					TH2D *h2;

					h2 = (TH2D*) key2->ReadObj ();
					ModifyHist (h2, cl,Lumi,lumiweights[cl-1],titles[cl-1],norm_);
					h2->SetStats(0);
					hh[cl] = h2;

					if (cl==2){
						allbkg  = (TH2D*) h2->Clone();
						allbkg->Reset();
						htt  = (TH2D*) h2->Clone();
						httx  = (TH2D*) h2->Clone();
						hstop  = (TH2D*) h2->Clone();
						hwj  = (TH2D*) h2->Clone();
						hdyj  = (TH2D*) h2->Clone();
						hrare  = (TH2D*) h2->Clone();
						hdib  = (TH2D*) h2->Clone();
						hqcd  = (TH2D*) h2->Clone();
						hrest  = (TH2D*) h2->Clone();

						htt->Reset();
						httx->Reset();
						hstop->Reset();
						hdyj->Reset();
						hwj->Reset();
						hrare->Reset();
						hdib->Reset();
						hqcd->Reset();
						hrest->Reset();


					}



					string  hn_ = obj->GetName();
					cout<<"  "<<fname<<endl;

						if (std::string::npos != fname.find(sn) || std::string::npos != fname.find(cc1) || std::string::npos != fname.find(sdata) ||  std::string::npos != fname.find(sdata2)  )
							flagg=true;


						string  title_ = fname;
							//cout<<"  "<<fname<<endl;

			int col = 0;
						if (!flagg)
						{
			if (std::string::npos != title_.find("wJets")|| std::string::npos != title_.find("WJetsToLNu"))  { col=mycolorwjet ; hwj->Add(h2); hwj->SetLineColor(col);}
			if (std::string::npos != title_.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8") || std::string::npos != title_.find("TTPow")) { col= mycolortt;htt->Add(h2); htt->SetLineColor(col) ;}
			if (std::string::npos != title_.find("QCD"))  {col= mycolorqcd;hqcd->Add(h2); hqcd->SetLineColor(col); cout<<" QCD ======================== "<<hqcd->GetSumOfWeights()<<endl;}
			if (std::string::npos != title_.find("DYJets"))  {col= mycolordyj;hdyj->Add(h2); hdyj->SetLineColor(col);}
			if (std::string::npos != title_.find("ST_") || std::string::npos != title_.find("channel") )  {col= mycolortt; hstop->Add(h2);hstop->SetLineColor(col);}
			if ( std::string::npos != title_.find("WW") || std::string::npos != title_.find("ZZ") ||  std::string::npos != title_.find("WZ") || std::string::npos != title_.find("WG") || std::string::npos != title_.find("ZG") ) {col=mycolorvv; hdib->Add(h2); hdib->SetLineColor(col);}
			if ( std::string::npos != title_.find("TTW") || std::string::npos != title_.find("TTZ") || std::string::npos != title_.find("tZq") || std::string::npos != title_.find("TG") || std::string::npos != title_.find("tG") || std::string::npos != title_.find("TTG") || std::string::npos != title_.find("ttW") || std::string::npos != title_.find("ttZ") || std::string::npos != title_.find("tZ") || std::string::npos != title_.find("TTT_") ) {col=mycolorttx ;httx->Add(h2);   httx->SetLineColor(col);}


							//cout<<" will add histogram "<<h2->GetName()<< " for "<<titles[cl-1]<<"  "<<fname<<"  cl  "<<cl<<endl;
							hs->Add(h2);
							allbkg->Add(h2);
						}

				}


				nextsource = (TFile *) sourcelist->After (nextsource);
			}				// while ( nextsource )
		}


		cout<<" Will now extract and save the 1D "<<endl;
		//const char *s1;
		TString s1;
		TH1D * histbkg;
		string n = obj->GetName();
		vector<float> scales;
		bool b_scale=false;
		if (std::string::npos != n.find(variable) && obj)
		{
			TH2D *hsum = ((TH2D*)(hs->GetStack()->Last())); // the "SUM"
			scales.push_back(15939.);
			for (int sc = 0 ; sc<scales.size();sc++){
	

				float scale = float(scales[sc])/float(Lumi);

				//stringstream ss (stringstream::in | stringstream::out);

				//ss << scales[sc];
	
				TString lumistring = "16invfb";
				TString smFilename = "Templates_"+variable+"_"+lumistring+"_mt_stau-stau.root";
				//TString smFilename = "Templates_mt.root";
				TFile *smFile = TFile::Open (smFilename, "recreate");


				smFile->cd();
				//if (sc==0) smFile->mkdir("mt");
				smFile->cd();
				allbkg = hsum;
				allbkg->SetMinimum(0.1);
				allbkg->SetLineColor(kRed);
				allbkg->SetMarkerColor(kRed);
				allbkg->SetFillColor(kRed);
				allbkg->SetFillStyle(3011);
				allbkg->SetMinimum(0.1);
				s1 = "data_obsMC_" + variable;
				//s1 = "data_obsMC";
				allbkg->SetName(s1);
				cout<<" Total Integral before scaling up - "<<allbkg->GetSumOfWeights()<<"  "<<variable<<" data "<<hh[1]->GetSumOfWeights()<<" tt "<<
					htt->GetSumOfWeights()<<" dy "<<hdyj->GetSumOfWeights()<<" wj "<<hwj->GetSumOfWeights()<<" stop "<<hstop->GetSumOfWeights()<<
					" vv "<<hdib->GetSumOfWeights()<<" qcd "<<hqcd->GetSumOfWeights()<<endl;
	if(b_scale)		allbkg->Scale(scale);
				cout<<" Total Integral after scaling up - "<<allbkg->GetSumOfWeights()<<"  "<<variable<<endl;
				allbkg->Write();


				histbkg = Unroll(allbkg,s1,1.);
				s1 = "1D_allbkg_" + variable;
				//s1 = "1D_data_obsMC";//_" + variable+"_"+lumistring;
				histbkg->SetName(s1);
				histbkg->SetMarkerColor(kRed);
				histbkg->SetLineColor(kRed);
				histbkg->Write();
	 if(b_scale){
				htt->Scale(scale);
				hwj->Scale(scale);
				hdyj->Scale(scale);
				hstop->Scale(scale);
				hdib->Scale(scale);
				hqcd->Scale(scale);
				httx->Scale(scale);
				hrest->Scale(scale);
	 }

				htt->SetMinimum(0.1);
				hdyj->SetMinimum(0.1);
				hwj->SetMinimum(0.1);
				httx->SetMinimum(0.1);
				hstop->SetMinimum(0.1);
				hdib->SetMinimum(0.1);
				hqcd->SetMinimum(0.1);

				s1 = "tt_"+variable;
				htt->SetName(s1);
				htt->Write();
				histbkg = Unroll(htt,s1,1.);
				s1 = "1D_tt_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "wj_"+variable;
				hwj->SetName(s1);
				hwj->Write();
				histbkg =Unroll(hwj,s1,1.);
				s1 = "1D_wj_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "dyj_"+variable;
				hdyj->SetName(s1);
				hdyj->Write();
				histbkg =Unroll(hdyj,s1,1.);
				s1 = "1D_dyj_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "sT_"+variable;
				hstop->SetName(s1);
				hstop->Write();
				histbkg =Unroll(hstop,s1,1.);
				s1 = "1D_sT_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "dib_"+variable;
				hdib->SetName(s1);
				hdib->Write();
				histbkg =Unroll(hdib,s1,1.);
				s1 = "1D_dib_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "ttx_"+variable;
				httx->SetName(s1);
				httx->Write();
				histbkg =Unroll(httx,s1,1.);
				s1 = "1D_ttx_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				s1 = "qcd_"+variable;
				hqcd->SetName(s1);
				hqcd->Write();
				histbkg =Unroll(hqcd,s1,1.);
				s1 = "1D_qcd_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();

				//rest of bkg hold all but ttjets and wjets background	
				s1 = "rest_bkg_"+variable;
				//hrest->Add(hdyj,1);
				//hrest->Add(httj,1);
				//hrest->Add(hwj,1);
				hrest->Add(hstop,1);
				hrest->Add(hdib,1);
				hrest->Add(hqcd,1);
				hrest->Add(httx,1);
				hrest->SetName(s1);
				hrest->Write();
				histbkg =Unroll(hrest,s1,1.);
				s1 = "1D_rest_bkg_"+variable;
				histbkg->SetName(s1);
				histbkg->Write();
				

				hh[1]->SetMinimum(0.1);
				hh[1]->SetLineColor(kBlack);
				hh[1]->SetFillColor(kBlack);
				hh[1]->SetMarkerColor(kBlack);
				s1 = "data_obs_"+variable;
				hh[1]->SetName(s1);

				//hh[1]->Scale(scale);
				hh[1]->Write();
				histbkg =Unroll(hh[1],s1,1.);
				//histbkg =Unroll(allbkg,s1,1.);
				s1 = "1D_data_obs_";
				s1 = "data_obs";
				histbkg->SetName(s1);
				histbkg->Write();

				for (unsigned int ij=0;ij<signal_names.size();++ij){
					//cout<<" again  "<<signal_names[ij]<<endl;
					//cout<<" again  "<<hh[ij+2]->GetName()<<"  "<<ij<<"  "<<signal_names[ij].c_str()<<endl;
					TString ss1 = signal_names[ij].c_str();
					s1 = "1D_"+ss1 +"_"+variable;
	if(b_scale)			hh[ij+2]->Scale(scale);
					histbkg = Unroll(hh[ij+2],s1,1.);
					smFile->cd();
					histbkg->SetName(s1);
					histbkg->SetMarkerColor(kBlue);
					histbkg->SetLineColor(kBlue);
					histbkg->Write();

				}
				smFile->SaveSelf (kTRUE);
				smFile->Close();
				delete smFile;
			}//loop over different predictions

		}///find variable





			}

	//delete c1;
	// save modifications to target file
	//target->SaveSelf (kTRUE);
	TH1::AddDirectory (status);
	cout << "	" << "========================================================" << endl;
	cout<< " Ended SuperImpose of files.... " <<endl;




}

void
//ModifyHist (TH1D* &h, int cl_ ,float & lumi,float & weight,string & title_, bool norm_=false,TLegend *& legend)
ModifyHist (TH2D* &h, int cl_ ,float & lumi,float & weight,string & title_, bool norm_=false)
{


	//cout<<"  will now reweight  "<<title_<<"  "<<h->GetName()<<"  "<<lumi<<"  "<<weight<<endl;
	h->Scale(weight);
	h->GetXaxis()->SetTitleOffset(0.8);
	h->GetYaxis()->SetTitleOffset(1.5);


}

TH1D* Unroll(TH2D *& hist2D,char *&histName, float & norm){
	int nBinsX = hist2D->GetNbinsX();
	int nBinsY = hist2D->GetNbinsY();

	int nBins = nBinsX * nBinsY;

	//hist = (TH1D*)unrolled(hist2D,histName);

	//hist2D->Scale(norm);

	hist1D = new TH1D(histName,hist2D->GetTitle(),nBins,double(0.),double(nBins));
	for (int ix=0; ix<nBinsX; ++ix) {
		for (int iy=0; iy<nBinsY; ++iy) {
			int bin = iy*nBinsX + ix + 1;
			double x = hist2D->GetBinContent(ix+1,iy+1);
			double ex = hist2D->GetBinError(ix+1,iy+1);
			hist1D->SetBinContent(bin,x);
			hist1D->SetBinError(bin,ex);
		}
	}

	return hist1D;
	/*
	   for (int iB=1; iB<=nBins; ++iB) {
	   double x = hist->GetBinContent(iB);
	   double e = hist->GetBinError(iB);
	   hist->SetBinContent(iB,norm*x);
	   hist->SetBinError(iB,norm*e);

	   }*/

}

Unroll(TH1D *& hist1D, TH2D *& hist2D,TString histName, float & norm){
	int nBinsX = hist2D->GetNbinsX();
	int nBinsY = hist2D->GetNbinsY();

	int nBins = nBinsX * nBinsY;

	//hist = (TH1D*)unrolled(hist2D,histName);

//	hist2D->Scale(norm);

	cout<<" For unrolling    you have in total  "<<nBins<<endl;


	hist1D = new TH1D(histName,hist2D->GetTitle(),nBins,double(0.),int(nBins));
	for (int ix=0; ix<nBinsX; ++ix) {
		for (int iy=0; iy<nBinsY; ++iy) {
			int bin = iy*nBinsX + ix + 1;
			double x = hist2D->GetBinContent(ix+1,iy+1);
			double ex = hist2D->GetBinError(ix+1,iy+1);
			hist1D->SetBinContent(bin,x);
			hist1D->SetBinError(bin,ex);
		}
	}


	/*
	   for (int iB=1; iB<=nBins; ++iB) {
	   double x = hist->GetBinContent(iB);
	   double e = hist->GetBinError(iB);
	   hist->SetBinContent(iB,norm*x);
	   hist->SetBinError(iB,norm*e);

	   }*/

}


TH1D * unrolled(TH2D * hist2D,	TString name) {


	int nBinsX = hist2D->GetNbinsX();
	int nBinsY = hist2D->GetNbinsY();

	int nBins = nBinsX * nBinsY;

	TH1D * hist1D = new TH1D(name,"1D",nBins,double(0.),double(nBins));
	for (int ix=0; ix<nBinsX; ++ix) {
		for (int iy=0; iy<nBinsY; ++iy) {
			int bin = iy*nBinsX + ix + 1;
			double x = hist2D->GetBinContent(ix+1,iy+1);
			double ex = hist2D->GetBinError(ix+1,iy+1);
			hist1D->SetBinContent(bin,x);
			hist1D->SetBinError(bin,ex);
		}
	}


	return hist1D;

}

void OverFlow(TH1D *& h, int &last_bin){

	int nb = h->GetNbinsX();
	float over_ = h->GetBinContent(last_bin);
	float contlast = 0.;//h->GetBinContent(last_bin);
	for (int b=last_bin; b <= nb+1; b++) {contlast +=h->GetBinContent(b);h->SetBinContent(b,0.);}

	h->SetBinContent(last_bin,0);
	h->SetBinContent(last_bin,contlast);
}

