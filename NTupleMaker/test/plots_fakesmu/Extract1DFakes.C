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
#include "TColor.h"
#include "TPie.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
//#include "tdrstyle.C"
#include "THStack.h"
#include <sstream>

using namespace std;

vector <string > signal_names;
TString signalnames[100];
TString variable;
TString syst_;

const int mycolor=TColor::GetColor("#ffcc66");
const int mycolorvv=TColor::GetColor("#6F2D35");
//int mycolorvv=TColor::GetColor("#FF6633");
int mycolorqcd=TColor::GetColor("#ffccff");
int mycolortt=TColor::GetColor("#9999cc");
//int mycolorttx=TColor::GetColor("#bbccdd");
int mycolorttx=TColor::GetColor("#33CCFF");
int mycolorwj=TColor::GetColor("#de5a6a");
//int mycolorwj=TColor::GetColor("#66CC66");
int mycolordyj=TColor::GetColor("#ffcc66");
int mycolorztt=TColor::GetColor("#58d885");
int mycolorww=TColor::GetColor("#6F2D35");

int col = -1;

void Impose( TDirectory *ttarget, TList *ssourcelist, string &np_legend , vector<string> titles_ ,vector<float> xsecs);
void Impose (TList * sourcel, string & np_title, vector<string> title,vector<float> xsecs, TString &variable, string &syst, string &region);
void ModifyHist (TH1D* &h1, int cl ,float & lumi_,float & weight_,string & title, bool norm=false);
void CheckHist (TH1D* &h1);
void CheckHistZero (TH1D* &h1);
void OverFlow (TH1D* &h, int bin);
void Unroll(TH1D *&hist,TH2D *& hist2D,char *& histName);



void Extract1DFakes(string syst="Nominal", string reg="SR")
{

	gROOT->SetStyle ("Plain");
	gStyle->SetPalette (1);
	gStyle->SetTextFont(22) ;
	gStyle->SetTitleFont(22,"xyz") ;
	gStyle->SetLabelFont(22,"xyz") ;
	//setTDRStyle();

	vector <string> titles;
	//TTrees
	TList *FileList;
	TFile *Target;
	titles.clear();
	int np=1;
	Float_t value=0;
	vector<float> xsecs_;
	signal_names.clear();
	ifstream ifs("datasets2D_C1N2_"+reg);
	string channel="mutau";
	string dirr="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/plots_"+channel+"/";
	dirr="";
	string line;
	while(std::getline(ifs, line)) // read one line from ifs
	{
		istringstream iss(line); // access line as a stream
		string dataname;
		float XSec;	
		float xs,fact,fact2,fact3;
		xs=0;fact=1;fact2=1;fact3=1;
		iss >> dataname >> xs >> fact >> fact2 >> fact3;
		if (std::string::npos != dataname.find("Single") || std::string::npos != dataname.find("MuonEG") || syst=="Nominal" ) {titles.push_back(dataname+"_B.root");}
		else{ titles.push_back(dataname+"_"+syst+"_B.root"); 
		}
		
		//titles.push_back(dirr+dataname+"_"+syst+"_B.root"); 
		//cout<<dataname+"_"+syst+"_B.root"<<endl;

		//titles.push_back(dataname+"_B.root");
		XSec= xs*fact*fact2*fact3;
		//cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dataname<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
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


	string np_title = titles[0];
	vector <string>  variables;
				
	TString textfilenameIn ="bins_"+reg;
			 
	 ifstream ifss(textfilenameIn);
				
	 string linee; string blabel;
	variables.push_back("CutFlowUnWFakeRate_17_0");
	variables.push_back("CutFlowUnWFakeRateJet_17_0");
	 while(std::getline(ifss, linee)) // read one line from ifs
	       	{
             istringstream iss(linee); // access line as a stream
	     int bin_;
	     string bl_;
	     iss >> bin_>> bl_;
					
	     stringstream ss;
	     ss << bin_;
	cout<<" ============================================ will try to push  CutFlowUnWFakeRate_17_"<<ss.str()<<endl;
//	variables.push_back("CutFlowUnWFakeRate_17_"+ss.str());
		}




	syst_= syst.c_str();

	for (int vr=0;vr<variables.size();vr++){

	variable=variables[vr].c_str();

	cout<< " the variable name  "<<variable<<" syst "<<syst<<endl;
	Impose (FileList, np_title,titles,xsecs_, variable,syst,reg);
	}
//	delete FileList;
//	delete Target;
}


void
Impose (TList * sourcelist, string & np_title_, vector<string> titles,vector<float> xsecs, TString &variable, string &syst, string &region)
{
	cout << "	" << "========================================================" << endl;
	cout << "	" << "This is a macro to superimpose plots of different root files." << endl;


	float Lumi = 35864.;
	bool norm_=false;
	bool blind = true;
        int MaxEventsBin = 100;
	//cout<<titles[0]<<"   "<<titles.size()<<endl;

	vector <float > lumiweights;
	lumiweights.clear();

	TH1D* allbkg, *htt,*hstop,*hwj,*hdyj,*hztt,*hrare,*hdib,*hww,*hqcd,*httx;
	THStack *hs;
	TFile *first_source = (TFile *) sourcelist->First ();
	first_source->cd ("mutau");

	TH1D* eventCount = (TH1D*)first_source->Get("mutau/histWeightsH");
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
			TH1D* hh[1500];
			TH1D* hsignal[1500];
	while ((key = (TKey *) nextkey ())) {

	int count=0;
		count++;
		//keep only the highest cycle number for each key
		//        if (oldkey && !strcmp (oldkey->GetName (), key->GetName ()))
		//            continue;

		// read object from first source file and create a canvas
		// first_source->cd (path);
		first_source->cd ("mutau");
		TObject *obj = key->ReadObj ();



		string nn = obj->GetName();
		bool flcont = true;
		bool NormTT = false;
		NormTT=true;
		//if (string::npos == nn.find(variable) ) continue;
		if ( nn != variable ) continue;
		//		if (!flcont) continue;
		if (obj->IsA ()->InheritsFrom ("TTree") ) continue;

		if (obj->IsA ()->InheritsFrom ("TH2") ) continue;
		if (obj->IsA ()->InheritsFrom ("TH3") ) continue;
		if (obj->IsA ()->InheritsFrom ("TH1") ) {

			cout<<"=================================================== OK for variable "<<variable<<" syst "<<syst<<endl;


			TH1D* h1 = (TH1D*) obj;

			ModifyHist (h1,1,Lumi,lumiweights[0],titles[0],norm_);
			TFile *nextsource = (TFile *) sourcelist->After (first_source);

			int cl, countsignal;
			h1->SetStats(000000);
			cl=1;
			countsignal=1;

			hh[cl]=h1;

			hs = new THStack(h1->GetName(),h1->GetTitle());



			string sn="stau";
			string sdata="Single";
			string sdata2="MuonEG";
			string cc1="C1";



			while (nextsource) {

			string fname= nextsource->GetName();


				bool flagg= false;

				if (std::string::npos != fname.find(sn) || std::string::npos != fname.find(cc1) || std::string::npos != fname.find(sdata)    ) 	flagg=true;

				nextsource->cd("mutau");
				TH1D* eventCountt ;
				

				  eventCountt = (TH1D*)nextsource->Get("mutau/histWeightsH");
                                  if (std::string::npos != fname.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") && syst != "TopPtUp" && syst != "TopPtDown") eventCountt = (TH1D*)nextsource->Get("mutau/histTopPt");
                                  if (std::string::npos != fname.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") && syst == "TopPtDown" ) eventCountt = (TH1D*)nextsource->Get("mutau/histWeightsH");
                                  if (std::string::npos != fname.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") &&  syst == "TopPtUp") eventCountt = (TH1D*)nextsource->Get("mutau/histTopPtSq");
				
				float xsecc = xsecs[cl];

				float nGenn = eventCountt->GetSumOfWeights();
				float normm = float(xsecc*Lumi) / float(nGenn)  ;
				if (std::string::npos != fname.find("DataDriven")) normm=1.;

				lumiweights.push_back(normm);

				TKey *key2 = (TKey *) gDirectory->GetListOfKeys ()->FindObject (h1->GetName ());

				if (key2) {
					cl++;
					countsignal++;
					TH1D *h2;

					h2 = (TH1D*) key2->ReadObj ();
					ModifyHist (h2, cl,Lumi,lumiweights[cl-1],titles[cl-1],norm_);
					h2->SetStats(0);
					hh[cl] = h2;

					if (cl==2){
						allbkg  = (TH1D*) h2->Clone("allbkg");
						allbkg->Reset();
						htt  = (TH1D*) h2->Clone();
						httx  = (TH1D*) h2->Clone();
						hstop  = (TH1D*) h2->Clone();
						hwj  = (TH1D*) h2->Clone();
						hdyj  = (TH1D*) h2->Clone();
						hztt  = (TH1D*) h2->Clone();
						hrare  = (TH1D*) h2->Clone();
						hdib  = (TH1D*) h2->Clone();
						hww  = (TH1D*) h2->Clone();
						hqcd  = (TH1D*) h2->Clone();

						htt->Reset();
						httx->Reset();
						hstop->Reset();
						hdyj->Reset();
						hztt->Reset();
						hwj->Reset();
						hrare->Reset();
						hdib->Reset();
						hww->Reset();
						hqcd->Reset();


					}



					string  hn_ = obj->GetName();
					cout<<"  "<<fname<<endl;

						if (std::string::npos != fname.find(sn) || std::string::npos != fname.find(cc1) || std::string::npos != fname.find(sdata) ||  std::string::npos != fname.find(sdata2)  )
							flagg=true;


						string  title_ = fname;

						if (!flagg)
						{
			if (std::string::npos != title_.find("JetsToLNu"))  { col=mycolorwj ; hwj->Add(h2,1); hwj->SetLineColor(col);}
			if (std::string::npos != title_.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") || std::string::npos != title_.find("TTPow")) { col=kBlue;col= mycolortt;htt->Add(h2); htt->SetLineColor(col) ;}
			if (std::string::npos != title_.find("QCD"))  {col= mycolorqcd;hqcd->Add(h2,1); hqcd->SetLineColor(col); cout<<" QCD ======================== "<<hqcd->GetSumOfWeights()<<endl;}
			if (std::string::npos != title_.find("JetsToLL") && std::string::npos == title_.find("isZTT"))  {col= mycolordyj;hdyj->Add(h2); hdyj->SetLineColor(col);}
			if (std::string::npos != title_.find("isZTT"))  {col= mycolorztt;hztt->Add(h2); hztt->SetLineColor(col);}
			if (std::string::npos != title_.find("ST_") || std::string::npos != title_.find("channel") )  {col= mycolortt; hstop->Add(h2);hstop->SetLineColor(col);}
			if (  ( std::string::npos != title_.find("WW") || std::string::npos != title_.find("ZZ") ||  std::string::npos != title_.find("WZ") || std::string::npos != title_.find("WG") || std::string::npos != title_.find("ZG")) &&  std::string::npos == title_.find("WWTo") )  {col=mycolorvv; hdib->Add(h2) ; hdib->SetLineColor(col);}
			if ( std::string::npos != title_.find("WWTo"))  {col=mycolorww; hww->Add(h2); hww->SetLineColor(col);}
			if ( std::string::npos != title_.find("TTW") || std::string::npos != title_.find("TTZ") || std::string::npos != title_.find("tZ") || std::string::npos != title_.find("TG") || std::string::npos != title_.find("tG")  || std::string::npos != title_.find("ttW")  || std::string::npos != title_.find("TTT_") ) {col=mycolorttx ;httx->Add(h2);   httx->SetLineColor(col);}


							hs->Add(h2);
							//allbkg->Add(h2,1);
						}

				}


				nextsource = (TFile *) sourcelist->After (nextsource);
			}				// while ( nextsource )
		}


		cout<<" Will now extract and save the 1D "<<endl;
		//const char *s1;
		TString s1;
		TH1D * histbkg;
		TH1D* htest_tt;
		TH1D* htest_wj;
		TH1D* htest_dyj;
		TH1D* htest_ztt;
		TH1D* htest_dib;
		TH1D* htest_ww;
		TH1D* htest_stop;
		TH1D* htest_qcd;
		TH1D* htest_ttx;
		TH1D* htest_signal;
		TH1D* htest_data;

		string n = obj->GetName();
		vector<float> scales;
		if (std::string::npos != n.find(variable) && obj)
		{
			TH1D *hsum = ((TH1D*)(hs->GetStack()->Last())); // the "SUM"
			scales.push_back(40000.);
			for (int sc = 0 ; sc<scales.size();sc++){
	

				float scale = float(scales[sc])/float(Lumi);

				//stringstream ss (stringstream::in | stringstream::out);
				string var ;
				TString lumistring = "35invfb";
				TString smFilename ; 
				stringstream ssB;
				ssB << MaxEventsBin;
				string str = ssB.str();
				
				//if (MaxEventsBin<10000) lumistring = lumistring+"_"+str+"MaxEvntsBin_";
				if (syst!="Nominal") smFilename =  region+"/Templ_"+variable+"_"+lumistring+"_mt_C1N2_"+region+"_"+syst+".root";
					else smFilename =  region+"/Templ_"+variable+"_"+lumistring+"_mt_C1N2_"+region+".root";

				if (syst!="Nominal") 
					variable=variable+"_"+syst;

				TFile *smFile = TFile::Open (smFilename, "recreate");
				
				TH1D* all_;
						
				all_ = (TH1D*) hsum->Clone("all_");
				smFile->cd();
				smFile->cd();
				allbkg = hsum;
				allbkg->SetMinimum(0.1);
				allbkg->SetLineColor(kRed);
				allbkg->SetMarkerColor(kRed);
				allbkg->SetFillColor(kRed);
				allbkg->SetFillStyle(3011);
				s1 = "data_obsMC_" + variable;
				//s1 = "data_obsMC";
				allbkg->SetName(s1);
				cout<<" Total Integral before scaling up - "<<allbkg->GetSumOfWeights()<<"  "<<variable<<" data "<<hh[1]->GetSumOfWeights()<<" tt "<<
					htt->GetSumOfWeights()<<" dy "<<hdyj->GetSumOfWeights()<<" ztt "<<hztt->GetSumOfWeights()<<" wj "<<hwj->GetSumOfWeights()<<" stop "<<hstop->GetSumOfWeights()<<
					" vv "<<hdib->GetSumOfWeights()+hww->GetSumOfWeights()<<" qcd "<<hqcd->GetSumOfWeights()<<" ttx "<<httx->GetSumOfWeights()<<endl;
				cout<<" Total Integral after scaling up - "<<allbkg->GetSumOfWeights()<<"  "<<variable<<"  all_  "<<all_->GetSumOfWeights()<<endl;
				allbkg->Write();
			
				CheckHist(all_);
				ofstream tfile;
				TString textfilename ="bins_"+region;
				vector <int> keep_bin;keep_bin.clear();
				vector <string> label_bin;label_bin.clear();
			 if (syst=="Nominall"){
				tfile.open(textfilename);
				//tfile.open("bins_"+region);
				for (int nb=1;nb<=all_->GetNbinsX();++nb){

				float bc_ = all_->GetBinContent(nb);
				
				bool SelEvents=false;


				if (region =="SR" && bc_<MaxEventsBin) SelEvents = true;
				else if (region =="SR_Fit" && bc_<MaxEventsBin) SelEvents = true;
				else if (region =="CRA" && bc_>MaxEventsBin) SelEvents = true;
				else if (region =="CRB" && bc_>MaxEventsBin && bc_<1000) SelEvents = true;
				else if (region =="CRC" && bc_>1000) SelEvents = true;
				else if (region =="SR_CR" ) SelEvents = true;
				else if (region =="SR_CR1" && bc_<100) SelEvents = true;

				else SelEvents=false;

				if (SelEvents && bc_>0.) {
			//		cout<<" will keep this bin "<<nb<<"  "<<all_->GetBinContent(nb)<<"  "<<keep_bin.size()<<" for Region  "<<region<<endl;
					keep_bin.push_back(nb);
					stringstream ss;
					ss << nb;
					string str = ss.str();
					label_bin.push_back(str);
					tfile <<nb<<"  "<<nb<<endl;
				}
					}
				tfile.close();
			 }

				TString textfilenameIn ="bins_"+region;
			 if (syst!="NominalL"){
				 ifstream ifs(textfilenameIn);
				
				 string line; string blabel;
			        while(std::getline(ifs, line)) // read one line from ifs
        				{
			                istringstream iss(line); // access line as a stream
			                int bin_;
					string bl_;
			                iss >> bin_>> bl_;
					keep_bin.push_back(bin_);
					stringstream ss;
					ss << bin_;
					string str = ss.str();
					label_bin.push_back(str);
					cout<<" Read in for systematic "<<syst<<" nbin "<<bin_<<"  from  "<<textfilenameIn<<endl;
        				}


				 }

				const int sb_ = keep_bin.size();
				cout<<" in total "<<sb_<<" bins "<<"  "<<htt->GetNbinsX()<<"  "<<htt->Integral()<<endl;

				htest_tt = new TH1D (htt->GetName(),htt->GetTitle(),sb_,1,sb_+1);
				htest_wj = new TH1D (hwj->GetName(),hwj->GetTitle(),sb_,1,sb_+1);
				htest_dyj = new TH1D (hdyj->GetName(),hdyj->GetTitle(),sb_,1,sb_+1);
				htest_ztt = new TH1D (hztt->GetName(),hztt->GetTitle(),sb_,1,sb_+1);
				htest_stop = new TH1D (hstop->GetName(),hstop->GetTitle(),sb_,1,sb_+1);
				htest_dib = new TH1D (hdib->GetName(),hdib->GetTitle(),sb_,1,sb_+1);
				htest_ww = new TH1D (hww->GetName(),hww->GetTitle(),sb_,1,sb_+1);
				htest_qcd = new TH1D (hqcd->GetName(),hqcd->GetTitle(),sb_,1,sb_+1);
				htest_ttx = new TH1D (httx->GetName(),httx->GetTitle(),sb_,1,sb_+1);
				htest_data = new TH1D (hh[1]->GetName(),hh[1]->GetTitle(),sb_,1,sb_+1);

				for (int nbb=0;nbb<keep_bin.size();++nbb){
		/*		htest_tt->SetBinContent(nbb+1,0.);
				htest_wj->SetBinContent(nbb+1,0.);
				htest_dyj->SetBinContent(nbb+1,0.);
				htest_stop->SetBinContent(nbb+1,0.);
				htest_dib->SetBinContent(nbb+1,0.);
				htest_ww->SetBinContent(nbb+1,0.);
				htest_qcd->SetBinContent(nbb+1,0.);
				htest_ttx->SetBinContent(nbb+1,0.);
				htest_data->SetBinContent(nbb+1,0.);
						
				htest_tt->SetBinContent(nbb+1,htt->GetBinContent(keep_bin[nbb]));
				htest_wj->SetBinContent(nbb+1,hwj->GetBinContent(keep_bin[nbb]));
				htest_dyj->SetBinContent(nbb+1,hdyj->GetBinContent(keep_bin[nbb]));
				htest_stop->SetBinContent(nbb+1,hstop->GetBinContent(keep_bin[nbb]));
				htest_dib->SetBinContent(nbb+1,hdib->GetBinContent(keep_bin[nbb]));
				htest_ww->SetBinContent(nbb+1,hww->GetBinContent(keep_bin[nbb]));
				htest_qcd->SetBinContent(nbb+1,hqcd->GetBinContent(keep_bin[nbb]));
				htest_ttx->SetBinContent(nbb+1,httx->GetBinContent(keep_bin[nbb]));
				htest_data->SetBinContent(nbb+1,hh[1]->GetBinContent(keep_bin[nbb]));
		*/	
				TString lab_ = label_bin[nbb].c_str();
				cout<<" will set the label  "<<lab_<<endl;
				htest_tt->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_wj->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_dyj->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_ztt->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_stop->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_dib->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_ww->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_qcd->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_ttx->GetXaxis()->SetBinLabel(nbb+1,lab_);
				htest_data->GetXaxis()->SetBinLabel(nbb+1,lab_);
				
				}

/*
				htt->Reset();
				hwj->Reset();
				hdyj->Reset();
				hstop->Reset();
				hdib->Reset();
				hww->Reset();
				hqcd->Reset();
				httx->Reset();
				hh[1]->Reset();
				htt = (TH1D* )htest_tt->Clone();
				hdyj = (TH1D* )htest_dyj->Clone();
				hwj = (TH1D* )htest_wj->Clone();
				hstop = (TH1D* )htest_stop->Clone();
				hdib = (TH1D* )htest_dib->Clone();
				hww = (TH1D* )htest_ww->Clone();
				hqcd = (TH1D* )htest_qcd->Clone();
				httx = (TH1D* )htest_ttx->Clone();
				hh[1] = (TH1D* )htest_data->Clone();
*/				
				CheckHistZero(htt);
				CheckHistZero(hdyj);
				CheckHistZero(hztt);
				CheckHistZero(hwj);
				CheckHistZero(hstop);
				CheckHistZero(hdib);
				CheckHistZero(hww);
				CheckHistZero(hqcd);
				CheckHistZero(httx);
				//
				htt->SetLineColor(mycolortt);
				hwj->SetLineColor(mycolorwj);
				hdyj->SetLineColor(mycolordyj);
				hztt->SetLineColor(mycolorztt);
				hstop->SetLineColor(mycolortt);
				hdib->SetLineColor(mycolorvv);
				hww->SetLineColor(mycolorww);
				hqcd->SetLineColor(mycolorqcd);
				httx->SetLineColor(mycolortt);
				hh[1]->SetLineColor(kBlack);

				TH1D* allnew_ = (TH1D* )htt->Clone("allnew_");
				allnew_ ->Add(hwj);
				allnew_ ->Add(hdyj);
				allnew_ ->Add(hztt);
				allnew_ ->Add(hstop);
				allnew_ ->Add(hdib);
				allnew_ ->Add(hww);
				allnew_ ->Add(hqcd);
				allnew_ ->Add(httx);

	
				histbkg = allnew_;
				
				//hh[1] = allbkg;
				s1 = "1D_allbkg_" + variable;
				//s1 = "1D_data_obsMC";//_" + variable+"_"+lumistring;
				histbkg->SetName(s1);
				histbkg->SetMarkerColor(kRed);
				histbkg->SetLineColor(kRed);
				histbkg->Write();

				htt->SetMinimum(0.1);
				hdyj->SetMinimum(0.1);
				hztt->SetMinimum(0.1);
				hwj->SetMinimum(0.1);
				httx->SetMinimum(0.1);
				hstop->SetMinimum(0.1);
				hdib->SetMinimum(0.1);
				hww->SetMinimum(0.1);
				hqcd->SetMinimum(0.1);

				s1 = "tt_"+variable;
				htt->SetLineColor(mycolortt);
				htt->SetName(s1);
				htt->Write();
				s1 = "1D_tt_"+variable;
				htt->SetName(s1);
				htt->Write();

				s1 = "wj_"+variable;
				hwj->SetLineColor(mycolorwj);
				hwj->SetName(s1);
				hwj->Write();
				s1 = "1D_wj_"+variable;
				hwj->SetName(s1);
				hwj->Write();

				s1 = "dyj_"+variable;
				hdyj->SetLineColor(mycolordyj);
				hdyj->SetName(s1);
				hdyj->Write();
				s1 = "1D_dyj_"+variable;
				hdyj->SetName(s1);
				hdyj->Write();

				s1 = "ztt_"+variable;
				hztt->SetLineColor(mycolorztt);
				hztt->SetName(s1);
				hztt->Write();
				s1 = "1D_ztt_"+variable;
				hztt->SetName(s1);
				hztt->Write();

				s1 = "sT_"+variable;
				//CheckHistZero(hstop);
				hstop->SetLineColor(mycolortt);
				hstop->SetName(s1);
				hstop->Write();
				s1 = "1D_sT_"+variable;
				hstop->SetName(s1);
				hstop->Write();

				s1 = "dib_"+variable;
				//CheckHistZero(hdib);
				hdib->SetLineColor(mycolorvv);
				hdib->SetName(s1);
				hdib->Write();
				s1 = "1D_dib_"+variable;
				hdib->SetName(s1);
				hdib->Write();

				s1 = "ww_"+variable;
				//CheckHistZero(hww);
				hww->SetLineColor(mycolorvv);
				hww->SetName(s1);
				hww->Write();
				s1 = "1D_ww_"+variable;
				hww->SetName(s1);
				hww->Write();

				s1 = "ttx_"+variable;
				//CheckHistZero(httx);
				httx->SetLineColor(mycolorttx);
				httx->SetName(s1);
				httx->Write();
				s1 = "1D_ttx_"+variable;
				httx->SetName(s1);
				httx->Write();

				s1 = "qcd_"+variable;
				//CheckHistZero(hqcd);
				hqcd->SetLineColor(mycolorqcd);
				hqcd->SetName(s1);
				hqcd->Write();
				s1 = "1D_qcd_"+variable;
				hqcd->SetName(s1);
				hqcd->Write();
				cout<<" Adding qcd ================= "<<hqcd->GetSumOfWeights()<<endl;


				s1 = "data_obs_"+variable;
				hh[1]->SetMinimum(0.1);
				hh[1]->SetLineColor(kBlack);
				hh[1]->SetFillColor(kBlack);
				hh[1]->SetMarkerColor(kBlack);
				hh[1]->SetName(s1);
				hh[1]->Write();
				
	if (blind)		hh[1] = (TH1D*) allnew_->Clone("data_obs");
				s1 = "data_obs";
				hh[1]->SetName(s1);
				hh[1]->Write();

				cout<<" comparing MC to data  , MC "<<allnew_->GetSumOfWeights()<<"  data  "<<hh[1]->GetSumOfWeights()<<" tt "<<htt->GetSumOfWeights()<<"  "<<" wj "<<hwj->GetSumOfWeights()<<" dyj "<<hdyj->GetSumOfWeights()<<" qcd "<<hqcd->GetSumOfWeights()<<" sT "<<hstop->GetSumOfWeights()<<" dib "<<hdib->GetSumOfWeights()<<" ttx "<<httx->GetSumOfWeights()<<endl;

				for (unsigned int ij=0;ij<signal_names.size();++ij){
					//cout<<" again  "<<signal_names[ij]<<endl;
					//cout<<" again  "<<hh[ij+2]->GetName()<<"  "<<ij<<"  "<<signal_names[ij].c_str()<<endl;
					string str = signal_names[ij];
					string pattern = "_"+syst;
					 std::string::size_type i = str.find(pattern);
   					while (i != std::string::npos) {
					str.erase(i, pattern.length());
					     i = str.find(pattern, i);
  					 }



					TString ss1 = str.c_str();
					s1 = "1D_"+ss1 +"_"+variable;
					hh[ij+2]->SetLineColor(kBlue);
				//	histbkg->Reset();
					histbkg = hh[ij+2];
					smFile->cd();
					histbkg->SetName(s1);
					histbkg->SetMarkerColor(kBlue);
					histbkg->SetLineColor(kBlue);
					int nB = hh[ij+2]->GetNbinsX();
					cout<<" signal has  "<<nB<<" and Integral "<<hh[ij+2]->Integral()<<endl;
					//htest_signal = new TH1D (histbkg->GetName(),histbkg->GetTitle(),nB,1,nB+1);
					htest_signal  = (TH1D*) hh[ij+2]->Clone();

				//for (int nbb=0;nbb<keep_bin.size();++nbb){
				/*for (int nbb=0;nbb<htest_signal->GetNbinsX();++nbb){
					htest_signal->SetBinContent(nbb+1,histbkg->GetBinContent(keep_bin[nbb]));
					htest_signal->GetXaxis()->SetBinLabel(nbb+1,label_bin[nbb].c_str());
				cout<<" bin  "<<nbb<<"  "<<histbkg->GetBinContent(keep_bin[nbb])<<endl;
//			float sing =  htest_signal->GetSumOfWeights()/sqrt(hh[1]->GetSumOfWeights());
//				cout<<"  sign "<<sing<<" nb "<<nbb<<endl;
				}*/

				htest_signal->SetName(s1);
				htest_signal->SetMarkerColor(kBlue);
				htest_signal->SetLineColor(kBlue);
				htest_signal->Write();
				}///end in signals


				//smFile->SaveSelf (kTRUE);
				smFile->Write();
				smFile->Close();
				delete smFile;
			}//loop over different predictions

		}///find variable





			}

	// save modifications to target file
	//target->SaveSelf (kTRUE);
	TH1::AddDirectory (status);
	cout << "	" << "========================================================" << endl;
	cout<< " Ended SuperImpose of files.... " <<endl;




}

void
//ModifyHist (TH1D* &h, int cl_ ,float & lumi,float & weight,string & title_, bool norm_=false,TLegend *& legend)
ModifyHist (TH1D* &h, int cl_ ,float & lumi,float & weight,string & title_, bool norm_=false)
{


	h->GetXaxis()->SetTitleOffset(0.8);
	h->GetYaxis()->SetTitleOffset(1.5);
			
	h->Scale(weight);
	int col=cl_;

			if (std::string::npos != title_.find("Data") || std::string::npos != title_.find("Single") || std::string::npos != title_.find("MuonEG")  )  {


				col=kBlack ;
				h->SetLineStyle(1);
				h->SetMarkerStyle(20);
				h->SetMarkerSize(1.2);
				h->SetMarkerColor(col);
				h->SetLineColor(col);
			}


			if ( std::string::npos != title_.find("tau") || std::string::npos != title_.find("C1") || std::string::npos != title_.find("SMS")){

				if (cl_>9) col=30+cl_;
				col=cl_+1;
				h->SetLineStyle(1);
				h->SetMarkerStyle(32);
				h->SetMarkerSize(0.3);
				h->SetMarkerColor(col);
				h->SetLineColor(col);
				//h->SetLineStyle(6);
				h->SetFillColor(0);
			}


			if (std::string::npos != title_.find("JetsToLNu") )  { col=mycolorwj ;}
			if (std::string::npos != title_.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") || std::string::npos != title_.find("TTPow")) { col= mycolortt; }
			if (std::string::npos != title_.find("QCD"))  {col= mycolorqcd; }
			if (std::string::npos != title_.find("JetsToLL") && std::string::npos == title_.find("isZTT"))  {col= mycolordyj;}
			if (std::string::npos != title_.find("isZTT"))  {col= mycolorztt;}
			if (std::string::npos != title_.find("ST_") || std::string::npos != title_.find("channel") )  {col= mycolortt; }

			if (  ( std::string::npos != title_.find("WW") || std::string::npos != title_.find("ZZ") ||  std::string::npos != title_.find("WZ") || std::string::npos != title_.find("WG") || std::string::npos != title_.find("ZG")) &&  std::string::npos == title_.find("WWTo") )  {col=mycolorvv;}
			if ( std::string::npos != title_.find("WWTo"))  {col=mycolorww;}


			if ( std::string::npos != title_.find("TTW") || std::string::npos != title_.find("TTZ") || std::string::npos != title_.find("tZq") || std::string::npos != title_.find("TG") || std::string::npos != title_.find("tG") || std::string::npos != title_.find("TTG") || std::string::npos != title_.find("ttW") || std::string::npos != title_.find("ttZ") || std::string::npos != title_.find("tZ") || std::string::npos != title_.find("TTT_") ) {col=mycolorttx ;}
	
			h->SetLineColor(col);	
			h->SetMarkerColor(col);	
			h->SetFillColor(col);	


			for (int nb=1;nb<=h->GetNbinsX();++nb)
			{
				float bc_ = h->GetBinContent(nb);
	//		if (bc_ > 0. && bc_ <0.001) h->SetBinContent(nb,0.001);
			}


}

void
CheckHist (TH1D* &h)
{


			for (int nb=0;nb<=h->GetNbinsX();++nb)
			{
				float bc_ = h->GetBinContent(nb);
			if (bc_ <0.001) h->SetBinContent(nb,0.001);
			}


}

void
CheckHistZero (TH1D* &h)
{


			for (int nb=0;nb<=h->GetNbinsX();++nb)
			{
				float bc_ = h->GetBinContent(nb);
			if (bc_ <0) h->SetBinContent(nb,0.01);
			}


}



TH1D* Unroll(TH2D *& hist2D,char *&histName, float & norm){
	int nBinsX = hist2D->GetNbinsX();
	int nBinsY = hist2D->GetNbinsY();

	int nBins = nBinsX * nBinsY;

	int bin = 0;
	//hist = (TH1D*)unrolled(hist2D,histName);

	//hist2D->Scale(norm);
//		cout<<" "<<endl;
	TH1D* hist1D = new TH1D(histName,hist2D->GetTitle(),nBins,double(0.),double(nBins));
	for (int ix=0; ix<nBinsX; ++ix) {
		for (int iy=0; iy<nBinsY; ++iy) {
			bin++;
			//int bin = iy*nBinsX + ix + 1;
			double x = hist2D->GetBinContent(ix+1,iy+1);
			double ex = hist2D->GetBinError(ix+1,iy+1);
			hist1D->SetBinContent(bin,x);
			hist1D->SetBinError(bin,ex);
			if (hist1D->GetBinContent(bin)< 0.) {hist1D->SetBinContent(bin,float(0.1));hist1D->SetBinError(bin,float(sqrt(0.1)));}
//			cout<<" Filling now the bin "<<bin<<" x  "<<x<<"  "<<hist2D->GetName()<<endl;
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

