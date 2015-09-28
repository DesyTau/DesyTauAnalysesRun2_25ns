#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/bin/AnalysisMacro_mm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // Run range
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // **** end of configuration

  // file name and tree name
 
  char ff[100];

  sprintf(ff,"%s/%s",argv[3],argv[2]);

  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  string SelectionSign="mumu";

  TString era=argv[3];
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  //TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  TFile * file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());


  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);

  //initialize histograms
  InitAllHisto();

  int nFiles = 0;
  int nEvents = 0;
  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;

  //-------------------------------------
  //from Alexis
  //read the xsection
  Float_t XSec=-1;
  Float_t xs,fact,fact2;

  ifstream ifs("xsecs");
  string line;

  while(std::getline(ifs, line)) // read one line from ifs
    {

      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt;
      iss >> dt >> xs >> fact >> fact2;
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      //if (dt==argv[2]) {
      //if (std::string::npos != dt.find(argv[2])) {
      if (  dt == argv[2]) {
        XSec= xs*fact*fact2;
        cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
      }

    }

  if (XSec<0) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}
  //------------------------------------------
  


  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
 
  //----Attention----//
  //if(XSec!=1) nTotalFiles=20;
  //nTotalFiles=5;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      //------------------------------------------------
      //from Alexis
      //controll of correct weight for data or MD
      bool isData= false;
      bool lumi=false;

      if (XSec == 1)  isData = true;
      if (!isData && XSec !=1 ){ 
	weight *=analysisTree.genweight;
	lumi=true;
	histWeights->Fill(1,weight);
      }
    
      else histWeights->Fill(1);

      //-----------------------------------------------

      //cout << "after fill weight" << endl;

      /*
	if (analysisTree.event_run<RunRangeMin) continue;
	if (analysisTree.event_run>RunRangeMax) continue;
      */

      //cout << "ater range cut" << endl;

      if (isData){
	if (analysisTree.event_run<RunRangeMin) continue;
	if (analysisTree.event_run>RunRangeMax) continue;
      

	if (analysisTree.event_run<RunMin)
	  RunMin = analysisTree.event_run;
      
	if (analysisTree.event_run>RunMax)
	  RunMax = analysisTree.event_run;

	//std::cout << " Run : " << analysisTree.event_run << std::endl;

	bool isNewRun = true;
	if (allRuns.size()>0) {
	  for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	    if (analysisTree.event_run==allRuns.at(iR)) {
	      isNewRun = false;
	      break;
	    }
	  }
	}

	if (isNewRun) 
	  allRuns.push_back(analysisTree.event_run);

   
	std::vector<Period> periods;
    
	std::fstream inputFileStream("temp", std::ios::in);
	for(std::string s; std::getline(inputFileStream, s); )
	  {
	    periods.push_back(Period());
	    std::stringstream ss(s);
	    ss >> periods.back();
	  }
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;

	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {
        
	    if ( num.c_str() ==  a.name ) {
	      //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	      //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }
	
	  }
    
      }

      if (!lumi) continue;
      //if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
      //std::remove("myinputfile");


      /*
      //------------------------------------
      //trigger
      //from Illia
      bool isTriggerMuon = false;
      bool isTriggerElectron = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(MuonHLTName)) {
	  //	  std::cout << MuonHLTName << " : " << it->second << std::endl;
	  if (it->second==1)
	    isTriggerMuon = true;
	}
	if (trigName.Contains(ElectronHLTName)) {
	  //	  std::cout << ElectronHLTName << " : " << it->second << std::endl;
	  if (it->second==1)
            isTriggerElectron = true;
	}
      }

      bool acceptTrig = false;
      if (trigger==0) { 
	if (isTriggerMuon||isTriggerElectron) acceptTrig = true;
      }
      if (trigger==1) {
	if (isTriggerMuon) acceptTrig = true;
      }
      if (trigger==2) {
	if (isTriggerElectron) acceptTrig = true;
      }
      if (trigger==3) {
	if (isTriggerMuon&&!isTriggerElectron) acceptTrig = true;
      }
      if (trigger==4) {
	if (!isTriggerMuon&&isTriggerElectron) acceptTrig = true;
      }
      if (!acceptTrig) continue;
      weightsTriggerH->Fill(0.0,weight);
      //------------------------------------
      */



      //cuts start

      // vertex cuts
      Float_t Met = sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex+analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      FillAllHisto(0,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
      FillAllGeneralHisto(0,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);//no cuts
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      FillAllHisto(1,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
      FillAllGeneralHisto(1,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      FillAllHisto(2,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
      FillAllGeneralHisto(2,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      FillAllHisto(3,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
      FillAllGeneralHisto(3,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);

      
      // electron selection

      // muon selection

      vector<unsigned int> allMuons; allMuons.clear();
      vector<unsigned int> idMuons; idMuons.clear();
      vector<unsigned int> isoMuons; isoMuons.clear();

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	allMuons.push_back(im);
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
        FillAllHisto(4,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	FillAllGeneralHisto(4,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
        FillAllHisto(5,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	FillAllGeneralHisto(5,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
        FillAllHisto(6,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	FillAllGeneralHisto(6,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
        FillAllHisto(7,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	FillAllGeneralHisto(7,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
	if (!analysisTree.muon_isMedium[im]) continue;
        FillAllHisto(8,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	FillAllGeneralHisto(8,weight,analysisTree.muon_count,analysisTree.muon_pt,analysisTree.muon_eta,analysisTree.muon_dxy,analysisTree.muon_dz);
	idMuons.push_back(im);
	float absIso = analysisTree.muon_chargedHadIso[im];
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut) continue;
	isoMuons.push_back(im);
      }

      //std::cout << "allMuons : " << allMuons.size() << std::endl;
      //std::cout << "idMuons  : " << idMuons.size() << std::endl;
      //std::cout << "isoMuons : " << isoMuons.size() << std::endl;

      //      continue;

      bool isAllMuonsPair = false;
      if (allMuons.size()>1) {
	//	std::cout << "allMuons : " << allMuons.size() << std::endl;
	for (unsigned int im1=0; im1<allMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<allMuons.size(); ++im2) {
	    unsigned int index1 = allMuons[im1];
	    unsigned int index2 = allMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {//charge controll
	      isAllMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      h_tot_mass->Fill(mass,weight);
	      FillHistosMuon1(0,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(0,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      JPsiMassAllMuonsH->Fill(mass,weight);
	      YpsilonMassAllMuonsH->Fill(mass,weight);
	      ZMassAllMuonsH->Fill(mass,weight);
              FillAllHisto(9,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR<dRleptonsCut) continue;
	      JPsiMassAllMuonsDRCutH->Fill(mass,weight);
	      YpsilonMassAllMuonsDRCutH->Fill(mass,weight);
	      ZMassAllMuonsDRCutH->Fill(mass,weight);
              FillAllHisto(10,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	      FillHistosMuon1(1,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(1,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	    }
	  }
	}
	//	std::cout << "allMuons : " << allMuons.size() << "  :  OK" << std::endl;
      }

      //      continue;

      bool isIdMuonsPair = false;
      if (idMuons.size()>1) {
	for (unsigned int im1=0; im1<idMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<idMuons.size(); ++im2) {
	    unsigned int index1 = idMuons[im1];
	    unsigned int index2 = idMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {
	      isIdMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      //	      std::cout << "Mass = " << mass << std::endl;
	      FillHistosMuon1(2,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(2,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      JPsiMassIdMuonsH->Fill(mass,weight);
	      YpsilonMassIdMuonsH->Fill(mass,weight);
	      ZMassIdMuonsH->Fill(mass,weight);
              FillAllHisto(11,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR<dRleptonsCut) continue;
	      JPsiMassIdMuonsDRCutH->Fill(mass,weight);
	      YpsilonMassIdMuonsDRCutH->Fill(mass,weight);
	      ZMassIdMuonsDRCutH->Fill(mass,weight);
              FillAllHisto(12,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,Met);
	      FillHistosMuon1(3,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(3,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	    }
	  }
	}
      }

      double pt1 = 0.;
      double pt2 = 0.;
      double temp = 0.;
      int indx1 = -1;
      int indx2 = -1;
      int itemp = -1;

      //vector<double> pt_1muon; pt_1muon.clear();
      //vector<double> pt_2muon; pt_2muon.clear();
      bool isIsoMuonsPair = false;
      if (isoMuons.size()>1) {
	for (unsigned int im1=0; im1<isoMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
	    unsigned int index1 = isoMuons[im1];
	    unsigned int index2 = isoMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {
	      //pt_1muon.push_back(analysisTree.muon_pt[index1]);
	      //pt_2muon.push_back(analysisTree.muon_pt[index2]);
	      if(analysisTree.muon_pt[index1]>pt1 && analysisTree.muon_pt[index1]!=pt2){
		pt1=analysisTree.muon_pt[index1];
		indx1=index1;
	      }
	      if(analysisTree.muon_pt[index2]>pt2 && analysisTree.muon_pt[index2]!=pt1){
		pt2=analysisTree.muon_pt[index2];
		indx2=index2;
	      }
	      if(pt1<pt2){
		temp=pt1;
		pt1=pt2;
		pt2=temp;
		itemp=indx1;
		indx1=indx2;
		indx2=itemp;
	      }
	      isIsoMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      ZMassIsoMuonsH->Fill(mass,weight);
	      FillAllHisto(13,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,PtoPt(analysisTree.pfmet_ex,analysisTree.pfmet_ey));
	      FillHistosMuon1(4,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(4,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR<dRleptonsCut) continue;
	      ZMassIsoMuonsDRCutH->Fill(mass,weight);
	      FillAllHisto(14,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,PtoPt(analysisTree.pfmet_ex,analysisTree.pfmet_ey));
	      h_deltaR->Fill(dR,weight);
	      h_tot_pt1pt2->Fill(analysisTree.muon_pt[index1]+analysisTree.muon_pt[index2]);
	      h_tot_m1m2->Fill(mass,weight);
	      MassvsPt->Fill(mass,float(analysisTree.muon_pt[index1]+analysisTree.muon_pt[index2]));
	      double tot_ptjet=0;
	      for(UInt_t k=0;k<analysisTree.pfjet_count;k++){
		tot_ptjet=tot_ptjet+analysisTree.pfjet_pt[k];
	      }
	      h_Ht->Fill(tot_ptjet,weight);
	      FillHistosMuon1(5,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
	      FillHistosMuon2(5,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      //-------------------------
	      //Z mass window
	      if((70.<mass) && (mass<110.)){
		FillAllHisto(15,weight,analysisTree.primvertex_count,analysisTree.primvertex_z,analysisTree.pfjet_count,PtoPt(analysisTree.pfmet_ex,analysisTree.pfmet_ey));
		FillHistosMuon1(6,weight,analysisTree.muon_pt[index1],analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
		FillHistosMuon2(6,weight,analysisTree.muon_pt[index2],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
		h_tot_m1m2_Z->Fill(mass,weight);
		h_tot_pt1pt2_Z->Fill(analysisTree.muon_pt[index1]+analysisTree.muon_pt[index2]);
	      }
	      //-------------------------
	    }
	  }
	}
      }
      /*
      h_massTo12->Fill();
      h_massUpto12->Fill();
      */

      //h_coupleHPt->Fill();

      /*
      double temp1;
      double temp2;
      for(unsigned int j=0;j<pt_1muon.size();j++){
	for(unsigned int k=0;k<pt_1muon.size();k++){
	  if(pt_1muon.at(k)<pt_1muon.at(k+1)){
	    temp1=pt_1muon.at(k);
	    pt_1muon.at(k)=pt_1muon.at(k+1);
	    pt_1muon.at(k+1)=temp1;
	    temp2=pt_2muon.at(k);
	    pt_2muon.at(k)=pt_2muon.at(k+1);
	    pt_2muon.at(k+1)=temp2;
	  }
	  if(fabs(pt_1muon.at(k)-pt_1muon.at(k+1))<1.e-5){
	    if(pt_2muon.at(k)<pt_2muon.at(k+1)){
	      temp2=pt_2muon.at(k);
	      pt_2muon.at(k)=pt_2muon.at(k+1);
	      pt_2muon.at(k+1)=temp2;
	    }
	  }
	}
      }
      */


      /*
	      //jets cuts---------------------------
	      double ptjetsCut=30.;//GeV
	      double etajetsCut=4.7;
	      double ptjets_tot=0;
	      double etajets_tot=0;
	      for(UInt_t j=0;j<pfjet_count<++){
		ptjets_tot=ptjets_tot+pfjet_pt[j];
		etajets_tot=etajets_tot+pfjet_eta[j];
	      }
	      if(ptjets_tot>ptjetsCut && etajets_tot<etajetsCut)h_Njets[15]->Fill(pfjet_count,weight);
	      //------------------------------------
	      */

      if (isAllMuonsPair) selEventsAllMuons++;
      if (isIdMuonsPair)  selEventsIdMuons++;
      if (isIsoMuonsPair) selEventsIsoMuons++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                     = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events (muon pairs)     = " << selEventsAllMuons << std::endl;
  std::cout << "Total number of selected events (id muon pairs)  = " << selEventsIdMuons << std::endl;
  std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  //cout << "weight used:" << weight << std::endl;

  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;

  file->cd(SelectionSign.c_str());

  WriteAllHisto();
  h_tot_m1m2->Write();
  h_tot_pt1pt2->Write();
  //h_primvert->Write();
  h_Ht->Write();
  h_deltaR->Write();
  h_tot_mass->Write();

  h_tot_m1m2_Z->Write();
  h_tot_pt1pt2_Z->Write();

  //---------------------
  //from Alexis
  hxsec->Fill(XSec);
  hxsec->Write();
  histWeights->Write();
  //----------------------

  ZMassAllMuonsH->Write();
  ZMassAllMuonsDRCutH->Write();
  ZMassIdMuonsH->Write();
  ZMassIdMuonsDRCutH->Write();
  ZMassIsoMuonsH->Write();
  ZMassIsoMuonsDRCutH->Write();
  MassvsPt->Write();
  //  file->Write();
  file->Close();
  delete file;
  
}



