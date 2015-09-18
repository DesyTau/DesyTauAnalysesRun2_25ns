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
#include "TChain.h"
#include "TMath.h"
#include "Riostream.h"

#include "TRandom.h"
#include "TRandom.h"

#include "AnalysisMacro.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"

using namespace std;





int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);
  string SelectionSign="muel";

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");
  
  // kinematic cuts on Jets
  const float etaJetCut   = cfg.get<float>("etaJetCut");
  const float ptJetCut   = cfg.get<float>("ptJetCut");
  
  const Float_t dRleptonsCut   = cfg.get<Float_t>("dRleptonsCut");
  const Float_t dZetaCut       = cfg.get<Float_t>("dZetaCut");
  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
 





	// tau
  const Float_t taupt    = cfg.get<Float_t>("taupt");
  const Float_t taueta    = cfg.get<Float_t>("taueta");
  const Float_t decayModeFinding    = cfg.get<Float_t>("decayModeFinding");
  const Float_t   decayModeFindingNewDMs  = cfg.get<Float_t>("decayModeFindingNewDMs");
  const Float_t   againstElectronVLooseMVA5  = cfg.get<Float_t>("againstElectronVLooseMVA5");
  const Float_t   againstMuonTight3  = cfg.get<Float_t>("againstMuonTight3");
  const Float_t   vertexz =  cfg.get<Float_t>("vertexz");
  const Float_t   byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<Float_t>("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  

  const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
  const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");


  const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
  const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");


  // topological cuts
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
  
  const float Lumi   = cfg.get<float>("Lumi");
  const Float_t bTag   = cfg.get<Float_t>("bTag");
  const Float_t metcut = cfg.get<Float_t>("metcut");

  Bool_t          os;
  Bool_t          dilepton_veto;
  Bool_t          extraelec_veto;
  Bool_t          extramuon_veto;


  TString LowPtLegElectron(lowPtLegElectron);
  TString HighPtLegElectron(highPtLegElectron);
  
  TString LowPtLegMuon(lowPtLegMuon);
  TString HighPtLegMuon(highPtLegMuon);
  // **** end of configuration
 
   //TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   //dir.ReplaceAll("basic.C","");
   //dir.ReplaceAll("/./","/");
   //ifstream in;
 CutList.clear();
 CutList.push_back("No cut");
 CutList.push_back("2l dR > "+to_string(dRleptonsCut));
 CutList.push_back("Trigger");
 CutList.push_back("3rd lept-Veto");
 CutList.push_back("b-Veto ");
  CutList.push_back("Zmass cut");
 //CutList.push_back("lep SumpT> 0");
 CutList.push_back("MET $>$ 50");
 CutList.push_back("MET $>$ 100");
 CutList.push_back("dPhi > 1");



  int CutNumb = int(CutList.size());
  xs=1;fact=1;fact2=1;
 
  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;
        
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

  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
	char ff[100];

	
	sprintf(ff,"%s/%s",argv[3],argv[2]);

  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  
  TString era=argv[3];
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());



  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
  int iCut=0;
    double CFCounter[CutNumb];
    double statUnc[CutNumb];
    int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; i++){
          CFCounter[i] = 0;
         iCFCounter[i] = 0;
         statUnc[i] =0;
        }

  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
 
  SetupHists(CutNumb); 
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
    //NE=1000;

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
     //numberOfEntries = 1000;
    
     std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;

     
     numberOfEntries=100;

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
     
      analysisTree.GetEntry(iEntry);
      nEvents++;
 
      iCut = 0;
      
      
      Float_t weight = 1;
      isData= false;
      bool lumi=false;

      if (XSec == 1)  isData = true;
      if (!isData && XSec !=1 )  { 

     	      histWeights->Fill(1,weight); 
	      weight *=analysisTree.genweight;   
	      lumi=true;
      } 
     
     if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
     
         
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
      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();
	mu_index=-1;
	tau_index=-1;
	el_index=-1;
      Float_t MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      
      METV.SetPx(analysisTree.pfmet_ex);	      
      METV.SetPy(analysisTree.pfmet_ey);
 
 
      for (unsigned int ijj = 0; ijj<analysisTree.pfjet_count; ++ijj) {
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[ijj], analysisTree.pfjet_py[ijj], analysisTree.pfjet_pz[ijj], analysisTree.pfjet_e[ijj]);
	JetsMV.push_back(JetsV);
      } 


      for (unsigned int imm = 0; imm<analysisTree.muon_count; ++imm) {
	MuV.SetPtEtaPhiM(analysisTree.muon_pt[imm], analysisTree.muon_eta[imm], analysisTree.muon_phi[imm], muonMass);
	MuMV.push_back(MuV);
	//mu_index=0;
      }

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
	ElMV.push_back(ElV);
//	el_index=0;
      }
   
      for (unsigned int itt = 0; itt<analysisTree.tau_count; ++itt) {
	TauV.SetPtEtaPhiM(analysisTree.tau_pt[itt], analysisTree.tau_eta[itt], analysisTree.tau_phi[itt], tauMass);
	TauMV.push_back(TauV);
//	tau_index=0;
      }


      
     // vector <string> ss; ss.push_back(SelectionSign.c_str());
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;






     // hnJets[->Fill(pfjet_count);

      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      // int indexW  = -1;
      // int indexNu = -1; 
      // int indexMu = -1;
      // int indexE  = -1;
      // int nGenMuons = 0;
      // int nGenElectrons = 0;
      // for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      // 	float pxGen = analysisTree.genparticles_px[igen];
      // 	float pyGen = analysisTree.genparticles_py[igen];
      // 	float pzGen = analysisTree.genparticles_pz[igen];
      // 	float etaGen = PtoEta(pxGen,pyGen,pzGen);
      // 	float ptGen  = PtoPt(pxGen,pyGen);

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
      // 	  indexW = igen;
      // 	if ((fabs(analysisTree.genparticles_pdgid[igen])==12 
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==14
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==16) 
      // 	    && analysisTree.genparticles_info[igen]== (1<<1) )
      // 	  indexNu = igen;

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==13) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexMu = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenMuons++;
      // 	  }
      // 	}
      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==11) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexE = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenElectrons++;
      // 	  }
      // 	}
      // }

      // trigger selection
 
      //selecTable.Fill(1,0, weight );      
      bool trigAccept = false;

      unsigned int nLowPtLegElectron = 0;
      bool isLowPtLegElectron = false;
      
      unsigned int nHighPtLegElectron = 0;
      bool isHighPtLegElectron = false;

      unsigned int nLowPtLegMuon = 0;
      bool isLowPtLegMuon = false;
      
      unsigned int nHighPtLegMuon = 0;
      bool isHighPtLegMuon = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //      std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==LowPtLegElectron) {
	  nLowPtLegElectron = i;
	  isLowPtLegElectron = true;
	}
	if (HLTFilter==HighPtLegElectron) {
	  nHighPtLegElectron = i;
	  isHighPtLegElectron = true;
	}
	if (HLTFilter==LowPtLegMuon) {
	  nLowPtLegMuon = i;
	  isLowPtLegMuon = true;
	}
	if (HLTFilter==HighPtLegMuon) {
	  nHighPtLegMuon = i;
	  isHighPtLegMuon = true;
	}
      }
      if (!isLowPtLegElectron) {
	std::cout << "HLT filter " << LowPtLegElectron << " not found" << std::endl;
	return(-1);
      }
      if (!isHighPtLegElectron) {
	std::cout << "HLT filter " << HighPtLegElectron << " not found" << std::endl;
	return(-1);
      }
      if (!isLowPtLegMuon) {
	std::cout << "HLT filter " << LowPtLegMuon << " not found" << std::endl;
	return(-1);
      }
      if (!isHighPtLegMuon) {
	std::cout << "HLT filter " << HighPtLegMuon << " not found" << std::endl;
	return(-1);
      }
      //      std::cout << "LowPtE  : " << LowPtLegElectron << " : " << nLowPtLegElectron << std::endl;
      //      std::cout << "HighPtE : " << HighPtLegElectron << " : " << nHighPtLegElectron << std::endl;
      //      std::cout << "LowPtM  : " << LowPtLegMuon << " : " << nLowPtLegMuon << std::endl;
      //      std::cout << "HighPtM : " << HighPtLegMuon << " : " << nHighPtLegMuon << std::endl;
      //      std::cout << std::endl;
      //      continue;

/*
      MuMV.clear();
      ElMV.clear();
      TauMV.clear();
      LeptMV.clear();
  */
	// electron selection
      Float_t isoElMin = 9999;
      bool el_iso=false;
      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	electronPtAllH->Fill(analysisTree.electron_pt[ie],weight);
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ie];
	if (relIso>isoElectronHighCut) continue;
	if (relIso<isoElectronLowCut) continue;
	bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[ie],
						analysisTree.electron_mva_id_nontrigPhys14[ie]);
	if (!electronMvaId&&applyElectronId) continue;

	    if (relIso<isoElMin) {
	      isoElMin  = relIso;
	      el_index = ie;
	      el_iso=true;
	      electrons.push_back(ie);
	      ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.muon_phi[ie], electronMass);
	      LeptMV.push_back(ElV);
        //      ElMV.push_back(ElV);
	    }

	    if (relIso!=0 && relIso==isoElMin && ie != el_index) {
             analysisTree.electron_pt[ie] > analysisTree.electron_pt[el_index] ? el_index = ie : el_index = el_index;
	    cout<<" found a pair  " <<relIso <<"  "<<el_index<<"  "<<ie<<endl;
	  }

      }
      if (electrons.size()==0 || !el_iso) continue;

      Float_t isoMuMin = 9999;
      bool mu_iso=false;
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	muonPtAllH->Fill(analysisTree.muon_pt[im],weight);
	if (analysisTree.muon_pt[im]<ptMuonHighCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	Float_t neutralIso = 
	  analysisTree.muon_neutralHadIso[im] + 
	  analysisTree.muon_photonIso[im] - 
	  0.5*analysisTree.muon_puIso[im];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = analysisTree.muon_chargedHadIso[im] + neutralIso;
	Float_t relIso = absIso/analysisTree.muon_pt[im];
        hmu_relISO[1]->Fill(relIso,weight);
	if (relIso>isoMuonHighCut) continue;
	if (relIso<isoMuonLowCut) continue;

	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;

	    if (relIso<isoMuMin) {
	      isoMuMin  = relIso;
	      mu_index = im;
	      mu_iso=true;
	      muons.push_back(im);
	      MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
	      LeptMV.push_back(MuV);
           //   MuMV.push_back(MuV);
	    }

	    if (relIso!=0 && relIso==isoMuMin && im != mu_index) {
             analysisTree.muon_pt[im] > analysisTree.muon_pt[mu_index] ? mu_index = im : mu_index = mu_index;
	    cout<<" found a pair  " <<relIso <<"  "<<mu_index<<"  "<<im<<endl;
	  }
   
      }
      if (muons.size()==0 || !mu_iso) continue;




      ///sort leptons vector
      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      if (LeptMV.size() == 0 ) continue; 
      
      float dR = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			    analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);

      if (dR<dRleptonsCut) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

	//for (unsigned int j=0;j<LeptMV.size();++j) cout<<" j "<<j<<"  "<<LeptMV.at(j).Pt()<<endl;
	//cout<<""<<endl;
      ////////jets cleaning 
 
	bool isMu23 = false;
	bool isMu8 = false;
	bool isEle23 = false;
	bool isEle12 = false;

      for (unsigned int im=0; im<muons.size(); ++im) {
	isMu23 = false;
	isMu8 = false;
	      unsigned int mIndex  = muons.at(im);
	
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
	float photonIsoMu = analysisTree.muon_photonIso[mIndex];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
	float puIsoMu = analysisTree.muon_puIso[mIndex];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
	}
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];
	
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]) { // Mu23 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMu23 = true;
	    }
	  }
	  if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]) { // Mu8 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMu8 = true;
	    }
	  }
	}

	if ((!isMu23) && (!isMu8)) continue;



	for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	isEle23 = false;
	isEle12 = false ;

	  unsigned int eIndex = electrons.at(ie);

	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;


	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][nHighPtLegElectron]) { // Ele23 Leg
	      float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		isEle23 = true;
	      }
	    }
	  }

	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][nLowPtLegElectron]) { // Ele12 Leg
	      float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		isEle12 = true;
	      }
	    }
	  }

	}//ele
	  trigAccept = (isMu23&&isEle12) || (isMu8&&isEle23);
	  //	  std::cout << "Trigger match = " << trigMatch << std::endl;
      }//muons
        if (!trigAccept) continue;
	  
	  bool isKinematicMatch = false;
	  if (isMu23&&isEle12) {
	    if (analysisTree.muon_pt[mu_index]>ptMuonHighCut&&analysisTree.electron_pt[el_index]>ptElectronLowCut)
	      isKinematicMatch = true;
	  }
	  if (isMu8&&isEle23) {
            if (analysisTree.muon_pt[mu_index]>ptMuonLowCut&&analysisTree.electron_pt[el_index]>ptElectronHighCut)
              isKinematicMatch = true;
          }
	  if (!isKinematicMatch) continue;

          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;


      os = (analysisTree.muon_charge[mu_index]*analysisTree.electron_charge[el_index]) < 0;



      bool ThirdLeptVeto=false;

      if (doThirdLeptVeto){
  	if (analysisTree.electron_count>0) {
	  for (unsigned int iev = 0; iev<analysisTree.electron_count; ++iev) {


	    Float_t neutralIsoV = analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev] -  4*TMath::Pi()*(0.3*0.3)*analysisTree.rho;
	    Float_t IsoWithEA =  analysisTree.electron_r03_sumChargedHadronPt[iev] + TMath::Max(Float_t(0), neutralIsoV);
	    Float_t relIsoV = IsoWithEA/analysisTree.electron_pt[iev];

	     bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[iev], analysisTree.electron_superclusterEta[iev], analysisTree.electron_mva_id_nontrigPhys14[iev]);


	    if ( analysisTree.electron_pt[iev] > 10 &&  fabs(analysisTree.electron_eta[iev]) < 2.5 && fabs(analysisTree.electron_dxy[iev])<0.045
		&& fabs(analysisTree.electron_dz[iev]) < 0.2 && relIsoV< 0.3 && electronMvaId ) ThirdLeptVeto=true;

	  }
	}


     	if (analysisTree.muon_count>0){
	  for (unsigned int imvv = 0; imvv<analysisTree.muon_count; ++imvv) {

//       if ( imvv != mu_index  && analysisTree.muon_charge[imvv] != analysisTree.muon_charge[mu_index] ){
	    Float_t neutralIso = 
	      analysisTree.muon_neutralHadIso[imvv] + 
	      analysisTree.muon_photonIso[imvv] - 
	      0.5*analysisTree.muon_puIso[imvv];
	    neutralIso = TMath::Max(Float_t(0),neutralIso); 
	    Float_t absIso = analysisTree.muon_chargedHadIso[imvv] + neutralIso;
	    Float_t relIso = absIso/analysisTree.muon_pt[imvv];
	    if ( imvv != mu_index &&  analysisTree.muon_isMedium[imvv] &&  analysisTree.muon_pt[imvv]> 10 &&  fabs(analysisTree.muon_eta[imvv])< 2.4 && fabs(analysisTree.muon_dxy[imvv])<0.045 
		 && fabs(analysisTree.muon_dz[imvv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imvv]) ) ThirdLeptVeto=true;
	 	 }
		}
	}
  
    	if (ThirdLeptVeto) continue;



/*
      // selecting muon and electron pair (OS or SS);
      float ptScalarSum = -1;
      float dRleptons = -1;
      int electronIndex = -1;
      int muonIndex = -1;
      



      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	int eIndex = electrons[ie];
	    
	if ( analysisTree.electron_pt[eIndex]< ptElectronLowCut ) ElMV.erase(ElMV.begin()+ie);
	
			for (unsigned int im=0; im<muons.size(); ++im) {
	  int mIndex = muons[im];
	  float qProd = analysisTree.electron_charge[eIndex]*analysisTree.muon_charge[mIndex];
	
	  if ( analysisTree.muon_pt[mIndex]<ptMuonLowCut ) MuMV.erase(MuMV.begin()+im);

	  //// This is for the SS to be true
	  if (SelectionSign == "SS") {
	  if (oppositeSign || qProd<0) continue;
	  //if (!oppositeSign && qProd<0) continue;
	  
	  }
	  if (SelectionSign == "OS") {
	  if (!oppositeSign || qProd>0) continue;
	  //if (oppositeSign &&  qProd<0) continue;
	  }
           
	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;

		
           Float_t leadLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] > analysisTree.muon_pt[mIndex]	: analysisTree.muon_pt[mIndex];
           Float_t trailLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] < analysisTree.muon_pt[mIndex] : analysisTree.muon_pt[mIndex];

 
	   //cout<<"  Leading "<<leadLep<<" trailing "<<trailLep<<endl;
	   //cout<<"analysisTree.electron_pt[eIndex] "<<analysisTree.electron_pt[eIndex]<<"  nalysisTree.muon_pt[mIndex] "<<analysisTree.muon_pt[mIndex]<<endl;
           //cout<<""<<endl;

	  bool kinematics = 
            (analysisTree.electron_pt[eIndex]> ptElectronLowCut && analysisTree.muon_pt[mIndex]>ptMuonHighCut) || 
	    (analysisTree.electron_pt[eIndex]>ptElectronHighCut && analysisTree.muon_pt[mIndex]>ptMuonLowCut);

	  if (!kinematics) continue;

	  float sumPt = analysisTree.electron_pt[eIndex] + analysisTree.muon_pt[mIndex];
	  if (sumPt>ptScalarSum) {
	    ptScalarSum = sumPt;
	    dRleptons = dR;
	    electronIndex = ie;
	    muonIndex = im;
	  	}
	  
		}

      }//electrons

      */
          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	
      Float_t DRmax=0.4;
      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;
      
      
      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      for (unsigned int il = 0; il<LeptMV.size(); ++il) {
      
	 for (unsigned int ij = 0; ij<JetsMV.size(); ++ij) {
        
		 if(fabs(JetsMV.at(ij).Eta())>etaJetCut) continue;
                 if(fabs(JetsMV.at(ij).Pt())<ptJetCut) continue;
      
       Float_t Dr= deltaR(LeptMV.at(il).Eta(), LeptMV.at(il).Phi(),JetsMV.at(ij).Eta(),JetsMV.at(ij).Phi());

     if (  Dr  < DRmax) {
	     
	     JetsMV.erase (JetsMV.begin()+ij);
    		 }	
		       
	 }
      }

      bool btagged= false;
      bool JetsPt30C =false;
      

      if (JetsMV.size() >3) continue;

	int xj = -1;
      
	for (unsigned int ib = 0; ib <JetsMV.size();++ib){
	
		if (JetsMV.at(ib).Pt()>30) JetsPt30C = true;
		for (unsigned int il = 0; il < analysisTree.pfjet_count; ++il)
	    	{
	
    		if (float(JetsMV.at(ib).Pt()) == float(analysisTree.pfjet_pt[il])) {xj=il;
	      //cout<<" found Jets "<<analysisTree.pfjet_pt[il]<<"  "<<JetsMV.at(ib).Pt()<<"  "<<ib<<"  "<<xj<<"  "<<il<<endl;
		        }
		}

	//     cout<<" Jets "<<analysisTree.pfjet_pt[xj]<<"  "<<JetsMV.at(ib).Pt()<<"  "<<ib<<"  "<<xj<<endl;
      	      if (analysisTree.pfjet_btag[xj][8]  > bTag) btagged = true;
      }
      
      if (JetsPt30C || btagged ||  JetsMV.size() >3) continue;
	

          // Jets
	  FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
          // pt Scalar
  /*  	
	  if (ptScalarSum<0  ) continue;
          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
    */
	  // computations of kinematic variables

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[mu_index],
					    analysisTree.muon_py[mu_index],
					    analysisTree.muon_pz[mu_index],
					    muonMass);

      TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[el_index],
						    analysisTree.electron_py[el_index],
						    analysisTree.electron_pz[el_index],
						    electronMass);
      

      TLorentzVector dileptonLV = muonLV + electronLV;
      float dileptonMass = dileptonLV.M();
      float dileptonPt = dileptonLV.Pt();
      float dileptonEta = dileptonLV.Eta();

      TLorentzVector diL = MuMV.at(mu_index) + TauMV.at(tau_index);
      if ( diL.M() <80 && diL.M()>40 ) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



      float ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // bisector of electron and muon transverse momenta
      float electronUnitX = electronLV.Px()/electronLV.Pt();
      float electronUnitY = electronLV.Py()/electronLV.Pt();
	
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();

      float zetaX = electronUnitX + muonUnitX;
      float zetaY = electronUnitY + muonUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + electronLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + electronLV.Py();
      
      float vectorVisX = muonLV.Px() + electronLV.Px();
      float vectorVisY = muonLV.Py() + electronLV.Py();

      // computation of DZeta variable
      float PZeta = vectorX*zetaX + vectorY*zetaY;
      float PVisZeta = vectorVisX*zetaX + vectorVisY*zetaY;
      float DZeta = PZeta - 1.85*PVisZeta;

      // computation of MT variable
      float dPhi=-999; 
      float MT = -999;

	   	dPhi=dPhiFrom2P( dileptonLV.Px(), dileptonLV.Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      		//MT=TMath::Sqrt(2*dileptonPt*ETmiss*(1-TMath::Cos(dPhi)));
	




      // filling histograms after dilepton selection

      electronPtH->Fill(electronLV.Pt(),weight);
      electronEtaH->Fill(electronLV.Eta(),weight);
      
      muonPtH->Fill(muonLV.Pt(),weight);
      muonEtaH->Fill(muonLV.Eta(),weight);
      
      dileptonMassH->Fill(dileptonMass,weight);
      dileptonPtH->Fill(dileptonPt,weight);
      dileptonEtaH->Fill(dileptonEta,weight);
     // dileptondRH->Fill(dRleptons,weight);
      
      ETmissH->Fill(ETmiss,weight);
      MtH->Fill(MT,weight);
      DZetaH->Fill(DZeta,weight);

      if (ETmiss < metcut) continue;
          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
       
       if (ETmiss < 2*metcut) continue;
          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
      // topological cut
      //if (DZeta<dZetaCut) continue;
       if (dPhi<1) continue; 
       
          FillMainHists(iCut, weight, ElMV, MuMV, TauMV, JetsMV,METV, analysisTree, SelectionSign,mu_index,el_index,tau_index);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	  

      electronPtSelH->Fill(electronLV.Pt(),weight);
      electronEtaSelH->Fill(electronLV.Eta(),weight);
      
      muonPtSelH->Fill(muonLV.Pt(),weight);
      muonEtaSelH->Fill(muonLV.Eta(),weight);
      
      dileptonMassSelH->Fill(dileptonMass,weight);
      dileptonPtSelH->Fill(dileptonPt,weight);
      dileptonEtaSelH->Fill(dileptonEta,weight);
     // dileptondRSelH->Fill(dRleptons,weight);
      
      ETmissSelH->Fill(ETmiss,weight);
      MtSelH->Fill(MT,weight);
      DZetaSelH->Fill(DZeta,weight);


      //      std::cout << std::endl;
      
      selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    
    file_->Close();
    delete file_;
}
 

cout<<"done"<<endl;
	
cout<<" Will use weight  "<<histWeights->GetSumOfWeights()<<" Norm Factor "<<XSec*Lumi/( histWeights->GetSumOfWeights())<<endl;
/*
 for (int i=0;i<CutNumb;++i){
    CFCounter[i] *= Float_t(XSec*Lumi/( histWeights->GetSumOfWeights()));
    if (iCFCounter[i] <0.2) statUnc[i] =0;
    else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
  }
*/

for (int i=0;i<CutNumb;++i){
 if (!isData) { cout << " i "<<i<<" "<<iCFCounter[i]<<"  "<<XSec*Lumi/( histWeights->GetSumOfWeights())<<endl;  
	 CFCounter[i] *= Float_t(XSec*Lumi/( histWeights->GetSumOfWeights()));}
    if (iCFCounter[i] <0.2) statUnc[i] =0;
    //else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
    else statUnc[i] = sqrt(CFCounter[i]);
  }
  //ofstream tfile1;
  //TString textfile_Con = "CMG_cutflow_Con_Mu_"+outname+".txt";
  //tfile1.open(textfile_Con);
  //tfile1 << "########################################" << endl;
  //tfile1 << "RCS:" << endl;



  //write out cutflow
  ofstream tfile;
  // TString outname = argv[argc-1];
  TString outname=argv[2];
  TString textfilename = "cutflow_"+outname+"_"+SelectionSign+"_"+argv[3]+".txt";
  tfile.open(textfilename);
  tfile << "########################################" << endl;
  //tfile << "Cut efficiency numbers:" << endl;

  //    tfile << " Cut "<<"\t & \t"<<"#Evnts for "<<Lumi/1000<<" fb-1 & \t"<<" Uncertainty \t"<<" cnt\t"<<endl;
  for(int ci = 0; ci < CutNumb; ci++)
    {
      tfile << CutList[ci]<<"\t & \t"
	    << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
       CutFlow->SetBinContent(1+ci,CFCounter[ci]);
      //CutFlow->SetBinContent(1+ci,);
    }

  tfile.close();


    //ofstream tfile1;
    //TString textfile_Con = "CMG_cutflow_Con_Mu_"+outname+".txt";
    //tfile1.open(textfile_Con);
    //tfile1 << "########################################" << endl;
    //tfile1 << "RCS:" << endl;






  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd(SelectionSign.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeights->Write();

  muonPtAllH->Write();
  electronPtAllH->Write();

  // histograms (dilepton selection)
  electronPtH->Write();  
  electronEtaH ->Write();
  muonPtH ->Write();
  muonEtaH ->Write();

  dileptonMassH ->Write();
  dileptonPtH ->Write();
  dileptonEtaH ->Write();
  dileptondRH ->Write();
  ETmissH ->Write();
  MtH ->Write();
  DZetaH ->Write();

  // histograms (dilepton selection + DZeta cut DZeta)
  electronPtSelH ->Write();
  electronEtaSelH ->Write();
  muonPtSelH  ->Write();
  muonEtaSelH ->Write();

  dileptonMassSelH ->Write();
  dileptonPtSelH ->Write();
  dileptonEtaSelH ->Write();
  dileptondRSelH ->Write();
  ETmissSelH ->Write();
  MtSelH ->Write();
  DZetaSelH ->Write();


  file->Write();
  file->Close();
  
  delete file;
  
}



