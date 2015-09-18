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

#include "AnalysisMacro.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  //const int CutNumb = 8;
  //string CutList[CutNumb]={"No cut","Trigger","1l","lept-Veto","b-Veto","MET $>$ 50","MET $>$ 100","dPhi $>$ 1"};

  // **** configuration
  Config cfg(argv[1]);
  string SelectionSign="mutau";

  // kinematic cuts on electrons
  const Float_t ptElectronLowCut   = cfg.get<Float_t>("ptElectronLowCut");
  const Float_t ptElectronHighCut  = cfg.get<Float_t>("ptElectronHighCut");
  const Float_t etaElectronCut     = cfg.get<Float_t>("etaElectronCut");
  const Float_t dxyElectronCut     = cfg.get<Float_t>("dxyElectronCut");
  const Float_t dzElectronCut      = cfg.get<Float_t>("dzElectronCut");
  const Float_t isoElectronLowCut  = cfg.get<Float_t>("isoElectronLowCut");
  const Float_t isoElectronHighCut = cfg.get<Float_t>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // kinematic cuts on muons
  const Float_t ptMuonLowCut   = cfg.get<Float_t>("ptMuonLowCut");
  const Float_t ptMuonHighCut  = cfg.get<Float_t>("ptMuonHighCut");
  const Float_t etaMuonCut     = cfg.get<Float_t>("etaMuonCut");
  const Float_t dxyMuonCut     = cfg.get<Float_t>("dxyMuonCut");
  const Float_t dzMuonCut      = cfg.get<Float_t>("dzMuonCut");
  const Float_t isoMuonLowCut  = cfg.get<Float_t>("isoMuonLowCut");
  const Float_t isoMuonHighCut = cfg.get<Float_t>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  


  const string Mu24Leg  = cfg.get<string>("Mu24Leg");
  const string Mu27Leg  = cfg.get<string>("Mu27Leg");
  const string Mu17Tau20MuLegA  = cfg.get<string>("Mu17Tau20MuLegA");
  const string Mu17Tau20MuLegB  = cfg.get<string>("Mu17Tau20MuLegB");
  const string Mu17Tau20TauLegA  = cfg.get<string>("Mu17Tau20TauLegA");
  const string Mu17Tau20TauLegB  = cfg.get<string>("Mu17Tau20TauLegB");



  const Float_t leadchargedhadrcand_dz = cfg.get<Float_t>("leadchargedhadrcand_dz");
  const Float_t leadchargedhadrcand_dxy = cfg.get<Float_t>("leadchargedhadrcand_dxy");


//  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");  
 
//  const float zVertexCut     = cfg.get<float>("ZVertexCut");
//  const float dVertexCut     = cfg.get<float>("DVertexCut");


  // kinematic cuts on Jets
  const Float_t etaJetCut   = cfg.get<Float_t>("etaJetCut");
  const Float_t ptJetCut   = cfg.get<Float_t>("ptJetCut");
  
  
  // topological cuts
  const Float_t dRleptonsCutmutau   = cfg.get<Float_t>("dRleptonsCutmutau");
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
  

  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
  
  TString Muon24Leg(Mu24Leg);
  TString Muon27Leg(Mu27Leg);
  TString Muon17Tau20MuLegA (Mu17Tau20MuLegA );
  TString Muon17Tau20MuLegB (Mu17Tau20MuLegB );
  TString Muon17Tau20TauLegA (Mu17Tau20TauLegA );
  TString Muon17Tau20TauLegB (Mu17Tau20TauLegB );



  const Float_t Lumi   = cfg.get<Float_t>("Lumi");
  const Float_t bTag   = cfg.get<Float_t>("bTag");
  const Float_t metcut = cfg.get<Float_t>("metcut");
 
  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("$mu$");
  CutList.push_back("$tau_h$");
  CutList.push_back("$DeltaR>0.5$");
  CutList.push_back("Trigger");
  CutList.push_back("2nd lept-Veto");
  CutList.push_back("3rd lept-Veto");
  CutList.push_back("b-Veto ");
  CutList.push_back("ZMass cut ");
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
  bool doMuVeto=true;

  //CutList[CutNumb]=CutListt[CutNumb];
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

  Float_t Weight=0;
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
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
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
    //if ((TH1D*)file_->Get("histoWeights")) 
    // histWeights2= (TH1D*)file_->Get("makeroottree/histoWeights");
     //histWeights= (TH1D*)file_->Get("makeroottree/histoWeights");
    // Weight =  histWeights2->GetSumOfWeights();
     //histWeights->Fill(histWeights2->GetBinContent(1));
     //if (histWeights2 ==NULL) continue;
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    //NE=1000;

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    //numberOfEntries = 1000;
    
    std::cout << "      number of entries in Tree = " << numberOfEntries <<" Weight  "<<Weight<< std::endl;
   // numberOfEntries = 1000;
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
     
      analysisTree.GetEntry(iEntry);
      nEvents++;
    
      iCut = 0;
      
      
      Float_t weight = 1;
      isData= false;
      bool lumi=false;

      if (XSec == 1)  isData = true;
      if (!isData && XSec !=1 )  { 
	      
      histWeights->Fill(1,analysisTree.genweight);  
	     weight *= analysisTree.genweight;
	     //weight *= histWeights2->GetBinContent(1); 
	      lumi=true;
      		} 

    // cout<<"  "<<histWeights->GetBinContent(1)<<"  "<<histWeights2->GetSumOfWeights()<<"  "<<weight<<endl; //(TH1D*)file_->Get("makeroottree/histoWeights");
   	

     
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
//	mu_index=0;
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

      // 	Float_t pxGen = analysisTree.genparticles_px[igen];
      // 	Float_t pyGen = analysisTree.genparticles_py[igen];
      // 	Float_t pzGen = analysisTree.genparticles_pz[igen];
      // 	Float_t etaGen = PtoEta(pxGen,pyGen,pzGen);
      // 	Float_t ptGen  = PtoPt(pxGen,pyGen);

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

 
      //selecTable.Fill(1,0, weight );      
      bool trigAccept = false;
      
      unsigned int nMuon24Leg = 0;
      bool isMuon24Leg = false;
      
      unsigned int nMuon27Leg = 0;
      bool isMuon27Leg = false;
      
      unsigned int nMuon17Tau20MuLegA = 0;
      bool isMuon17Tau20MuLegA = false;
      
      unsigned int nMuon17Tau20MuLegB = 0;
      bool isMuon17Tau20MuLegB = false;
      
      unsigned int nMuon17Tau20TauLegA = 0;
      bool isMuon17Tau20TauLegA = false;

      unsigned int nMuon17Tau20TauLegB = 0;
      bool isMuon17Tau20TauLegB = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
          //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==Muon24Leg) {
	  nMuon24Leg = i;
	  isMuon24Leg = true;
	}
	if (HLTFilter==Muon27Leg) {
	  nMuon27Leg = i;
	  isMuon27Leg = true;
	}
      
	if (HLTFilter==Muon17Tau20MuLegA) {
	  nMuon17Tau20MuLegA = i;
	  isMuon17Tau20MuLegA = true;
	}
      
	if (HLTFilter==Muon17Tau20MuLegB) {
	  nMuon17Tau20MuLegB = i;
	  isMuon17Tau20MuLegB = true;
	}


	if (HLTFilter==Muon17Tau20MuLegA) {
	  nMuon17Tau20MuLegA = i;
	  isMuon17Tau20TauLegA = true;
	}
	if (HLTFilter==Muon17Tau20MuLegB) {
	  nMuon17Tau20MuLegB = i;
	  isMuon17Tau20TauLegB = true;
	}

      
      }



      if (!isMuon24Leg) {
	std::cout << "HLT filter " << Muon24Leg << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon27Leg) {
	std::cout << "HLT filter " << Muon27Leg << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20MuLegA) {
	std::cout << "HLT filter " << Muon17Tau20MuLegA << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20MuLegB) {
	std::cout << "HLT filter " << Muon17Tau20MuLegB << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20TauLegA) {
	std::cout << "HLT filter " << Muon17Tau20TauLegA << " not found" << std::endl;
	return(-1);
      }
     
      if (!isMuon17Tau20TauLegB) {
	std::cout << "HLT filter " << Muon17Tau20TauLegB << " not found" << std::endl;
	return(-1);
      }
/*
            std::cout << "LowPtE  : " << LowPtLegElectron << " : " << nLowPtLegElectron << std::endl;
            std::cout << "HighPtE : " << HighPtLegElectron << " : " << nHighPtLegElectron << std::endl;
            std::cout << "LowPtM  : " << LowPtLegMuon << " : " << nLowPtLegMuon << std::endl;
            std::cout << "HighPtM : " << HighPtLegMuon << " : " << nHighPtLegMuon << std::endl;
            std::cout << std::endl;
      //      continue;
      // vertex cuts
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
*/

      /////now clear the Mu.El.Jets again to fill them again after cleaning
    /*  MuMV.clear();
      ElMV.clear();
      TauMV.clear();
      LeptMV.clear();
*/
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
            //  MuMV.push_back(MuV);
	    }

	    if (relIso!=0 && relIso==isoMuMin && im != mu_index) {
             analysisTree.muon_pt[im] > analysisTree.muon_pt[mu_index] ? mu_index = im : mu_index = mu_index;
	    cout<<" found a pair  " <<relIso <<"  "<<mu_index<<"  "<<im<<endl;
	  }
   
      }
      if (muons.size()==0 || !mu_iso) continue;

      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      if (LeptMV.size() == 0 ) continue; 
      
      //mu_index=muons[0];
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      Float_t isoTauMin = 999;
      bool tau_iso = false;
      vector<int> tau; tau.clear();
      for (unsigned  int it = 0; it<analysisTree.tau_count; ++it) {

	tauPtAllH->Fill(analysisTree.tau_pt[it],weight);
	tauEtaAllH->Fill(analysisTree.tau_eta[it],weight);
	if (analysisTree.tau_pt[it] < 20 || fabs(analysisTree.tau_eta[it])> 2.3) continue;
	if (analysisTree.tau_decayModeFinding[it]<decayModeFinding && analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	if (analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
	
	if (analysisTree.tau_againstElectronVLooseMVA5[it]<againstElectronVLooseMVA5) continue;
	if (analysisTree.tau_againstMuonTight3[it]<againstMuonTight3) continue;
	//phys14 if ( fabs(analysisTree.tau_vertexz[it] - analysisTree.primvertex_z ) > vertexz ) continue;
	if (analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it] > byCombinedIsolationDeltaBetaCorrRaw3Hits ) continue;
	

	Float_t  tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];

	    if (tauIso<isoTauMin ) {
	//      cout<<"  there was a chenge  "<<tauIso<<"  "<<isoTauMin<<" it "<<it<<" tau_index "<<tau_index<<"  "<<analysisTree.tau_count<<endl;
	      isoTauMin  = tauIso;
	      tau_iso=true;
	  //   it > 0 ? tau_index= it -1: 
	      tau_index = it;
	      tau.push_back(it);
	      TauV.SetPtEtaPhiM(analysisTree.tau_pt[it], analysisTree.tau_eta[it], analysisTree.tau_phi[it], tauMass);
//	      TauMV.push_back(TauV);

	    }

	    if (tauIso!=0 && tauIso==isoTauMin && it != tau_index) {
             analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ? tau_index = it : tau_index = tau_index;
	      cout<<" found a pair  " <<tauIso <<"  "<<tau_index<<"  "<<it<<endl;
	  }
      }
      if (tau.size()==0 || !tau_iso) continue;


	//cout<< " Lets check  "<<mu_index <<"  "<<tau_index <<"  "<<endl;
	//cout<<"  "<<endl;

      //tau_index=tau[0];
      //if ( !tau_iso) continue;

      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      float dR = deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);

      if (dR<dRleptonsCutmutau) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
/*
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
*/
	bool isMu24 = false;
	bool isMu27 = false;
	bool isMuTau_MuLegA = false;
	bool isMuTau_MuLegB = false;
	bool isMuTau_TauLegA = false;
	bool isMuTau_TauLegB = false;
        for (unsigned int im=0; im<muons.size(); ++im) {
	 isMu24 = false;
	 isMu27 = false;
	 isMuTau_MuLegA = false;
	 isMuTau_MuLegB = false;
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
	  if (analysisTree.trigobject_filters[iT][nMuon24Leg]) { // Mu24 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMu24 = true;
	    }
	  }
	  if (analysisTree.trigobject_filters[iT][nMuon27Leg]) { // Mu27 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMu27 = true;
				}
			}

	  if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegA]) { // Mu27 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMuTau_MuLegA = true;
				}
			}

	  if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegB]) { // Mu27 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isMuTau_MuLegB = true;
				}
			}
		}
		
	
     	if ((!isMu24) && (!isMu27) && (!isMuTau_MuLegA || !isMuTau_MuLegB) ) continue;


        for (unsigned int it=0; it<tau.size(); ++it) {
	 isMuTau_TauLegA = false;
	 isMuTau_TauLegB = false;
        unsigned int tIndex  = tau.at(it);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegA] ) { 
	      float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		isMuTau_TauLegA = true;
	      }
	    }

	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegB]) { 
	      float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		isMuTau_TauLegB = true;
		}
     	     } 
	  }
	}

     	if (     ( (isMu24) || (isMu27) )  || ( (isMuTau_MuLegA && isMuTau_MuLegB) && (isMuTau_TauLegA && isMuTau_TauLegB) ) ) trigAccept=true;
	
//cout<<" mu_index "<<mu_index<<"  "<<isMu24<<"  "<<isMu27<<"  "<<isMuTau_MuLegA<<"  "<<isMuTau_MuLegB<<"  "<<isMuTau_TauLegA<<"  "<<isMuTau_TauLegB<<endl;
	}//muons 
        if (!trigAccept) continue;
	


      //Trigger
      //FillMainHists(iCut, weight, ElMV, MuMV, JetsMV,METV,analysisTree, SelectionSign, mu_index,el_index,tau_index);
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      // electron selection





	//Set this flag if there is an opposite-charge muon pair in the event with muons separated by DR>0.15 and both passing the loose selection: 
	


      bool MuVeto=false;

      if (doMuVeto){
     	if (muons.size()>1){
	  for (unsigned  int imv = 0; imv<analysisTree.muon_count; ++imv) {
       if ( imv != mu_index ){

	    Float_t neutralIso = 
	      analysisTree.muon_neutralHadIso[imv] + 
	      analysisTree.muon_photonIso[imv] - 
	      0.5*analysisTree.muon_puIso[imv];
	    neutralIso = TMath::Max(Float_t(0),neutralIso); 
	    Float_t absIso = analysisTree.muon_chargedHadIso[imv] + neutralIso;
	    Float_t relIso = absIso/analysisTree.muon_pt[imv];
      
	    float dRr = deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			    analysisTree.muon_eta[imv],analysisTree.muon_phi[imv]);

	    if ( analysisTree.muon_charge[imv] != analysisTree.muon_charge[mu_index] &&  analysisTree.muon_isGlobal[imv] && analysisTree.muon_isTracker[imv] && analysisTree.muon_isPF[imv]  
		 &&  analysisTree.muon_pt[imv]> 15 &&  fabs(analysisTree.muon_eta[imv])< 2.4 && fabs(analysisTree.muon_dxy[imv])<0.045 
		 && fabs(analysisTree.muon_dz[imv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imv]) && dRr > 0.15)
	
	      MuVeto=true;
	  }
	}
      }
}
      if (MuVeto) continue;

      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



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


      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      //	for (unsigned int j=0;j<LeptMV.size();++j) cout<<" j "<<j<<"  "<<LeptMV.at(j).Pt()<<endl;
      //	cout<<""<<endl;
      ////////jets cleaning 
      Float_t DRmax=0.4;
      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;
      
      //      continue;
      
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
      


      // selecting muon and electron pair (OS or SS);
      Float_t ptScalarSum = -1;


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
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      // pt Scalar
      // computations of kinematic variables

//cout<<"  "<<mu_index<<"  "<<tau_index<<"   "<<MuMV.at(mu_index).M()<<"  "<<TauMV.at(tau_index).M()<<endl;

      TLorentzVector diL = MuMV.at(mu_index) + TauMV.at(tau_index);
      if ( diL.M() <80 && diL.M()>40 ) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      Float_t ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // bisector of electron and muon transverse momenta

      // computation of MT variable
      Float_t dPhi=-999; 


      dPhi=dPhiFrom2P( LeptMV.at(0).Px(), LeptMV.at(0).Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      //MT = TMath::Sqrt(2*LeptMV.at(0).Pt()*ETmiss*(1-TMath::Cos(dPhi)));


      // filling histograms after dilepton selection

      
      // ETmissH->Fill(ETmiss,weight);
      // MtH->Fill(MT,weight);
      
      if (ETmiss < metcut) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
       
      if (ETmiss < 2*metcut) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      // topological cut
      //if (DZeta<dZetaCut) continue;
      if (dPhi>1) continue; 
       
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      
      //      std::cout << std::endl;
      
      selEvents++;
      //histWeights->SetBinContent(1,Weight+ histWeights->GetBinContent(1));  
    } // end of file processing (loop over events in one file)
	//histWeights->Fill(Weight);
  
    cout<< " Weight  "<<Weight<<  "   histWeigh Sum " <<histWeights->GetBinContent(1)<<"  "<<histWeights2->GetSumOfWeights()<<"   "<<" Events "<<inputEventsH->GetSum()<<"  "<<inputEventsH->GetSumOfWeights()<<endl;
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
  //histWeights2->Write();
  
  CutFlow->Write();

  /*
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
*/
  DZetaSelH ->Write();
  /*
    for(int cj = 0; cj < CutNumb; cj++)
    {
    file->cd("");
    //outf->mkdir(CutList[cj]);
    //outf->cd(CutList[cj]);
    h0JetpT[cj]->Write();
    hnJet[cj]->Write();
    hnOver[cj]->Write();
    hnBJet[cj]->Write();
    hnEl[cj]->Write();
    hElpt[cj]->Write();
    hnMu[cj]->Write();
    hMupt[cj]->Write();
    hLepeta[cj]->Write();
    hMET[cj]->Write();
    hHT[cj]->Write();
    hST[cj]->Write();
    hToppT[cj]->Write();
    hnTop[cj]->Write();
    hWTagpT[cj]->Write();
    hWTagMass[cj]->Write();
    hnW[cj]->Write();
    hWmassTagpT[cj]->Write();
    hWmassTagMass[cj]->Write();
    hnWmass[cj]->Write();
    hdPhiMETLep[cj]->Write();
    hdPhiJMET[cj]->Write();

    }
  */
  
  file->Write();
  file->Close();
  
  delete file;
  
}



