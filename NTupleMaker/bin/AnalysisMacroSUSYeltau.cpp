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
  string SelectionSign="eltau";

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
  

  const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
  const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");


  const string El32Leg  = cfg.get<string>("El32Leg");
  const string El22Tau20ElLegA  = cfg.get<string>("El22Tau20ElLegA");
  const string El22Tau20ElLegB  = cfg.get<string>("El22Tau20ElLegB");
  const string El22Tau20TauLegA  = cfg.get<string>("El22Tau20TauLegA");
  const string El22Tau20TauLegB  = cfg.get<string>("El22Tau20TauLegB");


  const Float_t leadchargedhadrcand_dz = cfg.get<Float_t>("leadchargedhadrcand_dz");
  const Float_t leadchargedhadrcand_dxy = cfg.get<Float_t>("leadchargedhadrcand_dxy");


//  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");  
 
//  const float zVertexCut     = cfg.get<float>("ZVertexCut");
//  const float dVertexCut     = cfg.get<float>("DVertexCut");


  // kinematic cuts on Jets
  const Float_t etaJetCut   = cfg.get<Float_t>("etaJetCut");
  const Float_t ptJetCut   = cfg.get<Float_t>("ptJetCut");
  
  
  // topological cuts
  const Float_t dRleptonsCuteltau   = cfg.get<Float_t>("dRleptonsCuteltau");
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
  
  const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
  const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");

  TString LowPtLegElectron(lowPtLegElectron);
  TString HighPtLegElectron(highPtLegElectron);
  
  TString LowPtLegMuon(lowPtLegMuon);
  TString HighPtLegMuon(highPtLegMuon);
  
  TString Elecs24Leg(El32Leg);
  TString Elecs22Tau20ElLegA(El22Tau20ElLegA);
  TString Elecs22Tau20ElLegB(El22Tau20ElLegB);
  TString Elecs22Tau20TauLegA(El22Tau20TauLegA);
  TString Elecs22Tau20TauLegB(El22Tau20TauLegB);


  const Float_t Lumi   = cfg.get<Float_t>("Lumi");
  const Float_t bTag 	     = cfg.get<Float_t>("bTag");
  const Float_t metcut         = cfg.get<Float_t>("metcut");
 
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
  bool doElVeto=true;


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
 // if (nTotalFiles>400) nTotalFiles=250;
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
    //numberOfEntries = 10000;
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
 
 

      for (unsigned int ij = 0; ij<analysisTree.pfjet_count; ++ij) {
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
	JetsMV.push_back(JetsV);
      } 



      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
	MuMV.push_back(MuV);
	//mu_index=0;

      }

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
	ElMV.push_back(ElV);
	//el_index=0;
      }
   
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	TauV.SetPtEtaPhiM(analysisTree.tau_pt[it], analysisTree.tau_eta[it], analysisTree.tau_phi[it], tauMass);
	TauMV.push_back(TauV);
	//tau_index=0;
      }

      
     // vector <string> ss; ss.push_back(SelectionSign.c_str());
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
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
      
      unsigned int nElecs24Leg = 0;
      bool isElecs24Leg = false;
 
      unsigned int nElecs22Tau20ElLegA = 0;
      bool isElecs22Tau20ElLegA = false;
   
      unsigned int nElecs22Tau20ElLegB = 0;
      bool isElecs22Tau20ElLegB = false;

      unsigned int nElecs22Tau20TauLegA = 0;
      bool isElecs22Tau20TauLegA = false;

      unsigned int nElecs22Tau20TauLegB = 0;
      bool isElecs22Tau20TauLegB = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
          //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));

	if (HLTFilter==Elecs24Leg) {
      	  nElecs24Leg = i;
	  isElecs24Leg = true;
	}
	if (HLTFilter==Elecs22Tau20ElLegA) {
      	  nElecs22Tau20ElLegA = i;
	  isElecs22Tau20ElLegA = true;
	}
	
	if (HLTFilter==Elecs22Tau20ElLegB) {
      	  nElecs22Tau20ElLegB = i;
	  isElecs22Tau20ElLegB = true;
	}


	if (HLTFilter==Elecs22Tau20ElLegB) {
      	  nElecs22Tau20ElLegB = i;
	  isElecs22Tau20ElLegB = true;
	}


	if (HLTFilter==Elecs22Tau20TauLegA) {
      	  nElecs22Tau20TauLegA = i;
	  isElecs22Tau20TauLegA = true;
	}

	if (HLTFilter==Elecs22Tau20TauLegB) {
      	  nElecs22Tau20TauLegB = i;
	  isElecs22Tau20TauLegB = true;
	}


 	}


      if (!isElecs24Leg) {
	std::cout << "HLT filter " << Elecs24Leg << " not found" << std::endl;
	return(-1);
      }
      if (!isElecs22Tau20ElLegA) {
	std::cout << "HLT filter " << Elecs22Tau20ElLegA << " not found" << std::endl;
	return(-1);
      }
      if (!isElecs22Tau20ElLegB) {
	std::cout << "HLT filter " << Elecs22Tau20ElLegB << " not found" << std::endl;
	return(-1);
      }
      if (!isElecs22Tau20TauLegA) {
	std::cout << "HLT filter " << Elecs22Tau20TauLegA << " not found" << std::endl;
	return(-1);
      }

      if (!isElecs22Tau20TauLegB) {
	std::cout << "HLT filter " << Elecs22Tau20TauLegB << " not found" << std::endl;
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
/*      MuMV.clear();
      ElMV.clear();
      TauMV.clear();
*/
      float isoElMin = 9999;
      bool el_iso=false;
	vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	electronPtAllH->Fill(analysisTree.electron_pt[ie],weight);
	if (analysisTree.electron_pt[ie]<ptElectronHighCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	Float_t neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	Float_t relIso = absIso/analysisTree.electron_pt[ie];
        hel_relISO[1]->Fill(relIso,weight);
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
              ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
	     LeptMV.push_back(ElV);
        //      ElMV.push_back(ElV);
	    //cout<<" RelIso " <<relIso <<"  "<<el_index<<endl;
	    }

	    if (relIso!=0 && relIso==isoElMin && ie != el_index) {
             analysisTree.electron_pt[ie] > analysisTree.electron_pt[el_index] ? el_index = ie : el_index = el_index;
	    cout<<" found a pair  " <<relIso <<"  "<<el_index<<"  "<<ie<<endl;
	  }


      }

      if (electrons.size()==0 || !el_iso) continue;
      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      if (LeptMV.size() == 0 ) continue; 
      // el_index=electrons[0];

      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      float isoTauMin = 999;
      bool tau_iso = false;
      vector<int> tau; tau.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
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
	      //TauMV.push_back(TauV);

	    }

	    if (tauIso!=0 && tauIso==isoTauMin && it != tau_index) {
             analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ? tau_index = it : tau_index = tau_index;
	      cout<<" found a pair  " <<tauIso <<"  "<<tau_index<<"  "<<it<<endl;
	  }
      }
      if (tau.size()==0 || !tau_iso) continue;


	
      //cout<<" tau, elects " <<tau_index<<"  "<<el_index<<"  "<<endl;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      float dR = deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index]);

      if (dR<dRleptonsCuteltau) continue;
   
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


	bool isEl32 = false;
	bool isEl22Tau20ElLegA = false;
	bool isEl22Tau20ElLegB = false;
	bool isEl22Tau20TauLegA = false;
	bool isEl22Tau20TauLegB = false;
 
        for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	 isEl32 = false;
	 isEl22Tau20ElLegA = false;
	 isEl22Tau20ElLegB = false;
        unsigned int eIndex  = electrons.at(ie);

/*
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
*/

	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nElecs24Leg]) { 
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isEl32 = true;
	    }
	  }
	
	  if (analysisTree.trigobject_filters[iT][nElecs22Tau20ElLegA]) { 
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isEl22Tau20ElLegA = true;
	    }
	  }

	  if (analysisTree.trigobject_filters[iT][nElecs22Tau20ElLegB ]) { 
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isEl22Tau20ElLegB  = true;
				}
			}


		}
	if ( (!isEl32) && (!isEl22Tau20ElLegA || !isEl22Tau20ElLegB) ) continue;


        for (unsigned int it=0; it<tau.size(); ++it) {
	 isEl22Tau20TauLegA = false;
	 isEl22Tau20TauLegB = false;
        unsigned int tIndex  = tau.at(it);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	 
	  if (analysisTree.trigobject_filters[iT][nElecs22Tau20TauLegA]) { 
	    float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isEl22Tau20TauLegA = true;
	    }
	  }

	  if (analysisTree.trigobject_filters[iT][nElecs22Tau20TauLegB]) { 
	    float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isEl22Tau20TauLegB = true;
	    }
	  }

	}
	}
     	
	if  ( isEl32 || (  ( isEl22Tau20ElLegA && isEl22Tau20ElLegB) && ( isEl22Tau20TauLegA && isEl22Tau20TauLegB) )) trigAccept=true;

	}//electrons
	

        if (!trigAccept) continue;



      //Trigger
      //FillMainHists(iCut, weight, ElMV, MuMV, JetsMV,METV,analysisTree, SelectionSign);
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      // electron selection



      bool ElVeto=false;
      if (doElVeto){

	  for (unsigned int ievv = 0; ievv<analysisTree.electron_count; ++ievv) {

       if ( ievv != el_index ){
	 //   Float_t neutralIsoV = analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev] -  4*TMath::Pi()*(0.3*0.3)*analysisTree.rho[iev];
	 //   Float_t IsoWithEA =  analysisTree.electron_r03_sumChargedHadronPt[iev] + TMath::Max(Float_t(0), neutralIsoV);
	 //   Float_t relIsoV = IsoWithEA/analysisTree.electron_pt[iev];

	Float_t neutralIso =   analysisTree.electron_neutralHadIso[ievv] + 	  analysisTree.electron_photonIso[ievv] -  0.5*analysisTree.electron_puIso[ievv];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = analysisTree.electron_chargedHadIso[ievv] + neutralIso;
	Float_t relIsoV = absIso/analysisTree.electron_pt[ievv];
        hel_relISO[1]->Fill(relIsoV,weight);
	if (relIsoV>isoElectronHighCut) continue;
	if (relIsoV<isoElectronLowCut) continue;

             bool ElVetoID = electronVetoTight(analysisTree.electron_superclusterEta[ievv], analysisTree.electron_eta[ievv],analysisTree.electron_phi[ievv],  analysisTree.electron_full5x5_sigmaietaieta[ievv], 
			     analysisTree.electron_ehcaloverecal[ievv],  analysisTree.electron_dxy[ievv], analysisTree.electron_dz[ievv], analysisTree.electron_ooemoop[ievv],
			     relIsoV,analysisTree.electron_nmissinginnerhits[ievv],analysisTree.electron_pass_conversion[ievv]);

	    if ( analysisTree.electron_charge[ievv] != analysisTree.electron_charge[electrons[0]] && analysisTree.electron_pt[ievv] > 15 &&  fabs(analysisTree.electron_eta[ievv]) < 2.5 && fabs(analysisTree.electron_dxy[ievv])<0.045
		&& fabs(analysisTree.electron_dz[ievv]) < 0.2 && relIsoV< 0.3 && ElVetoID) 
		    ElVeto=true;
		  	}
		}
 	}
      if (ElVeto) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      bool ThirdLeptVeto=false;
      if (doThirdLeptVeto){
  	if (analysisTree.electron_count>0) {
	  for (unsigned int iev = 0; iev<analysisTree.electron_count; ++iev) {
		if (iev !=el_index) {

	 //   Float_t neutralIsoV = analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev] -  4*TMath::Pi()*(0.3*0.3)*analysisTree.rho[iev];
	 //   Float_t IsoWithEA =  analysisTree.electron_r03_sumChargedHadronPt[iev] + TMath::Max(Float_t(0), neutralIsoV);
	 //   Float_t relIsoV = IsoWithEA/analysisTree.electron_pt[iev];

	Float_t neutralIso = 
	  analysisTree.electron_neutralHadIso[iev] + 	  analysisTree.electron_photonIso[iev] -  0.5*analysisTree.electron_puIso[iev];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = analysisTree.electron_chargedHadIso[iev] + neutralIso;
	Float_t relIsoV = absIso/analysisTree.electron_pt[iev];
        hel_relISO[1]->Fill(relIsoV,weight);
	if (relIsoV>isoElectronHighCut) continue;
	if (relIsoV<isoElectronLowCut) continue;
             bool ElVetoID = electronVetoTight(analysisTree.electron_superclusterEta[iev], analysisTree.electron_eta[iev],analysisTree.electron_phi[iev],  analysisTree.electron_full5x5_sigmaietaieta[iev], 
			     analysisTree.electron_ehcaloverecal[iev],  analysisTree.electron_dxy[iev], analysisTree.electron_dz[iev], analysisTree.electron_ooemoop[iev],
			     relIsoV,analysisTree.electron_nmissinginnerhits[iev],analysisTree.electron_pass_conversion[iev]);


	    if ( analysisTree.electron_pt[iev] > 10 &&  fabs(analysisTree.electron_eta[iev]) < 2.5 && fabs(analysisTree.electron_dxy[iev])<0.045
		&& fabs(analysisTree.electron_dz[iev]) < 0.2 && relIsoV< 0.3 && ElVetoID) ThirdLeptVeto=true;

	  }
	}

	}

     	if (analysisTree.muon_count>0){
	  for (unsigned int imvv = 0; imvv<analysisTree.muon_count; ++imvv) {
	    Float_t neutralIso = 
	      analysisTree.muon_neutralHadIso[imvv] + 
	      analysisTree.muon_photonIso[imvv] - 
	      0.5*analysisTree.muon_puIso[imvv];
	    neutralIso = TMath::Max(Float_t(0),neutralIso); 
	    Float_t absIso = analysisTree.muon_chargedHadIso[imvv] + neutralIso;
	    Float_t relIso = absIso/analysisTree.muon_pt[imvv];
	    if ( analysisTree.muon_isMedium[imvv] &&  analysisTree.muon_pt[imvv]> 10 &&  fabs(analysisTree.muon_eta[imvv])< 2.4 && fabs(analysisTree.muon_dxy[imvv])<0.045 
		 && fabs(analysisTree.muon_dz[imvv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imvv]) )
	
	      ThirdLeptVeto=true;
	}

      }
      }
      if (ThirdLeptVeto) continue;



      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
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

      TLorentzVector diL = MuMV.at(el_index) + TauMV.at(tau_index);
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
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
       
      if (ETmiss < 2*metcut) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      // topological cut
      //if (DZeta<dZetaCut) continue;
      if (dPhi>1) continue; 
       
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, analysisTree, SelectionSign);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      
      //      std::cout << std::endl;
      
      selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }


cout<< " " <<histWeights->GetSumOfWeights()<<"  "<<inputEventsH->GetSum()<<endl;

for (int i=0;i<CutNumb;++i){
    CFCounter[i] *= Float_t(XSec*Lumi/( histWeights->GetSumOfWeights()));
    if (iCFCounter[i] <0.2) statUnc[i] =0;
    else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
  }


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
      CutFlow->SetBinContent(1+ci,iCFCounter[ci]);
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



