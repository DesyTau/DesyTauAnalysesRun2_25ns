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
#include "TGraphAsymmErrors.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  //const int CutNumb = 8;
  //string CutList[CutNumb]={"No cut","Trigger","1l","lept-Veto","b-Veto","MET $>$ 50","MET $>$ 100","dPhi $>$ 1"};

  // **** configuration
  Config cfg(argv[1]);
  string Channel="mutau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
  const bool InvertTauIso = cfg.get<bool>("InvertTauIso");
  const bool InvertMuIso = cfg.get<bool>("InvertMuIso");
  const double ptElectronLowCut   = cfg.get<double>("ptElectronLowCut");
  const double ptElectronHighCut  = cfg.get<double>("ptElectronHighCut");
  const double etaElectronCut     = cfg.get<double>("etaElectronCut");
  const double dxyElectronCut     = cfg.get<double>("dxyElectronCut");
  const double dzElectronCut      = cfg.get<double>("dzElectronCut");
  const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCut");
  const double isoElectronHighCut = cfg.get<double>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");
  // kinematic cuts on muons
  const double ptMuonLowCuT   = cfg.get<double>("ptMuonLowCut");
  const double ptMuonHighCut  = cfg.get<double>("ptMuonHighCut");
  const double etaMuonCut     = cfg.get<double>("etaMuonCut");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonLowCut  = cfg.get<double>("isoMuonLowCut");
  const double isoMuonHighCut = cfg.get<double>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
 
  const double ptMuonLowCutmtau = cfg.get<double>("ptMuonLowCutmtau"); 
  const double etaMuonCutmtau = cfg.get<double>("etaMuonCutmtau"); 


  string TrigLeg  ;
   if (!isData) TrigLeg  = cfg.get<string>("Mu17LegMC");
   if (isData) TrigLeg  = cfg.get<string>("Mu18LegData");
  const string Mu17Tau20MuLegA  = cfg.get<string>("Mu17Tau20MuLegA");
  const string Mu17Tau20MuLegB  = cfg.get<string>("Mu17Tau20MuLegB");
  const string Mu17Tau20TauLegA  = cfg.get<string>("Mu17Tau20TauLegA");
  const string Mu17Tau20TauLegB  = cfg.get<string>("Mu17Tau20TauLegB");
  
//  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");
  string Region;


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");
  
  
  // topological cuts
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double dZetaCut       = cfg.get<double>("dZetaCut");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
 
  // tau
  const double taupt    = cfg.get<double>("taupt");
  const double taueta    = cfg.get<double>("taueta");
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double   decayModeFindingNewDMs  = cfg.get<double>("decayModeFindingNewDMs");
  const double   againstElectronVLooseMVA5  = cfg.get<double>("againstElectronVLooseMVA5");
  const double   againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double   vertexz =  cfg.get<double>("vertexz");
  const double   byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  

  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
  
  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  // lepton scale factors
  const string muonSfDataBarrel = cfg.get<string>("MuonSfDataBarrel");
  const string muonSfDataEndcap = cfg.get<string>("MuonSfDataEndcap");
  const string muonSfMcBarrel = cfg.get<string>("MuonSfMcBarrel");
  const string muonSfMcEndcap = cfg.get<string>("MuonSfMcEndcap");
  
  TString MainTrigger(TrigLeg);
  TString Muon17Tau20MuLegA (Mu17Tau20MuLegA );
  TString Muon17Tau20MuLegB (Mu17Tau20MuLegB );
  TString Muon17Tau20TauLegA (Mu17Tau20TauLegA );
  TString Muon17Tau20TauLegB (Mu17Tau20TauLegB );



  const double Lumi   = cfg.get<double>("Lumi");
  const double bTag   = cfg.get<double>("bTag");
  const double metcut = cfg.get<double>("metcut");
 
  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("No cut after PU");
  CutList.push_back("$\\mu$");
  CutList.push_back("$\\tau_h$");
  CutList.push_back("Trigger");
  CutList.push_back("2nd $\\ell$-Veto");
  CutList.push_back("3rd $\\ell$-Veto");
  CutList.push_back("Jets $<$3");
  CutList.push_back("b-Veto");
  CutList.push_back("$40<\\rm{Inv}_M<80");
  CutList.push_back("$ E_T^{\\rm miss}>$ 85");
  CutList.push_back("$1.5<\\Delta R<4$");
 

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
      string dt,st1,st2;st1="stau2_1";st2="stau5_2";
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
      /*
  	if ( argv[2] == st1) {ChiMass=100;mIntermediate=200;}
	else if (argv[2] == st2) {ChiMass=200;mIntermediate=500;}
	*/
      if (isData) XSec=1.;
	ChiMass=0.0;
    }

  if (XSec<0 && !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}

  std::vector<unsigned int> allRuns; allRuns.clear();

  cout<<" ChiMass is "<<ChiMass<<"  "<<mIntermediate<<endl;
  bool doThirdLeptVeto=true;
  bool doMuVeto=true;

  //CutList[CutNumb]=CutListt[CutNumb];
  char ff[100];

	
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // reading vertex weights
 // TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertDataFileName);
 // TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertMcFileName);
  TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/SingleMuon_Run2015D_Stau.root");
 TFile * fileMcNVert = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/MC.root");
//hardcoded
  //TH1D * vertexDataH = (TH1D*)fileDataNVert->Get(TString(vertHistName));
  //TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get(TString(vertHistName));
  TH1D * vertexDataH = (TH1D*)fileDataNVert->Get("mutau/npv_0");
  TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get("mutau/npv_0");

  double normVertexData = vertexDataH->GetSumOfWeights();
  double normVertexMc   = vertexMcH->GetSumOfWeights();

  vertexDataH->Scale(1/normVertexData);
  vertexMcH->Scale(1/normVertexMc);

  TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataBarrel);  // mu SF barrel data
  TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataEndcap); // mu SF endcap data
  TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcBarrel);  // mu SF barrel MC
  TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcEndcap); // mu SF endcap MC 
  
  TGraphAsymmErrors *hEffBarrelData = (TGraphAsymmErrors*)f10->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapData = (TGraphAsymmErrors*)f11->Get("ZMassEndcap");
  TGraphAsymmErrors *hEffBarrelMC = (TGraphAsymmErrors*)f12->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapMC = (TGraphAsymmErrors*)f13->Get("ZMassEndcap");
  
  double * dataEffBarrel = new double[10];
  double * dataEffEndcap = new double[10];
  double * mcEffBarrel = new double[10];
  double * mcEffEndcap = new double[10];
  
  dataEffBarrel = hEffBarrelData->GetY();
  dataEffEndcap = hEffEndcapData->GetY();
  mcEffBarrel = hEffBarrelMC->GetY();
  mcEffEndcap = hEffEndcapMC->GetY();
  double Weight=0;
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
  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  
  TString era=argv[3];
  TString invMu, invTau;

  /*
  if (InvertMuIso && !InvertTauIso && Sign=="OS") { invMu = "_InvMuIso_"; Region= "C";}
  if (!InvertMuIso && InvertTauIso && Sign=="OS") { invTau = "_InvTauIso_"; Region= "C";}
  if (InvertMuIso  && InvertTauIso && Sign=="OS") { invMu = "_InvMuIso_"; invTau ="_InvTauIso_"; Region= "C";}

  if (InvertMuIso && !InvertTauIso && Sign=="SS") { invMu = "_InvMuIso_"; Region= "D";}
  if (!InvertMuIso && InvertTauIso && Sign=="SS") { invTau = "_InvTauIso_"; Region= "D";}
  if (InvertMuIso  && InvertTauIso && Sign=="SS") { invMu = "_InvMuIso_"; invTau ="_InvTauIso_"; Region= "D";}


  if (!InvertMuIso  && !InvertTauIso && Sign=="SS") { Region= "A";}
  if (!InvertMuIso  && InvertTauIso && Sign=="SS") { invTau ="_InvTauIso_";Region= "A";}
  if (InvertMuIso  && !InvertTauIso && Sign=="SS") { invMu ="_InvMuIso_";Region= "A";}


  if (!InvertMuIso  && !InvertTauIso && Sign=="OS") { Region= "B";}
 */
  if(InvertMuIso) invMu = "_InvMuIso_";
  if(InvertTauIso) invMu = "_InvTauIso_";

  TString TStrName(rootFileName+invMu+invTau+Sign);
  std::cout <<TStrName <<std::endl;  
 
  // output fileName with histograms
  TFile * file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;
      bool lumi=false;
      bool isLowIsoMu=false;
      bool isHighIsoMu = false;
      bool isLowIsoTau=false;
      bool isHighIsoTau = false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
 
  SetupHists(CutNumb); 
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;


    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
 

    //if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("initroottree/nEvents");
	
    //histWeights2= (TH1D*)file_->Get("initroottree/nEvents"); 
    
    //if ((TH1D*)file_->Get("histoWeights")) 
    // histWeights2= (TH1D*)file_->Get("makeroottree/histoWeights");
    //histWeights= (TH1D*)file_->Get("makeroottree/histoWeights");
    // Weight =  histWeights2->GetSumOfWeights();
    //histWeights->Fill(histWeights2->GetBinContent(1));
    //if (histWeights2 ==NULL) continue;
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    ============>  " << NE << std::endl;
    //NE=1000;

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    double genweights=1;
    
    if(!isData) 
    {

    TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	     
    genweightsTree->SetBranchAddress("genweight",&genweights);
    Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
    for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; iEntryInit++) { 
      genweightsTree->GetEntry(iEntryInit);
	histWeightsH->Fill(double(0),genweights);
    }
 	}

    Long64_t numberOfEntries = analysisTree.GetEntries();
    //numberOfEntries = 1000;
    double weight=1;

    std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;
    // numberOfEntries = 1000;
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
     
      analysisTree.GetEntry(iEntry);
      nEvents++;
    
      iCut = 0;
      
      
      weight = 1.;
     
      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
     
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  

      //isData= false;
      lumi=false;
      isLowIsoMu=false;
      isHighIsoMu = false;
      isLowIsoTau=false;
      isHighIsoTau = false;

      if (!isData && XSec!=1)  { 

        //cout<<"  "<<genweights<<"  "<<numberOfEntries<<endl;
	//histWeights->Fill(1,analysisTree.genweight);  
	weight *= analysisTree.genweight;
	histWeights->Fill(1,weight);  
        //cout<<"  weight from init "<<genweights<< "  "<<analysisTree.genweight<<"  "<<weight<<endl;
	lumi=true;
      } 

      if (isData)  {
	XSec=1;
	histWeights->Fill(1);
	histWeightsH->Fill(double(1));
	histRuns->Fill(analysisTree.event_run);
      // cout<<"  "<<histWeights->GetBinContent(1)<<"  "<<histWeights2->GetSumOfWeights()<<"  "<<weight<<endl; //(TH1D*)file_->Get("makeroottree/histoWeights");

	//if (analysisTree.event_run<RunRangeMin || analysisTree.event_run>RunRangeMax)
	//      cout<< "  event_run "<<analysisTree.event_run<<"  Min  "<<RunRangeMin<<" Max "<<RunRangeMax<<endl;
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
	           //std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	       //   std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }
	
	  }
    
      }
      //lumi=true;
      if (!lumi) continue;
      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();
      mu_index=-1;
      tau_index=-1;
      el_index=-1;
      

      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
	
	if (!isData && applyPUreweighting) {
	  int binNvert = vertexDataH->FindBin(analysisTree.primvertex_count);
	  double_t dataNvert = vertexDataH->GetBinContent(binNvert);
	  double_t mcNvert = vertexMcH->GetBinContent(binNvert);
	  if (mcNvert < 1e-10){mcNvert=1e-10;}
	  double_t vertWeight = dataNvert/mcNvert;
	  weight *= vertWeight;
	//  	  cout << "NVert = " << analysisTree.primvertex_count << "   weight = " << vertWeight << endl;
	}



      double MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      
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


      
      // vector <string> ss; ss.push_back(Channel.c_str());
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


 
      //selecTable.Fill(1,0, weight );      
      bool trigAccept = false;
      
      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;
      
      
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
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}
      
	if (HLTFilter==Muon17Tau20MuLegA) {
	  nMuon17Tau20MuLegA = i;
	  isMuon17Tau20MuLegA = true;
	}
      
	if (HLTFilter==Muon17Tau20MuLegB) {
	  nMuon17Tau20MuLegB = i;
	  isMuon17Tau20MuLegB = true;
	}


	if (HLTFilter==Muon17Tau20TauLegA) {
	  nMuon17Tau20TauLegA = i;
	  isMuon17Tau20TauLegA = true;
	}
	if (HLTFilter==Muon17Tau20TauLegB) {
	  nMuon17Tau20TauLegB = i;
	  isMuon17Tau20TauLegB = true;
	}

      
      }


      if (!isMainTrigger) {
	std::cout << "HLT filter for Mu20 " << MainTrigger << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20MuLegA) {
	std::cout << "HLT filter for Mu17LegA " << Muon17Tau20MuLegA << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20MuLegB) {
	std::cout << "HLT filter fpr Mu17LegB " << Muon17Tau20MuLegB << " not found" << std::endl;
	return(-1);
      }
      if (!isMuon17Tau20TauLegA) {
	std::cout << "HLT filter TauLegA " << Muon17Tau20TauLegA << " not found" << std::endl;
	return(-1);
      }
     
      if (!isMuon17Tau20TauLegB) {
	std::cout << "HLT filter TauLegB " << Muon17Tau20TauLegB << " not found" << std::endl;
	return(-1);
      }
      /*
	std::cout << "LowPtE  : " << LowPtLegElectron << " : " << nLowPtLegElectron << std::endl;
	std::cout << "HighPtE : " << HighPtLegElectron << " : " << nHighPtLegElectron << std::endl;
	std::cout << "LowPtM  : " << LowPtLegMuon << " : " << nLowPtLegMuon << std::endl;
	std::cout << "HighPtM : " << HighPtLegMuon << " : " << nHighPtLegMuon << std::endl;
	std::cout << std::endl;
      */

      /////now clear the Mu.El.Jets again to fill them again after cleaning
      /*  MuMV.clear();
	  ElMV.clear();
	  TauMV.clear();
	  LeptMV.clear();
      */
      double isoMuMin = 9999;
      bool mu_iso=false;
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	muonPtAllH->Fill(analysisTree.muon_pt[im],weight);
	if (analysisTree.muon_pt[im]<ptMuonLowCuT) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;


         double absIso= analysisTree.muon_r03_sumChargedHadronPt[im]
	            + max(analysisTree.muon_r03_sumNeutralHadronEt[im] + analysisTree.muon_r03_sumPhotonEt[im]
	            - 0.5 * analysisTree.muon_r03_sumPUPt[im],0.0);


	double relIso = absIso/analysisTree.muon_pt[im];
 

	hmu_relISO[1]->Fill(relIso,weight);

	if (relIso<isoMuonLowCut) continue;

	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;

	if (double(relIso)<double(isoMuMin)) {
	  //cout<<" after muIso index "<<int(mu_index)<<" pT "<<analysisTree.muon_pt[im]<<" relIso "<<relIso<<" isoMuMin "<<isoMuMin<<" muon_count "<<analysisTree.muon_count<<" im "<<im<<" event "<<iEntry<<endl;
	  isoMuMin  = relIso;
	  mu_index = int(im);
	  mu_iso=true;
	  muons.push_back(im);
	  MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
	  MuMV.push_back(MuV);
          LeptMV.push_back(MuV);
	}

	//cout<<" Indexes here  "<<im<<"   "<<mu_index<<endl;
//	if (relIso!=0 && relIso==isoMuMin && im != mu_index) {

	if (relIso == isoMuMin && im != mu_index) {
	  //cout<<" found a pair  for muons " <<relIso <<" mu_index  "<<mu_index<<"  pT "<<analysisTree.muon_pt[int(mu_index)]<<" new  index "<<im<<"  pT  "<<analysisTree.muon_pt[int(im)]<<" event "<<iEntry<<endl;
	  analysisTree.muon_pt[im] > analysisTree.muon_pt[mu_index] ? mu_index = int(im) : mu_index = mu_index;
	}
   
      }
      if (muons.size()==0 || !mu_iso) continue;

        double absIso= analysisTree.muon_r03_sumChargedHadronPt[mu_index]
	            + max(analysisTree.muon_r03_sumNeutralHadronEt[mu_index] + analysisTree.muon_r03_sumPhotonEt[mu_index]
	            - 0.5 * analysisTree.muon_r03_sumPUPt[mu_index],0.0);

	double relIso = absIso/analysisTree.muon_pt[mu_index];

	if (relIso>isoMuonHighCut ) { isHighIsoMu=true ;isLowIsoMu=false;}
	else	{ isHighIsoMu = false;isLowIsoMu=true;}

      
       	sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
        if (LeptMV.size() == 0 ) continue; 
      
	if (InvertMuIso && !isHighIsoMu) continue;
	if (!InvertMuIso && isHighIsoMu) continue;
	if (InvertMuIso && isLowIsoMu) continue;


      //mu_index=muons[0];
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      double isoTauMin = 9999;
      bool tau_iso = false;
      vector<int> tau; tau.clear();
      for (unsigned  int it = 0; it<analysisTree.tau_count; ++it) {

	tauPtAllH->Fill(analysisTree.tau_pt[it],weight);
	tauEtaAllH->Fill(analysisTree.tau_eta[it],weight);
	if (analysisTree.tau_pt[it] < ptMuonLowCutmtau || fabs(analysisTree.tau_eta[it])> etaMuonCutmtau) continue;
	//if (analysisTree.tau_decayModeFinding[it]<decayModeFinding || analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	if (analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
	
	if (analysisTree.tau_againstElectronVLooseMVA5[it]<againstElectronVLooseMVA5) continue;
	if (analysisTree.tau_againstMuonTight3[it]<againstMuonTight3) continue;
	//phys14 if ( fabs(analysisTree.tau_vertexz[it] - analysisTree.primvertex_z ) > vertexz ) continue;
	

	double  tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];
      

 //     cout<< " Lets check   "<<mu_index<<" isLowIsoTau "<<isLowIsoTau<<"  isLowIsoMu "<<isLowIsoMu<<" isHighIsoMu  "<<isHighIsoMu<<"  isHighIsoTau  "<<isHighIsoTau<<" invertMu "<<InvertMuIso<<" InvertTauIso "<<InvertTauIso<<"  Charge "<<Sign<<" tauIso "<<tauIso<<" tauIsoCut "<<byCombinedIsolationDeltaBetaCorrRaw3Hits<< endl;


      if (tauIso<isoTauMin ) {
	  //      cout<<"  there was a chenge  "<<tauIso<<"  "<<isoTauMin<<" it "<<it<<" tau_index "<<tau_index<<"  "<<analysisTree.tau_count<<endl;
	  isoTauMin  = tauIso;
	  tau_iso=true;
	  //   it > 0 ? tau_index= it -1: 
	  tau_index = int(it);
	  tau.push_back(it);
	  TauV.SetPtEtaPhiM(analysisTree.tau_pt[it], analysisTree.tau_eta[it], analysisTree.tau_phi[it], tauMass);
	  //	      TauMV.push_back(TauV);

	}

	//if (tauIso!=0 && tauIso==isoTauMin && it != tau_index) {
	if (tauIso==isoTauMin && it != tau_index) {
	  analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ? tau_index = it : tau_index = tau_index;
	  //cout<<" found a pair  " <<tauIso <<"  "<<tau_index<<"  "<<it<<endl;
	}
      }
      if (tau.size()==0 || !tau_iso) continue;

      double tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
	
      if (tauIso > byCombinedIsolationDeltaBetaCorrRaw3Hits && !InvertTauIso) continue;

      if (tauIso > byCombinedIsolationDeltaBetaCorrRaw3Hits ) {isHighIsoTau =true ; isLowIsoTau=false;} 
      else {isHighIsoTau =false ; isLowIsoTau=true;} 

	if (InvertTauIso && !isHighIsoTau) continue;
	if (!InvertTauIso && isHighIsoTau) continue;
	if (InvertTauIso && isLowIsoTau) continue;

      //cout<< " Lets check  mu_index "<<mu_index <<" tau_index "<<tau_index <<" isLowIsoTau "<<isLowIsoTau<<"  isLowIsoMu "<<isLowIsoMu<<" isHighIsoMu  "<<isHighIsoMu<<"  isLowIsoTau  "<<isLowIsoTau<<endl;
      //cout<<"  "<<endl;

     double q = analysisTree.tau_charge[tau_index] * analysisTree.muon_charge[mu_index];
    
    // cout<< " Lets check first  here "<<mu_index <<" tau_index "<<tau_index <<" isLowIsoTau "<<isLowIsoTau<<"  isLowIsoMu "<<isLowIsoMu<<" isHighIsoMu  "<<isHighIsoMu<<"  isHighIsoTau  "<<isHighIsoTau<<" invertMu "<<InvertMuIso<<" InvertTauIso "<<InvertTauIso<<"  Charge "<<Sign<<" q "<<q<<" tauIso "<<analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index]<<" muIso "<<relIso<< "  "<<iEntry<<endl; 

    
      //if (q>0 &&  (Sign=="OS" || Region !="B" || Region !="C")  ) continue;
      //if (q<0 && (Sign=="SS" || Region !="A" || Region !="D") ) continue;
      
      if (q>0 &&  Sign=="OS" ) continue;
      if (q<0 && Sign=="SS" ) continue;

      bool regionB = (q<0 && isLowIsoMu);
      bool regionA = (q>0 && isLowIsoMu);
      bool regionC = (q<0 && isHighIsoMu);
      bool regionD = (q>0 && isHighIsoMu);

	


   //  cout<< " Lets check  here "<<mu_index <<" tau_index "<<tau_index <<" isLowIsoTau "<<isLowIsoTau<<"  isLowIsoMu "<<isLowIsoMu<<" isHighIsoMu  "<<isHighIsoMu<<"  isHighIsoTau  "<<isHighIsoTau<<"  A   "<< regionA<<"  B  "<<regionB<<" C "<<regionC<<"  D  "<<regionD<<" invertMu "<<InvertMuIso<<" InvertTauIso "<<InvertTauIso<<"  Charge "<<Sign<<" q "<<q<<" tauIso "<<analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index]<<" muIso "<<relIso<< "  invertMu ? "<<InvertMuIso<< " InverTau ? "<<InvertTauIso<<"  "<<iEntry<<" Region Category "<<Region<<"   "<<TStrName<<endl; 
 
     FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



	//cout<<"  HOW MANY MUONS DO I HAVE ??  "<<muons.size()<<endl;
	/*
      isMainTrigger = false;
      bool isMu27 = false;
      bool isMuTau_MuLegA = false;
      bool isMuTau_MuLegB = false;
      bool isMuTau_TauLegA = false;
      bool isMuTau_TauLegB = false;
      for (unsigned int im=0; im<muons.size(); ++im) {
	isMainTrigger = false;
	isMu27 = false;
	isMuTau_MuLegA = false;
	isMuTau_MuLegB = false;
        unsigned int mIndex  = muons.at(im);

	  double neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
	  double photonIsoMu = analysisTree.muon_photonIso[mIndex];
	  double chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
	  double puIsoMu = analysisTree.muon_puIso[mIndex];
	  if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
	  }
	  double neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	  neutralIsoMu = TMath::Max(double(0),neutralIsoMu); 
	  double absIsoMu = chargedHadIsoMu + neutralIsoMu;
	  double relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];
	*/

	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nMainTrigger]) { // Mu17 Leg
	    double dRtrig = deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>deltaRTrigMatch) continue;
	 //cout<< " "<<analysisTree.muon_pt[mu_index]<<"   "<<analysisTree.trigobject_pt[iT]<<endl;
	    //if (!isData && analysisTree.trigobject_filters[iT][nMainTrigger] && analysisTree.trigobject_pt[iT]>18) 
	      isMainTrigger = true;
	    
	  }
	}
	/*
	  if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegA]) { // Mu27 Leg
	    double dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>deltaRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegA] && analysisTree.trigobject_pt[iT]>17) 
	      isMuTau_MuLegA = true;
	    
	  }

	  if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegB]) { // Mu27 Leg
	    double dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>deltaRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20MuLegB] && analysisTree.trigobject_pt[iT]>17) 
	      isMuTau_MuLegB = true;
	    
	  }
	*/
		
	
	bool passedTrigger = isMainTrigger ;//|| (isMuTau_MuLegA && isMuTau_MuLegB);
	if (!passedTrigger) continue;

	double dR = deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);

      	if (dR<dRleptonsCutmutau) continue;

	//if ((!isMainTrigger && !isMu27) || (!isMuTau_MuLegA && !isMuTau_MuLegB) ) continue;
     	//if ((!isMainTrigger ) && (!isMuTau_MuLegA && !isMuTau_MuLegB) ) continue;

	/*
        for (unsigned int it=0; it<tau.size(); ++it) {
	  isMuTau_TauLegA = false;
	  isMuTau_TauLegB = false;
	  unsigned int tIndex  = tau.at(it);
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegA] ) { 
	      double dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig>deltaRTrigMatch) continue;
	      if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegA] && analysisTree.trigobject_pt[iT]>20) 
		isMuTau_TauLegA = true;
	      
	    }

	    if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegB]) { 
	      double dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig>deltaRTrigMatch) continue;
	      if (analysisTree.trigobject_filters[iT][nMuon17Tau20TauLegA] && analysisTree.trigobject_pt[iT]>20) 
		isMuTau_TauLegB = true;
	      
	    } 
	  }
	

	//if (   ( (isMainTrigger) || (isMu27) )  || ( (isMuTau_MuLegA && isMuTau_MuLegB) && (isMuTau_TauLegA && isMuTau_TauLegB) ) ) trigAccept=true;
     	//if ((!isMainTrigger ) && (!isMuTau_MuLegA || !isMuTau_MuLegB) ) continue;
  
   	bool passedMuTauTriggerMuLeg	= isMuTau_MuLegA && isMuTau_MuLegB;
   	bool passedMuTauTriggerTauLeg   = isMuTau_TauLegA && isMuTau_TauLegB;
	bool OnlyMuTrigger = isMainTrigger && (!passedMuTauTriggerMuLeg && !passedMuTauTriggerTauLeg);

      	//if ( (isMainTrigger || isMu27) && ( (!isMuTau_MuLegA && !isMuTau_MuLegB) && (!isMuTau_TauLegA && !isMuTau_TauLegB) )  && 
//         if ( (isMainTrigger ) && ( (!isMuTau_MuLegA && !isMuTau_MuLegB) && (!isMuTau_TauLegA && !isMuTau_TauLegB) )  && 
	if (OnlyMuTrigger  && ( analysisTree.muon_pt[mu_index] < ptMuonHighCut || fabs(analysisTree.muon_eta[mu_index]) > 2.1) ) continue;
	
	
//	trigAccept=true;

//	if (!trigAccept) continue;
	
	}//taus

      }//muons 
      	*/
	
	//cout<<" mu_index "<<mu_index<<"  "<<isMainTrigger<<"  "<<isMu27<<"  "<<isMuTau_MuLegA<<"  "<<isMuTau_MuLegB<<"  "<<isMuTau_TauLegA<<"  "<<isMuTau_TauLegB<<endl;
      

	  double dataMu1 = 1;
	  double mcMu1 = 1;
	if (!isData && applyLeptonSF) {
	  double ptMu1 = TMath::Min(float(59),analysisTree.muon_pt[mu_index]);
	  int ptBinMu1 = binNumber(ptMu1,nPtBins,ptBins);
	  double absEtaMu1 = fabs(analysisTree.muon_eta[mu_index]);
	  dataMu1 = 1;
	  mcMu1 = 1;
	  if(absEtaMu1 < 1.48){
	    dataMu1 = dataEffBarrel[ptBinMu1];
	    mcMu1 = mcEffBarrel[ptBinMu1];
	  }
	  else {
	    dataMu1 = dataEffEndcap[ptBinMu1];
	    mcMu1 = mcEffEndcap[ptBinMu1];
	  }
	  double wMu1 = dataMu1/mcMu1;
	  	//  cout << "Muons SF : Mu1 = " << wMu1 << endl;
	  weight = weight*wMu1;
	}


      //Trigger
      //FillMainHists(iCut, weight, ElMV, MuMV, JetsMV,METV,analysisTree, Channel, mu_index,el_index,tau_index);
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;





      //Set this flag if there is an opposite-charge muon pair in the event with muons separated by DR>0.15 and both passing the loose selection: 
	


      bool MuVeto=false;

      if (doMuVeto){

     	if (muons.size()>1){
	  for (unsigned  int imv = 0; imv<analysisTree.muon_count; ++imv) {
	    if ( imv != mu_index ){
/*
	      double neutralIso = 
		analysisTree.muon_neutralHadIso[imv] + 
		analysisTree.muon_photonIso[imv] - 
		0.5*analysisTree.muon_puIso[imv];
	      neutralIso = TMath::Max(double(0),neutralIso); 
	      double absIso = analysisTree.muon_chargedHadIso[imv] + neutralIso;
	      double relIso = absIso/analysisTree.muon_pt[imv];
  */    
              double absIso= analysisTree.muon_r03_sumChargedHadronPt[imv]
	            + max(analysisTree.muon_r03_sumNeutralHadronEt[imv] + analysisTree.muon_r03_sumPhotonEt[imv]
	            - 0.5 * analysisTree.muon_r03_sumPUPt[imv],0.0);


	     double relIso = absIso/analysisTree.muon_pt[imv];


	      double dRr = deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
				 analysisTree.muon_eta[imv],analysisTree.muon_phi[imv]);
	        bool OSCharge = false;  
		if ( imv != mu_index  && analysisTree.muon_charge[imv] != analysisTree.muon_charge[mu_index] ) OSCharge=true;

	      //if ( analysisTree.muon_charge[imv] != analysisTree.muon_charge[mu_index] &&  analysisTree.muon_isGlobal[imv] && analysisTree.muon_isTracker[imv] && analysisTree.muon_isPF[imv]  
	      if ( analysisTree.muon_charge[imv] != analysisTree.muon_charge[mu_index] &&  analysisTree.muon_isGlobal[imv] && analysisTree.muon_isTracker[imv] && analysisTree.muon_isPF[imv]  
		   &&  analysisTree.muon_pt[imv]> 15 &&  fabs(analysisTree.muon_eta[imv])< 2.4 && fabs(analysisTree.muon_dxy[imv])<0.045 
		   && fabs(analysisTree.muon_dz[imv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imv]) && dRr > 0.15 && OSCharge) //removed from last recipe
		MuVeto=true;
	    }
	  }
	}
      }
      if (MuVeto) continue;

      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



      bool ThirdLeptVeto=false;

      if (doThirdLeptVeto){
  	if (analysisTree.electron_count>0) {
	  for (unsigned int iev = 0; iev<analysisTree.electron_count; ++iev) {

/*
	    double neutralIsoV = analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev] -  4*TMath::Pi()*(0.3*0.3)*analysisTree.rho;
	    double IsoWithEA =  analysisTree.electron_r03_sumChargedHadronPt[iev] + TMath::Max(double(0), neutralIsoV);
	

*/

           double IsoWithEA = analysisTree.electron_r03_sumChargedHadronPt[iev] 
			   + max(analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev]
		           - 0.5 * analysisTree.electron_r03_sumPUPt[iev], 0.0) ;

           double relIsoV = IsoWithEA/analysisTree.electron_pt[iev];


	    bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[iev], analysisTree.electron_superclusterEta[iev], analysisTree.electron_mva_id_nontrigPhys14[iev]);


	    if ( analysisTree.electron_pt[iev] > 10 &&  fabs(analysisTree.electron_eta[iev]) < 2.5 && fabs(analysisTree.electron_dxy[iev])<0.045
		 && fabs(analysisTree.electron_dz[iev]) < 0.2 && relIsoV< 0.3 && electronMvaId && analysisTree.electron_pass_conversion[iev] 
		 && analysisTree.electron_nmissinginnerhits[iev] <=1) ThirdLeptVeto=true;

	  }
	}


     	if (analysisTree.muon_count>0){
	  for (unsigned int imvv = 0; imvv<analysisTree.muon_count; ++imvv) {

	    //       if ( imvv != mu_index  && analysisTree.muon_charge[imvv] != analysisTree.muon_charge[mu_index] ){

         double absIso= analysisTree.muon_r03_sumChargedHadronPt[imvv]
	            + max(analysisTree.muon_r03_sumNeutralHadronEt[imvv] + analysisTree.muon_r03_sumPhotonEt[imvv]
	            - 0.5 * analysisTree.muon_r03_sumPUPt[imvv],0.0);

	double relIso = absIso/analysisTree.muon_pt[imvv];


	    if ( imvv != mu_index &&  analysisTree.muon_isMedium[imvv] &&  analysisTree.muon_pt[imvv]> 10 &&  fabs(analysisTree.muon_eta[imvv])< 2.4 && fabs(analysisTree.muon_dxy[imvv])<0.045 
		 && fabs(analysisTree.muon_dz[imvv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imvv]) ) ThirdLeptVeto=true;
	  }
	}
      }
  
      if (ThirdLeptVeto) continue;


      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      //	for (unsigned int j=0;j<LeptMV.size();++j) cout<<" j "<<j<<"  "<<LeptMV.at(j).Pt()<<endl;
      //	cout<<""<<endl;
      ////////jets cleaning 
      double DRmax=0.4;
      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;
      
      //      continue;
      
      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      for (unsigned int il = 0; il<LeptMV.size(); ++il) {
      
	for (unsigned int ij = 0; ij<JetsMV.size(); ++ij) {
        
	//  if(fabs(JetsMV.at(ij).Eta())>etaJetCut) continue;
       //   if(fabs(JetsMV.at(ij).Pt())<ptJetCut) continue;
      

	  double Dr= deltaR(LeptMV.at(il).Eta(), LeptMV.at(il).Phi(),JetsMV.at(ij).Eta(),JetsMV.at(ij).Phi());

	  if (  Dr  < DRmax) {
	     
	    JetsMV.erase (JetsMV.begin()+ij);
	  }	
		       
	}
      }
      


      double ptScalarSum = -1;


      bool btagged= false;
      bool JetsPt30C =false;
      

      if (JetsMV.size() >3) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      int xj = -1;
      
      for (unsigned int im = 0; im <JetsMV.size();++im)
	
	if (JetsMV.at(im).Pt()>30 ) JetsPt30C = true;
    
      //if (JetsPt30C ) continue;
      //if (JetsMV.size() >3) continue;
    /*
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
*/
      for (unsigned int ib = 0; ib <JetsMV.size();++ib){
	
	//	if (JetsMV.at(ib).Pt()>30 ) JetsPt30C = true;
	for (unsigned int il = 0; il < analysisTree.pfjet_count; ++il)
	  {
	
	    if (double(JetsMV.at(ib).Pt()) == double(analysisTree.pfjet_pt[il])) {xj=il;
	      //cout<<" found Jets "<<analysisTree.pfjet_pt[il]<<"  "<<JetsMV.at(ib).Pt()<<"  "<<ib<<"  "<<xj<<"  "<<il<<endl;
	    }
	  }

	//     cout<<" Jets "<<analysisTree.pfjet_pt[xj]<<"  "<<JetsMV.at(ib).Pt()<<"  "<<ib<<"  "<<xj<<endl;
	if (analysisTree.pfjet_btag[xj][8]  > bTag) btagged = true;
      }
      
      if (btagged ) continue;
	
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
  
  
      // pt Scalar
      // computations of kinematic variables

      //cout<<"  "<<mu_index<<"  "<<tau_index<<"   "<<MuMV.at(mu_index).M()<<"  "<<TauMV.at(tau_index).M()<<endl;

      TLorentzVector diL = MuMV.at(mu_index) + TauMV.at(tau_index);
      if ( diL.M() <80 && diL.M()>40 ) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      double ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // bisector of electron and muon transverse momenta

      // computation of MT variable
      double dPhi=-999; 


      dPhi=dPhiFrom2P( LeptMV.at(0).Px(), LeptMV.at(0).Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      //MT = TMath::Sqrt(2*LeptMV.at(0).Pt()*ETmiss*(1-TMath::Cos(dPhi)));


      // filling histograms after dilepton selection

      
      // ETmissH->Fill(ETmiss,weight);
      // MtH->Fill(MT,weight);
      
      if (ETmiss < 85) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
       
       double MTv = mT(diL,METV);
     /* if (ETmiss < 100) continue;
      //if (MTv>80 && MTv<120) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      if (ETmiss < 120) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      // topological cut
      //if (DZeta<dZetaCut) continue;
       */

       TLorentzVector muVc ; MuV.SetPtEtaPhiM(analysisTree.muon_pt[mu_index], analysisTree.muon_eta[mu_index], analysisTree.muon_phi[mu_index], muonMass);
       TLorentzVector tauVc; TauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_index], analysisTree.tau_eta[tau_index], analysisTree.tau_phi[tau_index], tauMass);

       TLorentzVector Dil =  muVc + tauVc;
      double dRr = deltaR(Dil.Eta(), Dil.Phi(), METV.Eta(), METV.Phi());
      
      if (dRr<1.5 || dRr>4) continue;
      
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      /*
      if (dPhi<1) continue; 
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      //      std::cout << std::endl;
      */
      selEvents++;
      //histWeights->SetBinContent(1,Weight+ histWeights->GetBinContent(1));  
    } // end of file processing (loop over events in one file)
    //histWeights->Fill(Weight);
  
    //cout<< " Weight  "<<Weight<<  "   histWeigh Sum " <<histWeights->GetBinContent(1)<<"  "<<histWeights2->GetSumOfWeights()<<"   "<<" Events "<<inputEventsH->GetSum()<<"  "<<inputEventsH->GetSumOfWeights()<<endl;
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;
	
  cout<<" Total events  "<<nEvents<<"  Will use weight  "<<histWeights->GetSumOfWeights()<<" Norm Factor "<<XSec*Lumi/( histWeights->GetSumOfWeights())<<endl;
  cout<<" First content "<<CFCounter[0]<<endl;
  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;
  /*
    for (int i=0;i<CutNumb;++i){
    CFCounter[i] *= double(XSec*Lumi/( histWeights->GetSumOfWeights()));
    if (iCFCounter[i] <0.2) statUnc[i] =0;
    else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
    }
  */
  //write out cutflow
  ofstream tfile;
  // TString outname = argv[argc-1];
  TString outname=argv[2];
  TString textfilename = "cutflow_"+outname+"_"+Channel+"_"+argv[3]+".txt";
//  tfile.open(textfilename);
//  tfile << "########################################" << endl;
  for(int ci = 0; ci < CutNumb; ci++)
    {
     // tfile << CutList[ci]<<"\t & \t"
//	    << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
//      if (!isData) { 
    cout << " i "<<ci<<" "<<iCFCounter[ci]<<"  "<<XSec*Lumi/( histWeights->GetSumOfWeights())<<endl;  
   
      CutFlowUnW->SetBinContent(1+ci,CFCounter[ci] );
      CFCounter[ci] *= double(XSec*Lumi/( histWeights->GetSumOfWeights()));
    //}
      CutFlow->SetBinContent(1+ci,CFCounter[ci]);
      
      if (iCFCounter[ci] <0.2) statUnc[ci] =0;
    //else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
    else statUnc[ci] = sqrt(CFCounter[ci]);
  }
  //ofstream tfile1;
  //TString textfile_Con = "CMG_cutflow_Con_Mu_"+outname+".txt";
  //tfile1.open(textfile_Con);
  //tfile1 << "########################################" << endl;
  //tfile1 << "RCS:" << endl;



  //tfile << "Cut efficiency numbers:" << endl;

  //    tfile << " Cut "<<"\t & \t"<<"#Evnts for "<<Lumi/1000<<" fb-1 & \t"<<" Uncertainty \t"<<" cnt\t"<<endl;

 // tfile.close();



  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd(Channel.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeights->Write();
  histWeightsH->Write();
  histRuns->Write();
  //histWeights2->Write();
  
  CutFlow->Write();
  CutFlowUnW->Write();
   /*
  CutFlowUnWNorm = (TH1D*) CutFlowUnW->Clone();
  double normCF = CutFlowUnW->GetSumOfWeights();
  CutFlowUnWNorm->Scale(1/normCF);
  CutFlowUnWNorm->Write();
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
  
  //WriteHists(CutNumb,file,"mutau");
  file->Write();
  file->Close();
   
  delete file;
  
}



