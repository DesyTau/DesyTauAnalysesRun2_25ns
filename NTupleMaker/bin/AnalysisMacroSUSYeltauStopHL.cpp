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
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
int main(int argc, char * argv[]) {



  // **** configuration
  Config cfg(argv[1]);
  string Channel="eltau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  
  

////////////muons

  const float  ptElectronCut       = cfg.get<float>("ptElectronCuteltau");
  const double etaElectronCut     = cfg.get<double>("etaElectronCuteltau");
  const double dxyElectronCut     = cfg.get<double>("dxyElectronCuteltau");
  const double dzElectronCut      = cfg.get<double>("dzElectronCuteltau");
  const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCuteltau");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

////// tau
  const double ptTauCut = cfg.get<double>("ptTauCut"); 
  const double etaTauCut = cfg.get<double>("etaTauCut"); 
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double  againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
  const double  againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double  vertexz =  cfg.get<double>("vertexz");
  const double  byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");


  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // dilemuon veto 
  const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
  const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
  const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
  const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
  const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");


//veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");



  const string dataBaseDir = cfg.get<string>("DataBaseDir");



  //const string MuonidIsoEffFileBCDEFGH = cfg.get<string>("MuonidIsoEffFileBCDEFGH");
  //const string MuontrigEffFileBCDEFGH = cfg.get<string>("MuontrigEffFileBCDEFGH");
  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");

  const double dRleptonsCuteltau   = cfg.get<double>("dRleptonsCuteltau");


  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");

  // topSingleMuonTriggerFile
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");

  const string TrigLeg  = cfg.get<string>("SingleMuonFilterName") ;
  const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");
  const float SingleElectronTriggerPtCut = cfg.get<float>("SingleElectronTriggerPtCut");


  // vertex distributions filenames and histname

  const double bTag   = cfg.get<double>("bTag");
  string cmsswBase = (getenv ("CMSSW_BASE"));
 

  xs=1;fact=1;fact2=1;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;



  bool doThirdLeptVeto=true;
  bool doMuVeto=true;
/*
// PU reweighting
  PileUp * PUofficial = new PileUp();
  //TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_Cert_271036-277148_13TeV_PromptReco_Collisions16_xsec69p2mb.root","read");
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_RunBCDEFGH_ReReco.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring16_PU25ns_V1.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);
*/
  char ff[100];
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  int nTotalFiles = 0;

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  //std::string ntupleName("Delphes");
  //std::string ntupleName("makeroottree/AC1B");
  std::string ntupleName("ntuple/AC1B");


  bool isDelphes = false;
  if (  string::npos != rootFileName.find("ntuple") ) isDelphes=true;
  if (isDelphes) ntupleName = "AC1B";
	
  std::string initNtupleName("ntuple/AC1B");

  TString era=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<std::endl;  
  // output fileName with histograms


  TFile * file;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
 /* if (SUSY)
	{  
	TString TStrNameS(rootFileName+invMuStr+invTauStr+invMETStr+"_"+Region+"_"+Sign+"_"+st1+"_"+st2);
	file = new TFile(era+"/"+TStrNameS+TString(".root"),"update");

	}*/
  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  //string treename = rootFileName+"_tree.root";
  Int_t CutNumb=10;
  SetupTree(); 
  //SetupHists(CutNumb); 
  //std::cout <<" Test 1  "<<std::endl;  
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
 
for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));


bool WithInit = false;

if (WithInit) cout << "With initroottree AAAAAAAA"<<endl;
if (!WithInit) cout << "Without initroottree"<<endl;


    TTree * _inittree = NULL;
if (!WithInit)  _inittree = (TTree*)file_->Get(TString(ntupleName));
if (WithInit)  _inittree = (TTree*)file_->Get(TString(initNtupleName));

    if (_inittree==NULL) continue;
    Float_t genweight;
    if (!isData)
      _inittree->SetBranchAddress("genweight",&genweight);
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
    std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; ++iEntry) {
      _inittree->GetEntry(iEntry);
      if (isData)
	histWeightsH->Fill(0.,1.);
      //else
      //histWeightsH->Fill(0.,genweight);
    }


    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);
	if (!isData && !WithInit)
	//if (!isData)
		{    
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
			{
			analysisTree.GetEntry(iEntry);
		/*	if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
			histWeightsH->Fill(0.,analysisTree.genweight);
			}
		}
  	float genweights=1.;

	//bool isWJ=false;
      //if (  string::npos != filen.find("WToLNu") ) isWJ=true;



	TFile * fileFakeJets = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/test/mutau/JetFake/WjetsFakesNewBinning.root","read");

	TH2D * FakeJets = (TH2D*)fileFakeJets->Get("eltau/FakeJets");
	TH2D * AllJets = (TH2D*)fileFakeJets->Get("eltau/AllJets");
	
	TFile * SFfile = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/ElectronTight_PTEta.root","read");
	TH2D * SFhisto = (TH2D*)SFfile->Get("FullSimOverDelphes");


    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

histCutFlow->Fill(1,analysisTree.genweight);

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
/*	if (SUSY)
		{
		if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
		&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
		}*/
      nEvents++;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%1000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
	if (!isDelphes)
	{	
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  
	}

      histCutFlow->Fill(2,analysisTree.genweight);
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;



	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;
					

/*
      if (!isData) 
	{
	    puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	//	puweight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
	    weight *=puweight; 
	    pu_weight = puweight;
	}
*/


/*
      if (isData){
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}

      }
	}//if isData check for filters

      if (!isData ) isMainTrigger = true;

      if (!isMainTrigger) {
	std::cout << "HLT filter for Mu Trigger " << MainTrigger << " not found" << std::endl;
	return(-1);
      }
*/
	
//      ptMuonCut = SingleMuonTriggerPtCut;

    vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<SingleElectronTriggerPtCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	//if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	//if (fabs(analysisTree.muon_dz[ie])>dzMuonCut) continue;
	if (!isDelphes && !analysisTree.electron_cutId_tight_Summer16[ie]) continue;
	if ( fabs(analysisTree.electron_charge[ie]) != 1) continue;

	electrons.push_back((int)ie);
std::cout <<analysisTree.electron_pt[ie]<<"  "<<fabs(analysisTree.electron_eta[ie])<<"  "<<analysisTree.electron_cutId_tight_Summer16[ie]<<"  "<<analysisTree.electron_charge[ie]<<"  " <<electrons.size()<<std::endl;

      }
      
      if (electrons.size()==0) continue;
      
histCutFlow->Fill(3,analysisTree.genweight);

      vector<int> taus; taus.clear();
      vector<int> tausDelphes; tausDelphes.clear();

      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (!isDelphes && analysisTree.tau_decayModeFinding[it]<-1) continue;
	if (!isDelphes && analysisTree.tau_decayModeFinding[it]<6.1 && analysisTree.tau_decayModeFinding[it]>5.9) continue;
	if (!isDelphes && analysisTree.tau_decayModeFinding[it]<5.1 && analysisTree.tau_decayModeFinding[it]>4.9) continue;
	if (!isDelphes && analysisTree.tau_chargedIsoPtSum[it]<-1) continue;
	if (!isDelphes && analysisTree.tau_chargedIsoPtSum[it]>2.5) continue;

	//if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
	//if (!isDelphes && analysisTree.tau_byIsolationMVArun2v1Medium[it]<0.25) continue;
        if ( fabs(analysisTree.tau_charge[it]) != 1 ) continue;
	  taus.push_back((int)it);
std::cout <<analysisTree.tau_pt[it]<<"  "<<fabs(analysisTree.tau_eta[it])<<"  "<<analysisTree.tau_charge[it]<<"  " <<taus.size()<<std::endl;

	}

      if (taus.size()==0)  continue;

histCutFlow->Fill(4,analysisTree.genweight);

      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float isoElMin  = 1e+10;
      float isoTauMin = 1.; 
      float isoTau = 1.; 
	bool CutBasedTauId = false;
      if (CutBasedTauId) isoTauMin = 1e+10;
      if (!CutBasedTauId) isoTauMin = -10;
      float ptEl = 0;
      float ptTau = 0;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;
      
	bool isLegMatch = false;
	for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	//isLegMatch = false;
	//	bool isMuonTauMuonLegMatch = false;
	//	bool isMuonTauOverlapMuonMatch = false;
	unsigned int eIndex  = electrons.at(ie);

	float relIsoElec = analysisTree.electron_relIso[eIndex];

	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex]);
cout <<"dR "<<dR <<endl;
	  if (dR<dRleptonsCuteltau) continue;
histCutFlow->Fill(5,analysisTree.genweight);
 

if (!CutBasedTauId){
   //isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
   //isoTau = analysisTree.tau_chargedIsoPtSum[tIndex];

          if ((int)eIndex!=(int)el_index) {
            if (relIsoElec==isoElMin) {
              if (analysisTree.electron_pt[eIndex]>ptEl) {
                isoElMin  = relIsoElec;
                ptEl = analysisTree.electron_pt[eIndex];
                el_index =(int)eIndex;
                //isoTauMin = isoTau;
                ptTau = analysisTree.tau_pt[tIndex];
                tau_index =(int)tIndex;
              }
            }
            else if (relIsoElec<isoElMin) {
              isoElMin  = relIsoElec;
              ptEl = analysisTree.electron_pt[eIndex];
              el_index =(int)eIndex;
              //isoTauMin = isoTau;
              ptTau = analysisTree.tau_pt[tIndex];
              tau_index =(int)tIndex;
            }
          }
          else {
            //if (isoTau==isoTauMin) {
              if (analysisTree.tau_pt[tIndex]>ptTau) {
                ptTau = analysisTree.tau_pt[tIndex];
                //isoTauMin = isoTau;
                tau_index =(int)tIndex;
              //}
            }
          }
	}

      }
 }




      //int tau_index = taus.at(0);//////////////// addeddddddddddddd
	if (isDelphes)
	{
	for (unsigned int it=0; it<taus.size(); ++it) 	
		{
				
	  	unsigned int tIndex = taus.at(it);

	  	float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index]);
	 	if (dR<dRleptonsCuteltau) continue;
		tausDelphes.push_back(tIndex);

		}
	}


						

      bool TauId = false;
cout <<"tau_index "<<tau_index <<endl;
//cout <<"el_index "<<el_index <<endl;
      if ((int)tau_index<0) continue;
//      if ((int)el_index<0) continue;
histCutFlow->Fill(7,analysisTree.genweight);
	


      el_relIso[0]=isoElMin;
      
      double q = analysisTree.tau_charge[tau_index] * analysisTree.electron_charge[el_index];
      event_sign  = q;

	double dReltau = deltaR(analysisTree.tau_eta[(int)tau_index],analysisTree.tau_phi[(int)tau_index],
				analysisTree.electron_eta[(int)el_index],analysisTree.electron_phi[(int)el_index]);
	if (dReltau < 0.5) continue;
histCutFlow->Fill(8,analysisTree.genweight);

	TFR_weight=1;



	//weight *= TFR_weight;


  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;

  event_secondLeptonVeto = false;
  event_thirdLeptonVeto = false;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<20) continue;
	if (fabs(analysisTree.electron_eta[ie])>2.4) continue;
	//if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	//if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	//bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
	if (!isDelphes && analysisTree.electron_cutId_loose_Summer16[ie]) continue;
	//if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
	//if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
	if (analysisTree.electron_relIso[ie]>0.5) continue;
	foundExtraElectron = true;
      }

      // looking for extra muon's (dimuon veto)
      bool foundExtraMuon = false;
      vector<int> mu_dimuons; mu_dimuons.clear(); 
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if ((int)im==(int)mu_index) continue;

	float relIsoMu = analysisTree.muon_relIso[im];

	 // mu_dimuons.push_back(im);

	if (analysisTree.muon_pt[im]<20) continue;
	if (fabs(analysisTree.muon_eta[im])>2.4) continue;
	//if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	//if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (!isData && !isDelphes && applyVetoMuonId && !analysisTree.muon_isLoose[im]) continue;
	if (relIsoMu>0.5) continue;
	foundExtraMuon = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;
    

   	event_secondLeptonVeto = dilepton_veto;
//	if (dilepton_veto)  continue;


//	if (extraelec_veto) continue;
//	if (extramuon_veto) continue;

      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;




///////////////Trigger weight
 /*
      double ptMu1 = (double)analysisTree.muon_pt[mu_index];
      double etaMu1 = (double)analysisTree.muon_eta[mu_index];
      float trigweight=1.;

      float EffFromData = 1.;
      float EffFromDataA = 1.;
      float EffFromDataB = 1.;
      float Lumi,LumiA,LumiB;
      Lumi=36590.;
      LumiA = 16357.;
      LumiB = 20233.;
      
      
      if (!isData) { 
	      
	//EffFromDataA = (float)SF_muonTriggerBCDEF->get_EfficiencyData(double(ptMu1),double(etaMu1));
	//EffFromDataB = (float)SF_muonTriggerGH->get_EfficiencyData(double(ptMu1),double(etaMu1));
	EffFromData = (float)SF_muonTrigger->get_EfficiencyData(double(ptMu1),double(etaMu1));



	
      //bool Signal = true;

	//if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) Signal=true;
 
	//if (!isData) trigweight = EffFromDataA * LumiA/Lumi + EffFromDataB * LumiB/Lumi;
      trigweight = EffFromData;
      weight *= trigweight;
      trig_weight = trigweight;
*/
	LSF_weight=1;
	if (isDelphes)
	{
	///LSF 

	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
	float IdIsoSF_el = 1;
	double ptEl1 = (double)analysisTree.electron_pt[el_index];
    double etaEl1 = (double)analysisTree.electron_eta[el_index];
    IdIsoSF_el = SFhisto->GetBinContent(FakeJets->GetXaxis()->FindBin(ptEl1),FakeJets->GetYaxis()->FindBin(etaEl1));

	//if (iEntry%2==0) IdIsoSF_mu1=  SF_muonIdIsoBCDEF->get_ScaleFactor(ptMu1, etaMu1);
	//if (iEntry%2!=0) IdIsoSF_mu1 = SF_muonIdIsoGH->get_ScaleFactor(ptMu1, etaMu1);
	// IdIsoSF_mu =  SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
	 //IdIsoSF_mu1 =  SF_muonIdIsoBCDEF->get_ScaleFactor(ptMu1, etaMu1);
	 //IdIsoSF_mu2 = SF_muonIdIsoGH->get_ScaleFactor(ptMu1, etaMu1);

	//LSF_weight = IdIsoSF_mu1 * LumiA/Lumi + IdIsoSF_mu2 * LumiB/Lumi;
	LSF_weight = IdIsoSF_el;
	weight *= LSF_weight;
      }
	


/*

      bool isTauMatched = false;
      bool isGenLeptonMatched = false;
      bool isGenLeptonMatchedMu = false;
      bool isGenLeptonMatchedEl = false;
      bool isGenTauDecayedElMatched = false;
      bool isGenTauDecayedMuMatched = false;
      if (!isData && !isDelphes){
	TLorentzVector genTauV;  
	TLorentzVector genLepV;  

	for (unsigned int gt = 0 ; gt < analysisTree.gentau_count; ++gt){

	 // genTauV.SetXYZT(0.,0.,0.,0.);
	  genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);


	  double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    genTauV.Eta(), genTauV.Phi());

//cout << "Drr  "<<Drr<<endl;  
	  if (Drr < 0.2 && analysisTree.gentau_isPrompt[gt] > 0.5  && genTauV.Pt() > 15. ) genTauMatched = true;

	}
      
      
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      		  if ( (abs(analysisTree.genparticles_pdgid[igen])==11 || abs(analysisTree.genparticles_pdgid[igen])==13) ){

	  genLepV.SetXYZT(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen], analysisTree.genparticles_e[igen]);

	  double Drm=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  genLepV.Eta(),genLepV.Phi());
//cout << "Drm  "<<Drm<<endl;  
//cout << "analysisTree.genparticles_pdgid[igen]  "<<analysisTree.genparticles_pdgid[igen]<<endl;  

		if (Drm < 0.2 && genLepV.Pt() > 8. ) genLeptonMatched = true;
		if (Drm < 0.2 && genLepV.Pt() > 8. && abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptEl = true;
		if (Drm < 0.2 && genLepV.Pt() > 8. && abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptMu = true;
		if (Drm < 0.2 && genLepV.Pt() > 8. && abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecEl = true;
		if (Drm < 0.2 && genLepV.Pt() > 8. && abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecMu = true;

		}
      
	  }
      
      
      }//!isData
*/

      //////////////////////////////////////////////
      //muon_index = (int)mu_index;
      //electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)analysisTree.muon_count;
      //cout<<" here ============================> "<<iEntry<<"  "<<mu_count<<"  "<<(int)analysisTree.muon_count<<"  "<<analysisTree.muon_count<<endl;
      for (unsigned int im=0;im<analysisTree.muon_count; ++im){
	mu_px[im]=analysisTree.muon_px[im];
	mu_py[im]=analysisTree.muon_py[im];
	mu_pz[im]=analysisTree.muon_pz[im];
	mu_eta[im]=analysisTree.muon_eta[im];
	mu_pt[im]=analysisTree.muon_pt[im];
	mu_phi[im]=analysisTree.muon_phi[im];
	mu_charge[im]=analysisTree.muon_charge[im];
	//mu_dxy[im]=analysisTree.muon_dxy[im];
	//mu_dz[im]=analysisTree.muon_dz[im];
	//mu_dxyerr[im]=analysisTree.muon_dxyerr[im];
	//mu_dzerr[im]=analysisTree.muon_dzerr[im];

	mu_relIsoMu[im]  = analysisTree.muon_relIso[im] ;
   
     	}



      el_count=(int)analysisTree.electron_count;
      for (unsigned int ie=0;ie<analysisTree.electron_count; ++ie){
	el_px[ie]=analysisTree.electron_px[ie];
	el_py[ie]=analysisTree.electron_py[ie];
	el_pz[ie]=analysisTree.electron_pz[ie];
	el_eta[ie]=analysisTree.electron_eta[ie];
	el_pt[ie]=analysisTree.electron_pt[ie];
	el_phi[ie]=analysisTree.electron_phi[ie];
	el_charge[ie]=analysisTree.electron_charge[ie];
	//el_dxy[ie]=analysisTree.electron_dxy[ie];
	//el_dz[ie]=analysisTree.electron_dz[ie];
	//el_dxyerr[ie]=analysisTree.electron_dxyerr[ie];
	//el_dzerr[ie]=analysisTree.electron_dzerr[ie];

	el_relIsoEl[ie]  = analysisTree.electron_relIso[ie] ;

      }

				
      ta_count=(int)analysisTree.tau_count;
      for (unsigned int it=0;it<analysisTree.tau_count; ++it){
	//ta_mass[it]=analysisTree.tau_mass[it];
	ta_px[it]=analysisTree.tau_px[it];
	ta_py[it]=analysisTree.tau_py[it];
	ta_pz[it]=analysisTree.tau_pz[it];
	ta_eta[it]=analysisTree.tau_eta[it];
	ta_pt[it]=analysisTree.tau_pt[it];
	ta_phi[it]=analysisTree.tau_phi[it];
	ta_charge[it]=analysisTree.tau_charge[it];
	ta_relIso[it]=analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];

	ta_isLoose[it]=analysisTree.tau_againstElectronLooseMVA6[it];
	ta_isMedium[it]=analysisTree.tau_againstElectronMediumMVA6[it];
	ta_isTight[it]=analysisTree.tau_againstElectronTightMVA6[it];


	//ta_dxy[it]=analysisTree.tau_dxy[it];
	//ta_dz[it]=analysisTree.tau_dz[it];
     	//ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[it];
     	ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[it];
     	//ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[it];



      }
      jet_count=(int)analysisTree.pfjet_count;
      for (unsigned int jj=0;jj<analysisTree.pfjet_count; ++jj){

	jet_e[jj] = analysisTree.pfjet_e[jj];
	jet_px[jj] = analysisTree.pfjet_px[jj];
	jet_py[jj] = analysisTree.pfjet_py[jj];
	jet_pz[jj] = analysisTree.pfjet_pz[jj];
	jet_pt[jj] = analysisTree.pfjet_pt[jj];
	jet_eta[jj] = analysisTree.pfjet_eta[jj];
	jet_phi[jj] = analysisTree.pfjet_phi[jj];
	jet_flavour[jj] = analysisTree.pfjet_flavour[jj];
	jet_btag[jj] = analysisTree.pfjet_btag[jj][0];
      }



      ////////jets cleaning 
      TLorentzVector leptonsV, muonJ, jetsLV;


      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);




      float jetEta = 2.4;
      float DRmax = 0.5;
      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;

	int counter_cleaned_jets = 0;


      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


	//bool isPFJetId = false ; 
      	bool btagged= false;
	//isPFJetId =looseJetiD(analysisTree,jet);
	//isPFJetId =tightLepVetoJetiD(analysisTree,jet);

	//if (!isPFJetId) continue;
	bool cleanedJet = true;

	double Dr=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	//if (  Dr  < DRmax)  cleanedJet=false;


	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	//if ( Drr < DRmax) cleanedJet=false;

	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

	if (analysisTree.pfjet_btag[jet][0]  > bTag) btagged = true;
	

	  if (btagged && cleanedJet) bjets.push_back(jet);
	}


	if (cleanedJet){
		
	//	cout<<"  will push to save now cleaned jet  "<<(int)jet<<"  for counter_cleaned_jet "<<(int)counter_cleaned_jets<<" event "<<iEntry<<endl;

	jets.push_back((int)jet);
	jets_cleaned[counter_cleaned_jets]=(int)jet;
	counter_cleaned_jets++;
	}


      }///loop in all jets

      njets = jets.size();
      jet_count = jets.size();
      //njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      //nbtag_nocleaned = bjets_nocleaned.size();

      //npv =  analysisTree.primvertex_count;
      //npu = analysisTree.numtruepileupinteractions;
	
      //SusyMother = SusyMotherMassF;
      //SusyLSP = SusyLSPMassF;



      met_ex = analysisTree.pfmet_pt*TMath::Cos(analysisTree.pfmet_phi);
      met_ey = analysisTree.pfmet_pt*TMath::Sin(analysisTree.pfmet_phi);
      met_pt = analysisTree.pfmet_pt;
      met_phi = analysisTree.pfmet_phi;

      all_weight = weight;
histCutFlow->Fill(9,analysisTree.genweight);






	if (isDelphes)
	{
	vector<float> TFR_weights; TFR_weights.clear();
		for (unsigned int it=0; it<tausDelphes.size(); ++it) 	
		{
		float NFakeJets=1;
		float NAllJets=1;
		float t_pt=analysisTree.tau_pt[tausDelphes.at(it)];
		float t_eta=fabs(analysisTree.tau_eta[tausDelphes.at(it)]);
		tau_index = tausDelphes.at(it);

		if (t_pt>1000)t_pt=900;
		if (t_eta>2.3)t_eta=2.2;

		NFakeJets = FakeJets->GetBinContent(FakeJets->GetXaxis()->FindBin(t_pt),FakeJets->GetYaxis()->FindBin(t_eta));
		NAllJets = AllJets->GetBinContent(AllJets->GetXaxis()->FindBin(t_pt),AllJets->GetYaxis()->FindBin(t_eta));
		TFR_weight = NFakeJets/NAllJets;
		TFR_weights.push_back(TFR_weight);
		if (it == 1) {TFR_weight = TFR_weight*(1-TFR_weights.at(0)); }
		if (it == 2) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1)); }
		if (it == 3) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1))*(1-TFR_weights.at(2)); }
		if (it == 4) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1))*(1-TFR_weights.at(2))*(1-TFR_weights.at(3)); }
		if (it > 4.5) continue;

		cout <<"FakeJets->GetSumOfWeights()  "<< FakeJets->GetSumOfWeights()<<"  FakeJets->GetXaxis()->FindBin(analysisTree.tau_pt[(int)tau_index])  "<< FakeJets->GetXaxis()->FindBin(t_pt)<<"  FakeJets->GetYaxis()->FindBin(analysisTree.tau_eta[(int)tau_index])  "<< FakeJets->GetYaxis()->FindBin(t_eta)<<endl;
		cout <<"SFFakeRate  "<< TFR_weight<<"  NFakeJets  "<< NFakeJets<<"  NAllJets  "<< NAllJets<<endl;


		T->Fill();
		selEvents++;
		}
	}






      	if (!isDelphes) T->Fill();
	
        if (!isDelphes) selEvents++;
      continue;
      /////////////////////////////////////////////////



    } // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }


  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;


  std::cout << std::endl;
  int allEvents = (int)inputEventsH->GetEntries();
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;


  file->cd(Channel.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histCutFlow->Write();
  histRuns->Write();
  //CutFlowUnW->Write();
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}
