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
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZMuMuTree.h"

#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

typedef std::vector<std::pair<int,int> > lumi_json;

struct compare_lumi { //accepts two pairs, return 1 if left.first < right.first or left.first = right.first e left.second < right.second
  bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
    if (left.first < right.first)
      return 1;
    else if (left.first > right.first)
      return 0;
    else
      return left.second < right.second;
  }
};

int read_json(std::string filename, lumi_json& json); 
bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json);
bool isGoodLumi(int run, int lumi, const lumi_json& json);
float abs_Iso(int Index, const AC1B * analysisTree); 
bool isICHEPmed(int Index, const AC1B * analysisTree); 
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData); // select medium id or ICHEP medium id for different runs
void SetTreeDefaultValues(ZMuMuTree * otree);


float topPtWeight_run1(float pt1,
		  float pt2) {

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter

  if (pt1>400) pt1 = 400;
  if (pt2>400) pt2 = 400;

  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  return TMath::Sqrt(w1*w2);

}


float topPtWeight_run2(float pt1,
		  float pt2) {
    
  if (pt1>400) pt1 = 400;
  if (pt2>400) pt2 = 400;
    
  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter
    
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
    
  return TMath::Sqrt(w1*w2);  
}


void computeRecoil(float metx, float mety,
		   float unitX,float unitY,
		   float perpUnitX, float perpUnitY,
		   float dimuonPt,
		   float & recoilParal,
		   float & recoilPerp,
		   float & responseHad) {

  recoilParal = metx*unitX + mety*unitY;
  recoilPerp = metx*perpUnitX + mety*perpUnitY;
  responseHad = 1 + recoilParal/dimuonPt;
}


int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist
  // third argument - optional index of the first file to be analyzed
  // fourth argument - optional index of the last file to be analyzed

  using namespace std;

  bool printDebugMessages = false;

  // **** configuration
  Config cfg(argv[1]);
  const string infiles = argv[2];

  const bool isData = cfg.get<bool>("isData");
  const string year = cfg.get<string>("year");

  string cmsswBase = (getenv ("CMSSW_BASE"));

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  // trigger
  const bool applyTrigger    = cfg.get<bool>("ApplyTrigger");
  const string hltFilterName = cfg.get<string>("HLTFilterName");
  TString HLTFilterName(hltFilterName);

  // recoil corrections
  const bool applyRecoilCorrections                 = cfg.get<bool>("ApplyRecoilCorrections");
  const bool applyRecoilOnGenerator                 = cfg.get<bool>("ApplyRecoilOnGenerator");
  const bool applyRecoilCorrectionsByMeanResolution = cfg.get<bool>("ApplyRecoilCorrectionsByMeanResolution");
  const string recoilFileName                       = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  // lepton scale factors
  const bool applyLeptonSF             = cfg.get<bool>("ApplyLeptonSF");
  const string MuonIdIsoFile           = cfg.get<string>("MuonIdIsoEffFile");
  const string MuonTrigFile            = cfg.get<string>("MuonTrigEffFile"); 
  const string correctionWorkspaceFile = cfg.get<string>("CorrectionWorkspaceFile");
  TString CorrectionWorkspaceFile(correctionWorkspaceFile);
  
  // pileup reweghting
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");
  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  // other weights 
  const bool applyTopPtReweighting      = cfg.get<bool>("ApplyTopPtReweighting");
  bool applyRun1TopPtWeights = false;
  if (applyTopPtReweighting){
    const bool ApplyRun1TopPtWeights    = cfg.get<bool>("ApplyRun1TopPtWeights");
    applyRun1TopPtWeights = ApplyRun1TopPtWeights; 
  }
  const bool applyZMassPtReweighting    = cfg.get<bool>("ApplyZMassPtReweighting");
  const bool interpolateZMassPtWeight   = cfg.get<bool>("InterpolateZMassPtWeight");
  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);


  // muon kinematic cuts
  const float leadingMuonPtCut   = cfg.get<float>("leadingMuonPtCut");
  const float leadingMuonEtaCut  = cfg.get<float>("leadingMuonEtaCut");
  const float leadingMuonDxyCut  = cfg.get<float>("leadingMuonDxyCut");
  const float leadingMuonDzCut   = cfg.get<float>("leadingMuonDzCut");
  const float leadingMuonIsoCut  = cfg.get<float>("leadingMuonIsoCut");
  const float subleadMuonPtCut   = cfg.get<float>("subleadMuonPtCut");
  const float subleadMuonEtaCut  = cfg.get<float>("subleadMuonEtaCut");
  const float subleadMuonDxyCut  = cfg.get<float>("subleadMuonDxyCut");
  const float subleadMuonDzCut   = cfg.get<float>("subleadMuonDzCut");

  // jet cuts
  const float jetEtaCut      = cfg.get<float>("JetEtaCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
  const float dimuonMassCut  = cfg.get<float>("DimuonMassCut");

  // **** end of configuration

  //file list reading 
  int ifile = 0;
  int jfile = -1;

  if (argc > 3)
    ifile = atoi(argv[3]);
  if (argc > 4)
    jfile = atoi(argv[4]);

  // create input files list
  std::vector<std::string> fileList;  
  if (infiles.find(".root") != std::string::npos){
    ifile = 0;
    jfile = 1;

    fileList.push_back(infiles);
  }
  else{
    ifstream input;
    std::string infile;
    
    input.open(infiles);

    while(true){
      input>>infile;
      if(!input.eof()){
	if (infile.length() > 0)
	  fileList.push_back(infile);
      }
      else
	break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }  

  //output root file inizialization
  const string sample = argv[2];
  TString rootFileName(sample);
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_MuMuTree.root";
  
  std::string ntupleName("makeroottree/AC1B");
  
  TH1::SetDefaultSumw2(true);

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);

  // reweighting official recipe 
  // initialize pile up object
  PileUp * PUdistributions = new PileUp();
  
  if (applyPUreweighting) {
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read"); 
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read"); 
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUdistributions->set_h_data(PU_data); 
    PUdistributions->set_h_MC(PU_mc);
  }

  // HTT Met recoil corrections

  RecoilCorrector recoilPFMetCorrector(RecoilFileName);

  // Lepton Scale Factors 
  ScaleFactor * SF_muonIdIso; 
  ScaleFactor * SF_muonTrig;
  if (applyLeptonSF) {
    SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
  }

  // tracking efficiency SF
  TFile * workspaceFile = new TFile(TString(cmsswBase)+"/src/"+CorrectionWorkspaceFile);
  RooWorkspace *correctionWS = (RooWorkspace*)workspaceFile->Get("w");

  // Z mass pt weights
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName); 
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "/src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName); 
  if (histZMassPtWeights==NULL) {
    std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName << std::endl;
    exit(-1);
  }


  file->cd("");  
  TTree * tree = new TTree("ZMuMu","ZMuMu");
  ZMuMuTree *otree = new ZMuMuTree(tree);

  int nTotalFiles = 0;
  int nFiles = 0;
  int nEvents = 0;
  int nSelEvents = 0;

  for (int iF=ifile; iF<jfile; ++iF) {

    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  

    if (_tree==NULL) continue; 
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");

    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();


    // EVENT LOOP //
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) cout << "      processed " << nEvents << " events" << endl; 

      if (printDebugMessages) std::cout << "      processed " << nEvents << " events" << std::endl; 

      otree->run = int(analysisTree.event_run);
      otree->lumi = int(analysisTree.event_luminosityblock);
      otree->evt = int(analysisTree.event_nr); 
      otree->npv = int(analysisTree.primvertex_count);

      SetTreeDefaultValues(otree);
  
      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;

      // Generator-level analysis
      TLorentzVector genZ; genZ.SetXYZT(0,0,0,0); 
      TLorentzVector genV; genV.SetXYZT(0,0,0,0);
      TLorentzVector genL; genL.SetXYZT(0,0,0,0);
      TLorentzVector genMuonsLV; genMuonsLV.SetXYZT(0,0,0,0);
      TLorentzVector genElectronsLV; genElectronsLV.SetXYZT(0,0,0,0);
      TLorentzVector genTausLV; genTausLV.SetXYZT(0,0,0,0);
      TLorentzVector genVisTausLV; genVisTausLV.SetXYZT(0,0,0,0);
      std::vector<unsigned int> genMuons; genMuons.clear();
      std::vector<unsigned int> genElectrons; genElectrons.clear();
      std::vector<unsigned int> genTaus; genTaus.clear();
      std::vector<unsigned int> genVisTaus; genVisTaus.clear();

      if (!isData) {

	    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	      TLorentzVector genPart; genPart.SetXYZT(analysisTree.genparticles_px[igen],
						  analysisTree.genparticles_py[igen],
						  analysisTree.genparticles_pz[igen],
						  analysisTree.genparticles_e[igen]);
	      if (analysisTree.genparticles_pdgid[igen]==23) {
	        if (analysisTree.genparticles_fromHardProcess[igen])
	           genZ.SetXYZT(analysisTree.genparticles_px[igen],
			      analysisTree.genparticles_py[igen],
			      analysisTree.genparticles_pz[igen],
			      analysisTree.genparticles_e[igen]);
	      }

	      bool isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
	      bool isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
	      bool isLepton = isMuon || isElectron;
	      bool isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	           fabs(analysisTree.genparticles_pdgid[igen])==14||
	           fabs(analysisTree.genparticles_pdgid[igen])==16;
	      bool isPrompt = analysisTree.genparticles_isPrompt[igen]||
	           analysisTree.genparticles_isPromptTauDecayProduct[igen];

	      if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	        if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	          genElectrons.push_back(igen);
	          genElectronsLV += genPart;
	        }
	      }
	  
	     if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	       if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	          genMuons.push_back(igen);
	          genMuonsLV += genPart;
	       }
	     }

	     if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
	       if (isLepton&& fabs(genPart.Eta())<2.4&& genPart.Pt()>10) {
	               genV += genPart;
	               genL += genPart;
	       }
	       if (isNeutrino) 
	        genV += genPart;
	     }
	  } // end of loop on genparticles

	  if (genV.Pt()<0.1) genV.SetXYZM(0.1,0.1,0.,0.);

	  float genZPt = -1.0;
	  float genZMass = -1.0;

	  if (genElectrons.size()==2) {
	    genZPt   = genElectronsLV.Pt();
	    genZMass = genElectronsLV.M();
        otree-> is_gen_Zee = 1;
	  }
	  else if (genMuons.size()==2) {
	    genZPt   = genMuonsLV.Pt();
        genZMass = genMuonsLV.M();
        otree-> is_gen_Zmm = 1;
	  }
	  else if (genTaus.size()==2) {
        genZPt   = genTausLV.Pt();
        genZMass = genTausLV.M();
        otree-> is_gen_Ztt = 1;
      }

      otree-> genZ_mass = genZMass;
      otree-> genZ_pt = genZPt;
      otree-> genV_px = genV.Px(); 
      otree-> genV_py = genV.Py(); 
      otree-> genV_pt = genV.Pt();
      otree-> genV_mass = genV.M();  
      otree-> genV_eta  = genV.Eta();
      otree-> genV_phi  = genV.Phi();
      otree->mc_weight = analysisTree.genweight;

      if (printDebugMessages) std::cout << "Filled gen quantities" << std::endl;
  
	  float dyWeight = 1;
      float Zmass = genZMass;
      float Zpt = genZPt;
	  if (applyZMassPtReweighting) {

	    if (Zmass>50.0&& Zpt>0.0) {

          if (Zmass > (histZMassPtWeights-> GetXaxis()->GetXmax()) ) 
            Zmass = (histZMassPtWeights-> GetXaxis()->GetXmax()) - 0.1;

          if (Zpt > (histZMassPtWeights-> GetYaxis()->GetXmax()) ) 
            Zpt = (histZMassPtWeights-> GetYaxis()->GetXmax()) - 0.1;

          if (Zpt < (histZMassPtWeights-> GetYaxis()->GetXmin()) )
            Zpt = histZMassPtWeights-> GetYaxis()-> GetXmin() + 0.1;

	      if (interpolateZMassPtWeight) 
	        dyWeight = histZMassPtWeights->Interpolate(Zmass,Zpt);
	      else 
	        dyWeight = histZMassPtWeights->GetBinContent(histZMassPtWeights->FindBin(Zmass,Zpt));
	    }
      }
	  otree->zpt_weight= dyWeight;

      if (applyPUreweighting) {
        otree->pu_weight = float(PUdistributions->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      }

	  if (applyTopPtReweighting) {
	    float topPt = -1;
	    float antitopPt = -1;
	    for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	      if (analysisTree.genparticles_pdgid[igen]==6)
	        topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	      if (analysisTree.genparticles_pdgid[igen]==-6)
	        antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);  
	    }

	    if (topPt>0&&antitopPt>0) {
          float topptweight =1;

          if (applyRun1TopPtWeights)
	        topptweight = topPtWeight_run1(topPt,antitopPt);
          else 
            topptweight = topPtWeight_run2(topPt,antitopPt);
 
          otree->toppt_weight = topptweight; 
	    }
	  }
      
    } // end of if !isData
    if (printDebugMessages)  std::cout << "End of !isData" << std::endl;
      
    //find the trigger leg to be matched 
    unsigned int nfilters = analysisTree.run_hltfilters->size();
    unsigned int nMuonFilter = 0;
    bool checkMuonFilter = false;
    if(isData || applyTrigger){  
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
        if (HLTFilter==HLTFilterName) {
          nMuonFilter = i;
          checkMuonFilter = true;
        }
      }
      if (!checkMuonFilter) {
        std::cout << "HLT filter " << HLTFilterName << " not found" << std::endl;
        exit(-1);
      }
    }
      
    // muon selection
    vector<unsigned int> leadingMuons; leadingMuons.clear();
    vector<unsigned int> subleadMuons; subleadMuons.clear();

    if (printDebugMessages)   std::cout << " Muons in the event: " <<  analysisTree.muon_count << std::endl;

    // leading muons, passing cuts, id, and matched to trigger. Isolation cut not applied. 
    for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
 
      // check requirements for leading muons
      if (analysisTree.muon_pt[im]<leadingMuonPtCut) continue;
	  if (fabs(analysisTree.muon_eta[im])>leadingMuonEtaCut) continue;
	  if (fabs(analysisTree.muon_dxy[im])>leadingMuonDxyCut) continue;
	  if (fabs(analysisTree.muon_dz[im])>leadingMuonDzCut) continue;
      if (!isIdentifiedMediumMuon(im, &analysisTree, isData)) continue;

	  bool isTriggerMatched = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMuonFilter])
	    isTriggerMatched = true;
	  }
	  if (!isTriggerMatched) continue;

      leadingMuons.push_back(im);
    }

    // subleading muons, passing cuts and id. Isolation cut not applied.
    for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
    
      // check requirements for subleading muons
      if (analysisTree.muon_pt[im]<subleadMuonPtCut) continue;
	  if (fabs(analysisTree.muon_eta[im])>subleadMuonEtaCut) continue;
	  if (fabs(analysisTree.muon_dxy[im])>subleadMuonDxyCut) continue;
	  if (fabs(analysisTree.muon_dz[im])>subleadMuonDzCut) continue;
      if (!isIdentifiedMediumMuon(im, &analysisTree, isData)) continue;

      subleadMuons.push_back(im);
    }

    // loop on leading and subleading muons to find the best pair. Requires also opposite sign and dR cut. 
    // Then, fill the otree variables for leading and subleading muon
    if (printDebugMessages) std::cout << "Leading and subleading muons size: " <<  leadingMuons.size() << " , " <<     subleadMuons.size() << std::endl;

    if (leadingMuons.size() <= 0 || subleadMuons.size() <=0 ) continue; 

    int index_lead =-1;
    int index_subl =-1;

    for (unsigned int iS =0; iS<subleadMuons.size(); ++iS){
      unsigned int index1 = subleadMuons.at(iS);
      for (unsigned int iL=0; iL<leadingMuons.size(); ++iL){
        unsigned int index2 = leadingMuons.at(iL);
		if (index1==index2) continue;
        float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	    if (dR<dRleptonsCut) continue;
        if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]){
          index_lead = index1;
          index_subl = index2;
		} 
	    else {
          index_lead = index2;
          index_subl = index1;
		}
      }
    }
    if (printDebugMessages) std::cout << "Found leading and subleading muons" << std::endl; 
    if (printDebugMessages) std::cout << "Index of leading and subleading muons : " << index_lead << " , " << index_subl << std::endl; 
    // Define dimuon LV of leading and subleading muons
    TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[index_lead],
					analysisTree.muon_py[index_lead],
					analysisTree.muon_pz[index_lead],
					muonMass);

	TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[index_subl],
					analysisTree.muon_py[index_subl],
					analysisTree.muon_pz[index_subl],
					muonMass);

	TLorentzVector dimuon = mu1 + mu2;

    float sumMuonPt = (mu1.Pt()+mu2.Pt());
	float ptRatio = 0.0;
	if (sumMuonPt != 0)
	  ptRatio = (dimuon.Pt()/sumMuonPt);

    // apply dimuon mass cut
    if (dimuon.M()<dimuonMassCut) continue;

    otree-> q_1 =  analysisTree.muon_charge[index_lead];
    otree-> q_2 =  analysisTree.muon_charge[index_subl];
    otree-> pt_1 = mu1.Pt();
    otree-> pt_2 = mu2.Pt();
    otree-> eta_1 = mu1.Eta();
    otree-> eta_2 = mu2.Eta();
    otree-> phi_1 = mu1.Phi();
    otree-> phi_2 = mu2.Phi();
    otree-> iso_1 = abs_Iso(index_lead, &analysisTree) / ( otree-> pt_1);
    otree-> iso_2 = abs_Iso(index_subl, &analysisTree) / ( otree-> pt_2);
    otree-> pt_ll = dimuon.Pt(); 
    otree-> px_ll = dimuon.Px();
    otree-> py_ll = dimuon.Py();
    otree-> eta_ll = dimuon.Eta();
    otree-> phi_ll = dimuon.Phi(); 
    otree-> m_ll = dimuon.M(); 
    otree-> pt_ratio_ll = ptRatio;
    otree-> delta_phi_ll = dPhiFrom2P(mu1.Px(), mu1.Py(), mu2.Px(), mu2.Py());
    // opposite sign or same sign muons
	float q1 = analysisTree.muon_charge[index_lead];
	float q2 = analysisTree.muon_charge[index_subl];
    if (q1 * q2 < 0)  otree-> os = 1;
    else if (q1 * q2 >0) otree-> os = 0; 

 
    // idiso and trigger scale facotrs, tracking weight 
	if (!isData && applyLeptonSF) {

	  //leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
      otree-> idiso_weight_1 =  SF_muonIdIso->get_ScaleFactor(mu1.Pt(), mu1.Eta());
      otree-> idiso_weight_2 =  SF_muonIdIso->get_ScaleFactor(mu2.Pt(), mu2.Eta());

      // tracking efficiency weights 
	  correctionWS->var("m_eta")->setVal(mu1.Eta());
	  correctionWS->var("m_pt")->setVal(mu1.Pt());
	  otree-> track_weight_1 = correctionWS->function("m_trk_ratio")->getVal();
	  correctionWS->var("m_eta")->setVal(mu2.Eta());
      correctionWS->var("m_pt")->setVal(mu2.Pt());
	  otree-> track_weight_2 = correctionWS->function("m_trk_ratio")->getVal();

      // trigger weight. Trigger match required only for leading muon. 
	  double effDataTrig1 = SF_muonTrig->get_EfficiencyData(mu1.Pt(), mu1.Eta());  
	  double weightTrig = 1;
	  if (applyTrigger) {
	    double effMcTrig1 = SF_muonTrig->get_EfficiencyMC(mu1.Pt(), mu1.Eta());
	    if (effDataTrig1>0&&effMcTrig1>0) 
	      weightTrig = effDataTrig1/effMcTrig1;
	  }
	  else weightTrig = effDataTrig1;

      otree->trig_weight = weightTrig;
      if (printDebugMessages) std::cout << "Applied lepton SF" << std::endl; 
	}

    // next:  compute jet variables

    int nJets30 =0;
    int nJets20 =0;
    float HT30 = 0;
    float HT20 = 0;

	int indexLeadingJet = -1;
	float ptLeadingJet = -1;
	int indexSubLeadingJet = -1;
	float ptSubLeadingJet = -1;

    for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	  float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	  if (absJetEta>jetEtaCut) continue;
      float jetPt = analysisTree.pfjet_pt[jet];

	  // pfJetId according to the year
      bool isPFJetId = false;
      bool noisyJet = false;
      if (year.find("2016") != std::string::npos) 
        isPFJetId = tightJetiD_2016(analysisTree,int(jet));
      else if (year.find("2017") != std::string::npos){
	    isPFJetId = tightJetiD_2017(analysisTree,int(jet));
     	noisyJet = analysisTree.pfjet_pt[jet]<50 && absJetEta > 2.65 && absJetEta < 3.139;
      }

	  if (!isPFJetId) continue;
      if (noisyJet) continue;

      // exclude jets overalpping to the leading and subleading muons
	  float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.muon_eta[index_lead],analysisTree.muon_phi[index_lead]);
	  if (dR1<dRJetLeptonCut) continue;
		
	  float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.muon_eta[index_subl],analysisTree.muon_phi[index_subl]);
	  if (dR2<dRJetLeptonCut) continue;      

		
	  if (analysisTree.pfjet_pt[jet]>30) {nJets30++; HT30 += analysisTree.pfjet_pt[jet];}
	  if (analysisTree.pfjet_pt[jet]>20) {nJets20++; HT20 += analysisTree.pfjet_pt[jet];}

      if (indexLeadingJet>=0) {
        if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
          indexSubLeadingJet = jet;
          ptSubLeadingJet = jetPt;
        }
      }
      if (jetPt>ptLeadingJet) {
        indexSubLeadingJet = indexLeadingJet;
        ptSubLeadingJet = ptLeadingJet;
        indexLeadingJet = jet;
        ptLeadingJet = jetPt;
      }

	}	 // end of loop on jets

    if (printDebugMessages) std::cout << " After jets loop " << std::endl;
    otree->njets30 = nJets30;
    otree->njets20 = nJets20;
    otree->HT20 = HT20;
    otree->HT30 = HT30;

	float mjj = -1;
	float etaLeadingJet = -999;
	float phiLeadingJet = -999;
    TLorentzVector jet1, jet2;

	if (indexLeadingJet>=0) {
      jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
					       analysisTree.pfjet_py[indexLeadingJet],
					       analysisTree.pfjet_pz[indexLeadingJet],
					       analysisTree.pfjet_e[indexLeadingJet]);

      otree-> pt_jet_1 = jet1.Pt();
      otree-> eta_jet_1 = jet1.Eta();
      otree-> phi_jet_1 = jet1.Phi();
	}

	if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {

	 if (printDebugMessages) std::cout << "Found leading and subleading jets" << std::endl;
 
     jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
                     analysisTree.pfjet_py[indexSubLeadingJet],
                     analysisTree.pfjet_pz[indexSubLeadingJet],
                     analysisTree.pfjet_e[indexSubLeadingJet]);

	  mjj = (jet1+jet2).M();
      otree-> m_jj = mjj;
      otree-> pt_jet_2 = jet2.Pt();
      otree-> eta_jet_2 = jet2.Eta();
      otree-> phi_jet_2 = jet2.Phi();
      otree->delta_phi_jj = dPhiFrom2P(jet1.Px(), jet1.Py(), jet2.Px(), jet2.Py());
	}

	float visiblePx = dimuon.Px();
	float visiblePy = dimuon.Py();
	if (applyRecoilOnGenerator) {
	  visiblePx = genL.Px();
	  visiblePy = genL.Py();
	}
     
    // MET variables, w/o recoil applied.
    float pfmet_ex = analysisTree.pfmetcorr_ex;
    float pfmet_ey = analysisTree.pfmetcorr_ey;
    float pfmet_phi = analysisTree.pfmetcorr_phi;
    float pfmet = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);

    float puppimet_ex = analysisTree.puppimet_ex;
    float puppimet_ey = analysisTree.puppimet_ey;
    float puppimet_phi = analysisTree.puppimet_phi;
    float puppimet = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);
    if (printDebugMessages){ 
    std::cout << "puppimet_ex" <<  puppimet_ex << std::endl;
    std::cout << "puppimet_ey" <<  puppimet_ey << std::endl;
    std::cout << "puppimet" <<  puppimet << std::endl;
    }

    // store met values w/o recoil corrections applied
    otree-> met = pfmet;
    otree-> met_phi = pfmet_phi;
    otree-> puppimet = puppimet;
    otree-> puppimet_phi = puppimet_phi;
    if (printDebugMessages)  std::cout << "Filled MET variables before corrections " << std::endl; 

    // calculate recoil and response w/o recoil corrections applied
    float dimuonPt = dimuon.Pt(); 
	float unitX = dimuon.Px()/dimuon.Pt();
	float unitY = dimuon.Py()/dimuon.Pt();
	float phiUnit = TMath::ATan2(unitY,unitX);
	float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
	float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
	  
	if (!isData) {
      // response wrt. to gen V pT
	  float Hparal = 0;
	  float Hperp  = 0;

	  ComputeHadRecoilFromMet(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	  float response = - Hparal/genV.Pt();
      otree->met_response_gen = response;

	  ComputeHadRecoilFromMet(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	  response = - Hparal/genV.Pt();
      otree->puppimet_response_gen = response;
	  }

	float recoilParal = 0;
	float recoilPerp  = 0;
	float responseHad = 0;
	computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilParal,recoilPerp,responseHad);

    otree-> met_recoil_paral = recoilParal;
    otree-> met_recoil_perp = recoilPerp;
    otree-> met_response_had = responseHad;

	float recoilPuppiParal = 0;
	float recoilPuppiPerp  = 0;
	float responsePuppiHad = 0;
	computeRecoil(puppimet_ex,puppimet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilPuppiParal,recoilPuppiPerp,responsePuppiHad);
    otree-> puppimet_recoil_paral = recoilPuppiParal;
    otree-> puppimet_recoil_perp = recoilPuppiPerp;
    otree-> puppimet_response_had = responsePuppiHad;

    // Fill the recoil corrected quantities with their original values w/o corrections, and update them later with the corrected values. 
    otree-> met_rcorr = otree->met;
    otree-> met_phi_rcorr = otree-> met_phi;
    otree-> met_recoil_paral_rcorr = otree-> met_recoil_paral;
    otree-> met_recoil_perp_rcorr  = otree-> met_recoil_perp;
    otree-> met_response_had_rcorr = otree-> met_response_had; 

    otree-> puppimet_rcorr = otree->puppimet;
    otree-> puppimet_phi_rcorr = otree-> puppimet_phi;
    otree-> puppimet_response_had_rcorr = otree-> puppimet_response_had;        
    otree-> puppimet_recoil_paral_rcorr = otree-> puppimet_recoil_paral;
    otree-> puppimet_recoil_perp_rcorr  = otree-> puppimet_recoil_perp;
    

	if (!isData && applyRecoilCorrections) {
      if (printDebugMessages)  std::cout << "Recoil corrections going to be applied" << std::endl; 
      //correct PF MET
	  float pfmetcorr_ex = pfmet_ex;
	  float pfmetcorr_ey = pfmet_ey;

	  if (applyRecoilCorrectionsByMeanResolution) 
      // scale mean and reoslution
	    recoilPFMetCorrector.CorrectByMeanResolution(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);
	  else 
      // quantile mapping recoil correcitons 
	    recoilPFMetCorrector.Correct(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);

      // store corrected values in the tree
 	  otree-> met_phi_rcorr = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
      otree-> met_rcorr = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);

      // update pfmet_ex and pfmet_ey to the recoil corrected values
      pfmet_ex = pfmetcorr_ex;
      pfmet_ey = pfmetcorr_ey;

      // correct Puppi MET
	  float puppimetcorr_ex = puppimet_ex;
	  float puppimetcorr_ey = puppimet_ey;

	  if (applyRecoilCorrectionsByMeanResolution) 
	    recoilPFMetCorrector.CorrectByMeanResolution(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,puppimetcorr_ex,puppimetcorr_ey);
	  else 
	    recoilPFMetCorrector.Correct(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,puppimetcorr_ex,puppimetcorr_ey);

 	  otree-> puppimet_phi_rcorr = TMath::ATan2(puppimetcorr_ey,puppimetcorr_ex);
      otree-> puppimet_rcorr = TMath::Sqrt(puppimetcorr_ex*puppimetcorr_ex+puppimetcorr_ey*puppimetcorr_ey);
 
      // update pfmet_ex and pfmet_ey to the recoil corrected values
      puppimet_ex = puppimetcorr_ex;
      puppimet_ey = puppimetcorr_ey;
      if (printDebugMessages)  std::cout << "Recoil corrections have been applied " << std::endl;

      // response wrt. to gen V pT
	  float Hparal = 0;
	  float Hperp  = 0;
	  ComputeHadRecoilFromMet(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	  float response = - Hparal/genV.Pt();
      otree->met_response_gen_rcorr = response;

      Hparal = 0;
      Hperp = 0;
	  ComputeHadRecoilFromMet(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	  response = - Hparal/genV.Pt();
      otree->puppimet_response_gen_rcorr = response;

	  float recoilParal = 0;
	  float recoilPerp  = 0;
	  float responseHad = 0;
	  computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilParal,recoilPerp,responseHad);

      otree-> met_recoil_paral_rcorr = recoilParal;
      otree-> met_recoil_perp_rcorr = recoilPerp;
      otree-> met_response_had_rcorr = responseHad;

	  recoilPuppiParal = 0;
	  recoilPuppiPerp  = 0;
	  responsePuppiHad = 0;
	  computeRecoil(puppimet_ex,puppimet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilPuppiParal,recoilPuppiPerp,responsePuppiHad);
      otree-> puppimet_recoil_paral_rcorr = recoilPuppiParal;
      otree-> puppimet_recoil_perp_rcorr = recoilPuppiPerp;
      otree-> puppimet_response_had_rcorr = responsePuppiHad;
      
	}

    otree->Fill();

    nSelEvents++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Number of selected events                        = " << nSelEvents << std::endl;
  std::cout << std::endl;
  
  file->Write();
  file->Close();
  delete file;
  
  }

////// functions implementations

int read_json(std::string filename, lumi_json& json){

  std::pair <int,int> lumi;

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);

  BOOST_FOREACH(boost::property_tree::ptree::value_type &json_run, pt.get_child("")){
    int irun = atoi(json_run.first.data());
    BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_ranges, json_run.second.get_child("")){
      int ilumi[2] = {};

      int count = 0;
      BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_boundaries, lumi_ranges.second.get_child("")){
  ilumi[count] = atoi(lumi_boundaries.second.data().data());
  count++;
      }
      
      for (;ilumi[0] <= ilumi[1]; ilumi[0]++){
  lumi = std::make_pair(irun, ilumi[0]);
  json.push_back(lumi);
      }
    }
  }

  sort( json.begin(), json.end(),  compare_lumi());
  json.erase( unique( json.begin(), json.end() ), json.end() );
  
  return 0;
}


bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json){
  static compare_lumi compare;
  static std::pair<int,int> oldlumi = lumi;
  static bool old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  
  if(lumi.first != oldlumi.first || lumi.second != oldlumi.second){
    oldlumi = lumi;
    old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  }

  return old;
}

bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}

bool isICHEPmed(int Index, const AC1B * analysisTree) {
        bool goodGlob = analysisTree->muon_isGlobal[Index] && analysisTree->muon_normChi2[Index] < 3 && analysisTree->muon_combQ_chi2LocalPosition[Index] < 12
                                   && analysisTree->muon_combQ_trkKink[Index] < 20;

        bool isICHEPmedium  = analysisTree->muon_isLoose[Index] &&
                                          analysisTree->muon_validFraction[Index] >0.49 &&
                                          analysisTree->muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
        return isICHEPmedium;
}

// select medium id or ICHEP medium id for different runs, in data
// use ICHEP medium ID fro runs up to 278808 (end of Run2016F)
// on MC, always apply medium ID
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData){
  bool isGoodMuon;
  if (isData){
	if (analysisTree->event_run<=278808) isGoodMuon= isICHEPmed(Index, analysisTree);
    else isGoodMuon=analysisTree->muon_isMedium[Index]; 
	}
  else isGoodMuon=analysisTree->muon_isMedium[Index]; 
  return isGoodMuon;
}

float abs_Iso (int Index, const AC1B * analysisTree){
  float neutralHadIso, photonIso, chargedHadIso, puIso;
 
  neutralHadIso =     analysisTree->muon_r04_sumNeutralHadronEt[Index];
  photonIso =         analysisTree->muon_r04_sumPhotonEt[Index];
  chargedHadIso =     analysisTree->muon_r04_sumChargedHadronPt[Index];
  puIso =             analysisTree->muon_r04_sumPUPt[Index];
 
  float neutralIso = neutralHadIso + photonIso -0.5*puIso;
  neutralIso = TMath::Max(float(0), neutralIso);
  return(chargedHadIso + neutralIso);
}

void SetTreeDefaultValues (ZMuMuTree * otree){

  otree-> pt_1 = -999;
  otree-> eta_1 = -999;
  otree-> phi_1 = -999;
  otree-> iso_1 = -999;
  otree-> q_1 = -999;

  otree-> pt_2 = -999;
  otree-> eta_2 = -999;
  otree-> phi_2 = -999;
  otree-> iso_2 = -999;
  otree-> q_2 = -999; 

  otree-> os = -999;
  otree-> pt_ll = -999;
  otree-> px_ll = -999;
  otree-> py_ll = -999;
  otree-> m_ll = -999;
  otree-> eta_ll = -999;
  otree-> phi_ll = -999;
  otree-> delta_phi_ll = -999;

  otree-> genV_px = -999;
  otree-> genV_py = -999;
  otree-> genV_pt = -999;
  otree-> genV_mass = -999;
  otree-> genV_eta = -999;
  otree-> genV_phi = -999;
  otree-> genZ_pt = -999;
  otree-> genZ_mass = -999;

  otree-> is_gen_Zee = 0;
  otree-> is_gen_Zmm = 0;
  otree-> is_gen_Ztt = 0;

  otree-> pt_jet_1 = -999;
  otree-> eta_jet_1 = -999;
  otree-> phi_jet_1 = -999;

  otree-> pt_jet_2 = -999;
  otree-> eta_jet_2 = -999;
  otree-> phi_jet_2 = -999;

  otree-> m_jj= -999;
  otree-> delta_phi_jj = -999;
  otree-> njets20 = -999;
  otree-> njets30 = -999;
  otree-> HT20 = -999;
  otree-> HT30 = -999;

  otree-> met = -999;
  otree-> met_phi = -999;
  otree-> puppimet = -999;
  otree-> puppimet_phi = -999;
  otree-> met_rcorr = -999;
  otree-> met_phi_rcorr = -999;
  otree-> puppimet_rcorr = -999;
  otree-> puppimet_phi_rcorr = -999;
  otree-> pt_ratio_ll = -999;
  otree-> met_response_gen = -999;
  otree-> puppimet_response_gen = -999;
  otree-> met_recoil_paral = -999;
  otree-> met_recoil_perp = -999;
  otree-> met_response_had = -999;

  otree-> puppimet_recoil_paral = -999;
  otree-> puppimet_recoil_perp = -999;
  otree-> puppimet_response_had = -999;

  otree-> mc_weight = 1;
  otree-> pu_weight = 1; 
  otree-> toppt_weight = 1;
  otree-> zpt_weight = 1;
  otree-> track_weight_1 = 1;
  otree-> track_weight_2 = 1;
  otree-> idiso_weight_1 = 1;
  otree-> idiso_weight_2 = 1;
  otree-> trig_weight = 1;

}


