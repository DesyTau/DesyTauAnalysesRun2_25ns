#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>
#include "HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h"
#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
//#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys_WIP.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h"
//#include "HTT-utilities/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"
#include "Math/GenVector/Boost.h"
#define pi 3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass  0.000511
#define muonMass    0.105658
#define tauMass    1.77682
#define pionMass 0.1396
#define expectedtauspinnerweights 5

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;

struct DiTauInfo 
{ 
  DiTauInfo(){}; 
  int diTauCharge_; 
  double sumPt_; 
  double sumIso_;
  int index1_;
  int index2_;
}; 

struct SortDiTauPairs 
{ 
  bool operator() (const DiTauInfo t1, const DiTauInfo t2) 
  { 
   
    if(t1.sumIso_ > t2.sumIso_ ) return true;
    if(t1.sumIso_ < t2.sumIso_ ) return false;
   
    return (t1.sumPt_ > t2.sumPt_);  
  } 
};
float muon_id_sf(float eta){
  float sf;
  float Aeta = fabs(eta);
  if(Aeta>0&&Aeta<=0.4)
    sf=1.17;
  else if(Aeta>0.4&&Aeta<=0.8)
    sf=1.29;
  else if(Aeta>0.8&&Aeta<=1.2)
    sf=1.14;
  else if(Aeta>1.2&&Aeta<=1.7)
    sf=0.93;  
  else if(Aeta>1.7&&Aeta<=2.3)
    sf=1.61;
  else
    sf=1;
  return sf;
} 
float elec_id_sf(float eta){
  float sf;
  float Aeta = fabs(eta);
  if(Aeta<1.460)
    sf=1.4;
  if(Aeta>1.558)
    sf=1.21;
  return sf;
}


  

int main(int argc, char * argv[]) {

  cout<<endl<<"ok1"<<endl;

  using namespace std;


  string cmsswBase = (getenv ("CMSSW_BASE"));
  // **** configuration
  Config cfg(argv[1]);
  // configuration process
  const string sample = argv[2];
  const int era = cfg.get<int>("era");
  const bool isData = cfg.get<bool>("isData");
  const bool isQCD = cfg.get<bool>("isQCD");
  const bool isDY = cfg.get<bool>("IsDY");
  const bool isWJets = cfg.get<bool>("IsW");
  const bool isEmbedded = cfg.get<bool>("isEmbedded");
  const string infiles = argv[2];
  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyLeptonSFfromKIT  = cfg.get<bool>("ApplyLeptonSFfromKIT"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters = cfg.get<bool>("ApplyMetFilters");
  const bool applyTauSpinnerWeights= cfg.get<bool>("applyTauSpinnerWeights");
  const bool usePuppiMET      = cfg.get<bool>("UsePuppiMET");  
  const bool ApplyIpCorrection = cfg.get<bool>("ApplyIpCorrection");
  cout<<endl<<"ok2"<<endl;

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");
  const string ipCorrFileName = cfg.get<string>("IpCorrFileName");
  TString IpCorrFileName(ipCorrFileName);    
  
  IpCorrection *ip = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileName);

  // tau trigger efficiency
  TauTriggerSFs2017 * tauTriggerSF = new TauTriggerSFs2017(cmsswBase+"/src/TauAnalysisTools/TauTriggerSFs/data/tauTriggerEfficiencies2017.root","ditau","2017","tight","MVAv2");
  //svfit
  
  cout<<endl<<"ok3"<<endl;

  const string svFitPtResFile = TString(TString(cmsswBase)+"/src/"+TString(cfg.get<string>("svFitPtResFile"))).Data();
  TFile* inputFile_visPtResolution = new TFile(svFitPtResFile.data());
  //zptweight file 
  const string ZptweightFile = cfg.get<string>("ZptweightFile");
  //b-tag scale factors
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  TString pathToBtagScaleFactors = (TString) cmsswBase+"/src/"+BtagSfFile;

  BTagCalibration calib("csvv2", (string) pathToBtagScaleFactors );
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central");
  if( ApplyBTagScaling ){
    reader_B.load(calib,BTagEntry::FLAV_B,"comb");
    reader_C.load(calib,BTagEntry::FLAV_C,"comb");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");
  }
  const string TaggingEfficienciesFile = cfg.get<string>("BtagMCeffFile");
  TString pathToTaggingEfficiencies = (TString) cmsswBase+"/src/"+TaggingEfficienciesFile;
  TFile *fileTagging = new TFile( pathToTaggingEfficiencies );
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  
  cout<<endl<<"ok4"<<endl;

  if( ApplyBTagScaling ){
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  }
  TRandom3 *rand = new TRandom3();
  const struct btag_scaling_inputs inputs_btag_scaling_medium = { reader_B, reader_C, reader_Light, tagEff_B, tagEff_C, tagEff_Light, rand };

  const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  

  RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
  if(!isData && applyRecoilCorrections && (isDY || isWJets)){
    TString RecoilDir("HTT-utilities/RecoilCorrections/data/");
    
    TString RecoilFileName = RecoilDir; RecoilFileName += "Type1_PFMET_2017.root";
    std::cout<<RecoilFileName<<std::endl;
    recoilPFMetCorrector = new RecoilCorrector( RecoilFileName);
  }
  // Read in HLT filter
  vector<string> filterDiTau;
  filterDiTau.push_back(cfg.get<string>("filterDiTau1"));
  filterDiTau.push_back(cfg.get<string>("filterDiTau2"));
  filterDiTau.push_back(cfg.get<string>("filterDiTau3"));

   cout<<endl<<"ok5"<<endl;

  // tau cuts
  const float ptTauCut    = cfg.get<float>("ptTauCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut       = cfg.get<float>("dzTauCut");
  const bool applyTauId = cfg.get<bool>("ApplyTauId");
  // for tau Id efficiency
  const float tau_id_sf = cfg.get<float>("TauIdSF");

  // pair selection
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  // extra electron veto
  const float ptVetoElectronCut  = cfg.get<float>("ptVetoElectronCut");  
  const float etaVetoElectronCut = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut = cfg.get<float>("dxyVetoElectronCut");  
  const float dzVetoElectronCut  = cfg.get<float>("dzVetoElectronCut"); 
  const bool applyVetoElectronId = cfg.get<bool>("applyVetoElectronId");
  const float isoVetoElectronCut = cfg.get<float>("isoVetoElectronCut");  
  
  // extra muon veto
  const float ptVetoMuonCut  = cfg.get<float>("ptVetoMuonCut");  
  const float etaVetoMuonCut = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut = cfg.get<float>("dxyVetoMuonCut");  
  const float dzVetoMuonCut  = cfg.get<float>("dzVetoMuonCut"); 
  const bool  applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
  const float isoVetoMuonCut = cfg.get<float>("isoVetoMuonCut");

  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
  const float dRiso = cfg.get<float>("dRiso");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");
  
  // tau energy scale corrections
  const float shift_tes_1prong = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  const float shift_tes_3p1p0 = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0");

  const float shift_tes_1prong_efake = cfg.get<float>("TauEnergyScaleShift_OneProng_efake");
  const float shift_tes_1p1p0_efake  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_efake");

  const float shift_tes_1prong_mufake = cfg.get<float>("TauEnergyScaleShift_OneProng_mufake");
  const float shift_tes_1p1p0_mufake  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_mufake");

      // Workspace with corrections
      
      const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");

      TString workspace_filename = TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/" + CorrectionWorkspaceFileName;
      cout << "Taking correction workspace from " << workspace_filename << endl;
      TFile *f_workspace = new TFile(workspace_filename, "read");
      if (f_workspace->IsZombie()) {
	std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
	exit(-1);
      }
      
      cout<<endl<<"before "<<endl;
      RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");
     
      f_workspace->Close();

   cout<<endl<<"ok6"<<endl;

  // **** end of configuration analysis
  // create input files list
  int ifile = 0;
  int jfile = -1;

  std::vector<std::string> fileList; 
  int NumberOfFiles = 0;
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
	if (infile.length() > 0){
	  cout<<endl<<"push"<<endl;
	  fileList.push_back(infile);
	  NumberOfFiles +=1 ;
	}
      }
      else
	break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

   cout<<endl<<"ok7"<<endl;

  if(NumberOfFiles < jfile) jfile = NumberOfFiles;
  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }
  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");

  std::string ntupleNameinitroottree("initroottree/AC1B"); //Since update of big tupler, need to fetch all genweights from here normalisation

  //Fix the name for the tauspinner tree
  std::string TauSpinnerWeightTreeName("icTauSpinnerProducer/TauSpinnerWeightTree");

     cout<<endl<<"ok8"<<endl;

  // PU reweighting - initialization
  PileUp * PUofficial = new PileUp();
  if(ApplyPUweight && isData==false){
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(pileUpInDataFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(pileUpInMCFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    //std::cout << filePUdistribution_data << std::endl;
    //std::cout << filePUdistribution_MC << std::endl;
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(TString(pileUpforMC));
    if (PU_mc==NULL) {
      std::cout << "Histogram " << pileUpforMC << " is not present in pileup file" << std::endl;
      exit(-1);
    }
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }

       cout<<endl<<"ok9"<<endl;

  unsigned int lhc_run_era = 2;
  bool applyRun1topPtWeights = true;
  if (applyRun1topPtWeights) lhc_run_era =1;
  // Zpt reweighting for LO DY samples 
  TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+ZptweightFile,"read");
  //TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+"DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016_BtoH.root","read");
  //Merijn: the file will point now for 2017  instead to DesyTauAnalyses/NTupleMaker/data/zpt_weights_2017.root
  TH2D * h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");
  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_tt_Sync.root";
  
     cout<<endl<<"ok10"<<endl;

  std::cout <<rootFileName <<std::endl;  

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5,0.5);
  TH1D * nWeightedEventsHMiniAOD = new TH1D("nWeightedEventsHMiniAOD", "nWeightedEventsHMiniAOD", 1, -0.5,0.5);  

  TTree * tree = new TTree("TauCheck","TauCheck");
  TTree * gtree = new TTree("GenTauCheck","GenTauCheck");

  Synch17Tree *otree = new Synch17Tree(tree);
  Synch17GenTree *gentree = new Synch17GenTree(gtree);
  Synch17GenTree *gentreeForGoodRecoEvtsOnly = new Synch17GenTree(tree);
  //initializeCPvar(otree);
  int nTotalFiles = 0;

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  std::vector<TString> met_filters_list ;
  met_filters_list.push_back("Flag_HBHENoiseFilter");
  met_filters_list.push_back("Flag_HBHENoiseIsoFilter");
  met_filters_list.push_back("Flag_globalSuperTightHalo2016Filter");
  met_filters_list.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  met_filters_list.push_back("Flag_goodVertices");
  if (isData)
    met_filters_list.push_back("Flag_eeBadScFilter");
  met_filters_list.push_back("Flag_BadPFMuonFilter");
  //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
  met_filters_list.push_back("ecalBadCalibReducedMINIAODFilter"); //WAS NOT IN LIST ABOVE
  
  ofstream check;
  check.open("check.txt");
  check<<"checking...";

  int counter[20];

     cout<<endl<<"ok11"<<endl;
      check<<endl<<jfile<<endl;


  ///////////////FILE LOOP///////////////
  for (int iF=ifile; iF<1; ++iF) {
    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
     cout<<endl<<"ok12"<<endl;

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;



    TTree * _treenorm = NULL;
    _treenorm = (TTree*)file_->Get(TString(ntupleNameinitroottree));        
    if (_treenorm==NULL) continue;
    AC1B analysisTreenorm(_treenorm, isData);

    Long64_t numberOfEntriesNorm = analysisTreenorm.GetEntries();
    for (Long64_t iEntry=0; iEntry<numberOfEntriesNorm; iEntry++) {
      analysisTreenorm.GetEntry(iEntry);
      nWeightedEventsHMiniAOD->Fill(0., analysisTreenorm.genweight);
    }

   

    double * TSweight=new double[expectedtauspinnerweights];
    TTree * _treeTauSpinnerWeights = NULL;

    if(applyTauSpinnerWeights){ 
      _treeTauSpinnerWeights = (TTree*)file_->Get(TString(TauSpinnerWeightTreeName));
      _treeTauSpinnerWeights->SetBranchAddress("TauSpinnerWeights",TSweight);
    }  

    
    TH1D * histoInputEvents = NULL;

    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    Long64_t numberOfEntries = analysisTree.GetEntries();



    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;

 
    ///////////////EVENT LOOP///////////////
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {
      // cout<<"iEntry "<<iEntry<<endl;
      
      counter[0]++;
      analysisTree.GetEntry(iEntry);
      nEvents++;
   

      if(!applyTauSpinnerWeights){
	for(int tsindex = 0; tsindex < expectedtauspinnerweights; tsindex++) 
          TSweight[tsindex] = 1;

	otree->TauSpinnerWeightsEven = TSweight[0];
	gentree->sm_htt125 = TSweight[0];
	gentreeForGoodRecoEvtsOnly->sm_htt125 = TSweight[0];

	otree->TauSpinnerWeightsMaxMix = TSweight[1];
	gentree->mm_htt125 = TSweight[1];
	gentreeForGoodRecoEvtsOnly->mm_htt125 = TSweight[1];

	otree->TauSpinnerWeightsOdd = TSweight[2];
	gentree->ps_htt125 = TSweight[2];
	gentreeForGoodRecoEvtsOnly->ps_htt125 = TSweight[2];

	otree->TauSpinnerWeightsMinusMaxMix = TSweight[3];
	gentree->minusmm_htt125 = TSweight[3];
	gentreeForGoodRecoEvtsOnly->minusmm_htt125 = TSweight[3];

	otree->TauSpinnerWeightsMix0p375 = TSweight[4];
	gentree->mix0p375_htt125 = TSweight[4];
	gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = TSweight[4];
      
      }
      else if(era==2017){
        _treeTauSpinnerWeights->GetEntry(iEntry);

	otree->TauSpinnerWeightsEven = TSweight[0];
	gentree->sm_htt125 = TSweight[0];
	gentreeForGoodRecoEvtsOnly->sm_htt125 = TSweight[0];

	otree->TauSpinnerWeightsMaxMix = TSweight[1];
	gentree->mm_htt125 = TSweight[1];
	gentreeForGoodRecoEvtsOnly->mm_htt125 = TSweight[1];

	otree->TauSpinnerWeightsOdd = TSweight[2];
	gentree->ps_htt125 = TSweight[2];
	gentreeForGoodRecoEvtsOnly->ps_htt125 = TSweight[2];

	otree->TauSpinnerWeightsMinusMaxMix = TSweight[3];
	gentree->minusmm_htt125 = TSweight[3];
	gentreeForGoodRecoEvtsOnly->minusmm_htt125 = TSweight[3];

	otree->TauSpinnerWeightsMix0p375 = TSweight[4];
	gentree->mix0p375_htt125 = TSweight[4];
	gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = TSweight[4];
      }
      else{
	otree->TauSpinnerWeightsEven = analysisTree.TauSpinnerWeight[0];
	gentree->sm_htt125 = analysisTree.TauSpinnerWeight[0];
	gentreeForGoodRecoEvtsOnly->sm_htt125 = analysisTree.TauSpinnerWeight[0];

	otree->TauSpinnerWeightsMaxMix = analysisTree.TauSpinnerWeight[1];
	gentree->mm_htt125 = analysisTree.TauSpinnerWeight[1];
	gentreeForGoodRecoEvtsOnly->mm_htt125 = analysisTree.TauSpinnerWeight[1];

	otree->TauSpinnerWeightsOdd = analysisTree.TauSpinnerWeight[2];
	gentree->ps_htt125 = analysisTree.TauSpinnerWeight[2];
	gentreeForGoodRecoEvtsOnly->ps_htt125 = analysisTree.TauSpinnerWeight[2];

	if(analysisTree.TauSpinAngles_count>=5){
	  otree->TauSpinnerWeightsMinusMaxMix = analysisTree.TauSpinnerWeight[3];
	  gentree->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
	  gentreeForGoodRecoEvtsOnly->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
	    
	  otree->TauSpinnerWeightsMix0p375 = analysisTree.TauSpinnerWeight[4];
	  gentree->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
	  gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
	}
      }


      if (isData)
     	nWeightedEventsH->Fill(0., 1.);
      else {
     	nWeightedEventsH->Fill(0., analysisTree.genweight);
	
      }

    

      //Skip events not passing the MET filters, if applied

      
      
      //      if (ApplyMetFilters && !passedAllMetFilters(&analysisTree, met_filters_list, isData)) continue;

     

      counter[1]++;
      vector<int> nDiTauTrig(filterDiTau.size(),-1);

      vector<bool> checkFilterDiTauTrig(filterDiTau.size(), false);
      /*
      for(unsigned int i=0; i<analysisTree.run_hltfilters->size(); i++){
	cout<<"analysisTree.run_hltfilters->at(i)="<<analysisTree.run_hltfilters->at(i)<<endl;
      }
      */
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
     	TString HLTFilter(analysisTree.run_hltfilters->at(i));
     	for(unsigned int i_trig=0; i_trig<filterDiTau.size(); i_trig++){
     	  if (HLTFilter==filterDiTau.at(i_trig)){ nDiTauTrig.at(i_trig) = i; checkFilterDiTauTrig.at(i_trig) = true;}
     	}
      }
      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt = analysisTree.event_nr;
    
      check<<endl<<"ok1"<<endl;

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
     	continue;

      
      check<<endl<<"ok2"<<endl;   

            if(ApplyPUweight) fill_weight(&analysisTree, otree, PUofficial, isData);
      
          check<<endl<<"ok3"<<endl;

      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;

      // tau selection
      std::vector<DiTauInfo> sortDiTauInfos;  sortDiTauInfos.clear();
      vector<int> taus; taus.clear();

      check<<endl<<"ok4"<<endl;
      /*
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) { 
	cout<<endl<<"tau_pt[it]="<<analysisTree.tau_pt[it]<<endl<<"tau_eta[it]="<<analysisTree.tau_eta[it]<<endl<<"tau_phi[it]="<<analysisTree.tau_phi[it]<<endl;
      }

      for (unsigned int it = 0; it<analysisTree.l1tau_count; ++it) { 
	cout<<endl<<"l1tau_pt[it]="<<analysisTree.l1tau_pt[it]<<endl<<"l1tau_eta[it]="<<analysisTree.l1tau_eta[it]<<endl<<"l1tau_phi[it]="<<analysisTree.l1tau_phi[it]<<endl;
      }
      */


      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) { 
     	check<<endl<<"ok4.1"<<endl;

	float ptTauCut_corr = ptTauCut;

	if (analysisTree.tau_decayMode[it]==0 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/shift_tes_1prong;
	if (analysisTree.tau_decayMode[it]==1 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/shift_tes_1p1p0;
	if (analysisTree.tau_decayMode[it]==10 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/shift_tes_3prong;
	if (analysisTree.tau_decayMode[it]==11 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/shift_tes_3p1p0;
	if (analysisTree.tau_decayMode[it]==0 && (analysisTree.tau_genmatch[it]==1 || analysisTree.tau_genmatch[it]==3)) ptTauCut_corr = ptTauCut/shift_tes_1prong_efake;
	if (analysisTree.tau_decayMode[it]==1 && (analysisTree.tau_genmatch[it]==1 || analysisTree.tau_genmatch[it]==3)) ptTauCut_corr = ptTauCut/shift_tes_1p1p0_efake;
	if (analysisTree.tau_decayMode[it]==0 && (analysisTree.tau_genmatch[it]==2 || analysisTree.tau_genmatch[it]==4)) ptTauCut_corr = ptTauCut/shift_tes_1prong_mufake;
	if (analysisTree.tau_decayMode[it]==1 && (analysisTree.tau_genmatch[it]==2 || analysisTree.tau_genmatch[it]==4)) ptTauCut_corr = ptTauCut/shift_tes_1p1p0_mufake;

      
	//cout<<endl<<"analysisTree.tau_decayMode[it]="<<analysisTree.tau_decayMode[it]<<endl<<"ptTauCut_corr="<<ptTauCut_corr<<endl;

        if (analysisTree.tau_pt[it]<=ptTauCut_corr) continue;
      check<<endl<<"ok4.2"<<endl;
        if (fabs(analysisTree.tau_eta[it])>=etaTauCut) continue;
      check<<endl<<"ok4.3"<<endl;
        if (fabs(fabs(analysisTree.tau_charge[it])-1)>0.001) continue;
      check<<endl<<"ok4.4"<<endl;
     	if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>=dzTauCut) continue;
      check<<endl<<"ok4.5"<<endl;
        //Added AN
        if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] < 0.5) continue;//VVVLoose tau mva against jet
      check<<endl<<"ok4.6"<<endl;
     	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[it] < 0.5) continue;//VVVLoose mva aginst e
      check<<endl<<"ok4.7"<<endl;
     	if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] < 0.5) continue;//very loose mva agaist muon
      check<<endl<<"ok4.8"<<endl;
        if (applyTauId && (analysisTree.tau_decayModeFindingNewDMs[it] < 0.5 || analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6)) continue;
      check<<endl<<"ok4.9"<<endl;
        //-----------------------------------
        taus.push_back(it);

      }


      check<<endl<<"ok5"<<endl;
      if(taus.size()<2)continue;
       check<<endl<<"ok6"<<endl;

      

     //loop over pair of taus
      int tauIndex_1 = -1;
      int tauIndex_2 = -1;
      for (unsigned int it=0; it<taus.size(); ++it) {
        for (unsigned int itn=it+1; itn<taus.size(); ++itn) {
     	  unsigned int tIndex1 = taus.at(it);
          unsigned int tIndex2 = taus.at(itn);
     	  float dR = deltaR(analysisTree.tau_eta[tIndex1],analysisTree.tau_phi[tIndex1],analysisTree.tau_eta[tIndex2],analysisTree.tau_phi[tIndex2]);
	  

      check<<endl<<"ok7"<<endl;
          if (dR<dRleptonsCut) continue;
      check<<endl<<"ok8"<<endl;


     	  //Selecting isolated minimum tau pair                                                                                                                                                                 
     	  float sumPt = analysisTree.tau_pt[tIndex1]+analysisTree.tau_pt[tIndex2];
     	  //	  float isoTau1 = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex1];
     	  //          float isoTau2 = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex2];
     	  float isoTau1 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex1];  //by AN
          float isoTau2 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex2]; //by AN
     	  float sumIso = isoTau1+isoTau2;
          DiTauInfo sortDiTauInfo;
          sortDiTauInfo.index1_ = tIndex1;
          sortDiTauInfo.index2_ = tIndex2;
          sortDiTauInfo.sumPt_ = sumPt;
          sortDiTauInfo.sumIso_ = sumIso;
          
          sortDiTauInfos.push_back(sortDiTauInfo);
     	}//tau-
      }//tau+

      std::sort(sortDiTauInfos.begin(), sortDiTauInfos.end(), SortDiTauPairs());
      int diTauCounter = -1;

      for(std::vector<DiTauInfo>::iterator iter = sortDiTauInfos.begin(); iter != sortDiTauInfos.end() ; iter++){
     	if(diTauCounter >= 0) continue;

    

     	uint tIndex1 = iter->index1_;
     	uint tIndex2 = iter->index2_;
     	LV temp_Leg1_(analysisTree.tau_px[tIndex1], analysisTree.tau_py[tIndex1], analysisTree.tau_pz[tIndex1], analysisTree.tau_e[tIndex1]);
     	LV temp_Leg2_(analysisTree.tau_px[tIndex2], analysisTree.tau_py[tIndex2], analysisTree.tau_pz[tIndex2], analysisTree.tau_e[tIndex2]);
     	LV Leg1P4_, Leg2P4_;
     	// if(analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex1] > analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex2]){
     	//   Leg1P4_ = temp_Leg1_;
     	//   Leg2P4_ = temp_Leg2_;
     	// }
     	// else if(analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex1] < analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex2]){
     	//   Leg1P4_ = temp_Leg2_;tIndex1 = iter->index2_;
     	//   Leg2P4_ = temp_Leg1_;tIndex2 = iter->index1_;
     	// }
     	//else 
     	if(temp_Leg1_.pt() > temp_Leg2_.pt()){
     	  Leg1P4_ = temp_Leg1_;tIndex1 = iter->index1_;
     	  Leg2P4_ = temp_Leg2_;tIndex2 = iter->index2_;
     	}
     	else {
     	  Leg1P4_ = temp_Leg2_; tIndex1 = iter->index2_;
     	  Leg2P4_ = temp_Leg1_; tIndex2 = iter->index1_;
     	}
     	++diTauCounter;
     	tauIndex_1=tIndex1;
     	tauIndex_2=tIndex2;
      }
    
      if(tauIndex_1==-1 || tauIndex_2==-1) continue; 
  
    

      //make isTrigger to false                                                                                                                                                                             
      bool isDiTauTrig = false;
      bool isTauTrig = false;
      bool isTau2Trig = false;

     

      for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
     	float dRtrigTau1 = deltaR(analysisTree.tau_eta[tauIndex_1], analysisTree.tau_phi[tauIndex_1],
     				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
     	float dRtrigTau2 = deltaR(analysisTree.tau_eta[tauIndex_2], analysisTree.tau_phi[tauIndex_2],
     				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    
	    
     	if(dRtrigTau1<deltaRTrigMatch){
	      
     	  for(unsigned int i_trig = 0; i_trig<filterDiTau.size(); i_trig++){
     	    if (nDiTauTrig.at(i_trig) == -1) continue;
     	    if (analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)]) isTauTrig = true;
     	  }
     	}
     	if(dRtrigTau2<deltaRTrigMatch){
     	  for(unsigned int i_trig = 0; i_trig<filterDiTau.size(); i_trig++){
     	    //cout<<nDiTauTrig.at(i_trig)<<endl;                                                                                                                                                               
     	    if (nDiTauTrig.at(i_trig) == -1) continue;
	    
     	    if (analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)]) isTau2Trig = true;
	    
     	  }
     	}
	
      }//trigger obj match
      /*
      for (unsigned int it = 0; it<analysisTree.l1tau_count; ++it) { 
	TLorentzVector l1tau_LV(analysisTree.l1tau_px[it],analysisTree.l1tau_py[it],analysisTree.l1tau_pz[it],sqrt(pow(analysisTree.l1tau_px[it],2)+pow(analysisTree.l1tau_py[it],2)+pow(analysisTree.l1tau_pz[it],2)));
	//cout<<endl<<"analysisTree.l1tau_px[it]="<<analysisTree.l1tau_px[it]<<endl<<"analysisTree.l1tau_py[it]="<<analysisTree.l1tau_py[it]<<endl<<"analysisTree.l1tau_pz[it]="<<analysisTree.l1tau_pz[it]<<endl;
	//cout<<endl<<"l1tau_pt[it]="<<analysisTree.l1tau_pt[it]<<endl<<"l1tau_eta[it]="<<analysisTree.l1tau_eta[it]<<endl<<"l1tau_phi[it]="<<analysisTree.l1tau_phi[it]<<endl;
	cout<<endl<<"l1tau_LV.Pt()="<<l1tau_LV.Pt()<<endl<<"l1tau_LV.Eta()="<<l1tau_LV.Eta()<<endl<<"l1tau_LV.Phi()="<<l1tau_LV.Phi()<<endl;
      }

      cout<<endl<<"analysisTree.tau_pt[tauIndex_1]="<<analysisTree.tau_pt[tauIndex_1]<<endl<<"analysisTree.tau_eta[tauIndex_1]="<<analysisTree.tau_eta[tauIndex_1]<<endl<<"analysisTree.tau_phi[tauIndex_1]="<<analysisTree.tau_phi[tauIndex_1]<<endl;

	cout<<endl<<"analysisTree.tau_pt[tauIndex_2]="<<analysisTree.tau_pt[tauIndex_2]<<endl<<"analysisTree.tau_eta[tauIndex_2]="<<analysisTree.tau_eta[tauIndex_2]<<endl<<"analysisTree.tau_phi[tauIndex_2]="<<analysisTree.tau_phi[tauIndex_2]<<endl;
      */
 
      bool is_l1tau_match_1=false, is_l1tau_match_2=false;

      for (unsigned int iT=0; iT<analysisTree.l1tau_count; ++iT) {
	
	TLorentzVector l1tau_LV(analysisTree.l1tau_px[iT],analysisTree.l1tau_py[iT],analysisTree.l1tau_pz[iT],sqrt(pow(analysisTree.l1tau_px[iT],2)+pow(analysisTree.l1tau_py[iT],2)+pow(analysisTree.l1tau_pz[iT],2)));

        float dRl1tau1 = deltaR(analysisTree.tau_eta[tauIndex_1], analysisTree.tau_phi[tauIndex_1],
				l1tau_LV.Eta(),l1tau_LV.Phi());
        float dRl1tau2 = deltaR(analysisTree.tau_eta[tauIndex_2], analysisTree.tau_phi[tauIndex_2],
				l1tau_LV.Eta(),l1tau_LV.Phi());

        if(dRl1tau1<deltaRTrigMatch){
	  if(analysisTree.l1tau_pt[iT]>32.0) is_l1tau_match_1 = true;
	}

        if(dRl1tau2<deltaRTrigMatch){
	  if(analysisTree.l1tau_pt[iT]>32.0) is_l1tau_match_2 = true;
	}
	if(is_l1tau_match_1&&is_l1tau_match_2) break;
      }


      bool is_l1tau_match_double = is_l1tau_match_1 && is_l1tau_match_2;
     


      check<<endl<<"ok9"<<endl;                                                                                                          
      isDiTauTrig =isTauTrig && isTau2Trig;
      otree->trg_doubletau = (isDiTauTrig);// && is_l1tau_match_double);
      
      ///////////CALCULATING WEIGHTS/////////////

      
      otree->trigweight_1 = 1;
      otree->trigweight_2 = 1;
      otree->idisoweight_1 = 1;
      otree->idisoweight_2 = 1;

      if ((!isData || isEmbedded)) { 
	if (!isEmbedded){
	  if(analysisTree.tau_genmatch[tauIndex_1] == 5){
	    w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
	    otree->idisoweight_1 = w->function("t_deeptauid_dm_medium")->getVal();
	  }
	  if(analysisTree.tau_genmatch[tauIndex_2] == 5){
	    w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
	    otree->idisoweight_2 = w->function("t_deeptauid_dm_medium")->getVal();
	  }
	}
	else{
	  double t_dm;
	  if(analysisTree.tau_genmatch[tauIndex_1] == 5){
	    t_dm = analysisTree.tau_decayMode[tauIndex_1];
	    otree->idisoweight_1 = 0.99*((t_dm==0)*0.975 + (t_dm==1)*0.975*1.051 + (t_dm==10)*pow(0.975,3) + (t_dm==11)*pow(0.975,3)*1.051);
	  }
	  if(analysisTree.tau_genmatch[tauIndex_2] == 5){
	    t_dm = analysisTree.tau_decayMode[tauIndex_2];
	    otree->idisoweight_2 = 0.99*((t_dm==0)*0.975 + (t_dm==1)*0.975*1.051 + (t_dm==10)*pow(0.975,3) + (t_dm==11)*pow(0.975,3)*1.051);
	  }
	}
      }

      if ((!isData || isEmbedded) && ApplyLepSF) {
	TString suffix = "mc";
	TString suffixRatio = "ratio";
	if (isEmbedded) {suffix = "embed"; suffixRatio = "embed_ratio";}
	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex_1]);
	w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex_1]);
	w->var("t_phi")->setVal(analysisTree.tau_phi[tauIndex_1]);
	w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
	otree->trigweight_1 = w->function("t_trg_mediumDeepTau_ditau_"+suffixRatio)->getVal();
	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex_2]);
	w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex_2]);
	w->var("t_phi")->setVal(analysisTree.tau_phi[tauIndex_2]);
	w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
	otree->trigweight_2 = w->function("t_trg_mediumDeepTau_ditau_"+suffixRatio)->getVal();
      }

      otree->mcweight = analysisTree.genweight;
      otree->effweight = otree->idisoweight_1 * otree->idisoweight_2 * otree->trigweight_1 * otree->trigweight_2;
      otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
      


     

      /////////////////////////////////////////////
      

      //      otree->is_l1tau_match_double = is_l1tau_match_double;
      int tauIndex=-1;
      //      if(!isDiTauTrig) continue;
      //extra letpn veto
           check<<endl<<"ok10"<<endl;


      otree->extraelec_veto = extra_electron_veto(tauIndex, "tt", &cfg, &analysisTree);
      otree->extramuon_veto = extra_muon_veto(tauIndex, "tt", &cfg, &analysisTree, isData);
      
                check<<endl<<"ok11"<<endl;

      /*
      if(!isData){
     	otree->trigweight_1=tauTriggerSF->getTriggerScaleFactor(analysisTree.tau_pt[tauIndex_1],analysisTree.tau_eta[tauIndex_1],analysisTree.tau_phi[tauIndex_1],analysisTree.tau_decayMode[tauIndex_1]);
     	otree->trigweight_2=tauTriggerSF->getTriggerScaleFactor(analysisTree.tau_pt[tauIndex_2],analysisTree.tau_eta[tauIndex_2],analysisTree.tau_phi[tauIndex_2],analysisTree.tau_decayMode[tauIndex_2]);
	
     	otree->trigweight =  otree->trigweight_1* otree->trigweight_2;
     	otree->idisoweight_1 = tau_id_sf;
     	otree->idisoweight_2 = tau_id_sf;
     	
     	otree->weight = otree->idisoweight_1 * otree->idisoweight_2 * otree->trigweight * otree->puweight * otree->mcweight;
      }
      */        

      //MET
      //Merijn 2019 6 20: overloaded the function, it takes the era as arugment now, to take pfmetcorr for 2016 and 2017..
      //fillMET("tt", tauIndex_1, tauIndex_2, &analysisTree, otree);

          check<<endl<<"ok12"<<endl;
     
      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);

      // Zpt weight
      otree->zptweight = 1.;
      if (!isData && (isDY  || isWJets)){
     	genV = genTools::genV(analysisTree); // gen Z boson ?
     	float bosonMass = genV.M();
     	float bosonPt = genV.Pt();

     	//Merijn determine here some min and max values:
     	double massxmin=h_zptweight->GetXaxis()->GetXmin();//Merijn
     	double massxmax=h_zptweight->GetXaxis()->GetXmax();

     	double ptxmin=h_zptweight->GetYaxis()->GetXmin();
     	double ptxmax=h_zptweight->GetYaxis()->GetXmax();

          check<<endl<<"ok13"<<endl;

     	//Merijn 2019 6 13: adjust to T/M functions, to get boundaries right. Otherwise, for 2017 data we get few outliers that screw up the weight histogram dramatically.
     	Float_t zptmassweight = 1;
     	if (bosonMass>50.0) {
     	  float bosonMassX = bosonMass;
     	  float bosonPtX = bosonPt;
     	  if (bosonMassX>massxmax) bosonMassX = massxmax-h_zptweight->GetXaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;//Merijn: if doesn't work, lower by 1/2 binwidth..
     	  if (bosonPtX<ptxmin)      bosonPtX = ptxmin+h_zptweight->GetYaxis()->GetBinWidth(1)*0.5;
     	  if (bosonPtX>ptxmax)   bosonPtX = ptxmax-h_zptweight->GetYaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;
	  
     	  zptmassweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(bosonMassX),
     						     h_zptweight->GetYaxis()->FindBin(bosonPtX));
     	}
     	otree->zptweight =zptmassweight;
      }

          check<<endl<<"ok14"<<endl;	
	
      // topPt weight
      otree->topptweight =1.;
      if(!isData){
     	otree->topptweight = genTools::topPtWeight(analysisTree, lhc_run_era);
      }

          check<<endl<<"ok15"<<endl;
	
      TLorentzVector tauLV_1; tauLV_1.SetXYZT(analysisTree.tau_px[tauIndex_1],
     					      analysisTree.tau_py[tauIndex_1],
     					      analysisTree.tau_pz[tauIndex_1],
     					      analysisTree.tau_e[tauIndex_1]);
      TLorentzVector tauLV_2; tauLV_2.SetXYZT(analysisTree.tau_px[tauIndex_2],
     					      analysisTree.tau_py[tauIndex_2],
     					      analysisTree.tau_pz[tauIndex_2],
     					      analysisTree.tau_e[tauIndex_2]);
      TLorentzVector diTauLV = tauLV_1 + tauLV_2;
      // visible mass
      otree->m_vis = diTauLV.M();


      /////for sync
      // bisector of lepton and tau transverse momenta
    
      float tauUnitX_1 = tauLV_1.Px() / tauLV_1.Pt();
      float tauUnitY_1 = tauLV_1.Py() / tauLV_1.Pt();
    
      float tauUnitX_2 = tauLV_2.Px() / tauLV_2.Pt();
      float tauUnitY_2 = tauLV_2.Py() / tauLV_2.Pt();
    
      float zetaX = tauUnitX_1 + tauUnitX_2;
      float zetaY = tauUnitY_1 + tauUnitY_2;
    
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
    
      zetaX = zetaX / normZeta;
      zetaY = zetaY / normZeta;
    
      float vectorVisX = tauLV_1.Px() + tauLV_2.Px();
      float vectorVisY = tauLV_1.Py() + tauLV_2.Py();
    
      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;


      ////////

      IpCorrection * ipCorrector = NULL;
      if (ApplyIpCorrection && (!isData)) ipCorrector = ip; 
      acott_Impr(&analysisTree, otree, tauIndex_1, tauIndex_2, "tt",ipCorrector);
      //acott_Impr(&analysisTree,otree,tauIndex_1,tauIndex_2,"tt");
      // phiCP_pipi(&analysisTree, otree, tauIndex_1,tauIndex_2); 
      //const AC1B * analysisTree;      
  
      //tau1 
      otree->pt_1 = analysisTree.tau_pt[tauIndex_1];
      otree->eta_1 = analysisTree.tau_eta[tauIndex_1];
      otree->phi_1 = analysisTree.tau_phi[tauIndex_1];
      otree->q_1 = analysisTree.tau_charge[tauIndex_1];
      otree->gen_match_1 = analysisTree.tau_genmatch[tauIndex_1];
      otree->mva_1 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_1];
      otree->mva17_1= analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->d0_1 = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex_1];
      otree->dZ_1 = analysisTree.tau_leadchargedhadrcand_dz[tauIndex_1];      
      otree->iso_1 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tauIndex_1];
      otree->m_1 = analysisTree.tau_mass[tauIndex_1];
      //CP measurement
      otree->dm_1 = analysisTree.tau_decayMode[tauIndex_1];
      otree->dmMVA_1 = analysisTree.tau_MVADM2017v1[tauIndex_1];

      otree->tau_pca2D_x_1 = analysisTree.tau_pca2D_x[tauIndex_1];
      otree->tau_pca2D_y_1 = analysisTree.tau_pca2D_y[tauIndex_1];
      otree->tau_pca2D_z_1 = analysisTree.tau_pca2D_z[tauIndex_1];
      otree->tau_pca3D_x_1 = analysisTree.tau_pca3D_x[tauIndex_1];
      otree->tau_pca3D_y_1 = analysisTree.tau_pca3D_y[tauIndex_1];
      otree->tau_pca3D_z_1 = analysisTree.tau_pca3D_z[tauIndex_1];
      otree->tau_SV_x_1 = analysisTree.tau_SV_x[tauIndex_1];
      otree->tau_SV_y_1 = analysisTree.tau_SV_y[tauIndex_1];
      otree->tau_SV_z_1 = analysisTree.tau_SV_z[tauIndex_1];
      otree->tau_SV_covxx_1 = analysisTree.tau_SV_cov[tauIndex_1][0];
      otree->tau_SV_covyx_1 = analysisTree.tau_SV_cov[tauIndex_1][1];
      otree->tau_SV_covzx_1 = analysisTree.tau_SV_cov[tauIndex_1][2];
      otree->tau_SV_covyy_1 = analysisTree.tau_SV_cov[tauIndex_1][3];
      otree->tau_SV_covzy_1 = analysisTree.tau_SV_cov[tauIndex_1][4];
      otree->tau_SV_covzz_1 = analysisTree.tau_SV_cov[tauIndex_1][5];

      otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_1];
      otree-> byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = analysisTree.tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex_1];
      otree->againstMuonLoose3_1 = analysisTree.tau_againstMuonLoose3[tauIndex_1];
      otree->againstMuonTight3_1 = analysisTree.tau_againstMuonTight3[tauIndex_1];
      otree->againstElectronVLooseMVA6_1 = analysisTree.tau_againstElectronVLooseMVA6[tauIndex_1];
      otree->againstElectronTightMVA6_1 = analysisTree.tau_againstElectronTightMVA6[tauIndex_1];
      ////for sync
      otree->againstElectronLooseMVA6_1 = analysisTree.tau_againstElectronLooseMVA6[tauIndex_1];
      otree->againstElectronMediumMVA6_1 = analysisTree.tau_againstElectronMediumMVA6[tauIndex_1];
      otree->againstElectronVTightMVA6_1 = analysisTree.tau_againstElectronVTightMVA6[tauIndex_1];
      ////

     ///AN///
      otree->deepTauVsJetRaw_1 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tauIndex_1];
      otree->deepTauVsEleRaw_1 = analysisTree.tau_byDeepTau2017v2p1VSeraw[tauIndex_1];
      otree->deepTauVsMuRaw_1 = analysisTree.tau_byDeepTau2017v2p1VSmuraw[tauIndex_1];
      //////


      //tau2 
      otree->pt_2 = analysisTree.tau_pt[tauIndex_2];
      otree->eta_2 = analysisTree.tau_eta[tauIndex_2];
      otree->phi_2 = analysisTree.tau_phi[tauIndex_2];
      otree->q_2 = analysisTree.tau_charge[tauIndex_2];
      otree->gen_match_2 = analysisTree.tau_genmatch[tauIndex_2];
      otree->mva_2 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex_2];
      otree->mva17_2= analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      
      otree->d0_2 = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex_2];
      otree->dZ_2 = analysisTree.tau_leadchargedhadrcand_dz[tauIndex_2];      
      otree->iso_2 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tauIndex_2];
      otree->m_2 = analysisTree.tau_mass[tauIndex_2];
      //CP measurement

      

      otree->dm_2 = analysisTree.tau_decayMode[tauIndex_2];
      otree->dmMVA_2 = analysisTree.tau_MVADM2017v1[tauIndex_2];

      otree->tau_pca2D_x_2 = analysisTree.tau_pca2D_x[tauIndex_2];
      otree->tau_pca2D_y_2 = analysisTree.tau_pca2D_y[tauIndex_2];
      otree->tau_pca2D_z_2 = analysisTree.tau_pca2D_z[tauIndex_2];
      otree->tau_pca3D_x_2 = analysisTree.tau_pca3D_x[tauIndex_2];
      otree->tau_pca3D_y_2 = analysisTree.tau_pca3D_y[tauIndex_2];
      otree->tau_pca3D_z_2 = analysisTree.tau_pca3D_z[tauIndex_2];
      otree->tau_SV_x_2 = analysisTree.tau_SV_x[tauIndex_2];
      otree->tau_SV_y_2 = analysisTree.tau_SV_y[tauIndex_2];
      otree->tau_SV_z_2 = analysisTree.tau_SV_z[tauIndex_2];
      otree->tau_SV_covxx_2 = analysisTree.tau_SV_cov[tauIndex_2][0];
      otree->tau_SV_covyx_2 = analysisTree.tau_SV_cov[tauIndex_2][1];
      otree->tau_SV_covzx_2 = analysisTree.tau_SV_cov[tauIndex_2][2];
      otree->tau_SV_covyy_2 = analysisTree.tau_SV_cov[tauIndex_2][3];
      otree->tau_SV_covzy_2 = analysisTree.tau_SV_cov[tauIndex_2][4];
      otree->tau_SV_covzz_2 = analysisTree.tau_SV_cov[tauIndex_2][5];

      otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex_2];
      otree-> byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = analysisTree.tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex_2];
      otree->againstMuonLoose3_2 = analysisTree.tau_againstMuonLoose3[tauIndex_2];
      otree->againstMuonTight3_2 = analysisTree.tau_againstMuonTight3[tauIndex_2];
      otree->againstElectronVLooseMVA6_2 = analysisTree.tau_againstElectronVLooseMVA6[tauIndex_2];
      otree->againstElectronTightMVA6_2 = analysisTree.tau_againstElectronTightMVA6[tauIndex_2];
      ////for sync
      otree->againstElectronLooseMVA6_2 = analysisTree.tau_againstElectronLooseMVA6[tauIndex_2];
      otree->againstElectronMediumMVA6_2 = analysisTree.tau_againstElectronMediumMVA6[tauIndex_2];
      otree->againstElectronVTightMVA6_2 = analysisTree.tau_againstElectronVTightMVA6[tauIndex_2];
      ////

      ///AN///
      otree->deepTauVsJetRaw_2 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tauIndex_2];
      otree->deepTauVsEleRaw_2 = analysisTree.tau_byDeepTau2017v2p1VSeraw[tauIndex_2];
      otree->deepTauVsMuRaw_2 = analysisTree.tau_byDeepTau2017v2p1VSmuraw[tauIndex_2];
      /////

      
      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;
      // svfit variables
      otree->m_sv   = -9999;//Merijn updated for the DNN
      otree->pt_sv  = -9999;
      otree->eta_sv = -9999;
      otree->phi_sv = -9999;
      otree->met_sv = -9999;
      otree->mt_sv = -9999;


      

	/*
	cout<<endl<<"tauLV_2.Px()="<<tauLV_2.Px()<<endl<<"tauLV_2.Py()="<<tauLV_2.Py()<<endl<<"tauLV_2.Pz()="<<tauLV_2.Pz()<<endl<<"tauLV_2.Pt()="<<tauLV_2.Pt()<<endl;
	*/

          check<<endl<<"ok16"<<endl;
      if(ApplySVFit)
     	if (otree->njetspt20>0) svfit_variables("tt", &analysisTree, otree, &cfg, inputFile_visPtResolution);
      cout<<endl<<"det0 check 1"<<endl;

      if(isDY||isWJets){
     	otree->gen_noutgoing=analysisTree.genparticles_noutgoing;
     	//gen_noutgoing_NLO=analysisTree.genparticles_noutgoing_NLO;
      }
          check<<endl<<"ok17"<<endl;
      
      cout<<endl<<"det0 check 2"<<endl;
      
      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
      ////////////////////////////////////////////////////////////
      // MET Recoil Corrections
      ////////////////////////////////////////////////////////////
      cout<<endl<<"det0 check 3"<<endl;
          check<<endl<<"ok18"<<endl;

      otree->njetshad = otree->njets;

      otree->puppimet = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex + analysisTree.puppimet_ey*analysisTree.puppimet_ey);
      otree->puppimetphi = TMath::ATan2(analysisTree.puppimet_ey,analysisTree.puppimet_ex);
      otree->puppimetcov00 = analysisTree.puppimet_sigxx;
      otree->puppimetcov01 = analysisTree.puppimet_sigxy;
      otree->puppimetcov10 = analysisTree.puppimet_sigyx;
      otree->puppimetcov11 = analysisTree.puppimet_sigyy;

      if (!isData && applyRecoilCorrections && (isDY || isWJets ) ){
     	genV = genTools::genV(analysisTree);
     	genL = genTools::genL(analysisTree);
     	if(isWJets) otree->njetshad += 1;
      }

      cout<<endl<<"puppimet="<<otree->puppimet<<endl<<"otree->puppimetphi="<<otree->puppimetphi<<endl;
      
      // applying recoil correc. to puppi met
      genTools::RecoilCorrections( *recoilPFMetCorrector, 
     				   (!isData && applyRecoilCorrections && (isDY || isWJets )) * genTools::MeanResolution,
     				   otree->puppimet, otree->puppimetphi,
     				   genV.Px(), genV.Py(),
     				   genL.Px(), genL.Py(),
     				   otree->njetshad,
     				   otree->met_rcmr, otree->metphi_rcmr
     				   );
      
          check<<endl<<"ok19"<<endl;

      // overwriting with recoil-corrected values 
      otree->puppimet = otree->met_rcmr;
      otree->puppimetphi = otree->metphi_rcmr; 

      cout<<endl<<"puppimet="<<otree->puppimet<<endl<<"otree->puppimetphi="<<otree->puppimetphi<<endl;
      
      TLorentzVector metLV; metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi),
     					  otree->met*TMath::Sin(otree->metphi),
     					  0,
     					  TMath::Sqrt( otree->met*TMath::Sin(otree->metphi)*otree->met*TMath::Sin(otree->metphi) +
     						       otree->met*TMath::Cos(otree->metphi)*otree->met*TMath::Cos(otree->metphi)));

      TLorentzVector puppimetLV; puppimetLV.SetXYZT(otree->puppimet*TMath::Cos(otree->puppimetphi),
     					  otree->puppimet*TMath::Sin(otree->puppimetphi),
     					  0,
     					  TMath::Sqrt( otree->puppimet*TMath::Sin(otree->puppimetphi)*otree->puppimet*TMath::Sin(otree->puppimetphi) + otree->puppimet*TMath::Cos(otree->puppimetphi)*otree->puppimet*TMath::Cos(otree->puppimetphi)));
      


      ///////////////////////////////////////////
      ////////TES and propagation to MET/////////
      //////////////////////////////////////////

      check<<endl<<"befroe tes:"<<endl<<"analysisTree.tau_pt[tauIndex_1]="<<analysisTree.tau_pt[tauIndex_1]<<endl<<"analysisTree.tau_pt[tauIndex_2]="<<analysisTree.tau_pt[tauIndex_2]<<endl<<"tauLV_1.Pt()="<<tauLV_1.Pt()<<endl<<"tauLV_2.Pt()="<<tauLV_2.Pt()<<endl;

      float shift_tes_1 = 1.0, shift_tes_2 = 1.0; 
	if (otree->gen_match_1 == 5){
	  if (otree->dm_1 == 0) shift_tes_1 = shift_tes_1prong; 
          else if (otree->dm_1 == 1) shift_tes_1 = shift_tes_1p1p0; 
	  else if (otree->dm_1 == 10) shift_tes_1 = shift_tes_3prong;
	  else if (otree->dm_1 == 11) shift_tes_1 = shift_tes_3p1p0;
	}
	else if(otree->gen_match_1 == 1 || otree->gen_match_1 == 3){
	  if (otree->dm_1 == 0) shift_tes_1 = shift_tes_1prong_efake; 
          else if (otree->dm_1 == 1) shift_tes_1 = shift_tes_1p1p0_efake;
	} 
	else if(otree->gen_match_1 == 2 || otree->gen_match_1 == 4){
	  if (otree->dm_1 == 0) shift_tes_1 = shift_tes_1prong_mufake; 
          else if (otree->dm_1 == 1) shift_tes_1 = shift_tes_1p1p0_mufake; 
	}

	if (otree->gen_match_2 == 5){
	  if (otree->dm_2 == 0) shift_tes_2 = shift_tes_1prong; 
          else if (otree->dm_2 == 1) shift_tes_2 = shift_tes_1p1p0; 
	  else if (otree->dm_2 == 10) shift_tes_2 = shift_tes_3prong;
	  else if (otree->dm_2 == 11) shift_tes_2 = shift_tes_3p1p0;
	}
	else if(otree->gen_match_2 == 1 || otree->gen_match_2 == 3){
	  if (otree->dm_2 == 0) shift_tes_2 = shift_tes_1prong_efake; 
          else if (otree->dm_2 == 1) shift_tes_2 = shift_tes_1p1p0_efake;
	} 
	else if(otree->gen_match_2 == 2 || otree->gen_match_2 == 4){
	  if (otree->dm_2 == 0) shift_tes_2 = shift_tes_1prong_mufake; 
          else if (otree->dm_2 == 1) shift_tes_2 = shift_tes_1p1p0_mufake; 
	}      
	
	      
	tauLV_1 = tauLV_1*shift_tes_1;
	tauLV_2 = tauLV_2*shift_tes_2;
	metLV = metLV - (1-shift_tes_1)*tauLV_1 - (1-shift_tes_2)*tauLV_2;
	puppimetLV = puppimetLV - (1-shift_tes_1)*tauLV_1 - (1-shift_tes_2)*tauLV_2;
	otree->pt_1 = tauLV_1.Pt();
        otree->m_1 = tauLV_1.M();
	otree->pt_2 = tauLV_2.Pt();
        otree->m_2 = tauLV_2.M();
	otree->puppimet = puppimetLV.Pt();
	otree->puppimetphi = puppimetLV.Phi();
	otree->met = metLV.Pt();
	otree->metphi = metLV.Phi();


	check<<endl<<"after tes:"<<endl<<"analysisTree.tau_pt[tauIndex_1]="<<analysisTree.tau_pt[tauIndex_1]<<endl<<"analysisTree.tau_pt[tauIndex_2]="<<analysisTree.tau_pt[tauIndex_2]<<endl<<"tauLV_1.Pt()="<<tauLV_1.Pt()<<endl<<"tauLV_2.Pt()="<<tauLV_2.Pt()<<endl;



	//////////////////////////////END TES//////////////////////

	otree->pt_tt = (diTauLV+puppimetLV).Pt();




      // mt TOT

      float mtTOT = 2*(otree->pt_1)*metLV.Pt()*(1-cos(DeltaPhi(tauLV_1,metLV)));
      mtTOT += 2*(otree->pt_2)*metLV.Pt()*(1-cos(DeltaPhi(tauLV_2,metLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(tauLV_1,tauLV_2))); 
      otree->mt_tot = TMath::Sqrt(mtTOT);
      //QCD supresser parameter                                                                                                                                                                                    
      double per_px=(analysisTree.tau_px[tauIndex_1]+analysisTree.tau_px[tauIndex_2]+otree->met*TMath::Cos(otree->metphi));
      double per_py=(analysisTree.tau_py[tauIndex_1]+analysisTree.tau_py[tauIndex_2]+otree->met*TMath::Sin(otree->metphi));
      otree->Prompt_pT=sqrt(per_px*per_px+per_py*per_py);
      
      //
      otree->mt_1=mT(tauLV_1,metLV);
      otree->mt_2=mT(tauLV_2,metLV);
      otree->pzetamiss = otree->met*TMath::Cos(otree->metphi)*zetaX + otree->met*TMath::Sin(otree->metphi)*zetaY;
      otree->puppipzetamiss = otree->puppimet*TMath::Cos(otree->puppimetphi)*zetaX + otree->puppimet*TMath::Sin(otree->puppimetphi)*zetaY;
      
      //additional Sync variables //by AN
      //find correct refitted vertex
      bool useRefittedVtx = true, withBS=true;
      otree->pvx=analysisTree.primvertex_x;
      otree->pvy=analysisTree.primvertex_y;
      otree->pvz=analysisTree.primvertex_z; //use nominal PV by default

          check<<endl<<"ok20"<<endl;

      //choose vtx
      if(useRefittedVtx && !withBS){
        int vtxIndex = -1;
        for(unsigned int indx = 0; indx < analysisTree.refitvertex_count; indx++){
          if((analysisTree.refitvertex_tauIndex[indx][0] == tauIndex_1 && analysisTree.refitvertex_tauIndex[indx][1] == tauIndex_2) ||
             (analysisTree.refitvertex_tauIndex[indx][1] == tauIndex_1 && analysisTree.refitvertex_tauIndex[indx][0] == tauIndex_2)){
            vtxIndex = indx;
            break;
          }
        }

        if(vtxIndex>= 0){
          otree->pvx=analysisTree.refitvertex_x[vtxIndex];
          otree->pvy=analysisTree.refitvertex_y[vtxIndex];
          otree->pvz=analysisTree.refitvertex_z[vtxIndex];
        }
      }
      else if(useRefittedVtx && withBS){
        int vtxIndex = -1;
        for(unsigned int indx = 0; indx < analysisTree.refitvertexwithbs_count; indx++){
          if((analysisTree.refitvertexwithbs_tauIndex[indx][0] == tauIndex_1 && analysisTree.refitvertexwithbs_tauIndex[indx][1] == tauIndex_2) || (analysisTree.refitvertexwithbs_tauIndex[indx][1] == tauIndex_1 && analysisTree.refitvertexwithbs_tauIndex[indx][0] == tauIndex_2)){vtxIndex = indx;
            break;
          }
        }

        if(vtxIndex>= 0){
          otree->pvx=analysisTree.refitvertexwithbs_x[vtxIndex];
          otree->pvy=analysisTree.refitvertexwithbs_y[vtxIndex];
          otree->pvz=analysisTree.refitvertexwithbs_z[vtxIndex];
        }
      }

          check<<endl<<"ok21"<<endl;
	  /*
      // Measurement of CP observable
      TVector3 IpvH_reco_noBS_1 = Helical_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,false,false);
      TVector3 IpvH_reco_BS_1 = Helical_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,false,true);
      TVector3 IpvH_refit_BS_1 = Helical_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,true,true);

      TVector3 IpvL_reco_noBS_1 = Linear_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,false,false);
      TVector3 IpvL_reco_BS_1 = Linear_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,false,true);
      TVector3 IpvL_refit_BS_1 = Linear_IPV(&analysisTree,otree,tauIndex_1,tauIndex_2,true,true);
      
      SMatrixStd3F IP_cov_pv_1 = ipVecCovHelic(&analysisTree, tauIndex_1, tauIndex_2, false, false);
      SMatrixStd3F IP_cov_pv_withBS_1 = ipVecCovHelic(&analysisTree, tauIndex_1, tauIndex_2, false, true);
      SMatrixStd3F IP_cov_Refitpv_withBS_1 = ipVecCovHelic(&analysisTree, tauIndex_1, tauIndex_2, true, true);

      TVector3 IpvH_reco_noBS_2 = Helical_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,false,false);
      TVector3 IpvH_reco_BS_2 = Helical_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,false,true);
      TVector3 IpvH_refit_BS_2 = Helical_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,true,true);

      TVector3 IpvL_reco_noBS_2 = Linear_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,false,false);
      TVector3 IpvL_reco_BS_2 = Linear_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,false,true);
      TVector3 IpvL_refit_BS_2 = Linear_IPV(&analysisTree,otree,tauIndex_2,tauIndex_1,true,true);
      
      SMatrixStd3F IP_cov_pv_2 = ipVecCovHelic(&analysisTree, tauIndex_2, tauIndex_1, false, false);
      SMatrixStd3F IP_cov_pv_withBS_2 = ipVecCovHelic(&analysisTree, tauIndex_2, tauIndex_1, false, true);
      SMatrixStd3F IP_cov_Refitpv_withBS_2 = ipVecCovHelic(&analysisTree, tauIndex_2, tauIndex_1, true, true);


      ImpactParameter* ip = new ImpactParameter();
      double IPsigHelical_pv_1 = ip->CalculateIPSignificanceHelical(IpvH_reco_noBS_1,IP_cov_pv_1);
      double IPsigHelical_PVwithBS_1 = ip->CalculateIPSignificanceHelical(IpvH_refit_BS_1,IP_cov_pv_withBS_1);
      double IPsigHelical_refitPVwithBS_1 = ip->CalculateIPSignificanceHelical(IpvH_refit_BS_1,IP_cov_Refitpv_withBS_1);
      
      double IPsigTangent_pv_1 = ipVecCovSigmaTangent(&analysisTree,tauIndex_1, tauIndex_2, false,false, IpvL_reco_noBS_1);
      double IPsigTangent_PVwithBS_1 = ipVecCovSigmaTangent(&analysisTree,tauIndex_1, tauIndex_2, false,true, IpvL_reco_BS_1);
      double IPsigTangent_refitPVwithBS_1 = ipVecCovSigmaTangent(&analysisTree,tauIndex_1, tauIndex_2, true,true, IpvL_refit_BS_1);

      double IPsigHelical_pv_2 = ip->CalculateIPSignificanceHelical(IpvH_reco_noBS_2,IP_cov_pv_2);
      double IPsigHelical_PVwithBS_2 = ip->CalculateIPSignificanceHelical(IpvH_refit_BS_2,IP_cov_pv_withBS_2);
      double IPsigHelical_refitPVwithBS_2 = ip->CalculateIPSignificanceHelical(IpvH_refit_BS_2,IP_cov_Refitpv_withBS_2);
      
      double IPsigTangent_pv_2 = ipVecCovSigmaTangent(&analysisTree,tauIndex_2, tauIndex_1, false,false, IpvL_reco_noBS_2);
      double IPsigTangent_PVwithBS_2 = ipVecCovSigmaTangent(&analysisTree,tauIndex_2, tauIndex_1, false,true, IpvL_reco_BS_2);
      double IPsigTangent_refitPVwithBS_2 = ipVecCovSigmaTangent(&analysisTree,tauIndex_2, tauIndex_1, true,true, IpvL_refit_BS_2);
      //IP branches
     
     	otree->ipvecHelic_pv_x_1 = IpvH_reco_noBS_1.X();
     	otree->ipvecHelic_pv_y_1 = IpvH_reco_noBS_1.Y();
     	otree->ipvecHelic_pv_z_1 = IpvH_reco_noBS_1.Z();

     	otree->ipvecsigHelic_pv_1 = IPsigHelical_pv_1;
	
     
     
     	otree->ipvecHelic_pvwBS_x_1 = IpvH_reco_BS_1.X();
     	otree->ipvecHelic_pvwBS_y_1 = IpvH_reco_BS_1.Y();
     	otree->ipvecHelic_pvwBS_z_1 = IpvH_reco_BS_1.Z();

     	otree->ipvecsigHelic_pvwBS_1 = IPsigHelical_PVwithBS_1;

     	otree->ipvecHelic_refpvBS_x_1 = IpvH_refit_BS_1.X();
     	otree->ipvecHelic_refpvBS_y_1 = IpvH_refit_BS_1.Y();
     	otree->ipvecHelic_refpvBS_z_1 = IpvH_refit_BS_1.Z();

     	otree->ipvecsigHelic_refpvBS_1 = IPsigHelical_refitPVwithBS_1;

     	otree->ipvecTangent_pv_x_1 = IpvL_reco_noBS_1.X();
     	otree->ipvecTangent_pv_y_1 = IpvL_reco_noBS_1.Y();
     	otree->ipvecTangent_pv_z_1 = IpvL_reco_noBS_1.Z();

     	otree->ipvecsigTangent_pv_1 = IPsigTangent_pv_1;

     	otree->ipvecTangent_pvwBS_x_1 = IpvL_reco_BS_1.X();
     	otree->ipvecTangent_pvwBS_y_1 = IpvL_reco_BS_1.Y();
     	otree->ipvecTangent_pvwBS_z_1 = IpvL_reco_BS_1.Z();

     	otree->ipvecsigTangent_pvwBS_1 = IPsigTangent_PVwithBS_1;
	
     	otree->ipvecTangent_refpvBS_x_1 = IpvL_refit_BS_1.X();
     	otree->ipvecTangent_refpvBS_y_1 = IpvL_refit_BS_1.Y();
     	otree->ipvecTangent_refpvBS_z_1 = IpvL_refit_BS_1.Z();

     	otree->ipvecsigTangent_refpvBS_1 = IPsigTangent_refitPVwithBS_1;



     	otree->ipvecHelic_pv_x_2 = IpvH_reco_noBS_2.X();
     	otree->ipvecHelic_pv_y_2 = IpvH_reco_noBS_2.Y();
     	otree->ipvecHelic_pv_z_2 = IpvH_reco_noBS_2.Z();

     	otree->ipvecsigHelic_pv_2 = IPsigHelical_pv_2;
	
     
     
     	otree->ipvecHelic_pvwBS_x_2 = IpvH_reco_BS_2.X();
     	otree->ipvecHelic_pvwBS_y_2 = IpvH_reco_BS_2.Y();
     	otree->ipvecHelic_pvwBS_z_2 = IpvH_reco_BS_2.Z();

     	otree->ipvecsigHelic_pvwBS_2 = IPsigHelical_PVwithBS_2;

     	otree->ipvecHelic_refpvBS_x_2 = IpvH_refit_BS_2.X();
     	otree->ipvecHelic_refpvBS_y_2 = IpvH_refit_BS_2.Y();
     	otree->ipvecHelic_refpvBS_z_2 = IpvH_refit_BS_2.Z();

     	otree->ipvecsigHelic_refpvBS_2 = IPsigHelical_refitPVwithBS_2;

     	otree->ipvecTangent_pv_x_2 = IpvL_reco_noBS_2.X();
     	otree->ipvecTangent_pv_y_2 = IpvL_reco_noBS_2.Y();
     	otree->ipvecTangent_pv_z_2 = IpvL_reco_noBS_2.Z();

     	otree->ipvecsigTangent_pv_2 = IPsigTangent_pv_2;

     	otree->ipvecTangent_pvwBS_x_2 = IpvL_reco_BS_2.X();
     	otree->ipvecTangent_pvwBS_y_2 = IpvL_reco_BS_2.Y();
     	otree->ipvecTangent_pvwBS_z_2 = IpvL_reco_BS_2.Z();

     	otree->ipvecsigTangent_pvwBS_2 = IPsigTangent_PVwithBS_2;
	
     	otree->ipvecTangent_refpvBS_x_2 = IpvL_refit_BS_2.X();
     	otree->ipvecTangent_refpvBS_y_2 = IpvL_refit_BS_2.Y();
     	otree->ipvecTangent_refpvBS_z_2 = IpvL_refit_BS_2.Z();

     	otree->ipvecsigTangent_refpvBS_2 = IPsigTangent_refitPVwithBS_2;
	  */  
	check<<endl<<"ok22"<<endl;		
            
      otree->Fill();
      selEvents++;

          check<<endl<<"ok23"<<endl;

    }// event loop
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;

  }// close each file
   file->cd("");
  file->Write();
  file->Close();
  delete file;
 
}//close of main
