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
  string Channel="mutau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  
  

////////////muons

  const double ptMuonCut   = cfg.get<double>("ptMuonCut");
  const double etaMuonCut     = cfg.get<double>("etaMuonCut");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonLowCut  = cfg.get<double>("isoMuonLowCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

////// tau
  const double ptTauCut = cfg.get<double>("ptTauCut"); 
  const double etaTauCut = cfg.get<double>("etaTauCut"); 


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


  const string dataBaseDir = cfg.get<string>("DataBaseDir");



  //const string MuonidIsoEffFileBCDEFGH = cfg.get<string>("MuonidIsoEffFileBCDEFGH");
  //const string MuontrigEffFileBCDEFGH = cfg.get<string>("MuontrigEffFileBCDEFGH");
  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");

  // topSingleMuonTriggerFile
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");

  const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");


  // vertex distributions filenames and histname

  const double bTag   = cfg.get<double>("bTag");
  string cmsswBase = (getenv ("CMSSW_BASE"));
 

  xs=1;fact=1;fact2=1;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;



  bool doThirdLeptVeto=true;
  bool doMuVeto=true;
  char ff[100];
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  int nTotalFiles = 0;

  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  std::ifstream fileList(ff);
  std::ifstream fileList0(ff);
  std::string ntupleName("ntuple/AC1B");


  bool isDelphes = false;
  if (  string::npos != rootFileName.find("ntuple") ) isDelphes=true;
  isDelphes=true;
  if (isDelphes) ntupleName = "AC1B";
	
  //std::string initNtupleName("initroottree/AC1B");
  std::string initNtupleName("ntuple/AC1B");


  TString SaveDir=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<std::endl;  
  // output fileName with histograms


  TFile * file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"update");
  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  SetupTree(); 
  
  Float_t         sumPt;
  T->Branch("sumPt", &sumPt, "sumPt/F");
  
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
    }


    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);
	if (!isData && !WithInit)
		{    
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
			{
			analysisTree.GetEntry(iEntry);
			histWeightsH->Fill(0.,analysisTree.genweight);
			}
		}
  	float genweights=1.;



    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 


      Float_t weight = 1.;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%5000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      gen_weight = 1.;


	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;
					


      vector<int> muons; muons.clear();
      sumPt = 0;
           float  sumPtx = 0;
           float  sumPty = 0;

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
		  
	sumPtx+=analysisTree.muon_px[im]; sumPty+=analysisTree.muon_py[im];
	
	if (analysisTree.muon_pt[im]<SingleMuonTriggerPtCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if ( fabs(analysisTree.muon_charge[im]) != 1) continue;

	muons.push_back((int)im);
//std::cout <<analysisTree.muon_pt[im]<<"  "<<fabs(analysisTree.muon_eta[im])<<"  "<<analysisTree.muon_isTight[im]<<"  "<<analysisTree.muon_charge[im]<<"  " <<muons.size()<<std::endl;

      }
      
      if (muons.size()==0) continue;
      

      vector<int> taus; taus.clear();
      vector<int> tausDelphes; tausDelphes.clear();

      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	 sumPtx+=analysisTree.tau_px[it]; sumPty+=analysisTree.tau_py[it];

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	
	bool isTauMatched = analysisTree.tau_genmatch[it];
	
	  taus.push_back((int)it);
//std::cout <<analysisTree.tau_pt[it]<<"  "<<fabs(analysisTree.tau_eta[it])<<"  "<<analysisTree.tau_charge[it]<<"  " <<taus.size()<<std::endl;

	}

      if (taus.size()==0)  continue;


      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float MuPtMin  = 0;
      float ptMu = 0;
      float ptTau = 0;
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;
      
	for (unsigned int im=0; im<muons.size(); ++im) {
	unsigned int mIndex  = muons.at(im);

	float MuPt = analysisTree.muon_pt[mIndex];

	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);
	  if (dR<dRleptonsCutmutau) continue;
 

          if ((int)mIndex!=(int)mu_index) {
            if (MuPt==MuPtMin) {
              if (analysisTree.muon_pt[mIndex]>ptMu) {
                MuPtMin  = MuPt;
                ptMu = analysisTree.muon_pt[mIndex];
                mu_index =(int)mIndex;
                ptTau = analysisTree.tau_pt[tIndex];
                tau_index =(int)tIndex;
              }
            }
            else if (MuPt<MuPtMin) {
              MuPtMin  = MuPt;
              ptMu = analysisTree.muon_pt[mIndex];
              mu_index =(int)mIndex;
              ptTau = analysisTree.tau_pt[tIndex];
              tau_index =(int)tIndex;
            }
          }
          else {
              if (analysisTree.tau_pt[tIndex]>ptTau) {
                ptTau = analysisTree.tau_pt[tIndex];
                tau_index =(int)tIndex;
              //}
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
			    analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);
	 	if (dR<dRleptonsCutmutau) continue;
		tausDelphes.push_back(tIndex);

		}
	}



      if ((int)tau_index<0) continue;
	

      
      double q = analysisTree.tau_charge[tau_index] * analysisTree.muon_charge[mu_index];
      event_sign  = q;

	double dRmutau = deltaR(analysisTree.tau_eta[(int)tau_index],analysisTree.tau_phi[(int)tau_index],
				analysisTree.muon_eta[(int)mu_index],analysisTree.muon_phi[(int)mu_index]);
	if (dRmutau < 0.5) continue;

	TFR_weight=1;


  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;

  event_secondLeptonVeto = false;
  event_thirdLeptonVeto = false;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	sumPtx+=analysisTree.electron_px[ie]; sumPty+=analysisTree.electron_py[ie];
	if (analysisTree.electron_pt[ie]<20) continue;
	if (fabs(analysisTree.electron_eta[ie])>2.4) continue;
	if (analysisTree.electron_relIso[ie]>0.5) continue;
	foundExtraElectron = true;
      }

      // looking for extra muon's (dimuon veto)
      bool foundExtraMuon = false;
      vector<int> mu_dimuons; mu_dimuons.clear(); 
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if ((int)im==(int)mu_index) continue;


	if (analysisTree.muon_pt[im]>20&&
	    fabs(analysisTree.muon_eta[im])<2.4&&
	    fabs(analysisTree.muon_charge[im]) ==1)  
	{
	    float dRmuons = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				   analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);
	    if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[im]*analysisTree.muon_charge[mu_index]<0)) dilepton_veto = true;
 	  }

	if (analysisTree.muon_pt[im]<20) continue;
	if (fabs(analysisTree.muon_eta[im])>2.4) continue;
	foundExtraMuon = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;
    

   	event_secondLeptonVeto = dilepton_veto;

     if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;




      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)muons.size();

      for (unsigned int im=0;im<muons.size(); ++im){
	unsigned int mIndex = muons[im];
	mu_px[im]=analysisTree.muon_px[mIndex];
	mu_py[im]=analysisTree.muon_py[mIndex];
	mu_pz[im]=analysisTree.muon_pz[mIndex];
	mu_eta[im]=analysisTree.muon_eta[mIndex];
	mu_pt[im]=analysisTree.muon_pt[mIndex];
	mu_phi[im]=analysisTree.muon_phi[mIndex];
	mu_charge[im]=analysisTree.muon_charge[mIndex];
	mu_dxy[im]=analysisTree.muon_dxy[mIndex];
	mu_dz[im]=analysisTree.muon_dz[mIndex];
	mu_dxyerr[im]=analysisTree.muon_dxyerr[mIndex];
	mu_dzerr[im]=analysisTree.muon_dzerr[mIndex];
	mu_relIso[im]  = analysisTree.muon_relIso[mIndex] ;

      }



      ta_count=(int)taus.size();

      for (unsigned int it=0;it<taus.size(); ++it){
	unsigned int itt = taus[it];
	ta_px[it]=analysisTree.tau_px[itt];
	ta_py[it]=analysisTree.tau_py[itt];
	ta_pz[it]=analysisTree.tau_pz[itt];
	ta_eta[it]=analysisTree.tau_eta[itt];
	ta_pt[it]=analysisTree.tau_pt[itt];
	ta_phi[it]=analysisTree.tau_phi[itt];
	ta_charge[it]=analysisTree.tau_charge[itt];
	ta_dxy[it]=analysisTree.tau_dxy[itt];
	ta_dz[it]=analysisTree.tau_dz[itt];
	//

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




      float jetEta = 2.4;
      float DRmax = 0.5;
      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> bjets; bjets.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;

	int counter_cleaned_jets = 0;

      for (int n=0;n<30;n++){
	jets_cleaned[n]=-1;
      }



      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


      	bool btagged= false;

	bool cleanedJet = true;

	double Dr=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
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
      nbtag = bjets.size();


      met_ex = analysisTree.pfmet_pt*TMath::Cos(analysisTree.pfmet_phi);
      met_ey = analysisTree.pfmet_pt*TMath::Sin(analysisTree.pfmet_phi);
      met_pt = analysisTree.pfmet_pt;
      met_phi = analysisTree.pfmet_phi;

      all_weight = weight;





//////////////// TFR for Delphes!!!!1


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

		//NFakeJets = FakeJets->GetBinContent(FakeJets->GetXaxis()->FindBin(t_pt),FakeJets->GetYaxis()->FindBin(t_eta));
		//NAllJets = AllJets->GetBinContent(AllJets->GetXaxis()->FindBin(t_pt),AllJets->GetYaxis()->FindBin(t_eta));
		
		//Teff_weight = 0.97*0.77*((0.32*(ta_pt[taus_index]*0 + 1) + 0.01*ta_pt[taus_index] - 0.000054*ta_pt[taus_index]*ta_pt[taus_index])*(ta_pt[taus_index]<100)+0.78*(ta_pt[taus_index]>100)); 
		//TFR_weight = ((( -0.00621816*(ta_pt[taus_index]*0+1)+0.00130097*ta_pt[taus_index]-0.0000219642*ta_pt[taus_index]*ta_pt[taus_index]+0.000000149393*pow(ta_pt[taus_index],3)-0.000000000458972e*pow(ta_pt[taus_index],4)+0.000000000000527983e-13*pow(ta_pt[taus_index],5)))*(ta_pt[taus_index]<250) + 0.0032*(ta_pt[taus_index]>250)); 

		TFR_weight = ((( -0.00621816*(t_pt*0+1)+0.00130097*t_pt-0.0000219642*t_pt*t_pt+0.000000149393*pow(t_pt,3)-0.000000000458972*pow(t_pt,4)+0.000000000000527983*pow(t_pt,5)))*(t_pt<250) + 0.0032*(t_pt>250)); 
		if (analysisTree.tau_genmatch[tausDelphes.at(it)]) TFR_weight = 0.97*0.77*((0.32*(t_pt*0 + 1) + 0.01*t_pt - 0.000054*t_pt*t_pt)*(t_pt<100)+0.78*(t_pt>100));
		TFR_weights.push_back(TFR_weight);
		if (it == 1) {TFR_weight = TFR_weight*(1-TFR_weights.at(0)); }
		if (it == 2) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1)); }
		if (it == 3) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1))*(1-TFR_weights.at(2)); }
		if (it == 4) {TFR_weight = TFR_weight*(1-TFR_weights.at(0))*(1-TFR_weights.at(1))*(1-TFR_weights.at(2))*(1-TFR_weights.at(3)); }
		if (it > 4.5) continue;

		//cout <<"FakeJets->GetSumOfWeights()  "<< FakeJets->GetSumOfWeights()<<"  FakeJets->GetXaxis()->FindBin(analysisTree.tau_pt[(int)tau_index])  "<< FakeJets->GetXaxis()->FindBin(t_pt)<<"  FakeJets->GetYaxis()->FindBin(analysisTree.tau_eta[(int)tau_index])  "<< FakeJets->GetYaxis()->FindBin(t_eta)<<endl;
		//cout <<"SFFakeRate  "<< TFR_weight<<"  NFakeJets  "<< NFakeJets<<"  NAllJets  "<< NAllJets<<endl;


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
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}
