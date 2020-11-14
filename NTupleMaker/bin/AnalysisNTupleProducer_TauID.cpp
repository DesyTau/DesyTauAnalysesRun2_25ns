#include "NtupleMaker_TauID_all_eras_definitions.h"
#include "NtupleMaker_TauID_all_eras_functions.h"

#define pi   3.14159265358979312
#define d2r  1.74532925199432955e-02
#define r2d  57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 	 0.105658
#define tauMass 	 1.77682
#define pionMass 	 0.1396

TLorentzVector leadingTrackLV(AC1B* analysisTree, unsigned int tauIndex) {

  TLorentzVector leadLV; leadLV.SetPxPyPzE(0,0,0,0);
  int decayMode = analysisTree->gentau_decayMode[tauIndex];
  double dRcut=0.3;
  TLorentzVector Tau;
  Tau.SetPxPyPzE(analysisTree->gentau_visible_px[tauIndex],
		 analysisTree->gentau_visible_py[tauIndex],
		 analysisTree->gentau_visible_pz[tauIndex],
		 analysisTree->gentau_visible_e[tauIndex]);
  
  int partId = 211;
  if (decayMode==8) partId = 13;
  if (decayMode==9) partId = 11;
  double maxPt = -1;
  int npart = analysisTree->genparticles_count;
  for(int i=0;i<npart;i++){
    if(fabs(analysisTree->genparticles_pdgid[i])==partId&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5||analysisTree->genparticles_info[i]==6)&&analysisTree->genparticles_isLastCopy[i]){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      double Pt = lvector.Pt();
      double dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut){
	if (Pt>maxPt) {
	  maxPt = Pt;
	  leadLV = lvector;
	}
      }
    }
  }
  return leadLV;

}
float getEmbeddedWeight(const AC1B * analysisTree, RooWorkspace* wEm) {

  std::vector<TLorentzVector> taus; taus.clear();
  float emWeight = 1;
  for (unsigned int igentau = 0; igentau < analysisTree->gentau_count; ++igentau) {
    TLorentzVector tauLV; tauLV.SetXYZT(analysisTree->gentau_px[igentau], 
					analysisTree->gentau_py[igentau],
					analysisTree->gentau_pz[igentau],
					analysisTree->gentau_e[igentau]);
    if (analysisTree->gentau_isPrompt[igentau]&&analysisTree->gentau_isFirstCopy[igentau]) {
      taus.push_back(tauLV);
    }
  }

  //  std::cout << "n taus = " << taus.size() << "  :  wEm = " << wEm << std::endl;

  if (taus.size() == 2) {
    double gt1_pt  = taus[0].Pt();
    double gt1_eta = taus[0].Eta();
    double gt2_pt  = taus[1].Pt();
    double gt2_eta = taus[1].Eta();
    wEm->var("gt_pt")->setVal(gt1_pt);
    wEm->var("gt_eta")->setVal(gt1_eta);
    double id1_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt_pt")->setVal(gt2_pt);
    wEm->var("gt_eta")->setVal(gt2_eta);
    double id2_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt1_pt")->setVal(gt1_pt);
    wEm->var("gt2_pt")->setVal(gt2_pt);
    wEm->var("gt1_eta")->setVal(gt1_eta);
    wEm->var("gt2_eta")->setVal(gt2_eta);
    double trg_emb = wEm->function("m_sel_trg_ic_ratio")->getVal();
    emWeight = id1_embed * id2_embed * trg_emb;
  }

  return emWeight;

}

int main(int argc, char * argv[]) {

   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   

   // ------------------- read config file --------------------
   // **** configuration
   Config cfg(argv[1]);

   string cmsswBase = (getenv("CMSSW_BASE"));
   
   // general settings
   const string era = cfg.get<string>("Era");
   const bool debug = cfg.get<bool>("Debug");
   const bool isData = cfg.get<bool>("IsData");
   const bool isEmbedded = cfg.get<bool>("IsEmbedded");
   const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
   const string jsonFile = cfg.get<string>("jsonFile");
   
   const string singleMuonHLTFilter = cfg.get<string>("SingleMuonHLTFilter");
   const string singleMuonHLTFilter1 = cfg.get<string>("SingleMuonHLTFilter1");

   const string singlePFTau180Trk50Name = cfg.get<string>("SinglePFTau180Trk50Name");
   const string singlePFTau180Trk50Name2 = cfg.get<string>("SinglePFTau180Trk50Name2");
   const string singlePFTau180Trk50oneprongName = cfg.get<string>("SinglePFTau180Trk50oneprongName");
   const string singlePFTau180Trk50oneprongName2 = cfg.get<string>("SinglePFTau180Trk50oneprongName2");

   const string doubleTauFilter1Name = cfg.get<string>("DoubleTauFilter1Name");
   const string doubleTauFilter2Name = cfg.get<string>("DoubleTauFilter2Name");
   const string doubleTauFilter3Name = cfg.get<string>("DoubleTauFilter3Name");
   const string doubleTauFilter4Name = cfg.get<string>("DoubleTauFilter4Name");

   const string muTauFilter1Name = cfg.get<string>("MuTauFilter1Name");
   const string muTauFilter2Name = cfg.get<string>("MuTauFilter2Name");
   const string muTauFilter3Name = cfg.get<string>("MuTauFilter3Name");
   const string muTauFilter4Name = cfg.get<string>("MuTauFilter4Name");


   const string metHLTName = cfg.get<string>("MetHLTName");

   TString SingleMuonHLTFilter(singleMuonHLTFilter);
   TString SingleMuonHLTFilter1(singleMuonHLTFilter1);
   TString SinglePFTau180Trk50Name(singlePFTau180Trk50Name);
   TString SinglePFTau180Trk50Name2(singlePFTau180Trk50Name2);
   TString SinglePFTau180Trk50oneprongName(singlePFTau180Trk50oneprongName);
   TString SinglePFTau180Trk50oneprongName2(singlePFTau180Trk50oneprongName2);

   TString DoubleTauFilter1Name(doubleTauFilter1Name);
   TString DoubleTauFilter2Name(doubleTauFilter2Name);
   TString DoubleTauFilter3Name(doubleTauFilter3Name);
   TString DoubleTauFilter4Name(doubleTauFilter4Name);

   TString MuTauFilter1Name(muTauFilter1Name);
   TString MuTauFilter2Name(muTauFilter2Name);
   TString MuTauFilter3Name(muTauFilter3Name);
   TString MuTauFilter4Name(muTauFilter4Name);

   TString MetHLTName(metHLTName); 
   
   // tau cuts
   const float ptTauCut  = cfg.get<float>("PtTauCut");
   const float etaTauCut = cfg.get<float>("EtaTauCut");
   const float ptTauTriggerCut = cfg.get<float>("PtTauTriggerCut");

   // muon selection
   const float ptMuCut       = cfg.get<float>("PtMuCut");
   const float etaMuCut      = cfg.get<float>("EtaMuCut");
   const float isoMuCut      = cfg.get<float>("IsoMuCut");
   const float dxyMuCut      = cfg.get<float>("dxyMuCut");  
   const float dzMuCut       = cfg.get<float>("dzMuCut");

   // electron selection
   const float ptEleCut   = cfg.get<float>("PtEleCut");
   const float etaEleCut  = cfg.get<float>("EtaEleCut");
   const float isoEleCut  = cfg.get<float>("IsoEleCut");
   const float dxyEleCut  = cfg.get<float>("dxyEleCut");
   const float dzEleCut   = cfg.get<float>("dzEleCut");
   
   const string puDataFile = cfg.get<string>("PileUpDataFile");
   const string puMCFile = cfg.get<string>("PileUpMCFile");
   const string samplenameForPUHist = cfg.get<string>("SampleNameForPUHist");
   
   // L1Tau cut   

   TString PUDataFile(puDataFile);
   TString PUMCFile(puMCFile);

   // trigger eff filename
   const string trigEffFileName = cfg.get<string>("TrigEffFileName");
   const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");

   TString workspace_filename = TString(cmsswBase) + "/src/" + TString(CorrectionWorkspaceFileName);
   TFile * workspaceFile = new TFile(workspace_filename);
   RooWorkspace * correctionWS = (RooWorkspace*)workspaceFile->Get("w");

   // **** end of configuration
   
   // ------------------- open files and trees --------------------------------
   std::string rootFileName(argv[2]);
   std::ifstream fileList(argv[2]);
   std::ifstream fileList0(argv[2]);
   std::string ntupleName("makeroottree/AC1B");
   std::string initNtupleName("initroottree/AC1B");
   std::string eventHistoName("eventCount/EventCount");
   std::string eventHistoNameData("makeroottree/nEvents");
   std::string weightsHistoName("eventCount/EventWeights");
   
   TString TStrName(rootFileName);
   std::cout <<TStrName <<std::endl;  

   TFile * file = new TFile(TStrName+TString(".root"),"recreate");
   file->cd("");
   
   ntuple_ = new TTree("NTuple","NTuple");
   
   TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
   TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
   
   SetupTauTree();

   // ------------------- initialize good run selection --------------------
   std::vector<Period> periods;
   string fullPathToJsonFile = cmsswBase + "/src/" + jsonFile;
   if (isData) ReadJson(periods,fullPathToJsonFile);
  
  
   // ------------------- initialize PU reweighting ------------------------
   PileUp * PUofficial = new PileUp(); 
   TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/"+PUDataFile,"read");
   TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/"+PUMCFile, "read");
   if (!isData) iniPU(PUofficial,filePUOfficial_data,filePUOfficial_MC, samplenameForPUHist);
   
   // MET trigger efficiency -----------
   map<int,TGraphAsymmErrors*>  map_trigEffData;
   map<int,TGraphAsymmErrors*>  map_trigEffMC;
   TFile * trigEffFile = new TFile(TString(cmsswBase)+"/src/"+trigEffFileName,"read");
   iniTriggerEfficiencies(trigEffFile,map_trigEffData,map_trigEffMC);
  
   // --------------------------------------------------------------------------
   // ------------------- open files --------------------------------------------
   // --------------------------------------------------------------------------
   
   // count number of files --->
   std::string dummy;
   while (fileList0 >> dummy) nTotalFiles++;
   
   for (int iF=0; iF<nTotalFiles; ++iF) { // loop over files
     
     std::string filen;
     fileList >> filen;
     
     // opening file
     std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
     TFile * file_ = TFile::Open(TString(filen));
     if (file_->IsZombie()) {
       cout << "Problems opening file : quitting program" << endl;
       exit(-1);
    }
     
     // accessing initroot tree
     TTree * _inittree = NULL;
     _inittree = (TTree*)file_->Get(TString(initNtupleName));
     FillHistGenWeights( _inittree, histWeightsH, isData);
     
     // accessing tree
     TTree * tree_ = (TTree*)file_->Get(TString(ntupleName));
     if (tree_==NULL) { 
        cout << "NTuple " << ntupleName << " is not found in file : quitting program" << endl;
        exit(-1);
     }
     AC1B analysisTree(tree_);
     Long64_t numberOfEntries = analysisTree.GetEntries();
     std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
     
     // ------------------ print era specific information -------------------------
     if( era == "2017") std::cout<<"EE noise cleaning applied!"<<std::endl;
     else std::cout<<"EE noise cleaning NOT applied!"<<std::endl;

     // ----------------------------------------- ---------------------------------
     // ------------------- event loop --------------------------------------------
     // ----------------------------------------- ---------------------------------
     for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
           
        analysisTree.GetEntry(iEntry);
        nEvents++;
        if (nEvents%10000==0) cout << "Processed " << nEvents << endl;
       
        // ------------------- initialize variables ------------------------------------
        SetDefaultValues();
        
        run_ = analysisTree.event_run;
        lumi_ = analysisTree.event_luminosityblock;
        event_ = analysisTree.event_nr;
        nVert_ = analysisTree.primvertex_count;
        
	// ------- Selecting events with 1 or 2 generator taus --------------------------
	ngentaus = 0;
        if (!isData || isEmbedded) {

	  // initial loop over gentaus
	  for (unsigned int igentau=0; igentau<analysisTree.gentau_count;++igentau) {
	    
	    if (analysisTree.gentau_decayMode[igentau]<0) continue;
	    if (analysisTree.gentau_decayMode[igentau]>9) continue;
	    if (!analysisTree.gentau_isPrompt[igentau]) continue;
	    if (!analysisTree.gentau_isLastCopy[igentau]) continue;
	    TLorentzVector GenVisTau; 
	    GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
			      analysisTree.gentau_visible_py[igentau],
			      analysisTree.gentau_visible_pz[igentau],
			      analysisTree.gentau_visible_e[igentau]);	    
	    //	    if (GenVisTau.Pt()<ptTauCut) continue;
	    //	    if (fabs(GenVisTau.Eta())>etaTauCut) continue;
	    
	    gentau_pt[ngentaus] = GenVisTau.Pt();
	    gentau_eta[ngentaus] = GenVisTau.Eta();
	    gentau_phi[ngentaus] = GenVisTau.Phi();
	    gentau_decay[ngentaus] = analysisTree.gentau_decayMode[igentau];
	    gentau_index[ngentaus] = igentau;
	    ngentaus++;
	    if (ngentaus>=2) break;
	  }

	}
	/*
	std::cout << "ngentaus = " << ngentaus << std::endl;
	for (unsigned int i=0; i<ngentaus; ++i) {
	  std::cout << i << "  pt=" << gentau_pt[i] 
		    << "  eta=" << gentau_eta[i]
		    << "  phi=" << gentau_phi[i]
		    << "  decay=" << gentau_decay[i] << std::endl;
	}
	
	std::cout << std::endl;
	*/

	if (ngentaus>2) continue;
	if (ngentaus==0) continue;
	bool isTwoLeptons = (ngentaus == 2) && (gentau_decay[0]>7) && (gentau_decay[1]>7);
	if (isTwoLeptons) continue;
	bool isOneTauLep = (ngentaus == 1) && (gentau_decay[0]>7);
	if (isOneTauLep) continue;
	

        if (event_<=0) {
           std::cout << "Event : " << event_ << std::endl;
           std::cout << "From NTuple = " << analysisTree.event_nr << std::endl;
        }
        
        if (debug) {
           cout << "Run                       = " << analysisTree.event_nr           << endl;
           cout << "Event                     = " << analysisTree.event_run          << endl; 
           cout << "Number of gen particles   = " << analysisTree.genparticles_count << endl;
           cout << "Number of taus            = " << analysisTree.gentau_count       << endl;
           cout << "Number of gen taus        = " << analysisTree.tau_count          << endl;
           cout << "Number of jets            = " << analysisTree.pfjet_count        << endl;
           cout << "Number of muons           = " << analysisTree.muon_count         << endl;
           cout << "Number of electrons       = " << analysisTree.electron_count     << endl;
           cout << "Number of tracks          = " << analysisTree.track_count        << endl;
           cout << "Number of PF jets         = " << analysisTree.pfjet_count        << endl;
           cout << "Number of dimuons         = " << analysisTree.dimuon_count       << endl;
           cout << "Number of l1 taus         = " << analysisTree.l1tau_count        << endl;
           cout << "Number of l1 muons        = " << analysisTree.l1muon_count       << endl;
           cout << "Number of l1 egamma       = " << analysisTree.l1egamma_count     << endl;
           cout << "Number of trigger objects = " << analysisTree.l1egamma_count     << endl;
        }	
        
        // Set MC relevant variables
        if (!isData || isEmbedded) {
           if (analysisTree.genweight<0) genWeight_ = -1;
           else                          genWeight_ = 1;
           
           weight_ *= genWeight_;
           
           genHt_ = analysisTree.genparticles_lheHt;
           npartons_ = analysisTree.genparticles_noutgoing;
           npartonsNLO_ = analysisTree.genparticles_noutgoing_NLO;
           lheWPt_   = analysisTree.genparticles_lheWPt;
        }
	bool isFlatPolarisation = TStrName.Contains("Uncorr");
	if (isFlatPolarisation) {
	  evenWeight_ = analysisTree.TauSpinnerWeight[0];
	  oddWeight_ = analysisTree.TauSpinnerWeight[2];
	}
	else {
	  evenWeight_ = 1.0;
	  oddWeight_ = 1.0;
	}
	if (isEmbedded) {
	  embWeight_ = getEmbeddedWeight(&analysisTree,correctionWS);	  
	}

        // *********************************
        // ***** accessing trigger info ****
        // *********************************        
	isSingleMuonHLTFilter = AccessTriggerInfo(analysisTree,SingleMuonHLTFilter,nSingleMuonHLTFilter);
	isSingleMuonHLTFilter1 = AccessTriggerInfo(analysisTree,SingleMuonHLTFilter1,nSingleMuonHLTFilter1);

        isSinglePFTau180Trk50Filter = AccessTriggerInfo(analysisTree,SinglePFTau180Trk50Name,nSinglePFTau180Trk50Filter);
        isSinglePFTau180Trk50Filter2 = AccessTriggerInfo(analysisTree,SinglePFTau180Trk50Name2,nSinglePFTau180Trk50Filter2);
        isSinglePFTau180Trk50oneprongFilter = AccessTriggerInfo(analysisTree,SinglePFTau180Trk50oneprongName,nSinglePFTau180Trk50oneprongFilter);
        isSinglePFTau180Trk50oneprongFilter2 = AccessTriggerInfo(analysisTree,SinglePFTau180Trk50oneprongName2,nSinglePFTau180Trk50oneprongFilter2);

	isDoubleTauFilter1 = AccessTriggerInfo(analysisTree,DoubleTauFilter1Name,nDoubleTauFilter1);
	isDoubleTauFilter2 = AccessTriggerInfo(analysisTree,DoubleTauFilter2Name,nDoubleTauFilter2);
	isDoubleTauFilter3 = AccessTriggerInfo(analysisTree,DoubleTauFilter3Name,nDoubleTauFilter3);
	isDoubleTauFilter4 = AccessTriggerInfo(analysisTree,DoubleTauFilter4Name,nDoubleTauFilter4);

	isMuTauFilter1 = AccessTriggerInfo(analysisTree,MuTauFilter1Name,nMuTauFilter1);
	isMuTauFilter2 = AccessTriggerInfo(analysisTree,MuTauFilter2Name,nMuTauFilter2);
	isMuTauFilter3 = AccessTriggerInfo(analysisTree,MuTauFilter3Name,nMuTauFilter3);
	isMuTauFilter4 = AccessTriggerInfo(analysisTree,MuTauFilter4Name,nMuTauFilter4);



	//	std::cout << "isSinglePFTau180Trk50Filter2 : " << isSinglePFTau180Trk50Filter2
	//		  << "   nSinglePFTau180Trk50Filter2 : " << nSinglePFTau180Trk50Filter2 << std::endl;

        // *****************************************
        // ***** primary vertex selection, *********
        // ***** PU weights, good run selection ****
        // *****************************************
        
        if (analysisTree.primvertex_count==0) continue; // at least one good primary vertex
        
        if (!isData && !isEmbedded) { 
           puWeight_ =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
           weight_ *= puWeight_; 
           if (debug)
              cout << "nPU = " << analysisTree.numtruepileupinteractions << " --> puweight = " << puWeight_ << endl;
        }
        
        if (isData && applyGoodRunSelection){
           int n=analysisTree.event_run;
           int lum = analysisTree.event_luminosityblock;      
           if (!GoodRunSelection(n,lum,periods)) continue;
        }        
        
	bool isMetHLT = false;
	for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	  TString trigName(it->first);
	  if (trigName.Contains(MetHLTName)) {
	    if (it->second==1) isMetHLT = true;
	  }
	}
	trigger_ = isMetHLT;
	  
        // *****************************
        // accessing PF MET ************
        // *****************************
        float pfmetcorr_ex = analysisTree.pfmetcorr_ex;
        float pfmetcorr_ey = analysisTree.pfmetcorr_ey;

        met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
        if (met_<1e-4) met_ = 1e-4;
        metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
        TLorentzVector lorentzVectorMet; lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);

        // **************************************************************************
        // **** accessing muons to be used in computation of metNoMu and mhtNoMu ****
        // **************************************************************************
        TLorentzVector lorentzVectorAllMuons; lorentzVectorAllMuons.SetXYZT(0,0,0,0);        
        // ----------------------- store muons that pass selection ----------------------- 
        for (unsigned int imuon=0; imuon<analysisTree.muon_count; ++imuon) {

	  float relIso=10;
	  // muon selection
	  bool isDRIso03 = false;
	  if (!PassesMuonSelection(analysisTree,imuon,ptMuCut,etaMuCut,dxyMuCut,dzMuCut,isDRIso03,isoMuCut, relIso, era)) continue;
	  TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[imuon],
						analysisTree.muon_py[imuon],
						analysisTree.muon_pz[imuon],
						muonMass);
	  lorentzVectorAllMuons += muonLV;
        }
        metNoMu_    =  (lorentzVectorMet+lorentzVectorAllMuons).Pt();

	// electron selection 
	vector<int> electrons; electrons.clear();
	vector<int> isoElectrons; isoElectrons.clear();
	vector<float> isoVarElectrons; isoVarElectrons.clear();
	for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {

	  if (analysisTree.electron_pt[ie]<ptEleCut) continue;
	  if (fabs(analysisTree.electron_eta[ie])>etaEleCut) continue;
	  if (fabs(analysisTree.electron_dxy[ie])>dxyEleCut) continue;
	  if (fabs(analysisTree.electron_dz[ie])>dzEleCut) continue;
	  bool electronMvaId = true;
	  electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie]>0.5;
	  if (!electronMvaId) continue;
	  if (!analysisTree.electron_pass_conversion[ie]) continue;
	  if (analysisTree.electron_nmissinginnerhits[ie]>1) continue;
	  electrons.push_back(ie);
         float rhoNeutral = analysisTree.rho;
         float  eA = getEffectiveArea( fabs(analysisTree.electron_superclusterEta[ie]) );
         float absIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie] +
            TMath::Max(0.0f,analysisTree.electron_r03_sumNeutralHadronEt[ie]+analysisTree.electron_r03_sumPhotonEt[ie]-eA*rhoNeutral);
	 float relIso = absIsoEle/analysisTree.electron_pt[ie];
	 if (relIso<isoEleCut) isoElectrons.push_back(ie);
	 isoVarElectrons.push_back(relIso);
	}
	nElec_ = electrons.size();
	nSelElec_ = isoElectrons.size();

	// muon selection 
	vector<int> muons; muons.clear();
	vector<int> isoMuons; isoMuons.clear();
	vector<float> isoVarMuons; isoVarMuons.clear();
	for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	  if (analysisTree.muon_pt[im]<ptMuCut) continue;
	  if (fabs(analysisTree.muon_eta[im])>etaMuCut) continue;
	  if (fabs(analysisTree.muon_dxy[im])>dxyMuCut) continue;
	  if (fabs(analysisTree.muon_dz[im])>dzMuCut) continue;
	  bool muonId = analysisTree.muon_isMedium[im];
	  if (!muonId) continue;
	  muons.push_back(im);
	  float absIso = analysisTree.muon_chargedHadIso[im];
	  float neutralIso = analysisTree.muon_neutralHadIso[im] +
	    analysisTree.muon_photonIso[im] -
	    0.5*analysisTree.muon_puIso[im];
	  neutralIso = TMath::Max(float(0),neutralIso);
	  absIso += neutralIso;
	  float relIso = absIso/analysisTree.muon_pt[im];
	  if (relIso<isoMuCut) isoMuons.push_back(im);
	  isoVarMuons.push_back(relIso);
	}
	nSelMuon_ = isoMuons.size();
	nMuon_ = muons.size();

	vector<unsigned int> indexSingleTau; indexSingleTau.clear();

	if (!isData || isEmbedded) {

	  unsigned int indexGenTau = 0;
	  unsigned int indexGenTau1 = 1;

	  if (ngentaus==1) {
	    indexGenTau = 0;
	  }
	  else {
	    if (gentau_decay[0]>7) {
	      indexGenTau1 = 0;
	      indexGenTau  = 1;
	    }
	    if (gentau_decay[1]>7) {
	      indexGenTau1 = 1;
	      indexGenTau = 0;
	    }
	    if (gentau_decay[0]<8&&gentau_decay[1]<8) {
	      float recopt1 = -1;
	      float recopt  = -1;
	      float dRmin1 = 0.2;
	      float dRmin  = 0.2;
	      for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) {
		if (!LooseTauSelection(analysisTree,itau)) continue;
		float dR = deltaR(gentau_eta[0],gentau_phi[0],
				  analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
		if (dR<dRmin) {
		  dRmin = dR;
		  recopt = analysisTree.tau_pt[itau];
		}
		float dR1 = deltaR(gentau_eta[1],gentau_phi[1],
				   analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
		if (dR1<dRmin1) {
		  dRmin1 = dR1;
		  recopt1 = analysisTree.tau_pt[itau];
		}
	      }
	      if (recopt1>recopt) {
		indexGenTau1 = 1;
		indexGenTau = 0;
	      }
	      else {
		indexGenTau1 = 0;
                indexGenTau = 1;
	      }
	    }
	  }

	  // *****************************************************
	  // *******  tau (trailing or single) *******************
	  // *****************************************************
	  unsigned int igentau = gentau_index[indexGenTau];
	  
	  TLorentzVector GenVisTau; GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
						      analysisTree.gentau_visible_py[igentau],
						      analysisTree.gentau_visible_pz[igentau],
						      analysisTree.gentau_visible_e[igentau]);
	  TLorentzVector GenTau; GenTau.SetXYZT(analysisTree.gentau_px[igentau],
						analysisTree.gentau_py[igentau],
						analysisTree.gentau_pz[igentau],
						analysisTree.gentau_e[igentau]);
	    
	  //	  if (GenVisTau.Pt()<ptTauCut) continue;
	  //	  if (fabs(GenVisTau.Eta())>etaTauCut) continue;
	  genTauPt_ = GenTau.Pt();
	  genTauEta_ = GenTau.Eta();
	  genTauPhi_ = GenTau.Phi();
	  genTauVisPt_ = GenVisTau.Pt();
	  genTauVisEta_ = GenVisTau.Eta();
	  genTauVisPhi_ = GenVisTau.Phi();
	  genTauVisMass_ = GenVisTau.M();
	  genTauDecay_ = analysisTree.gentau_decayMode[igentau];
	  TLorentzVector leadTrackLV = leadingTrackLV(&analysisTree,igentau);
	  genTauLeadingTrackPt_ = leadTrackLV.Pt();
	  genTauLeadingTrackEta_ = leadTrackLV.Eta();
	  genTauLeadingTrackPhi_ = leadTrackLV.Phi();
	  float dRmin = 0.2;
	  unsigned int indexTau = 0;
	  genTauFoundReco_ = false;
	  
	  for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) { // loop over taus
	    if (!LooseTauSelection(analysisTree,itau)) continue;
	    if (analysisTree.tau_pt[itau]<ptTauCut) continue;
	    if (fabs(analysisTree.tau_eta[itau])>etaTauCut) continue;
	    float dR = deltaR(GenVisTau.Eta(),GenVisTau.Phi(),
			      analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
	    if (dR<dRmin) {
	      dRmin = dR;
	      genTauFoundReco_ = true;
	      indexTau = itau;
	    }
	  }
	  
	  if (genTauFoundReco_) {
	    tauPt_ = analysisTree.tau_pt[indexTau];
	    tauPz_ = analysisTree.tau_pz[indexTau];
	    tauEta_ = analysisTree.tau_eta[indexTau];
	    tauPhi_ = analysisTree.tau_phi[indexTau];
	    tauMass_ = analysisTree.tau_mass[indexTau];
	    tauQ_ = int(analysisTree.tau_charge[indexTau]);
	    
	    TLorentzVector tauLeadLV;
	    tauLeadLV.SetXYZM(analysisTree.tau_leadchargedhadrcand_px[indexTau],
			      analysisTree.tau_leadchargedhadrcand_py[indexTau],
			      analysisTree.tau_leadchargedhadrcand_pz[indexTau],
			      analysisTree.tau_leadchargedhadrcand_mass[indexTau]);

	    tauLeadingTrackPt_ = tauLeadLV.Pt();
	    tauLeadingTrackEta_ = tauLeadLV.Eta();
	    tauLeadingTrackPhi_ = tauLeadLV.Phi();
	    tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];
	    tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dz[indexTau];

	    TLorentzVector lorentzVectorTau; 
	    lorentzVectorTau.SetXYZT(analysisTree.tau_px[indexTau],
				     analysisTree.tau_py[indexTau],
				     analysisTree.tau_pz[indexTau],
				     analysisTree.tau_e[indexTau]);
	    
	    recoilDPhi_  = dPhiFromLV(lorentzVectorTau,lorentzVectorMet);
	    
	    tauDecay_ = analysisTree.tau_decayMode[indexTau];
	    
	    tauGenMatch_ = analysisTree.tau_genmatch[indexTau];
	    
	    tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
	    tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;
	    
	    taubyDeepTau2017v2p1VSeraw_ = analysisTree.tau_byDeepTau2017v2p1VSeraw[indexTau] > 0.5;
	    taubyDeepTau2017v2p1VSjetraw_ = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexTau] > 0.5;
	    taubyDeepTau2017v2p1VSmuraw_ = analysisTree.tau_byDeepTau2017v2p1VSmuraw[indexTau] > 0.5;
	    taubyLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexTau] > 0.5;
	    taubyMediumDeepTau2017v2p1VSe_ = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyMediumDeepTau2017v2p1VSjet_ = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyMediumDeepTau2017v2p1VSmu_ = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexTau] > 0.5;
	    taubyTightDeepTau2017v2p1VSe_ = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyTightDeepTau2017v2p1VSmu_ = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTau] > 0.5;
	    taubyVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyVLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexTau] > 0.5;
	    taubyVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyVVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyVVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    taubyVVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
	    taubyVVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
	    
	    tauSinglePFTau180Trk50_ = false;
	    tauSinglePFTau180Trk50_2_ = false;
	    tauSinglePFTau180Trk50oneprong_ = false;
	    tauSinglePFTau180Trk50oneprong_2_ = false;
	    tauDoubleTauTrigger1_ = false;
	    tauDoubleTauTrigger2_ = false;
	    tauDoubleTauTrigger3_ = false;
	    tauDoubleTauTrigger4_ = false;
	    tauMuTauTrigger1_ = false;
	    tauMuTauTrigger2_ = false;
	    tauMuTauTrigger3_ = false;
	    tauMuTauTrigger4_ = false;

	    for (unsigned int iT=0; iT<analysisTree.trigobject_count;++iT) {
	      double dR = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau], 
				 analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dR>0.5) continue;

	      if (isSinglePFTau180Trk50Filter) {
		if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50Filter]) { 
		  tauSinglePFTau180Trk50_ = true;
		  indexSingleTau.push_back(iT);
		}
	      }
	      if (isSinglePFTau180Trk50Filter2) {
		if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50Filter2]) tauSinglePFTau180Trk50_2_ = true;
	      }
	      if (isSinglePFTau180Trk50oneprongFilter) {
		if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50oneprongFilter]) tauSinglePFTau180Trk50oneprong_ = true;
	      }
	      if (isSinglePFTau180Trk50oneprongFilter2) {
		if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50oneprongFilter2]) tauSinglePFTau180Trk50oneprong_2_ = true;
	      }

	      if (isDoubleTauFilter1) {
		if (analysisTree.trigobject_filters[iT][nDoubleTauFilter1]) tauDoubleTauTrigger1_ = true;
	      }
	      if (isDoubleTauFilter2) {
		if (analysisTree.trigobject_filters[iT][nDoubleTauFilter2]) tauDoubleTauTrigger2_ = true;
	      }
	      if (isDoubleTauFilter3) {
		if (analysisTree.trigobject_filters[iT][nDoubleTauFilter3]) tauDoubleTauTrigger3_ = true;
	      }
	      if (isDoubleTauFilter4) {
		if (analysisTree.trigobject_filters[iT][nDoubleTauFilter4]) tauDoubleTauTrigger4_ = true;
	      }

	      if (analysisTree.trigobject_pt[iT]>ptTauTriggerCut) {
		if (isMuTauFilter1) {
		  if (analysisTree.trigobject_filters[iT][nMuTauFilter1]) tauMuTauTrigger1_ = true;
		}
		if (isMuTauFilter2) {
		  if (analysisTree.trigobject_filters[iT][nMuTauFilter2]) tauMuTauTrigger2_ = true;
		}
		if (isMuTauFilter3) {
		  if (analysisTree.trigobject_filters[iT][nMuTauFilter3]) tauMuTauTrigger3_ = true;
		}
		if (isMuTauFilter4) {
		  if (analysisTree.trigobject_filters[iT][nMuTauFilter4]) tauMuTauTrigger4_ = true;
		}
	      }
	      
	    }
	    tauMatchedL1Tau80_ = false;
	    tauMatchedL1Tau140_ = false;
	    for (unsigned int iT=0; iT<analysisTree.l1tau_count; ++iT) {
	      
		if (analysisTree.l1tau_bx[iT]!=0) continue;
		TLorentzVector l1LV; l1LV.SetXYZM(analysisTree.l1tau_px[iT],
						  analysisTree.l1tau_py[iT],
						  analysisTree.l1tau_pz[iT],
						  pionMass);
		if (fabs(l1LV.Eta())>2.3) continue;
		float dR = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau],
				  l1LV.Eta(),l1LV.Phi());
		
		if (dR>0.5) continue;
		if (l1LV.Pt()>80.) tauMatchedL1Tau80_ = true;
		if (l1LV.Pt()>140.) tauMatchedL1Tau140_ = true;
	    }


	    /*
	    if (tauSinglePFTau180Trk50_) {
	      std::cout << "SinglePFTau180Trk50 = " << tauSinglePFTau180Trk50_ 
			<< "  SinglePFTau180Trk50L1Matched = " << tauSinglePFTau180Trk50_2_
			<< "  tauMatchedL1Tau80_" << tauMatchedL1Tau80_
			<< "  tauMatchedL1Tau140_ " << tauMatchedL1Tau140_
			<< std::endl;

	    }
	    */

	    
	  }
	  // *****************************************************
	  // *******  end tau (trailing or single) ***************
	  // *****************************************************
	  //	  std::cout << std::endl;

	  // *****************************************************
	  // ****** tau (leading) or lepton **********************
	  // *****************************************************
	  if (ngentaus==2) {
	    unsigned int igentau = gentau_index[indexGenTau1];
	    
	    TLorentzVector GenVisTau; GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
							analysisTree.gentau_visible_py[igentau],
							analysisTree.gentau_visible_pz[igentau],
							analysisTree.gentau_visible_e[igentau]);
	    TLorentzVector GenTau; GenTau.SetXYZT(analysisTree.gentau_px[igentau],
						  analysisTree.gentau_py[igentau],
						  analysisTree.gentau_pz[igentau],
						  analysisTree.gentau_e[igentau]);
	    
	    genTau1Pt_ = GenTau.Pt();
	    genTau1Eta_ = GenTau.Eta();
	    genTau1Phi_ = GenTau.Phi();
	    genTau1VisPt_ = GenVisTau.Pt();
	    genTau1VisEta_ = GenVisTau.Eta();
	    genTau1VisPhi_ = GenVisTau.Phi();
	    genTau1VisMass_ = GenVisTau.M();
	    genTau1Decay_ = analysisTree.gentau_decayMode[igentau];
	    TLorentzVector lead1TrackLV = leadingTrackLV(&analysisTree,igentau);
	    genTau1LeadingTrackPt_ = lead1TrackLV.Pt();
	    genTau1LeadingTrackEta_ = lead1TrackLV.Eta();
	    genTau1LeadingTrackPhi_ = lead1TrackLV.Phi();

	    /*
	    std::cout << "Tau decay  = " << genTau1Decay_ << std::endl;
	    std::cout << "TauVis pT  = " << genTau1VisPt_ << std::endl;
	    std::cout << "LeadTrk pT = " << genTau1LeadingTrackPt_ << std::endl;
	    */

	    // *************
	    // Tau->Mu Decay
	    // *************
	    if (analysisTree.gentau_decayMode[igentau]==8) {
	      genTau1FoundReco_ = false;
	      float dRmin = 0.2;
	      unsigned int index = 0;
	      float isolation = -1;
	      for (unsigned int im=0; im<muons.size(); ++im) { // loop over muons
		unsigned int indexMu = muons.at(im);
		float dR = deltaR(GenVisTau.Eta(),GenVisTau.Phi(),
				  analysisTree.muon_eta[indexMu],analysisTree.muon_phi[indexMu]);
		if (dR<dRmin) {
		  dRmin = dR;
		  genTau1FoundReco_ = true;
		  index = indexMu;
		  isolation = isoVarMuons.at(im);
		}
	      }
	      if (genTau1FoundReco_) {
		tau1Pt_ = analysisTree.muon_pt[index];
		tau1Pz_ = analysisTree.muon_pz[index];
		tau1Eta_ = analysisTree.muon_eta[index];
		tau1Phi_ = analysisTree.muon_phi[index];
		tau1Mass_ = muonMass;
		tau1Q_ = int(analysisTree.muon_charge[index]);
		tau1Iso_ = isolation;
	     
		tau1SingleMuon_ = false;
		for (unsigned int iT=0; iT<analysisTree.trigobject_count;++iT) {
		  double dR = deltaR(analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT], 
				     tau1Eta_,tau1Phi_);
		  if (dR>0.5) continue;
		  if (isSingleMuonHLTFilter) {
		    if (analysisTree.trigobject_filters[iT][nSingleMuonHLTFilter]) tau1SingleMuon_ = true;
		  }
		  if (isSingleMuonHLTFilter1) {
		    if (analysisTree.trigobject_filters[iT][nSingleMuonHLTFilter1]) tau1SingleMuon_ = true;
		  }

		}
	      }
	    }
	    // ****************
	    // Tau->E Decay
	    // ****************
	    if (analysisTree.gentau_decayMode[igentau]==9) {
	      genTau1FoundReco_ = false;
	      tau1SingleMuon_ = false;
	      float dRmin = 0.2;
	      unsigned int index = 0;
	      float isolation = -1;
	      for (unsigned int ie=0; ie<electrons.size(); ++ie) { // loop over electrons
		unsigned int indexEle = electrons.at(ie);
		float dR = deltaR(GenVisTau.Eta(),GenVisTau.Phi(),
				  analysisTree.electron_eta[indexEle],analysisTree.electron_phi[indexEle]);
		if (dR<dRmin) {
		  dRmin = dR;
		  genTau1FoundReco_ = true;
		  index = indexEle;
		  isolation = isoVarElectrons.at(ie);
		}
	      }
	      if (genTau1FoundReco_) {
		tau1Pt_ = analysisTree.electron_pt[index];
		tau1Pz_ = analysisTree.electron_pz[index];
		tau1Eta_ = analysisTree.electron_eta[index];
		tau1Phi_ = analysisTree.electron_phi[index];
		tau1Mass_ = electronMass;
		tau1Q_ = int(analysisTree.electron_charge[index]);
		tau1Iso_ = isolation;
	      }
	    }

	    // ****************
	    // Tau->hadrons
	    // ****************
	    if (analysisTree.gentau_decayMode[igentau]<8) {
	      float dRmin = 0.2;
	      unsigned int indexTau = 0;
	      genTau1FoundReco_ = false;
	  
	      for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) { // loop over taus
		if (!LooseTauSelection(analysisTree,itau)) continue;
		if (analysisTree.tau_pt[itau]<ptTauCut) continue;
		if (abs(analysisTree.tau_eta[itau])>etaTauCut) continue;
		float dR = deltaR(GenVisTau.Eta(),GenVisTau.Phi(),
				  analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
		if (dR<dRmin) {
		  dRmin = dR;
		  genTau1FoundReco_ = true;
		  indexTau = itau;
		}
	      }
	    

	      if (genTau1FoundReco_) {
		tau1Pt_ = analysisTree.tau_pt[indexTau];
		tau1Pz_ = analysisTree.tau_pz[indexTau];
		tau1Eta_ = analysisTree.tau_eta[indexTau];
		tau1Phi_ = analysisTree.tau_phi[indexTau];
		tau1Mass_ = analysisTree.tau_mass[indexTau];
		tau1Q_ = int(analysisTree.tau_charge[indexTau]);
		
		TLorentzVector tauLeadLV;
		tauLeadLV.SetXYZM(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				  analysisTree.tau_leadchargedhadrcand_py[indexTau],
				  analysisTree.tau_leadchargedhadrcand_pz[indexTau],
				  analysisTree.tau_leadchargedhadrcand_mass[indexTau]);
		
		tau1LeadingTrackPt_ = tauLeadLV.Pt();
		tau1LeadingTrackEta_ = tauLeadLV.Eta();
		tau1LeadingTrackPhi_ = tauLeadLV.Phi();
		tau1LeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];
		tau1LeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dz[indexTau];

		tau1Decay_ = analysisTree.tau_decayMode[indexTau];
		
		tau1DM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
		tau1NewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;
		
		tau1byDeepTau2017v2p1VSeraw_ = analysisTree.tau_byDeepTau2017v2p1VSeraw[indexTau] > 0.5;
		tau1byDeepTau2017v2p1VSjetraw_ = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexTau] > 0.5;
		tau1byDeepTau2017v2p1VSmuraw_ = analysisTree.tau_byDeepTau2017v2p1VSmuraw[indexTau] > 0.5;
		tau1byLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexTau] > 0.5;
		tau1byMediumDeepTau2017v2p1VSe_ = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byMediumDeepTau2017v2p1VSjet_ = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byMediumDeepTau2017v2p1VSmu_ = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexTau] > 0.5;
		tau1byTightDeepTau2017v2p1VSe_ = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byTightDeepTau2017v2p1VSmu_ = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTau] > 0.5;
		tau1byVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byVLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexTau] > 0.5;
		tau1byVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byVVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byVVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexTau] > 0.5;
		tau1byVVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexTau] > 0.5;
		tau1byVVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexTau] > 0.5;
		
		tau1SinglePFTau180Trk50_ = false;
		tau1SinglePFTau180Trk50_2_ = false;
		tau1SinglePFTau180Trk50oneprong_ = false;
		tau1SinglePFTau180Trk50oneprong_2_ = false;
		tau1DoubleTauTrigger1_ = false;
		tau1DoubleTauTrigger2_ = false;
		tau1DoubleTauTrigger3_ = false;
		tau1DoubleTauTrigger4_ = false;
		
		for (unsigned int iT=0; iT<analysisTree.trigobject_count;++iT) {
		  double dR = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau], 
				     analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
		  if (dR>0.5) continue;
		  if (isSinglePFTau180Trk50Filter) {
		    if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50Filter]) tau1SinglePFTau180Trk50_ = true;
		  }
		  if (isSinglePFTau180Trk50Filter2) {
		    if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50Filter2]) tau1SinglePFTau180Trk50_2_ = true;
		  }
		  if (isSinglePFTau180Trk50oneprongFilter) {
		    if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50oneprongFilter]) tau1SinglePFTau180Trk50oneprong_ = true;
		  }
		  if (isSinglePFTau180Trk50oneprongFilter2) {
		    if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50oneprongFilter2]) tau1SinglePFTau180Trk50oneprong_2_ = true;
		  }
		  if (isDoubleTauFilter1) {
		    if (analysisTree.trigobject_filters[iT][nDoubleTauFilter1]) tau1DoubleTauTrigger1_ = true;
		  }
		  if (isDoubleTauFilter2) {
		    if (analysisTree.trigobject_filters[iT][nDoubleTauFilter2]) tau1DoubleTauTrigger2_ = true;
		  }
		  if (isDoubleTauFilter3) {
		    if (analysisTree.trigobject_filters[iT][nDoubleTauFilter3]) tau1DoubleTauTrigger3_ = true;
		  }
		  if (isDoubleTauFilter4) {
		    if (analysisTree.trigobject_filters[iT][nDoubleTauFilter4]) tau1DoubleTauTrigger4_ = true;
		  }
		}
		tau1MatchedL1Tau80_ = false;
		tau1MatchedL1Tau140_ = false;
		for (unsigned int iT=0; iT<analysisTree.l1tau_count; ++iT) {
	      
		  if (analysisTree.l1tau_bx[iT]!=0) continue;
		  TLorentzVector l1LV; l1LV.SetXYZM(analysisTree.l1tau_px[iT],
						    analysisTree.l1tau_py[iT],
						    analysisTree.l1tau_pz[iT],
						    pionMass);
		  if (fabs(l1LV.Eta())>2.3) continue;
		  float dR = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau],
				    l1LV.Eta(),l1LV.Phi());
		  
		  if (dR>0.5) continue;
		  if (l1LV.Pt()>80.) tau1MatchedL1Tau80_ = true;
		  if (l1LV.Pt()>140.) tau1MatchedL1Tau140_ = true;
		}
	      }
	    }	    
	  }
	}

	// ************************
	// **** accessing jets ****
	// ************************
	
	TLorentzVector lorentzVectorAllJetsForMht; lorentzVectorAllJetsForMht.SetXYZT(0,0,0,0);
	std::vector<unsigned int> centralJets20Indexes; centralJets20Indexes.clear();
	std::vector<unsigned int> forwardJets20Indexes; forwardJets20Indexes.clear();
	std::vector<unsigned int> centralJets30Indexes; centralJets30Indexes.clear();
	std::vector<unsigned int> forwardJets30Indexes; forwardJets30Indexes.clear();

	int indexTauJet = -1;
	int indexTau1Jet = -1;
	float dRMinTau = 0.4;
	float dRMinTau1 = 0.4;
	for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
	    
	  // ----------------------- clean for EEnoise jets ----------------------- 
	  if( era == "2017" && analysisTree.pfjet_pt[ijet] < 50 && fabs(analysisTree.pfjet_eta[ijet]) < 3.139 && fabs(analysisTree.pfjet_eta[ijet]) > 2.65) continue;
          

	  float dr = deltaR(tauEta_,
			    tauPhi_,
			    analysisTree.pfjet_eta[ijet],
			    analysisTree.pfjet_phi[ijet]);

	  if (dr<dRMinTau) {
	    dRMinTau = dr;
	    indexTauJet = ijet;
	  }

	  dr = deltaR(tau1Eta_,
		      tau1Phi_,
		      analysisTree.pfjet_eta[ijet],
		      analysisTree.pfjet_phi[ijet]);

	  if (dr<dRMinTau1) {
	    dRMinTau1 = dr;
	    indexTau1Jet = ijet;
	  }


	  float absJetEta = fabs(analysisTree.pfjet_eta[ijet]);
	  if (absJetEta>5.2) continue;
	  if (analysisTree.pfjet_pt[ijet]<20.0) continue; 
          
	  // ----------------------- jet ID ----------------------- 
	  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
	  bool isPFJetId = false;
	  if (era == "2017") isPFJetId = tightJetiD_2017(analysisTree,int(ijet));
	  else if (era == "2018") isPFJetId = tightJetiD_2018(analysisTree,int(ijet));
	  else if (era == "2016") isPFJetId = looseJetiD_2016(analysisTree,int(ijet)); // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	  else {
	    std::cout<<"No Jet Id specified for era "<<era<<std::endl;
	    exit(-1);
	  }
	  if (!isPFJetId) continue;
          
	  // jet four-vector
	  TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[ijet],
					      analysisTree.pfjet_py[ijet],
					      analysisTree.pfjet_pz[ijet],
					      analysisTree.pfjet_e[ijet]);
	  
	  // counting jets for Mht
	  lorentzVectorAllJetsForMht += jetLV;
          
	  bool overlap = false;
	  if (genTauFoundReco_) {
	    float dR = deltaR(tauEta_,tauPhi_,jetLV.Eta(),jetLV.Phi());
	    if (dR<0.5) overlap = true;
	  }
	  if (genTau1FoundReco_) {
	    float dR = deltaR(tau1Eta_,tau1Phi_,jetLV.Eta(),jetLV.Phi());
	    if (dR<0.5) overlap = true;
	  }
	  if (overlap) continue; 

           // pt > 20 GeV
           if (analysisTree.pfjet_pt[ijet]>20.0) {
              if (absJetEta<2.4) {
                 centralJets20Indexes.push_back(ijet);
              }
              else if (absJetEta<4.7) {
                 forwardJets20Indexes.push_back(ijet);
              }
           }
           
           // pt > 30 GeV
           if (analysisTree.pfjet_pt[ijet]>30) {
	     if (absJetEta<2.4) {
		centralJets30Indexes.push_back(ijet);
		for (unsigned int iT=0; iT<indexSingleTau.size(); ++iT) {
		  unsigned int iTr = indexSingleTau.at(iT);
		  float dRjet = deltaR(jetLV.Eta(),jetLV.Phi(),
				       analysisTree.trigobject_eta[iTr],analysisTree.trigobject_phi[iTr]);
		  if (dRjet<0.5) {
		    float dRtau = deltaR(tauEta_,tauPhi_,
					 analysisTree.trigobject_eta[iTr],analysisTree.trigobject_phi[iTr]);
		    float dRtaujet = deltaR(tauEta_,tauPhi_,jetLV.Eta(),jetLV.Phi());
		    std::cout << "Trigger : eta " << analysisTree.trigobject_eta[iTr] << "   phi = " << analysisTree.trigobject_phi[iTr] << std::endl;
		    std::cout << "Tau     : eta " << tauEta_ << "   phi = " << tauPhi_ << std::endl;
		    std::cout << "Jet     : eta " << jetLV.Eta() << "   phi = " << jetLV.Phi() << std::endl;
		    std::cout << "dR(tau,trig) = " << dRtau << "   dR(jet,trig) = " << dRjet << "  dR(tau,jet) = " << dRtaujet << std::endl;
		  }
		  
		}

              }
              else if (absJetEta<4.7) {
		forwardJets30Indexes.push_back(ijet);
              }
           }
	}
        nJetsCentral20_ = centralJets20Indexes.size();
        nJetsCentral30_ = centralJets30Indexes.size();
        nJetsForward20_ = forwardJets20Indexes.size();
        nJetsForward30_ = forwardJets30Indexes.size();
        mht_        = (lorentzVectorAllJetsForMht).Pt();
        mhtNoMu_    = (lorentzVectorAllJetsForMht - lorentzVectorAllMuons).Pt();

	deltaR_ = 9999;
	if (genTauFoundReco_&&genTau1FoundReco_) 
	  deltaR_ = deltaR(tauEta_,tauPhi_,tau1Eta_,tau1Phi_);

	if (indexTauJet>=0) {
	  tauJetPt_ = analysisTree.pfjet_pt[indexTauJet];
	  tauJetEta_ = analysisTree.pfjet_eta[indexTauJet];
	  tauJetPhi_ = analysisTree.pfjet_phi[indexTauJet];
	}

	if (indexTau1Jet>=0) {
	  tau1JetPt_ = analysisTree.pfjet_pt[indexTau1Jet];
	  tau1JetEta_ = analysisTree.pfjet_eta[indexTau1Jet];
	  tau1JetPhi_ = analysisTree.pfjet_phi[indexTau1Jet];
	}

	vtx_trk = analysisTree.primvertexwithbs_ntracks;
	vtx_ex = TMath::Sqrt(analysisTree.primvertex_cov[0]);
	vtx_ey = TMath::Sqrt(analysisTree.primvertex_cov[3]);
	vtx_ez = TMath::Sqrt(analysisTree.primvertex_cov[5]);
	

        // ****************************
        // ****** trigger weight ******
        // ****************************
        float trigEffData = 1.0;
        float trigEffMC   = 1.0;
        
        trigWeight_ = 1;
        
        for (auto const& it : map_trigEffMC)
           {
              if( mhtNoMu_ < it.first)
                 { 
                    trigEffMC    =  it.second->Eval(metNoMu_);
                    trigEffData  =  map_trigEffData[it.first]->Eval(metNoMu_);
                    break;
                 }
           }
        
        if(trigEffMC !=0 ) trigWeight_ = trigEffData / trigEffMC;
        if(trigWeight_ < 0) trigWeight_=0;
        
        if (debug) {
           cout << "MetNoMu  = " << metNoMu_
                << "  MhtNoMu  = " << mhtNoMu_
                << "  trigWeight = " << trigWeight_ << endl;
        }
	//	std::cout << "tauMatchedL1Tau80 " << tauMatchedL1Tau80_ << std::endl;
	ntuple_->Fill();

     } // end of file processing (loop over events in one file)
     nFiles++;
     delete tree_;
     file_->Close();
     delete file_;
   }
   
   int allEvents   = int(inputEventsH->GetSumOfWeights());
   double sumOfWeights = histWeightsH->GetSumOfWeights();
   std::cout << "Total number of input events      = " << allEvents << std::endl;
   std::cout << "Total weight sum                  = " << sumOfWeights << std::endl;
   std::cout << "Total number of events in Tree    = " << nEvents << std::endl;
   std::cout << "Total number of trigger events    = " << TrigEvents << std::endl;
   std::cout << "Total number of W+Jet events      = " << WJetEvents << std::endl;
   std::cout << "Total number of W->muv events     = " << WMuNuEvents << std::endl;
   std::cout << "Total number of W->tauv events    = " << WTauNuEvents << std::endl;
   std::cout << "Total number of single jet events = " << SingleJetEvents << std::endl;
   std::cout << "Total number of jet+tau events    = " << JetTauEvents << std::endl;
   std::cout << "Total number of dijet events      = " << DiJetEvents << std::endl;
   file->cd("");
   file->Write();
   file->Close();
   delete file;
   
}



