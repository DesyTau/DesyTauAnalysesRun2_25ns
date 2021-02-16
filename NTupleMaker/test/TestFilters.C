#include <map>
#include <string>
#include <vector>
void TestFilters(int numberOfEvents = 10) {


  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016B-17Jul2018_ver2-v1/MuonEG_Run2016B-17Jul2018_ver2-v1_699.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016C-17Jul2018-v1/MuonEG_Run2016C-17Jul2018-v1_99.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016D-17Jul2018-v1/MuonEG_Run2016D-17Jul2018-v1_407.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016E-17Jul2018-v2/MuonEG_Run2016E-17Jul2018-v2_841.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016F-17Jul2018-v1/MuonEG_Run2016F-17Jul2018-v1_979.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016G-17Jul2018-v1/MuonEG_Run2016G-17Jul2018-v1_909.root");
  // TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG/MuonEG_Run2016H-17Jul2018-v1/MuonEG_Run2016H-17Jul2018-v1_998.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/MuonEG/MuonEG_Run2018D-PromptReco-v2/MuonEG_Run2018D-PromptReco-v2_437.root");

  //  TFile * file = new TFile("output_MC.root");
  //  TFile * file = new TFile("output_DATA.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc//WToTauNu_M-200_TuneCP5_13TeV-pythia8-tauola/WToTauNu_M-200_TuneCP5_13TeV-pythia8-tauola_173.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/mc_2/SUSYGluGluToHToTauTau_M-1600/SUSYGluGluToHToTauTau_M-1600_148.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc//WToTauNu_M-200_TuneCP5_13TeV-pythia8-tauola/WToTauNu_M-200_TuneCP5_13TeV-pythia8-tauola_305.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2017/embedded/EmbeddedEmu/EmbeddingRun2017B/EmbeddingRun2017B_105.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2016/emb/EmbeddingRun2016H_ElTau/EmbeddingRun2016H_ElTau_2354.root");
  // e-mu embedded ->
  // 2016
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/EmbeddedEmu/EmbeddingRun2016H/EmbeddingRun2016H_2717.root");
  // 2017
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2017/embedded/EmbeddedEmu/EmbeddingRun2017F/EmbeddingRun2017F_3385.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2017/emb/EmbeddingRun2017F_ElTau/EmbeddingRun2017F_ElTau_3283.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2017/mc/SUSYGluGluToHToTauTau_M-1200/SUSYGluGluToHToTauTau_M-1200_323.root");
  // 2018
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/embedded/Embedding_emu_v2//EmbeddingRun2018C/ElMuFinalState-v1/EmbeddingRun2018C/ElMuFinalState-v1_97.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/emb/EmbeddingRun2018B_ElTau/EmbeddingRun2018B_ElTau_276.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/mc/SUSYGluGluToHToTauTau_M-1200/SUSYGluGluToHToTauTau_M-1200_470.root");

  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/data/MuonEG_Run2017B/MuonEG_Run2017B_1.root");
  
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/data/MuonEG_Run2017B/MuonEG_Run2017B_1.root");
  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/data/MuonEG_Run2017C/MuonEG_Run2017C_1107.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/data/MuonEG_Run2017D/MuonEG_Run2017D_152.root");

  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");

  vector<string> * hltfilters = new vector<string>();
  ULong64_t event_nr;

  vector<TString> filterList = {
    "hltL3fL1sMu7EG23f0Filtered8",
    "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",
    /*
    "hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17",
    "hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8",
    "hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8",
    "hltDiMuonGlb17Glb8DzFiltered0p2",
    "hltDiMuonGlb17Glb8DzFiltered0p2SameSign"
    */
    // DZ filters
    //    "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter",
    //    "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter",
    // 2016 single-tau
    /*
    "hltPFTau120TrackPt50LooseAbsOrRelVLooseIso",
    "hltPFTau140TrackPt50LooseAbsOrRelVLooseIso",
    "hltL3crIsoL1sMu18erIsoTau26erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
    "hltOverlapFilterIsoMu19MediumIsoPFTau32Reg",
    "hltOverlapFilterIsoMu19MediumCombinedIsoPFTau32Reg",
    "hltPFTau32Reg"
    */
    // 2016 single leptons
    /*
    "hltEle25erWPTightGsfTrackIsoFilter",
    "hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09",
    "hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09",
    "hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09",
    "hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09",
    */
    // 2017 single leptons
    /*      
    "hltEle32WPTightGsfTrackIsoFilter",
    "hltEle27WPTightGsfTrackIsoFilter",
    "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter",
    "hltEGL1SingleEGOrFilter",
    "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
    "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07",
    */
    // 2018 
    /*    
    "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
    "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07",
    "hltEle32WPTightGsfTrackIsoFilter",
    "hltEle35noerWPTightGsfTrackIsoFilter",
    "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
    "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"
    //    "hltEle32WPTightGsfTrackIsoFilter",
    /*
    "hltEle35noerWPTightGsfTrackIsoFilter",
    "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", // IsoMu27
    "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09",
    "hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso",
    "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched",
    "hltSingleL2Tau80eta2p2",
    "hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIso1Prong",
    "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong",
    "hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg",
    "hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg",
    "hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg",
    "hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg",
    "hltDoubleL2IsoTau26eta2p2",
    "hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg",
    "hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau40MonitoringReg",
    "hltOverlapFilterIsoMu24TightChargedIsoPFTau40MonitoringReg",
    "hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg",
    "hltSingleL2IsoTau26eta2p2",
    */
  }; 

  tree->SetBranchAddress("event_nr",&event_nr);
  tree->SetBranchAddress("run_hltfilters",&hltfilters);

  int nEntries = tree->GetEntries();

  std::cout << "Entries = " << nEntries << std::endl;

  int nEvents = 0;
  for (int iE=0; iE<numberOfEvents; ++iE) {

    nEvents++;

    tree->GetEntry(iE);

    std::cout << std::endl;
    std::cout << "event number = " << nEvents << std::endl;
    std::cout << std::endl;

    unsigned int nfilters = hltfilters->size();
    std::cout << "nfiltres = " << nfilters << std::endl;
    for (unsigned int ifilter=0; ifilter<filterList.size(); ++ifilter) {
      bool found = false;
      for (unsigned int i=0; i<nfilters; ++i) {
	//      std::cout << "HLT Filter : " << hltfilters->at(i) << " i: " << i << std::endl;
	TString Filter(hltfilters->at(i));
	if (Filter==filterList[ifilter]) {
	  std::cout << i << "  " << Filter << std::endl;
	  found = true;
	  break;
	}
      }
      if (!found) 
	std::cout << "Filter : " << filterList.at(ifilter) << " not found" << std::endl;
    }
    
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    
  }


}
