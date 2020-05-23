#include <map>
#include <string>
#include <vector>
void TestTree_etau(int numberOfEvents = 10) {

  std::vector<TString> list_of_filters = {
    "hltEle27WPTightGsfTrackIsoFilter",
    "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter",
    "hltEGL1SingleEGOrFilter",
    "hltEle24erWPTightGsfTrackIsoFilterForTau",
    "hltSelectedPFTau30LooseChargedIsolationL1HLTMatched",
    "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30",
    "hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3"
  };

  // Data
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2017/data/SingleElectron_Run2017F-31Mar2018-v1/SingleElectron_Run2017F-31Mar2018-v1_1680.root");

  // MC
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2017/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_4119.root");

  // Embedding
  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2017/embedded/Embedding_eltau/EmbeddingRun2017F_ElTau/EmbeddingRun2017F_ElTau_2306.root");

  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");

  vector<string> * hltpaths = new vector<string>();
  vector<string> * hltfilters = new vector<string>();
  vector<string> * btagdiscriminators = new vector<string>();
  std::map<std::string, int>* hltriggerresults = new std::map<std::string, int>();
  std::map<std::string, int>* hltprescales = new std::map<std::string, int>();
  std::map<std::string, int>* flags = new std::map<std::string, int>();
  ULong64_t event_nr;

  tree->SetBranchAddress("event_nr",&event_nr);
  tree->SetBranchAddress("run_hltnames",&hltpaths);
  tree->SetBranchAddress("run_hltfilters",&hltfilters);
  tree->SetBranchAddress("run_btagdiscriminators",&btagdiscriminators);

  tree->SetBranchAddress("hltriggerresults",&hltriggerresults); 
  tree->SetBranchAddress("hltriggerprescales",&hltprescales); 
  tree->SetBranchAddress("flags",&flags);

  int nEntries = tree->GetEntries();

  std::cout << "Entries = " << nEntries << std::endl;

  int nEvents = 0;
  for (int iE=0; iE<numberOfEvents; ++iE) {

    nEvents++;

    tree->GetEntry(iE);

    std::cout << std::endl;
    std::cout << "event number = " << nEvents << std::endl;
    std::cout << std::endl;

    unsigned int nf = flags->size();
    std::cout << "nf = " << nf << std::endl;
    for (std::map<string,int>::iterator it=flags->begin(); it!=flags->end(); ++it)
      std::cout << it->first << "  :  "  << it->second << std::endl;
    std::cout << std::endl;

    unsigned int nhlt = hltpaths->size();
    std::cout << "nhlt = " << nhlt << std::endl;
    for (unsigned int i=0; i<nhlt; ++i)
      std::cout << "HLT Path : " << hltpaths->at(i) << std::endl;
    std::cout << std::endl;

    unsigned int ntrig = hltriggerresults->size();
    std::cout << "ntrig = " << ntrig << std::endl;
    for (std::map<string,int>::iterator it=hltriggerresults->begin(); it!=hltriggerresults->end(); ++it) 
      std::cout << it->first << "  :  "  << it->second << std::endl;
    std::cout << std::endl;

    unsigned int nprescales = hltprescales->size();
    std::cout << "npres = " << nprescales << std::endl;
    for (std::map<string,int>::iterator it=hltprescales->begin(); it!=hltprescales->end(); ++it) 
      std::cout << it->first << "  :  "  << it->second << std::endl;
    std::cout << std::endl;

    unsigned int nfilters = hltfilters->size();
    std::cout << "Selected filters from list_of_filters out of " << nfilters << " => " << std::endl;
    for (unsigned int i=0; i<nfilters; ++i) {
      TString HltFilter(hltfilters->at(i));
      for (unsigned int j=0; j<list_of_filters.size(); ++j) {
	if (HltFilter==list_of_filters.at(j)) {
	  std::cout << "HLT Filter : " << HltFilter << " i: " << i << std::endl;
	}
      }
    }
    std::cout << std::endl;

    unsigned int nbtag = btagdiscriminators->size();
    std::cout << "nbtags = " << nbtag << std::endl;
    for (unsigned int i=0; i<nbtag; ++i)
      std::cout << "BTag discriminator : " << btagdiscriminators->at(i) << std::endl;
    std::cout << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    
  }


}
