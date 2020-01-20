#include <map>
#include <string>
#include <vector>
void TestTree_mutau(int numberOfEvents = 1) {

  std::vector<TString> filtersToTest;

  //  2018 
  filtersToTest.push_back("hltL1sMu18erTau24erIorMu20erTau24er");
  filtersToTest.push_back("hltL1sBigORMu18erTauXXer2p1");
  filtersToTest.push_back("hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p07");
  //  2017
  //  filtersToTest.push_back("hltL1sMu18erTau24erIorMu20erTau24er");
  //  filtersToTest.push_back("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07");
  //  2016
  //  filtersToTest.push_back("hltL1sMu18erTau20er");

  // 2018
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau/EmbeddingRun2018A_MuTau/EmbeddingRun2018A_MuTau_1880.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau/EmbeddingRun2018B_MuTau/EmbeddingRun2018B_MuTau_802.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau/EmbeddingRun2018C_MuTau/EmbeddingRun2018C_MuTau_420.root");
  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau/EmbeddingRun2018D_MuTau/EmbeddingRun2018D_MuTau_1300.root");

  // 2017
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data_v2/Embedding_mutau//EmbeddingRun2017B_MuTau/EmbeddingRun2017B_MuTau_39.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data_v2/Embedding_mutau//EmbeddingRun2017C_MuTau/EmbeddingRun2017C_MuTau_190.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data_v2/Embedding_mutau//EmbeddingRun2017D_MuTau/EmbeddingRun2017D_MuTau_143.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data_v2/Embedding_mutau//EmbeddingRun2017E_MuTau/EmbeddingRun2017E_MuTau_364.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data_v2/Embedding_mutau//EmbeddingRun2017F_MuTau/EmbeddingRun2017F_MuTau_520.root");

  // 2016
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016B_MuTau/EmbeddingRun2016B_MuTau_609.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016C_MuTau/EmbeddingRun2016C_MuTau_799.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016D_MuTau/EmbeddingRun2016D_MuTau_31.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016E_MuTau/EmbeddingRun2016E_MuTau_933.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016F_MuTau/EmbeddingRun2016F_MuTau_1354.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016G_MuTau/EmbeddingRun2016G_MuTau_1941.root");
  //  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau/EmbeddingRun2016H_MuTau/EmbeddingRun2016H_MuTau_2614.root");

  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");

  vector<string> * hltpaths = new vector<string>();
  vector<string> * hltfilters = new vector<string>();
  vector<string> * btagdiscriminators = new vector<string>();
  std::map<std::string, int>* hltriggerresults = new std::map<std::string, int>();
  std::map<std::string, int>* hltprescales = new std::map<std::string, int>();
  std::map<std::string, int>* flags = new std::map<std::string, int>();
  ULong64_t event_nr;
  UInt_t event_run;

  tree->SetBranchAddress("event_nr",&event_nr);
  tree->SetBranchAddress("event_run",&event_run);
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
    std::cout << "event number = " << event_nr << std::endl;
    std::cout << "run number = " << event_run << std::endl;
    std::cout << std::endl;

    /*
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
    */
    unsigned int nfilters = hltfilters->size();
    std::cout << "nfiltres = " << nfilters << std::endl;
    for (unsigned int i=0; i<nfilters; ++i) {
      TString FILTER(hltfilters->at(i));
      for (unsigned int j=0; j<filtersToTest.size(); ++j) {
	if (FILTER==filtersToTest.at(j))
	  std::cout << "HLT Filter : " << hltfilters->at(i) << " i: " << i << std::endl;
      }
    }
    std::cout << std::endl;
    /*
    unsigned int nbtag = btagdiscriminators->size();
    std::cout << "nbtags = " << nbtag << std::endl;
    for (unsigned int i=0; i<nbtag; ++i)
      std::cout << "BTag discriminator : " << btagdiscriminators->at(i) << std::endl;
    std::cout << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    */
  }


}
