#include <map>
#include <string>
#include <vector>
void TestTree(int numberOfEvents = 10) {

 // //  TFile * file = new TFile("/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/Sync2015/HiggsSM/GluGluToHToTauTau_M-125_MC_TauTau_v1/ntuple_GGF125_13TeV_Phys14_v1.root");
  TFile * file = new TFile("/nfs/dust/cms/user/rasp/ntuples/Run2015B/SingleElectron_Run2015B/SingleE.root");

  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");

  vector<string> * hltpaths = new vector<string>();
  vector<string> * hltfilters = new vector<string>();
  vector<string> * btagdiscriminators = new vector<string>();
  std::map<std::string, int>* hltriggerresults = new std::map<std::string, int>();
  std::map<std::string, int>* hltprescales = new std::map<std::string, int>();

  tree->SetBranchAddress("run_hltnames",&hltpaths);
  tree->SetBranchAddress("run_hltfilters",&hltfilters);
  tree->SetBranchAddress("run_btagdiscriminators",&btagdiscriminators);

  tree->SetBranchAddress("hltriggerresults",&hltriggerresults); 
  tree->SetBranchAddress("hltriggerprescales",&hltprescales); 

  int nEntries = tree->GetEntries();

  std::cout << "Entries = " << nEntries << std::endl;

  int nEvents = 0;
  for (int iE=0; iE<numberOfEvents; ++iE) {

    nEvents++;

    tree->GetEntry(iE);

    std::cout << std::endl;
    std::cout << "event number = " << nEvents << std::endl;
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
    std::cout << "nfiltres = " << nfilters << std::endl;
    for (unsigned int i=0; i<nfilters; ++i)
      std::cout << "HLT Filter : " << hltfilters->at(i) << std::endl;
    std::cout << std::endl;

    unsigned int nbtag = btagdiscriminators->size();
    std::cout << "nbtags = " << nbtag << std::endl;
    for (unsigned int i=0; i<nbtag; ++i)
      std::cout << "BTag discriminator : " << btagdiscriminators->at(i) << std::endl;
    std::cout << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    
  }


}
