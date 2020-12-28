#define M_jetmaxcount 1000

void TestBTag() {

  TFile * file = new TFile("/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/data/SingleMuon_Run2018D/SingleMuon_Run2018D_25499.root");

  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");
  vector<string> * btagdiscriminators = new vector<string>();
  unsigned int pfjet_count;
  float pfjet_btag[M_jetmaxcount][10];
  float pfjet_pt[M_jetmaxcount];
  float pfjet_eta[M_jetmaxcount];
  tree->SetBranchAddress("run_btagdiscriminators",&btagdiscriminators);
  tree->SetBranchAddress("pfjet_count", &pfjet_count);
  tree->SetBranchAddress("pfjet_btag", pfjet_btag);
  tree->SetBranchAddress("pfjet_eta", pfjet_eta);
  tree->SetBranchAddress("pfjet_pt", pfjet_pt);
  
  for (int iE = 0; iE <= 20; ++iE) {
    tree->GetEntry(iE);
    std::cout << "Event " << iE << std::endl;
    std::cout << "number of jets = " << pfjet_count << std::endl;
    unsigned int nbtag = btagdiscriminators->size();
    for (unsigned int jet=0; jet<pfjet_count; ++jet) {
      if (pfjet_pt[jet]<20) continue;
      if (TMath::Abs(pfjet_eta[jet])>2.4) continue;
      std::cout << "Jet pT = " << pfjet_pt[jet] << " Eta = " << pfjet_eta[jet] << std::endl;
      std::cout << "   discriminants -> " << std::endl;
      for (unsigned int j=0; j<nbtag; ++j)
	std:: cout << "   " << btagdiscriminators->at(j) << " = " << pfjet_btag[jet][j] << std::endl;

    }
    std::cout << std::endl;

  }


}
