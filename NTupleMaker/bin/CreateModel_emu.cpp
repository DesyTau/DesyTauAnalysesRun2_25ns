#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"

int main(int argc, char * argv[]) {

  TString dir = "/nfs/dust/cms/user/rasp/Run/emu_MSSM";

  TString folder(argv[1]);
  TString era(argv[2]);
  TString inclusive(argv[3]);

  std::vector<TString> categories; 
  if (inclusive=="inclusive") {
    categories.push_back("ttbar");
    categories.push_back("inclusive_btag");
    categories.push_back("inclusive_nobtag");
  }
  else {
    categories.push_back("ttbar");
    categories.push_back("low_pzeta_btag");
    categories.push_back("low_pzeta_nobtag");
    categories.push_back("medium_pzeta_btag");
    categories.push_back("medium_pzeta_nobtag");
    categories.push_back("high_pzeta_btag");
    categories.push_back("high_pzeta_nobtag");
  }

  std::vector<TString> eras;
  if (era=="all") {
    eras.push_back("2016");
    eras.push_back("2017");
    eras.push_back("2018");
  }
  else {
    eras.push_back(era);
  }
  std::vector<TString> signals = {
    "ggH",
    "bbH"
  };
  std::vector<TString> bkgs = {
    "EmbedZTT",
    "TT",
    "ZLL",
    "VV",
    "W"
  };

  for (auto Era : eras) {
    for (auto cat : categories) {
      TString fileName = dir + "/" + folder + "/" + cat + "_" + Era + ".root";
      TFile * file = new TFile(fileName);
      if (file->IsZombie()) {
	std::cout << fileName << "  cannot be opened" << std::endl;
	exit(-1);
      }

      TH1D * data = (TH1D*)file->Get("data_obs");
      double data_yield = data->GetSumOfWeights();

      TH1D * EmbedZTT = (TH1D*)file->Get("EmbedZTT");
      double embed_yield = EmbedZTT->GetSumOfWeights();

      TH1D * TT = (TH1D*)file->Get("TT");
      double tt_yield = TT->GetSumOfWeights();

      TH1D * ZLL = (TH1D*)file->Get("ZLL");
      double zll_yield = ZLL->GetSumOfWeights();

      TH1D * W = (TH1D*)file->Get("W");
      double w_yield = W->GetSumOfWeights();

      TH1D * VV = (TH1D*)file->Get("VV");
      double vv_yield = VV->GetSumOfWeights();

      TH1D * QCD = (TH1D*)file->Get("QCD");
      double qcd_yield = QCD->GetSumOfWeights();

      TString channel = cat + "_" + Era;

      double totBkgd = embed_yield + tt_yield + zll_yield + w_yield + vv_yield + qcd_yield;

      std::cout << "making cards for channel " << channel << std::endl;
      std::cout << "EmbedZTT  : " << embed_yield << std::endl;
      std::cout << "TT        : " << tt_yield << std::endl;
      std::cout << "ZLL       : " << zll_yield << std::endl;
      std::cout << "W         : " << w_yield << std::endl;
      std::cout << "VV        : " << vv_yield << std::endl;
      std::cout << "QCD       : " << qcd_yield << std::endl;
      std::cout << "Total Bkg : " << totBkgd << std::endl;
      std::cout << "Data      : " << data_yield << std::endl;
      std::cout << "Data/Bkg  : " << data_yield/totBkgd << std::endl;
      std::cout << "signif    : " << (data_yield-totBkgd)/TMath::Sqrt(data_yield) << std::endl;
      std::cout << std::endl;

      for (auto mass : masses ) {
	for (auto sig : signals) {

	  TString sig_name = sig+mass;
	  TString BaseName = dir + "/" + folder + "/" + cat + "_" + Era + "_" + sig_name;

	  TH1D * SIG = (TH1D*)file->Get(sig_name);
	  double sig_yield = SIG->GetSumOfWeights();

	  ostringstream str;
	  str << BaseName << ".txt";
	  string nn = str.str();
	  const char * p = nn.c_str();

	  std::ofstream textFile(p,std::ofstream::out);
	  textFile << "imax 1   number of channels" << std::endl;
	  textFile << "jmax *   number of backgrounds" << std::endl;
	  textFile << "kmax *   number of nuisance parameters" << std::endl;
	  textFile << "-----------------" << std::endl;
	  textFile << "observation " << data_yield << std::endl;
	  textFile << "-----------------" << std::endl;
	  textFile << "shapes * * " << fileName << "  $PROCESS   $PROCESS_$SYSTEMATIC " << std::endl;
	  textFile << "-----------------" << std::endl;
	  textFile << "bin        "; 
	  for (unsigned int i=0; i<7; ++i)
	    textFile << "  " << channel;
	  textFile << std::endl;
	  textFile << "process      " << sig_name << "  EmbedZTT  TT  QCD  W  VV  ZLL" << std::endl;
	  textFile << "process                  0    1      2      3      4      5      6" << std::endl;
	  textFile << "rate     " 
		   << sig_yield <<  "  "
		   << embed_yield << "  " 
		   << tt_yield << "  "
		   << qcd_yield << "  " 
		   << w_yield << "  " 
		   << vv_yield << "  " 
		   << zll_yield << std::endl;
	  textFile << "-----------------------------" << std::endl;
	  textFile << "CMS_eff_m                   lnN   1.03   1.03   1.03      -   1.03   1.03   1.03" << std::endl;
	  textFile << "CMS_eff_e                   lnN   1.04   1.04   1.04      -   1.04   1.04   1.04" << std::endl;
	  textFile << "CMS_ztt_ttjNorm             lnN      -      -   1.07      -      -      -      -" << std::endl;
	  textFile << "CMS_ztt_wNorm               lnN      -      -      -      -   1.15      -      -" << std::endl;
	  textFile << "CMS_ztt_vvNorm              lnN      -      -      -      -      -   1.20      -" << std::endl;
	  textFile << "CMS_ztt_em_zllNorm          lnN      -      -      -      -      -      -   1.20" << std::endl;
	  textFile << "CMS_ztt_em_qcdNorm          lnN      -      -      -   1.15      -      -      -" << std::endl;
	  textFile << "lumi_13TeV_" << Era  << "                 lnN   1.03    1.03   1.03      -   1.03   1.03      -" << std::endl;
	  textFile << "CMS_ztt_boson_scale_met     lnN   1.04    1.02      -      -   1.02      -      -" << std::endl;
	  textFile << "CMS_ztt_boson_reso_met      lnN   1.04    1.02      -      -   1.02      -      -" << std::endl;
	  //	  textFile << "* autoMCStats 10 "<< std::endl;
	  textFile.close();


	}
      }
    }
  }

}
