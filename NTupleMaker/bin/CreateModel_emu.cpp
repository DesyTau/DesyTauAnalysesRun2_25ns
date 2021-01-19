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

std::vector<TString> sysToLnNames_Common = {
  "CMS_scale_j_FlavorQCD",
  "CMS_scale_j_RelativeBal",
  "CMS_scale_j_HF",
  "CMS_scale_j_BBEC1",
  "CMS_scale_j_EC2",
  "CMS_scale_j_Absolute",
};

std::vector<TString> sysToLnNames_Era = {
  "CMS_scale_met_unclustered",
  "CMS_htt_eff_b",
  "CMS_htt_mistag_b",
  "CMS_scale_j_Absolute",
  "CMS_scale_j_HF",
  "CMS_scale_j_EC2",
  "CMS_scale_j_RelativeSample",
  "CMS_scale_j_BBEC1",
  "CMS_res_j",
};

std::vector<TString> higgsMass = {
  //  "200",
  //  "400",
  //  "600",
  //  "800",
  "1200",
  //  "1400",
  //  "1600",
  //  "1800",
  //  "2000",
  //  "2600",
  //  "3200"
};

TString ExtractSys(TFile * file, 
		   TString cat, 
		   TString process, 
		   TString sysName) {
  
  TString output(" - ");
  TH1D * histNominal = (TH1D*)file->Get(cat+"/"+process);
  TH1D * histUp = (TH1D*)file->Get(cat+"/"+process+"_"+sysName+"Up");
  if (histUp!=NULL&&histNominal!=NULL) {
    double xNominal = histNominal->GetSumOfWeights();
    double xUp = histUp->GetSumOfWeights();
    double err = xUp/xNominal;
    if (err<0.999||err>1.001) {
      char number[10];
      sprintf(number,"%5.3f",err);
      output = " " + TString(number) + " ";
    }
  }

  return output;

}

int main(int argc, char * argv[]) {

  if (argc!=2) {
    std::cout << "Usage : CreateModel_emu [trigger=0,1,2]" << std::endl;
    exit(-1);
  }
  TString trigger(argv[1]);

  TString       dir = "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan10/datacards_"+trigger;
  TString dir_cards = "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan10/cards_"+trigger;

  std::vector<TString> categories = { 
    "em_ttbar_control",
    "em_ttbar_btag",
    "em_ttbar_nobtag",
    //    "em_btag_lowdzeta",
    //    "em_btag_mediumdzeta",
    //    "em_btag_highdzeta",
    //    "em_nobtag_lowdzeta",
    //    "em_nobtag_mediumdzeta",
    //    "em_nobtag_highdzeta"
  };

  std::vector<TString> eras = {
    "2016",
    "2017",
    "2018"
  };

  std::vector<TString> signals = {
    "ggH",
    "bbH"
  };

  for (auto Era : eras) {
    TString fileName = dir + "/" + Era + "/htt_em_mssm.root";
    TFile * file = new TFile(fileName);
    if (file->IsZombie()) {
      std::cout << fileName << "  cannot be opened" << std::endl;
      exit(-1);
    }
    for (auto cat : categories) {

      TH1D * data = (TH1D*)file->Get(cat+"/data_obs");
      TH1D * EMB = (TH1D*)file->Get(cat+"/EMB");
      TH1D * TTL = (TH1D*)file->Get(cat+"/TTL");
      TH1D * VVL = (TH1D*)file->Get(cat+"/VVL");
      TH1D * QCD = (TH1D*)file->Get(cat+"/QCD");
      TH1D * W   = (TH1D*)file->Get(cat+"/W");
      TH1D * ZLL = (TH1D*)file->Get(cat+"/ZLL");

      std::map<TString,TString> TTL_sys;
      std::map<TString,TString> VVL_sys;
      std::map<TString,TString> W_sys;
      std::map<TString,TString> ZLL_sys;

      double data_yield = data->GetSumOfWeights();
      double tt_yield = TTL->GetSumOfWeights();
      double zll_yield = ZLL->GetSumOfWeights();
      double w_yield = W->GetSumOfWeights();
      double vv_yield = VVL->GetSumOfWeights();
      double qcd_yield = QCD->GetSumOfWeights();
      double embed_yield = EMB->GetSumOfWeights();

      TString channel = cat + "_" + Era;

      double totBkgd = embed_yield + tt_yield + zll_yield + w_yield + vv_yield + qcd_yield;

      std::cout << "making cards for channel " << channel << std::endl;
      std::cout << "EMB       : " << embed_yield << std::endl;
      std::cout << "TT        : " << tt_yield << std::endl;
      std::cout << "ZLL       : " << zll_yield << std::endl;
      std::cout << "W         : " << w_yield << std::endl;
      std::cout << "VV        : " << vv_yield << std::endl;
      std::cout << "QCD       : " << qcd_yield << std::endl;
      std::cout << "Total Bkg : " << totBkgd << std::endl;
      std::cout << "Data      : " << data_yield << std::endl;
      std::cout << "Data/Bkg  : " << data_yield/totBkgd << std::endl;
      std::cout << std::endl;

      for (auto sysName : sysToLnNames_Common) {
	TTL_sys[sysName] = ExtractSys(file, cat, "TTL", sysName);
	VVL_sys[sysName] = ExtractSys(file, cat, "VVL", sysName);
	W_sys[sysName]   = ExtractSys(file, cat, "W", sysName);
	ZLL_sys[sysName] = ExtractSys(file, cat, "ZLL", sysName); 
      }

      for (auto sysNameEra : sysToLnNames_Era) {
	TString sysName = sysNameEra + "_" + Era;
	TTL_sys[sysName] = ExtractSys(file, cat, "TTL", sysName);
	VVL_sys[sysName] = ExtractSys(file, cat, "VVL", sysName);
	W_sys[sysName]   = ExtractSys(file, cat, "W", sysName);
	ZLL_sys[sysName] = ExtractSys(file, cat, "ZLL", sysName); 
      }

      if (Era=="2017"||Era=="2016") {	
	TString sysName = "CMS_Prefiring";
	TTL_sys[sysName] = ExtractSys(file, cat, "TTL", sysName);
	VVL_sys[sysName] = ExtractSys(file, cat, "VVL", sysName);
	W_sys[sysName]   = ExtractSys(file, cat, "W", sysName);
	ZLL_sys[sysName] = ExtractSys(file, cat, "ZLL", sysName);	
      }

      for (auto mass : higgsMass ) {
	for (auto sig : signals) {

	  TString sig_name = sig+mass;
	  TString BaseName = dir_cards + "/" + Era + "/" + cat + "_" + sig_name;

	  TH1D * SIG = (TH1D*)file->Get(cat+"/"+sig_name);
	  if (SIG==NULL) {
	    std::cout << "Signal " << sig_name << " is absent" << std::endl;
	    exit(-1);
	  }

	  double sig_yield = SIG->GetSumOfWeights();

	  std::map<TString,TString> SIG_sys;

	  for (auto sysName : sysToLnNames_Common) {
	    SIG_sys[sysName] = ExtractSys(file, cat, sig_name, sysName);
	  }
	  
	  for (auto sysNameEra : sysToLnNames_Era) {
	    TString sysName = sysNameEra + "_" + Era;
	    SIG_sys[sysName] = ExtractSys(file, cat, sig_name, sysName);
	  }

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
	  textFile << "shapes * * " << fileName << "  " << cat << "/$PROCESS " << cat << "/$PROCESS_$SYSTEMATIC " << std::endl;
	  textFile << "-----------------" << std::endl;
	  textFile << "bin        "; 
	  for (unsigned int i=0; i<7; ++i)
	    textFile << "  " << channel;
	  textFile << std::endl;
	  textFile << "process      " << sig_name << "  EMB  TTL  QCD  W  VVL  ZLL" << std::endl;
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
	  textFile << "CMS_htt_norm_QCD_" << Era << "             lnN      -     -     -     1.05  -  -  -  "<< std::endl;
	  textFile << "CMS_htt_norm_emb_" << Era << "             lnN      -  1.04     -       -      -      -      -" << std::endl;
	  textFile << "CMS_eff_m                        lnN   1.02     -   1.02      -   1.02   1.02   1.02" << std::endl;
	  textFile << "CMS_eff_e                        lnN   1.02     -   1.02      -   1.02   1.02   1.02" << std::endl;
	  textFile << "CMS_trigger_em_" << Era << "     lnN   1.02     -   1.02      -   1.02   1.02   1.02" << std::endl;
	  textFile << "CMS_trigger_emb_em_" << Era << "  lnN      -   1.02     -      -      -      -      -" << std::endl;
	  textFile << "CMS_eff_m_emb                    lnN      -   1.02      -      -     -     -   -" << std::endl;
	  textFile << "CMS_eff_e_emb                    lnN      -   1.02      -      -     -     -   -" << std::endl;

	  textFile << "CMS_htt_tjXsec                   lnN      -      -   1.06      -      -      -      -" << std::endl;
	  textFile << "CMS_htt_wjXsec                   lnN      -      -      -      -   1.06      -      -" << std::endl;
	  textFile << "CMS_htt_vvXsec                   lnN      -      -      -      -      -   1.06      -" << std::endl;
	  textFile << "CMS_htt_zjXsec                   lnN      -      -      -      -      -      -   1.02" << std::endl;
	  textFile << "CMS_htt_zjXsec                   lnN      -      -      -      -      -      -   1.02" << std::endl;
	  textFile << "lumi_" << Era  << "           lnN   1.015    1.015   1.015      -   1.015   1.015   1.015" << std::endl;

	  textFile << "lumi                             lnN   1.015    1.015   1.015      -   1.015   1.015   1.015" << std::endl;
	  textFile << "CMS_scale_e                    shape   1.00       -    1.00     -   1.00    1.00     -" << std::endl;
	  //	  textFile << "CMS_scale_m                    shape   1.00    1.00    1.00     -   1.00    1.00   1.00" << std::endl;
	  textFile << "CMS_scale_e_emb                shape      -    1.00       -     -      -       -      -" << std::endl;
	  textFile << "CMS_htt_ttbarShape             shape      -       -     1.00    -      -       -      -" << std::endl;
	  if (Era=="2016")
	    textFile << "CMS_htt_dyShape_2016           shape      -       -        -    -      -       -   0.10" << std::endl;
	  else
	    textFile << "CMS_htt_dyShape                shape      -       -        -    -      -       -   0.10" << std::endl;

	  for (auto sysName : sysToLnNames_Common) {
	    textFile << sysName << "   lnN  "
		     << SIG_sys[sysName] 
		     << " - " 
		     << TTL_sys[sysName] 
		     << " - "
		     << W_sys[sysName]
		     << VVL_sys[sysName]
		     << ZLL_sys[sysName] << std::endl;
	  }

	  for (auto sysNameEra : sysToLnNames_Era) {
	    TString sysName = sysNameEra + "_" + Era;
	    textFile << sysName << "   lnN  "
		     << SIG_sys[sysName] 
		     << " - " 
		     << TTL_sys[sysName] 
		     << " - "
		     << W_sys[sysName]
		     << VVL_sys[sysName]
		     << ZLL_sys[sysName] << std::endl;
	  }
	  
	  textFile << "CMS_htt_boson_res_met_" << Era << "  shape  1.0   -   -   -  -  -  1.0" << std::endl;
	  textFile << "CMS_htt_boson_scale_met_" << Era << "  shape 1.0  -   -   -  -  -  1.0" << std::endl;
	  
	  textFile << "* autoMCStats 10 "<< std::endl;
	  textFile.close();

	}
      }
    }
  }

}
