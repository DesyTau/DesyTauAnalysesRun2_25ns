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
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"

using namespace std;
vector<TString> categories = {
  "em_NbtagGt1_DZetaLtm35",
  "em_DZetaLtm35",
  "em_NbtagGt1_DZetam35Tom10",
  "em_NbtagGt1_DZetam10To30",
  "em_NbtagGt1_DZetaGt30",
  "em_Nbtag0_DZetam35Tom10",
  "em_Nbtag0_DZetam10To30",
  "em_Nbtag0_DZetaGt30",
  "em_Nbtag0_DZetam35Tom10_MHGt250",
  "em_Nbtag0_DZetam10To30_MHGt250",
  "em_Nbtag0_DZetaGt30_MHGt250",
};

vector<TString> masses = {
    "_60",
    "_80",
    "_95",
    "_100",
    "_120",
    "_125",
    "_130",
    "_140",
    "_160",
    "_180",
    "_200",
    "_250",
    "_300",
    "_350",
    "_400",
    "_450",
    "_500",
    "_600",
    "_700",
    "_800",
    "_900",
    "_1000",
    "_1200",
    "_1400",
    "_1600",
    "_1800",
    "_2000",
    "_2300",
    "_2600",
    "_2900",
    "_3200",
    "_3500"
};

vector<TString> systematics = {
  "CMS_scale_e",
  "CMS_scale_m",
  "CMS_scale_j_FlavorQCD",
  "CMS_scale_j_RelativeBal",
  "CMS_scale_j_HF",
  "CMS_scale_j_BBEC1",
  "CMS_scale_j_EC2",
  "CMS_scale_j_Absolute",
  "CMS_scale_met_unclustered_2016",
  "CMS_htt_boson_res_met_2016",
  "CMS_htt_boson_scale_met_2016",
  "CMS_htt_eff_b_2016",
  "CMS_htt_mistag_b_2016",
  "CMS_scale_j_Absolute_2016",
  "CMS_scale_j_HF_2016",
  "CMS_scale_j_EC2_2016",
  "CMS_scale_j_RelativeSample_2016",
  "CMS_scale_j_BBEC1_2016",
  "CMS_res_j_2016",
  "CMS_scale_met_unclustered_2017",
  "CMS_htt_boson_res_met_2017",
  "CMS_htt_boson_scale_met_2017",
  "CMS_htt_eff_b_2017",
  "CMS_htt_mistag_b_2017",
  "CMS_scale_j_Absolute_2017",
  "CMS_scale_j_HF_2017",
  "CMS_scale_j_EC2_2017",
  "CMS_scale_j_RelativeSample_2017",
  "CMS_scale_j_BBEC1_2017",
  "CMS_res_j_2017",
  "CMS_scale_met_unclustered_2018",
  "CMS_htt_boson_res_met_2018",
  "CMS_htt_boson_scale_met_2018",
  "CMS_htt_eff_b_2018",
  "CMS_htt_mistag_b_2018",
  "CMS_scale_j_Absolute_2018",
  "CMS_scale_j_HF_2018",
  "CMS_scale_j_EC2_2018",
  "CMS_scale_j_RelativeSample_2018",
  "CMS_scale_j_BBEC1_2018",
  "CMS_res_j_2018",
  "CMS_prefiring",
  "CMS_htt_dyShape_2016",
  "CMS_htt_dyShape",
  "CMS_htt_ttbarShape",
  "CMS_scale_ggH",
  "CMS_PS_ISR_ggH",
  "CMS_PS_FSR_ggH",
  "QCDscale_ggH_REWEIGHT",
  "Hdamp_ggH_t_REWEIGHT",
  "Hdamp_ggH_b_REWEIGHT",
  "Hdamp_ggH_i_REWEIGHT"
};

vector<TString> dummy = {""};

void WriteHisto(TFile * inputFile,
		TFile * outputFile,
		TString category,
		TString sampleToProcess) {
  TH1D * hist = (TH1D*)inputFile->Get(category+"/"+sampleToProcess);
  if (hist!=NULL) {
    outputFile->cd(category);
    hist->Write(sampleToProcess);
  }
  else {
    cout << "Histogram " << sampleToProcess << " is absent!!!!" << std::endl; 
  }
}

int main(int argc, char * argv[]) {

  if (argc<3) {
    cout << "Usage : CopyDatacards Sample={VVL,TTL,W,ZL,ggH_{t,b,i},ggA_{t,b,i},bbH,qqH95,ggH125,qqH125,WH125,ZH125,ggHWW125,qqHWW125,WHWW125,ZHWW125} Era={2016,2017,2018}" << endl;
    exit(-1);
  }

  TString sample(argv[1]);
  TString era(argv[2]);
  
  vector<TString> suffixes = dummy;

  if (sample=="bbH"||sample=="ggH_t"||sample=="ggH_b"||sample=="ggH_i"||sample=="ggA_t"||sample=="ggA_b"||sample=="ggA_i"||sample=="ggh_t"||sample=="ggh_b"||sample=="ggh_i") {
    suffixes = masses;
  }

  TString inputFileName = "/nfs/dust/cms/user/rasp/CMSSW/CMSSW_10_2_13/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes_Jun3/"+era+"/em/htt_em_mssm.root";

  TString outputFileName = "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/datacardsZTT/"+era+"/"+sample+".root";
  TFile * inputFile = new TFile(inputFileName);
  TFile * outputFile = new TFile(outputFileName,"recreate");
  
  for (auto category : categories) {
    outputFile->cd("");
    outputFile->mkdir(category);
    outputFile->cd(category);
    for (auto suffix : suffixes) {
      TString sampleToProcess = sample + suffix;
      WriteHisto(inputFile,outputFile,category,sampleToProcess);
      for (auto sys : systematics) {
	sampleToProcess = sample + suffix + "_" + sys + "Up";
	WriteHisto(inputFile,outputFile,category,sampleToProcess);
	sampleToProcess = sample + suffix + "_" + sys + "Down";
	WriteHisto(inputFile,outputFile,category,sampleToProcess);	
      }
    }
  }
  outputFile->Close();

}
