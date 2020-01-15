#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TLorentzVector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"


class JESUncertainties {

 public:

  JESUncertainties(std::string uncFileName);
  ~JESUncertainties();


  void runOnEvent(AC1B &analysisTree, float etalep1, float philep1, float etalep2, float philep2);
  std::vector<std::string> getUncertNames(); 
  int getNJets(std::string uncName, bool Up);
  int getNJets();
  int getNJetsPt20(std::string uncName, bool Up);
  int getNJetsPt20();
  float getMjj(std::string uncName, bool Up);
  float getMjj();
  float getJdeta(std::string uncName, bool Up);
  float getJdeta();
  float getUncertainty(std::string uncName, float pt, float eta);

 private:

  std::string uncertFile;

  const float etaJetMax = 4.7;
  const float etaBtagJetMax = 2.4;
  const float ptJetMin = 30;
  const float ptBtagJetMin = 20;
  const float btagDiscr = 0.8484;
  const float dRjetLep = 0.5;

  float mjj;
  float jdeta;
  int njets;
  int njetspt20;

  std::map<std::string, float> mjjUp;
  std::map<std::string, float> mjjDown;
  std::map<std::string, float> jdetaUp;
  std::map<std::string, float> jdetaDown;
  std::map<std::string, int> njetsUp;
  std::map<std::string, int> njetsDown;
  std::map<std::string, int> njetspt20Up;
  std::map<std::string, int> njetspt20Down;



  std::map<std::string, JetCorrectorParameters const *> JetCorParMap;
  std::map<std::string, JetCorrectionUncertainty* > JetUncMap;
  std::vector< std::string > uncertNames;

  // unc names 2016
  std::vector< std::string > uncertNames_2016 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_2016",
    "EC2_2016",
    "Absolute_2016",
    "HF_2016",
    "RelativeSample_2016",
  }; // end of 2016 uncertNames

  // unc names 2017
  std::vector< std::string > uncertNames_2017 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_2017",
    "RelativeSample_2017",
    "EC2_2017",
    "HF_2017",
    "Absolute_2017"
  }; // end of 2017 uncertNames

  // unc names 2018
  std::vector< std::string > uncertNames_2018 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "Absolute_2018",
    "HF_2018",
    "EC2_2018",
    "RelativeSample_2018",
    "BBEC1_2018"
  }; // end of 2018 uncertNames


};

