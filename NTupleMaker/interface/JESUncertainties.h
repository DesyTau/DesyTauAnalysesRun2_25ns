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

  // unc names
  std::vector< std::string > uncertNames = {
    "AbsoluteFlavMap",
    "AbsoluteMPFBias",
    "AbsoluteScale",
    "AbsoluteStat",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpEnvelope",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "RelativeBal",
    //    "RelativeSample",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta",
    "Total",
  }; // end uncertNames

  // groupped uncertainties
  std::vector< std::string > groupedUncertNames = {"Eta0To5",
						   "Eta0To3",
						   "Eta3To5",
						   "RelativeBal",
						   //						   "RelativeSample"
  };

  int nsrc_Eta0To5 = 13;
  std::vector<std::string> srcnames_Eta0To5 = {"SinglePionECAL",
					       "SinglePionHCAL",
					       "AbsoluteFlavMap",
					       "AbsoluteMPFBias",
					       "AbsoluteScale",
					       "AbsoluteStat",
					       "Fragmentation",
					       "FlavorQCD",
					       "TimePtEta",
					       "PileUpDataMC",
					       "RelativeFSR",
					       "RelativeStatFSR",
					       "PileUpPtRef"};
 int nsrc_Eta0To3 = 9;
  std::vector<std::string> srcnames_Eta0To3 = {"PileUpPtEC1",
					       "PileUpPtEC2",
					       "PileUpPtBB",
					       "RelativeJEREC1",
					       "RelativeJEREC2",
					       "RelativePtEC1",
					       "RelativePtEC2",
					       "RelativeStatEC",
					       "RelativePtBB"};
  int nsrc_Eta3To5 = 4;
  std::vector<std::string> srcnames_Eta3To5 = {"RelativeStatHF",
					       "RelativePtHF",
					       "PileUpPtHF",
					       "RelativeJERHF"};
  int nsrc_RelativeBal = 1;
  std::vector<std::string> srcnames_RelativeBal = {"RelativeBal"};
  
  int nsrc_RelativeSample = 1;
  std::vector<std::string> srcnames_RelativeSample = {"RelativeSample"};

  std::map<std::string, std::vector<std::string> > MapUncert = {
    {groupedUncertNames[0],srcnames_Eta0To5},
    {groupedUncertNames[1],srcnames_Eta0To3},
    {groupedUncertNames[2],srcnames_Eta3To5},
    {groupedUncertNames[3],srcnames_RelativeBal},
    //    {groupedUncertNames[4],srcnames_RelativeSample}
  };



};

