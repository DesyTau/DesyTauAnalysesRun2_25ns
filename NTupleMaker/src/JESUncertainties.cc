#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"

JESUncertainties::JESUncertainties(std::string uncFileName) {
  
  std::string cmsswBase( getenv ("CMSSW_BASE") );
  uncertFile = cmsswBase + "/src/" + uncFileName;

  TString UncertFile(uncertFile);
  if (UncertFile.Contains("Regrouped_Autumn18"))
    uncertNames = uncertNames_2018;
  else if (UncertFile.Contains("Regrouped_Fall17"))
    uncertNames = uncertNames_2017;
  else 
    uncertNames = uncertNames_2016;

  for (auto const& name : uncertNames) {

    JetCorrectorParameters const * JetCorPar = new JetCorrectorParameters(uncertFile, name);
    JetCorParMap[name] = JetCorPar;

    JetCorrectionUncertainty * jecUnc(new JetCorrectionUncertainty(*JetCorParMap[name]));
    JetUncMap[name] = jecUnc;

    JetUncMap[name]->setJetPt(50.);
    JetUncMap[name]->setJetEta(1.6);
      
  };
  

  std::cout << "++++++" << std::endl;

}

JESUncertainties::~JESUncertainties() {

}

std::vector<std::string> JESUncertainties::getUncertNames() {

  return uncertNames;

}

float JESUncertainties::getUncertainty(std::string name, float pt, float eta) {

  double unc = 0;
  if (JetUncMap[name]!=NULL) { 
    JetUncMap[name]->setJetEta(eta);
    JetUncMap[name]->setJetPt(pt);
    unc = JetUncMap[name]->getUncertainty(true);
  }
  //  std::cout << "OK" << std::endl;

  return unc;

}

int JESUncertainties::getNJets(std::string uncName, bool Up) {

  int result = njetsDown[uncName];
  if (Up) result = njetsUp[uncName]; 

  return result;

}

int JESUncertainties::getNJets() {
  return njets;
}

int JESUncertainties::getNJetsPt20(std::string uncName, bool Up) {

  int result = njetspt20Down[uncName];
  if (Up) result = njetspt20Up[uncName]; 

  return result;

}

int JESUncertainties::getNJetsPt20() {
  return njetspt20;
}

float JESUncertainties::getMjj(std::string uncName, bool Up) {

  float result = mjjDown[uncName];
  if (Up) result = mjjUp[uncName];

  return result;

}

float JESUncertainties::getMjj() {
  return mjj;
}

float JESUncertainties::getJdeta(std::string uncName, bool Up) {

  float result = jdetaDown[uncName];
  if (Up) result = jdetaUp[uncName];

  return result;

}

float JESUncertainties::getJdeta() {
  return jdeta;
}


void JESUncertainties::runOnEvent(AC1B &analysisTree, float etalep1, float philep1, float etalep2, float philep2) {

  //  std::map<std::string, int> indexLeadingUp;
  //  std::map<std::string, int> indexLeadingDown;
  //  std::map<std::string, int> indexSubLeadingUp;
  //  std::map<std::string, int> indexSubLeadingDown;

  int indexLeading = -1;
  int indexTrailing = -1;
  float ptLeading = 0;
  float ptTrailing = 0;
  mjj = -9999;
  jdeta = -9999;
  njets = 0;
  njetspt20 = 0;

  for (auto const& name : uncertNames) {
    mjjUp[name] = -9999;
    mjjDown[name] = -9999;
    jdetaUp[name] = -9999;
    jdetaDown[name] = -9999;
    njetsUp[name] = 0;
    njetsDown[name] = 0;
    njetspt20Up[name] = 0;
    njetspt20Down[name] = 0;
    
  }

  for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

    float absEta = fabs(analysisTree.pfjet_eta[jet]);
    if (absEta>etaJetMax) continue;
    float jetPt = analysisTree.pfjet_pt[jet];
    float jetEta = analysisTree.pfjet_eta[jet];
    float jetPhi = analysisTree.pfjet_phi[jet];

    float dR1 = deltaR(etalep1,philep1,jetEta,jetPhi);
    if (dR1<dRjetLep) continue;
    float dR2 = deltaR(etalep2,philep2,jetEta,jetPhi);
    if (dR2<dRjetLep) continue;


    bool jetId = looseJetiD(analysisTree,int(jet));
    if (!jetId) continue;    

    if (jetPt>ptBtagJetMin) 
      njetspt20++;

    if (jetPt>ptJetMin)
      njets++;

    if (indexLeading>=0) {
      if (jetPt<ptLeading&&jetPt>ptTrailing) {
	indexTrailing = jet;
	ptTrailing = jetPt;
      }
    }

    if (jetPt>ptLeading) {
      indexTrailing = indexLeading;
      ptTrailing = ptLeading;
      indexLeading = jet;
      ptLeading = jetPt;
    }

    for (auto const& name : uncertNames) {
      
      float unc = getUncertainty(name,jetPt,jetEta);

      float jetPtUp = (1+unc)*jetPt;
      float jetPtDown = (1-unc)*jetPt;

      if (jetPtUp>ptJetMin)
	njetsUp[name]++;

      if (jetPtDown>ptJetMin)
	njetsDown[name]++;

      if (jetPtUp>ptBtagJetMin)
	njetspt20Up[name]++;

      if (jetPtDown>ptBtagJetMin)
	njetspt20Down[name]++;

    }

  }
  

  if (indexLeading>=0&&indexTrailing>=0) {

    TLorentzVector jet1; jet1.SetXYZT(analysisTree.pfjet_px[indexLeading],
				      analysisTree.pfjet_py[indexLeading],
				      analysisTree.pfjet_pz[indexLeading],
				      analysisTree.pfjet_e[indexLeading]);

    TLorentzVector jet2; jet2.SetXYZT(analysisTree.pfjet_px[indexTrailing],
                                      analysisTree.pfjet_py[indexTrailing],
                                      analysisTree.pfjet_pz[indexTrailing],
                                      analysisTree.pfjet_e[indexTrailing]);

    if (njetspt20>=2) {
      mjj = (jet1+jet2).M();
      jdeta = fabs(jet1.Eta()-jet2.Eta());
    }

    for (auto const& name : uncertNames) {

      float jetEta1 = analysisTree.pfjet_eta[indexLeading];
      float jetPt1 = analysisTree.pfjet_pt[indexLeading];
      float unc1 = getUncertainty(name,jetPt1,jetEta1);

      float jetEta2 = analysisTree.pfjet_eta[indexTrailing];
      float jetPt2 = analysisTree.pfjet_pt[indexTrailing];
      float unc2 = getUncertainty(name,jetPt2,jetEta2);

      TLorentzVector jet1Up = (1+unc1)*jet1;
      TLorentzVector jet1Down = (1-unc1)*jet1;
      
      TLorentzVector jet2Up = (1+unc2)*jet2;
      TLorentzVector jet2Down = (1-unc2)*jet2;
      
      if (njetspt20Up[name]>=2) {
	mjjUp[name] = (jet1Up+jet2Up).M();
	jdetaUp[name] = fabs(jet1Up.Eta()-jet2Up.Eta());
      }

      if (njetspt20Down[name]>=2) {
	mjjDown[name] = (jet1Down+jet2Down).M();
	jdetaDown[name] = fabs(jet1Down.Eta()-jet2Down.Eta());
      }

    }

  }


}
