#include "TMath.h"
#include "TLorentzVector.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>


using namespace std;

typedef std::vector<std::pair<int,int> > lumi_json;

struct compare_lumi { //accepts two pairs, return 1 if left.first < right.first or left.first = right.first e left.second < right.second
  bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
    if (left.first < right.first)
      return 1;
    else if (left.first > right.first)
      return 0;
    else
      return left.second < right.second;
  }
};


int read_json(std::string filename, lumi_json& json);
bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json);
bool isGoodLumi(int run, int lumi, const lumi_json& json);
void fill_weight(const AC1B * analysisTree, Synch17Tree *otree, PileUp *PUofficial, bool isData);
float abs_Iso_mt(int Index, const AC1B * analysisTree, float dRiso);
float abs_Iso_et(int Index, const AC1B * analysisTree, float dRiso);
bool passedSummer16VetoId(const AC1B * analysisTree, int index);
bool isICHEPmed(int Index, const AC1B * analysisTree);
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData);
float getEffectiveArea(float eta);
int ZDecay(const AC1B * analysisTree);
float EmbedElectronES_SF(const AC1B * analysisTree, int era, int electronIndex );

///////////////////////////////////////////////
//////////////FUNCTION DEFINITION//////////////
///////////////////////////////////////////////

bool isICHEPmed(int Index, const AC1B * analysisTree) {

  bool goodGlob = analysisTree->muon_isGlobal[Index] && 
                  analysisTree->muon_normChi2[Index] < 3 && 
                  analysisTree->muon_combQ_chi2LocalPosition[Index] < 12 && 
                  analysisTree->muon_combQ_trkKink[Index] < 20;

  bool isICHEPmedium  = analysisTree->muon_isLoose[Index] &&
                        analysisTree->muon_validFraction[Index] >0.49 &&
                        analysisTree->muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
  
  return isICHEPmedium;

}

int read_json(std::string filename, lumi_json& json){

  std::pair <int,int> lumi;

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);

  BOOST_FOREACH(boost::property_tree::ptree::value_type &json_run, pt.get_child("")){
    int irun = atoi(json_run.first.data());
    BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_ranges, json_run.second.get_child("")){
      int ilumi[2] = {};

      int count = 0;
      BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_boundaries, lumi_ranges.second.get_child("")){
	ilumi[count] = atoi(lumi_boundaries.second.data().data());
	count++;
      }
      
      for (;ilumi[0] <= ilumi[1]; ilumi[0]++){
	lumi = std::make_pair(irun, ilumi[0]);
	json.push_back(lumi);
      }
    }
  }

  sort( json.begin(), json.end(),  compare_lumi());
  json.erase( unique( json.begin(), json.end() ), json.end() );
  
  return 0;
}


bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json){
  static compare_lumi compare;
  static std::pair<int,int> oldlumi = lumi;
  static bool old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  
  if(lumi.first != oldlumi.first || lumi.second != oldlumi.second){
    oldlumi = lumi;
    old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  }

  return old;
}

//accepts run number, lumi and json, make pair of (run,lumi) and starts isGoodLumi 
bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}

float getEffectiveArea(float eta) {
    float effArea =  0.1440;
    float absEta = fabs(eta);
    if (absEta<1.0) effArea = 0.1440;
    else if (absEta < 1.4790) effArea = 0.1562;
    else if (absEta < 2.0) effArea = 0.1032;
    else if (absEta < 2.2) effArea = 0.0859;
    else if (absEta < 2.3) effArea = 0.1116;
    else if (absEta < 2.4) effArea = 0.1321;
    else if (absEta < 5.0) effArea = 0.1654;
    return effArea;

  } 


float GenMatch(const AC1B * analysisTree, bool isData, TString particleType, int Index){

  float gen_match = 6;
  
  // isZTT = false;
  // isZL  = false;
  // isZJ  = false;
  
  float minDR = 0.2;
  float pt,eta,phi;
  if(particleType=="m"){
    pt =      analysisTree->muon_pt[Index]; 
    eta =     analysisTree->muon_eta[Index]; 
    phi =     analysisTree->muon_phi[Index];}
  if(particleType=="e"){         
    pt =      analysisTree->electron_pt[Index];
    eta =     analysisTree->electron_eta[Index]; 
    phi =     analysisTree->electron_phi[Index];}
  if(particleType=="t"){         
    pt =      analysisTree->tau_pt[Index];
    eta =     analysisTree->tau_eta[Index]; 
    phi =     analysisTree->tau_phi[Index];
  }
  
  if (!isData){
    for (unsigned int igen=0; igen < analysisTree->genparticles_count; ++igen) {
      
      TLorentzVector genLV; genLV.SetXYZT(analysisTree->genparticles_px[igen],
					  analysisTree->genparticles_py[igen],
					  analysisTree->genparticles_pz[igen],
					  analysisTree->genparticles_e[igen]);
      float ptGen = genLV.Pt();
      bool type1 = abs(analysisTree->genparticles_pdgid[igen])==11 && analysisTree->genparticles_isPrompt[igen] && ptGen>8;
      bool type2 = abs(analysisTree->genparticles_pdgid[igen])==13 && analysisTree->genparticles_isPrompt[igen] && ptGen>8;
      bool type3 = abs(analysisTree->genparticles_pdgid[igen])==11 && analysisTree->genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
      bool type4 = abs(analysisTree->genparticles_pdgid[igen])==13 && analysisTree->genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
      
      bool isAnyType = type1 || type2 || type3 || type4;
      if (isAnyType && analysisTree->genparticles_status[igen]==1) {
	float etaGen = genLV.Eta();
	float phigen = genLV.Phi();
	float dR = deltaR(eta,phi,etaGen,phigen);
	if (dR<minDR) {
	  minDR = dR;
	  if (type1) gen_match = 1;
	  else if (type2) gen_match = 2;
	  else if (type3) gen_match = 3;
	  else if (type4) gen_match = 4;
	}	
      }
    }
    if (gen_match<1||gen_match>4) {
      for (unsigned int igen=0; igen < analysisTree->gentau_count; ++igen) {
	
	if (analysisTree->gentau_visibleNoLep_pt[igen]>15.) {
	  TLorentzVector genTauLV; genTauLV.SetXYZT(analysisTree->gentau_px[igen],
						    analysisTree->gentau_py[igen],
						    analysisTree->gentau_pz[igen],
						    analysisTree->gentau_e[igen]);
	  float dR = deltaR(eta,phi,
			    genTauLV.Eta(),genTauLV.Phi());
	  if (dR<minDR) {
	    minDR = dR;
	    gen_match = 5;
	  }
	}
      }
    }
  } 

  return gen_match;  
}

//////////////////////////////////////////////
//            channel dependent             //
//////////////////////////////////////////////

// select medium id or ICHEP medium id for different runs, in data
// use ICHEP medium ID fro runs up to 278808 (end of Run2016F)
// on MC, always apply medium ID
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData){
  bool isGoodMuon;
  if (isData){
  	if (analysisTree->event_run <= 278808) 
      isGoodMuon = isICHEPmed(Index, analysisTree);
    else 
      isGoodMuon = analysisTree->muon_isMedium[Index]; 
	}
  else 
    isGoodMuon = analysisTree->muon_isMedium[Index]; 
  return isGoodMuon;
}

//compute the absolute isolation for a given lepton labeled by Index in channel ch
float abs_Iso_mt (int Index, const AC1B * analysisTree, float dRiso){

  float neutralHadIso = -9999.;
  float photonIso     = -9999.;
  float chargedHadIso = -9999.;
  float puIso         = -9999.;

  neutralHadIso = analysisTree->muon_neutralHadIso[Index];
  photonIso =     analysisTree->muon_photonIso[Index];
  chargedHadIso = analysisTree->muon_chargedHadIso[Index];
  puIso =         analysisTree->muon_puIso[Index];
  if (dRiso>0.29) {
    neutralHadIso =     analysisTree->muon_r03_sumNeutralHadronEt[Index];
    photonIso =         analysisTree->muon_r03_sumPhotonEt[Index];
    chargedHadIso =     analysisTree->muon_r03_sumChargedHadronPt[Index];
    puIso =             analysisTree->muon_r03_sumPUPt[Index];
  }
  if (dRiso>0.39) {
    neutralHadIso =     analysisTree->muon_r04_sumNeutralHadronEt[Index];
    photonIso =         analysisTree->muon_r04_sumPhotonEt[Index];
    chargedHadIso =     analysisTree->muon_r04_sumChargedHadronPt[Index];
    puIso =             analysisTree->muon_r04_sumPUPt[Index];
  }
  
  float neutralIso = neutralHadIso + photonIso -0.5*puIso;
  neutralIso = TMath::Max(float(0), neutralIso);
  return(chargedHadIso + neutralIso);
}

float abs_Iso_et (int Index, const AC1B * analysisTree, float dRiso){

  float absIso        = -9999.;
  float isoCone       = 0.3;
  float neutralIso    = -9999.;
  /*
  absIso  = analysisTree->electron_r03_sumChargedHadronPt[Index];
  float neutralIso =
    analysisTree->electron_r03_sumNeutralHadronEt[Index] +
    analysisTree->electron_r03_sumPhotonEt[Index] - 
    analysisTree->rho*TMath::Pi()*isoCone*isoCone;
  neutralIso = TMath::Max(float(0),neutralIso);
  return(absIso + neutralIso);
  */
 
  float  eA = getEffectiveArea( fabs(analysisTree->electron_superclusterEta[Index]) );
  absIso = analysisTree->electron_r03_sumChargedHadronPt[Index];
  neutralIso = analysisTree->electron_r03_sumNeutralHadronEt[Index] +
    analysisTree->electron_r03_sumPhotonEt[Index] -
    eA*analysisTree->rho;

  return (absIso+TMath::Max(float(0),neutralIso));


  /*

  neutralHadIso = analysisTree->electron_neutralHadIso[Index];
  photonIso =     analysisTree->electron_photonIso[Index];
  chargedHadIso = analysisTree->electron_chargedHadIso[Index];
  puIso =         analysisTree->electron_puIso[Index];
  if (dRiso>0.29) {
    neutralHadIso =     analysisTree->electron_r03_sumNeutralHadronEt[Index];
    photonIso =         analysisTree->electron_r03_sumPhotonEt[Index];
    chargedHadIso =     analysisTree->electron_r03_sumChargedHadronPt[Index];
    puIso =             analysisTree->electron_r03_sumPUPt[Index];
  }
  
  float neutralIso = neutralHadIso + photonIso -0.5*puIso;
  neutralIso = TMath::Max(float(0), neutralIso);
  return(chargedHadIso + neutralIso);
  */
}



//////DILEPTON FUNCTIONS

//returns the dilepton veto for mt channel
bool dilepton_veto_mt(const Config *cfg,const  AC1B *analysisTree){

  for (unsigned int im = 0; im < analysisTree->muon_count; ++im) {
    
    float ptDiMuonVeto = cfg->get<float>("ptDiMuonVeto");
    float etaDiMuonVeto = cfg->get<float>("etaDiMuonVeto");
    float dxyDiMuonVeto = cfg->get<float>("dxyDiMuonVeto");
    float dzDiMuonVeto = cfg->get<float>("dzDiMuonVeto");
    float isoDiMuonVeto = cfg->get<float>("isoDiMuonVeto");
    float dRiso = cfg->get<float>("dRisoDiMuonVeto");
    float drDiMuonVeto = cfg->get<float>("drDiMuonVeto");
    
    if (analysisTree->muon_pt[im] <= ptDiMuonVeto) continue;
    if (fabs(analysisTree->muon_eta[im]) >= etaDiMuonVeto) continue;	
    if (fabs(analysisTree->muon_dxy[im]) >= dxyDiMuonVeto) continue;
    if (fabs(analysisTree->muon_dz[im]) >= dzDiMuonVeto) continue;

    float relIsoMu = abs_Iso_mt(im, analysisTree, dRiso) / analysisTree->muon_pt[im];
    if(relIsoMu >= isoDiMuonVeto) continue;
		    if ( !(analysisTree->muon_isGlobal[im] && analysisTree->muon_isTracker[im] && analysisTree->muon_isPF[im]) ) continue;
    
    for (unsigned int je = im + 1; je < analysisTree->muon_count; ++je) {  
      if (analysisTree->muon_pt[je] <= ptDiMuonVeto) continue;
      if (fabs(analysisTree->muon_eta[je]) >= etaDiMuonVeto) continue;	
      if (fabs(analysisTree->muon_dxy[je]) >= dxyDiMuonVeto) continue;
      if (fabs(analysisTree->muon_dz[je]) >= dzDiMuonVeto) continue;
      if (analysisTree->muon_charge[im] * analysisTree->muon_charge[je] > 0.) continue;

      float relIsoMu = abs_Iso_mt(je, analysisTree, dRiso) / analysisTree->muon_pt[je];
      if(relIsoMu >= isoDiMuonVeto) continue;	
      if ( ! (analysisTree->muon_isGlobal[je] && analysisTree->muon_isTracker[je] && analysisTree->muon_isPF[je]) ) continue;
		  
      float dR = deltaR(analysisTree->muon_eta[im], analysisTree->muon_phi[im], analysisTree->muon_eta[je], analysisTree->muon_phi[je]);
      if(dR <= drDiMuonVeto) continue;

      return(1);
    }
  }
  return(0);
}


//returns the dilepton veto fot et channel
bool dilepton_veto_et(const Config *cfg,const  AC1B *analysisTree, int era, bool isEmbedded){
  float sf_eleES_i = 1.0;   
  float sf_eleES_j = 1.0;   
  for (unsigned int ie = 0; ie<analysisTree->electron_count; ++ie) {
    if (isEmbedded) sf_eleES_i = EmbedElectronES_SF(analysisTree, era, ie);
    if (sf_eleES_i*analysisTree->electron_pt[ie]<=cfg->get<float>("ptDiElectronVeto")) continue;
    if (fabs(analysisTree->electron_eta[ie])>=cfg->get<float>("etaDiElectronVeto")) continue;	
    if (fabs(analysisTree->electron_dxy[ie])>=cfg->get<float>("dxyDiElectronVeto")) continue;
    if (fabs(analysisTree->electron_dz[ie])>=cfg->get<float>("dzDiElectronVeto")) continue;

    float absIsoEle =   abs_Iso_et(ie, analysisTree, cfg->get<float>("dRisoDiElectronVeto"));
    float relIsoEle =   absIsoEle / (sf_eleES_i*analysisTree->electron_pt[ie]) ;
    if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;
    
    bool passedVetoId =  analysisTree->electron_cutId_veto_Fall17V2[ie];
    if (!passedVetoId) continue;
		
    for (unsigned int je = ie+1; je<analysisTree->electron_count; ++je) {
      if (isEmbedded) sf_eleES_j = EmbedElectronES_SF(analysisTree, era, je);
      if (sf_eleES_j*analysisTree->electron_pt[je]<=cfg->get<float>("ptDiElectronVeto")) continue;
      if (fabs(analysisTree->electron_eta[je])>=cfg->get<float>("etaDiElectronVeto")) continue;	
      if (fabs(analysisTree->electron_dxy[je])>=cfg->get<float>("dxyDiElectronVeto")) continue;
      if (fabs(analysisTree->electron_dz[je])>=cfg->get<float>("dzDiElectronVeto")) continue;
		  
      float absIsoEle =  abs_Iso_et(je, analysisTree, cfg->get<float>("dRiso"));
      float relIsoEle =  absIsoEle / (sf_eleES_j*analysisTree->electron_pt[je]);
      if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;	

      passedVetoId =  analysisTree->electron_cutId_veto_Fall17V2[je];
      if (!passedVetoId) continue;

      if (analysisTree->electron_charge[ie] * analysisTree->electron_charge[je] > 0.) continue;
		  
      float dr = deltaR(analysisTree->electron_eta[ie],analysisTree->electron_phi[ie],
                        analysisTree->electron_eta[je],analysisTree->electron_phi[je]);

      if(dr<=cfg->get<float>("drDiElectronVeto")) continue;

      return(1);
    }
  }
  return(0);
}

//////EXTRA LEPTON VETO FUNCTIONS

//returns the extra electron veto
bool extra_electron_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree, int era, bool isEmbedded){
  for (unsigned int ie = 0; ie < analysisTree->electron_count; ++ie) {
    float sf_eleES = 1.0;
    if (isEmbedded && ch == "et") sf_eleES = EmbedElectronES_SF(analysisTree, era, ie);
    float ptVetoElectronCut = cfg->get<float>("ptVetoElectronCut");
    float etaVetoElectronCut = cfg->get<float>("etaVetoElectronCut");
    float dxyVetoElectronCut = cfg->get<float>("dxyVetoElectronCut");
    float dzVetoElectronCut = cfg->get<float>("dzVetoElectronCut");
    float applyVetoElectronId = cfg->get<bool>("applyVetoElectronId");
    float dRisoExtraElecVeto = cfg->get<float>("dRisoExtraElecVeto");
    float isoVetoElectronCut = cfg->get<float>("isoVetoElectronCut");

    if (ch == "et" && int(ie) == leptonIndex) continue;
    if (sf_eleES*analysisTree->electron_pt[ie] <= ptVetoElectronCut) continue;
    if (fabs(analysisTree->electron_eta[ie]) >= etaVetoElectronCut) continue;
    if (fabs(analysisTree->electron_dxy[ie]) >= dxyVetoElectronCut) continue;
    if (fabs(analysisTree->electron_dz[ie]) >= dzVetoElectronCut) continue;

    bool electronMvaId = analysisTree->electron_mva_wp90_noIso_Fall17_v2[ie] > 0.5;
    if (!electronMvaId && applyVetoElectronId) continue;
    if (!analysisTree->electron_pass_conversion[ie] && applyVetoElectronId) continue;
    if (analysisTree->electron_nmissinginnerhits[ie] > 1 && applyVetoElectronId) continue;

    float relIsoEle = abs_Iso_et(ie, analysisTree, dRisoExtraElecVeto) / (sf_eleES*analysisTree->electron_pt[ie]);
    if (relIsoEle >= isoVetoElectronCut) continue;

    return(1);		
  }
  return(0);
}			

//returns the extra muon veto
bool extra_muon_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree, bool isData){
  for (unsigned int im = 0; im < analysisTree->muon_count; ++im) {
    float ptVetoMuonCut = cfg->get<float>("ptVetoMuonCut");
    float etaVetoMuonCut = cfg->get<float>("etaVetoMuonCut");
    float dxyVetoMuonCut = cfg->get<float>("dxyVetoMuonCut");
    float dzVetoMuonCut = cfg->get<float>("dzVetoMuonCut");
    float applyVetoMuonId = cfg->get<bool>("applyVetoMuonId");
    float dRisoExtraMuonVeto = cfg->get<float>("dRisoExtraMuonVeto");
    float isoVetoMuonCut = cfg->get<float>("isoVetoMuonCut");
    
    if (ch == "mt" && int(im) == leptonIndex) continue;
    if (analysisTree->muon_pt[im] <= ptVetoMuonCut) continue;
    if (fabs(analysisTree->muon_eta[im]) >= etaVetoMuonCut) continue;
    if (fabs(analysisTree->muon_dxy[im]) >= dxyVetoMuonCut) continue;
    if (fabs(analysisTree->muon_dz[im]) >= dzVetoMuonCut) continue;

    if (applyVetoMuonId && !(isIdentifiedMediumMuon(im, analysisTree, isData)) ) continue;
    float relIsoMu = abs_Iso_mt(im, analysisTree, dRisoExtraMuonVeto) / analysisTree->muon_pt[im];
    if (relIsoMu >= isoVetoMuonCut) continue;

    return(1);
  }
  return(0);
}


//////MET FUNCTIONS

//fill the otree with the met variables
void fillMET(TString ch, int leptonIndex, int tauIndex, const AC1B * analysisTree, Synch17Tree *otree){

   // pfmet variables
  
  otree->met = TMath::Sqrt(analysisTree->pfmet_ex*analysisTree->pfmet_ex + analysisTree->pfmet_ey*analysisTree->pfmet_ey);
  otree->metphi = TMath::ATan2(analysisTree->pfmet_ey,analysisTree->pfmet_ex);
  otree->metcov00 = analysisTree->pfmet_sigxx;
  otree->metcov01 = analysisTree->pfmet_sigxy;
  otree->metcov10 = analysisTree->pfmet_sigyx;
  otree->metcov11 = analysisTree->pfmet_sigyy;
  float met_x = analysisTree->pfmet_ex;
  float met_y = analysisTree->pfmet_ey;

  float met_x2 = met_x * met_x;
  float met_y2 = met_y * met_y;

}

//Merijn 2019 6 20: added overloaded function, takes era as argument. Will do MET correct for 2016 2017, will need later to extend to 2018
void fillMET(const AC1B * analysisTree, Synch17Tree *otree, int era){

   // pfmet variables
  if(era == 2018){
    otree->met = TMath::Sqrt(analysisTree->pfmet_ex*analysisTree->pfmet_ex + analysisTree->pfmet_ey*analysisTree->pfmet_ey);
    otree->metphi = TMath::ATan2(analysisTree->pfmet_ey,analysisTree->pfmet_ex);
    otree->metcov00 = analysisTree->pfmet_sigxx;
    otree->metcov01 = analysisTree->pfmet_sigxy;
    otree->metcov10 = analysisTree->pfmet_sigyx;
    otree->metcov11 = analysisTree->pfmet_sigyy;
  }
  else{
    otree->met = TMath::Sqrt(analysisTree->pfmetcorr_ex*analysisTree->pfmetcorr_ex + analysisTree->pfmetcorr_ey*analysisTree->pfmetcorr_ey);
    otree->metphi = TMath::ATan2(analysisTree->pfmetcorr_ey,analysisTree->pfmetcorr_ex);
    otree->metcov00 = analysisTree->pfmetcorr_sigxx;
    otree->metcov01 = analysisTree->pfmetcorr_sigxy;
    otree->metcov10 = analysisTree->pfmetcorr_sigyx;
    otree->metcov11 = analysisTree->pfmetcorr_sigyy;
  }

  otree->puppimet = TMath::Sqrt(analysisTree->puppimet_ex*analysisTree->puppimet_ex + analysisTree->puppimet_ey*analysisTree->puppimet_ey);
  otree->puppimetphi = TMath::ATan2(analysisTree->puppimet_ey,analysisTree->puppimet_ex);
  otree->puppimetcov00 = analysisTree->puppimet_sigxx;
  otree->puppimetcov01 = analysisTree->puppimet_sigxy;
  otree->puppimetcov10 = analysisTree->puppimet_sigyx;
  otree->puppimetcov11 = analysisTree->puppimet_sigyy;

}



///////////////////////////////////
// SV fit 
///////////////////////////////////

void svfit_variables(TString ch, const AC1B *analysisTree, Synch17Tree *otree, const Config *cfg, TFile * inputFile_visPtResolution){


  if (cfg->get<bool>("ApplySVFit")||cfg->get<bool>("ApplyFastMTT")) {
    
    double measuredMETx =  otree->met * cos(otree->metphi);
    double measuredMETy =  otree->met * sin(otree->metphi);
    
    bool isPuppiMET = cfg->get<bool>("UsePuppiMET");
    
    if (isPuppiMET) {
      measuredMETx = otree->puppimet * cos(otree->puppimetphi);
      measuredMETy = otree->puppimet * sin(otree->puppimetphi);
    }
    
    // define MET covariance
    TMatrixD covMET(2, 2);
    
    // using PF MET
    covMET[0][0] = otree->metcov00;
    covMET[1][0] = otree->metcov10;
    covMET[0][1] = otree->metcov01;
    covMET[1][1] = otree->metcov11;
    
    if (isPuppiMET) {
      covMET[0][0] = otree->puppimetcov00;
      covMET[1][0] = otree->puppimetcov10;
      covMET[0][1] = otree->puppimetcov01;
      covMET[1][1] = otree->puppimetcov11;
    }
    
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    classic_svFit::MeasuredTauLepton::kDecayType type_ = classic_svFit::MeasuredTauLepton::kUndefinedDecayType;
    
    if (ch == "mt") type_ = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    if (ch == "et") type_ = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    if (ch == "tt") type_ = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    if (ch == "em") type_ = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    // define lepton four vectors
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type_,
								  otree->pt_1,
								  otree->eta_1,
								  otree->phi_1,
								  otree->m_1));
    if (ch == "em") {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay,
								    otree->pt_2,
								    otree->eta_2,
								    otree->phi_2,
								    otree->m_2));
      
    }
    else {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,
								    otree->pt_2,
								    otree->eta_2,
								    otree->phi_2,
								    otree->m_2,
								    otree->tau_decay_mode_2));
    }
    if (cfg->get<bool>("ApplySVFit")) {
      int verbosity = 1;
      ClassicSVfit svFitAlgo(verbosity);
      double kappa = 4.; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
      
      svFitAlgo.addLogM_fixed(true, kappa);
      svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
      bool isValidSolution = svFitAlgo.isValidSolution();
      
      otree->m_sv   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
      otree->pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPt();
      otree->eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getEta();
      otree->phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPhi();      
      //otree->met_sv = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getfittedMET().Rho();
      otree->mt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
    }
    if (cfg->get<bool>("ApplyFastMTT")) {
      // FasMTT
      LorentzVector tau1P4;
      LorentzVector tau2P4;
      FastMTT aFastMTTAlgo;
      
      aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
      LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
      tau1P4 = aFastMTTAlgo.getTau1P4();
      tau2P4 = aFastMTTAlgo.getTau2P4();
      //    std::cout << "tau decay mode : " << otree->tau_decay_mode_2 << "   mvaDM = " << otree->dmMVA_2 << "  mass = " << otree->m_2 << "  DMFinding = " << otree->DM << std::endl;
      double dPhiTT = dPhiFrom2P( tau1P4.Px(), tau1P4.Py(), tau2P4.Px(), tau2P4.Py() );
      otree->mt_fast = TMath::Sqrt(2*tau1P4.Pt()*tau2P4.Pt()*(1 - TMath::Cos(dPhiTT)));
      otree->m_fast = ttP4.M();
      otree->pt_fast = ttP4.Pt();
      otree->eta_fast = ttP4.Eta();
      otree->phi_fast = ttP4.Phi();
    }
  }
}

float EmbedElectronES_SF(const AC1B * analysisTree, int era, int electronIndex ){
  float sf_ele=1.0;
  if (era == 2016){
     if (fabs(analysisTree->electron_eta[electronIndex]) < 1.479 ) sf_ele= (1-0.00243);
     else sf_ele = (1-0.007);
  }
  if (era == 2017){
     if (fabs(analysisTree->electron_eta[electronIndex]) < 1.479 ) sf_ele =  (1-0.00067);
     else sf_ele =  (1-0.01133);
  }
  if (era == 2018){
     if (fabs(analysisTree->electron_eta[electronIndex]) < 1.479 ) sf_ele = (1-0.00328);
     else sf_ele =  (1-0.00557);
  }
  return sf_ele;
}

/// shift tau energy scale and propagate it to the met. 
void correctTauES(TLorentzVector& Tau, TLorentzVector& Met, float relative_shift, bool tau_is_one_prong){
  Met.SetPx(Met.Px()- (Tau.Px()*relative_shift) ) ;
  Met.SetPy(Met.Py()- (Tau.Py()*relative_shift) );
  if (tau_is_one_prong){ // don't scale the mass
  	Tau.SetXYZM( Tau.Px()*(1+relative_shift), Tau.Py()*(1+relative_shift), Tau.Pz()*(1+relative_shift), Tau.M());
  }
  else {
    Tau.SetXYZM( Tau.Px()*(1+relative_shift), Tau.Py()*(1+relative_shift), Tau.Pz()*(1+relative_shift), Tau.M()*(1+relative_shift));
  }
}

// Implementation of cut-based veto ID Summer 16 (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Electron_ID_Working_Points_WP_de)
bool passedSummer16VetoId(const AC1B * analysisTree, int index){

  if(abs(analysisTree->electron_superclusterEta[index]) <= 1.479){

    if( analysisTree->electron_full5x5_sigmaietaieta[index] >= 0.0155 )                       return false;
    if( abs(analysisTree->electron_detaInSeed[index]) >= 0.00749 )                            return false;
    if( abs(analysisTree->electron_deltaphisuperclustertrack[index]) >= 0.228 )               return false;
    if( analysisTree->electron_he[index] >= 0.356 )                                           return false;
    if( analysisTree->electron_eaIsolation[index]/analysisTree->electron_pt[index] >= 0.175 ) return false;
    if( abs(analysisTree->electron_ooemoop[index]) >= 0.299 )                                 return false;
    if( analysisTree->electron_nmissinginnerhits[index] > 2 )                                 return false;
    if( analysisTree->electron_pass_conversion[index] == false )                              return false;
    
    return true;
  }
  else{

    if( analysisTree->electron_full5x5_sigmaietaieta[index] >= 0.037 )                        return false;
    if( abs(analysisTree->electron_detaInSeed[index]) >= 0.00895 )                            return false;
    if( abs(analysisTree->electron_deltaphisuperclustertrack[index]) >= 0.213 )               return false;
    if( analysisTree->electron_he[index] >= 0.211 )                                           return false;
    if( analysisTree->electron_eaIsolation[index]/analysisTree->electron_pt[index] >= 0.159 ) return false;
    if( abs(analysisTree->electron_ooemoop[index]) >= 0.15 )                                  return false;
    if( analysisTree->electron_nmissinginnerhits[index] > 3 )                                 return false;
    if( analysisTree->electron_pass_conversion[index] == false )                              return false;
    
    return true;
  }
}


bool SafeRatio(double denominator){

  if(denominator<0.001) return false;
  else                  return true;

}

bool passedAllMetFilters(const AC1B *analysisTree, std::vector<TString> met_filters){
  bool passed = true;
  unsigned int nfilters = met_filters.size();
  for (std::map<string,int>::iterator it = analysisTree->flags->begin(); it != analysisTree->flags->end(); ++it) {
    TString filter_name(it->first);
    for (unsigned int iFilter = 0; iFilter < nfilters; ++iFilter){
      if (filter_name.Contains(met_filters[iFilter]) && it->second == 0) {
        passed = false;
        break;
      }
    }      
  }
  return passed;
}


//fill the otree with the weights
void fill_weight(const AC1B * analysisTree, Synch17Tree *otree, PileUp *PUofficial, bool isData){
  
  otree->mcweight = 1.;
  otree->puweight = 0;

  if(!isData)
    otree->mcweight = analysisTree->genweight;

  otree->puweight = float(PUofficial->get_PUweight(double(analysisTree->numtruepileupinteractions)));

  otree->trigweight_1 = 0;
  otree->trigweight_2 = 1;
  otree->idisoweight_1 = 0;
  otree->idisoweight_2 = 1;
  otree->trigweight_antiiso_1 = 0;
  otree->idisoweight_antiiso_1 = 0;
  otree->trkeffweight=1;
  otree->effweight = 0;
  otree->weight = 1;
  otree->gen_noutgoing = analysisTree->genparticles_noutgoing;
}

// defines decay mode of the Z 
// 1 : Z->ee
// 2 : Z->mumu
// 3 : Z->tautau
int ZDecay(const AC1B * analysisTree) {
  int zdecay = 0;
  unsigned int ngen = analysisTree->genparticles_count;
  unsigned int nPromptMuons = 0;
  unsigned int nPromptElectrons = 0;
  for (unsigned int igen =0; igen<ngen; ++igen) {
    bool isMuon = TMath::Abs(analysisTree->genparticles_pdgid[igen])==13;
    bool isElectron = TMath::Abs(analysisTree->genparticles_pdgid[igen])==11;
    //    if (TMath::Abs(analysisTree.genparticles_pdgId[igen])==23) isZfound = true; 
    //    bool isDaughterOfZ = analysisTree.genparticles_info[igen] == 1;
    if (analysisTree->genparticles_fromHardProcess[igen]&&analysisTree->genparticles_status[igen]==1) {
      if (isMuon) nPromptMuons++;
      if (isElectron) nPromptElectrons++;
    }
  }
  if (nPromptElectrons==2)
    zdecay = 1;
  else if (nPromptMuons==2)
    zdecay = 2;    
  else
    zdecay = 3;

  return zdecay;

}
