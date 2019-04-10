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
#include "TF1.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "TLorentzVector.h"

#include "TRandom.h"


bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags, bool isData) {

  bool passed = true;
  unsigned int nFlags = metFlags.size();
  //  std::cout << "MEt filters : " << std::endl;
  for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
    TString flagName(it->first);
    //    std::cout << it->first << " : " << it->second << std::endl;
    for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
      if (flagName.Contains(metFlags[iFilter])) {
	if(flagName.Contains("eeBadScFilter") && !isData) continue;
	if (it->second==0) {
	  passed = false;
	  break;
	}
      }
    }
  }
  //  std::cout << "Passed : " << passed << std::endl;
  return passed;

}

double dPhiFromLV(TLorentzVector v1, TLorentzVector v2) {

  return dPhiFrom2P(v1.Px(),v1.Py(),v2.Px(),v2.Py());

}

void ReadJson(std::vector<Period> &periods, string fullPathToJsonFile)
{
   std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
   if (inputFileStream.fail()) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
   }
   for(std::string s; std::getline(inputFileStream, s); )
      {
         periods.push_back(Period());
         std::stringstream ss(s);
         ss >> periods.back();
      }
}

void iniPU(PileUp *PUofficial, TFile *filePUOfficial_data, TFile *filePUOfficial_MC, TString samplenameForPUHist)
{
   // Official PU reweighting
   TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
   TString NamePUHistMC = samplenameForPUHist + "pileup";
   TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get(NamePUHistMC);
   if( !PUOfficial_mc) {
      cout<<endl<<"MC pileup histogram "<<NamePUHistMC<<" does not exist in root file. Exiting."<<endl<<endl;
      exit(-1);
   }
   PUofficial->set_h_data(PUOfficial_data);
   PUofficial->set_h_MC(PUOfficial_mc);
}


void iniTriggerEfficiencies(TFile * trigEffFile,map<int,TGraphAsymmErrors*> &map_trigEffData,map<int,TGraphAsymmErrors*> &map_trigEffMC){
   
   if(trigEffFile->IsZombie()){
      cout<<"Trigger file does not exists. Exiting."<<endl;
      exit(-1);
   }
   TIter next(trigEffFile->GetListOfKeys());
   TKey *key = 0;
   
   while ((key = (TKey*)next()))
      {
         TClass *c = gROOT->GetClass(key->GetClassName());
         if (!c->InheritsFrom("TGraphAsymmErrors")) continue;
         TGraphAsymmErrors *g = (TGraphAsymmErrors*) key->ReadObj();
         TString gName = g->GetName();
         Ssiz_t pos    = gName.First("To");
         Int_t  len    = gName.Length();
         int upperBound   = atoi( (TString) gName(pos+2,len));
         if(gName.Contains("data"))    map_trigEffData.insert( std::make_pair( upperBound , (TGraphAsymmErrors*) g->Clone() ) );
         else if(gName.Contains("mc")) map_trigEffMC.insert(   std::make_pair( upperBound , (TGraphAsymmErrors*) g->Clone() ) );
         delete g;
      }
   delete key;
}
 

bool GoodRunSelection(int n, int lum,std::vector<Period> &periods){
   
   bool lumi = false;
   std::string num = std::to_string(n);
   std::string lnum = std::to_string(lum);
   for(const auto& a : periods)
      {
         if ( num.c_str() ==  a.name ) {
            for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
               if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
            }
            auto last = std::prev(a.ranges.end());
            if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
         }
      }
   return lumi;
}  

void AccessTriggerInfo(AC1B &analysisTree, TString HLTFilterName, unsigned int &nHLTFilter, bool &isHLTFilter)
{
   for (unsigned int i=0; i<analysisTree.run_hltfilters->size(); ++i) {
      TString HLTFilter(analysisTree.run_hltfilters->at(i));
      if (HLTFilter==HLTFilterName) {
         nHLTFilter = i;
         isHLTFilter = true;
      }
   }
   
   if (!isHLTFilter) {
      std::cout << "HLT filter " << HLTFilterName << " not found" << std::endl;
      exit(-1);
   }
}

bool PassesMuonSelection(AC1B &analysisTree, unsigned int imuon, const float ptMuCut, const float etaMuCut,const float dxyMuCut,const float dzMuCut,const bool isDRIso03, const float isoMuCut, float &relIso, const string era)
{
   bool passes =true;
   
   if (analysisTree.muon_pt[imuon]<ptMuCut) passes=false;
   if (fabs(analysisTree.muon_eta[imuon])>etaMuCut) passes=false;
   
   // source  : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
   bool passedId = false;
   if (era=="2017" || era=="2018" || era=="2016") passedId = analysisTree.muon_isMedium[imuon];
   else {
      std::cout<<"Muon Id not defined for era "<<era<<std::endl;
      exit(-1);
   }
   if (!passedId) passes=false;
   
   bool passedIpCuts = 
      fabs(analysisTree.muon_dxy[imuon]) < dxyMuCut &&
      fabs(analysisTree.muon_dz[imuon]) < dzMuCut;
   if (!passedIpCuts) passes=false;
   
   float absIso = 0; 
   if (isDRIso03) {
      absIso = analysisTree.muon_r03_sumChargedHadronPt[imuon];
      float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[imuon] + 
         analysisTree.muon_r03_sumPhotonEt[imuon] - 0.5*analysisTree.muon_r03_sumPUPt[imuon];
      neutralIso = TMath::Max(float(0),neutralIso);
      absIso += neutralIso;
   }
   else {
      absIso = analysisTree.muon_chargedHadIso[imuon];
      float neutralIso = analysisTree.muon_neutralHadIso[imuon] +
         analysisTree.muon_photonIso[imuon] -
         0.5*analysisTree.muon_puIso[imuon];
      neutralIso = TMath::Max(float(0),neutralIso);
      absIso += neutralIso;
   } 
   relIso = absIso/analysisTree.muon_pt[imuon];
   bool passedIso = relIso < isoMuCut;
   if (!passedIso) passes=false;

   return passes;
}

bool TriggerMatching(AC1B &analysisTree, Float_t eta, Float_t phi, unsigned int nFilter)
{
   bool trigMatch = false;
   for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
      float dRtrig = deltaR(eta,phi,analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
      if (dRtrig>0.5) continue;
      if (analysisTree.trigobject_filters[iT][nFilter]) trigMatch = true;
      
   }
   
   return trigMatch;
}

bool PassesElectronSelection(AC1B &analysisTree, unsigned int ielec , const float dxyEleCut, const float dzEleCut, const float isoEleCut, const float ptEleCut, const float etaEleCut, const string era)
{
   bool passes = false;
   
   // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Fall17v2
   bool passedId = false;
   if (era == "2017") passedId =
                         analysisTree.electron_cutId_veto_Summer16[ielec] &&
                         analysisTree.electron_pass_conversion[ielec] &&
                         analysisTree.electron_nmissinginnerhits[ielec] <= 1;           
   else if (era == "2018")  passedId =
                         analysisTree.electron_cutId_veto_Fall17[ielec] &&
                         analysisTree.electron_pass_conversion[ielec] &&
                         analysisTree.electron_nmissinginnerhits[ielec] <= 1;
   else if (era == "2016")  passedId =
                         analysisTree.electron_cutId_veto_Fall17V2[ielec] &&
                         analysisTree.electron_pass_conversion[ielec] &&
                         analysisTree.electron_nmissinginnerhits[ielec] <= 1;
   else {
      std::cout<<"Electron Id not defined for era "<<era<<std::endl;
      exit(-1);
   }

   bool passedIpCuts = 
      fabs(analysisTree.electron_dxy[ielec]) < dxyEleCut &&
      fabs(analysisTree.electron_dz[ielec]) < dzEleCut;
   
   float absIso = analysisTree.electron_r03_sumChargedHadronPt[ielec];
   float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[ielec] + analysisTree.electron_r03_sumPhotonEt[ielec] - 0.5*analysisTree.electron_r03_sumPUPt[ielec];
   neutralIso = TMath::Max(float(0),neutralIso);
   absIso += neutralIso;
   float relIso = absIso/analysisTree.electron_pt[ielec];
   bool passedIso = relIso < isoEleCut;
  
   if (analysisTree.electron_pt[ielec]>ptEleCut && 
       fabs(analysisTree.electron_eta[ielec])<etaEleCut &&
       passedId && passedIpCuts && passedIso) passes = true; 
      
   return passes;
}

bool OverlapWithLepton(AC1B &analysisTree, unsigned int ijet, std::vector<unsigned int> muonIndexes, std::vector<unsigned int> eleIndexes)
{
   // checking overlap with muons
   bool overlapWithLepton = false;
   
   for (unsigned int iMu=0; iMu<muonIndexes.size(); ++iMu) {
      unsigned int indexMu = muonIndexes.at(iMu);
      float dRJetLep = deltaR(analysisTree.muon_eta[indexMu],analysisTree.muon_phi[indexMu],
                              analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet]);
      if (dRJetLep<0.4) {
         overlapWithLepton = true;
         break;
      }
   }
   // checking overlap with electrons
   for (unsigned int iEle=0; iEle<eleIndexes.size(); ++iEle) {
      unsigned int indexEle = eleIndexes.at(iEle);
      float dRJetLep = deltaR(analysisTree.electron_eta[indexEle],analysisTree.electron_phi[indexEle],
                              analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet]);
      if (dRJetLep<0.4) {
         overlapWithLepton = true;
         break;
      }
   }
   return overlapWithLepton;
   
}

int MatchingMuon(AC1B &analysisTree, unsigned int itau, std::vector<unsigned int> muonIndexes)
{
   int matchedMuIndex = -1;
   float dRmin = 0.4;
   for (unsigned int im=0; im<muonIndexes.size(); im++) {
      unsigned int indexMu = muonIndexes.at(im);
      float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
                        analysisTree.muon_eta[indexMu],analysisTree.muon_phi[indexMu]);
      if (dR<dRmin) {
         dRmin = dR;
         matchedMuIndex = int(indexMu);
      }
   }
   return matchedMuIndex;
}

int MatchingElectron(AC1B &analysisTree, unsigned int itau, std::vector<unsigned int> eleIndexes)
{
   int matchedEleIndex = -1;
   float dRmin = 0.4;
   for (unsigned int ie=0; ie<eleIndexes.size(); ie++) {
      unsigned int indexEle = eleIndexes.at(ie);
      float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
                        analysisTree.electron_eta[indexEle],analysisTree.electron_phi[indexEle]);
      if (dR<dRmin) {
         dRmin = dR;
         matchedEleIndex = int(indexEle);
      }
   }
   return matchedEleIndex;
}

bool LooseTauSelection(AC1B &analysisTree, unsigned int itau){
   
   bool passes = true;
   
   if (fabs(analysisTree.tau_eta[itau])>2.4) passes = false; // loose eta cut
   if (analysisTree.tau_pt[itau]<20.) passes = false; // loose pt cut
   
   float dZ = fabs(analysisTree.tau_vertexz[itau]-analysisTree.primvertex_z);
   if (dZ>1e-4) passes = false; // dz criterion
   
   bool foundByDecayMode = analysisTree.tau_decayModeFindingNewDMs[itau]>0.5 || analysisTree.tau_decayModeFinding[itau]>0.5;
   if (!foundByDecayMode) passes = false; // DM finding

   return passes;

}

UInt_t FindTauGenMatch(AC1B &analysisTree, Float_t tauEta, Float_t tauPhi)
{
   UInt_t tauGenMatch = -1;
   float minDR = 0.2;
   for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
      TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                          analysisTree.genparticles_py[igen],
                                          analysisTree.genparticles_pz[igen],
                                          analysisTree.genparticles_e[igen]);
      float ptGen = genLV.Pt();
      bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
      bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
      bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
      bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
      bool isAnyType = type1 || type2 || type3 || type4;
      
      if (isAnyType && analysisTree.genparticles_status[igen]==1) {
         float etaGen = genLV.Eta();
         float phiGen = genLV.Phi();
         float dR = deltaR(tauEta,tauPhi,
                           etaGen,phiGen);
         if (dR<minDR) {
            minDR = dR;
            if (type1) tauGenMatch_ = 1;
            else if (type2) tauGenMatch = 2;
            else if (type3) tauGenMatch = 3;
            else if (type4) tauGenMatch = 4;
         }
      }
   }
   return tauGenMatch;
 }

int MatchingJet(AC1B &analysisTree, TLorentzVector lorentzVectorTau, bool &jetFound, TLorentzVector &lorentzVectorTauJet)
{
   float dRmin = 0.4;
   int indexMatchingJet  = -1;
   for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
      TLorentzVector lorentzVectorJ; lorentzVectorJ.SetXYZT(analysisTree.pfjet_px[ijet],
                                                            analysisTree.pfjet_py[ijet],
                                                            analysisTree.pfjet_pz[ijet],
                                                            analysisTree.pfjet_e[ijet]);
      float drJetTau = deltaR(lorentzVectorJ.Eta(),lorentzVectorJ.Phi(),
                              lorentzVectorTau.Eta(),lorentzVectorTau.Phi());
      
      if (drJetTau<dRmin) {
         dRmin = drJetTau;
         jetFound = true;
         indexMatchingJet = ijet;
         lorentzVectorTauJet = lorentzVectorJ;
      }
      
   }
   return indexMatchingJet;
}



                  
             
