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
#include "TSystem.h"
#include "TRandom.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

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

bool AccessTriggerInfo(AC1B &analysisTree, TString HLTFilterName, unsigned int &nHLTFilter)
{
   bool isHLTFilter = false;
   
   for (unsigned int i=0; i<analysisTree.run_hltfilters->size(); ++i) {
      TString HLTFilter(analysisTree.run_hltfilters->at(i));
      if (HLTFilter==HLTFilterName) {
         nHLTFilter = i;
         isHLTFilter = true;
      }
   }
   return isHLTFilter;
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

bool TriggerMatching(AC1B &analysisTree, Float_t eta, Float_t phi, unsigned int nFilter, float deltaRTrigMatch = 0.5)
{
   bool trigMatch = false;
   for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
      float dRtrig = deltaR(eta,phi,analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
      if (dRtrig> deltaRTrigMatch) continue;
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
   else if (era == "2016")  passedId =
                         analysisTree.electron_cutId_veto_Fall17V2[ielec] &&
                         analysisTree.electron_pass_conversion[ielec] &&
                         analysisTree.electron_nmissinginnerhits[ielec] <= 1;
   else if (era == "2018")  passedId =
                         analysisTree.electron_cutId_veto_Fall17[ielec] &&
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
            if (type1) tauGenMatch = 1;
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

float getEffectiveArea(float eta) {
    float effArea = 0.1566;
    float absEta = fabs(eta);
    if (absEta<1.0) effArea = 0.1566;
    else if (absEta < 1.4790) effArea = 0.1626;
    else if (absEta < 2.0) effArea = 0.1073;
    else if (absEta < 2.2) effArea = 0.0854;
    else if (absEta < 2.3) effArea = 0.1051;
    else if (absEta < 2.4) effArea = 0.1204;
    else if (absEta < 5.0) effArea = 0.1524;
    return effArea;

  } 

void LoadCrystalBallEfficiencyClass( TString pathToCrystalLib){
   int openSuccessful = gSystem->Load( pathToCrystalLib );
   if (openSuccessful !=0 ) {
      cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
      exit( -1 );
   }
}                  
             

void LoadZpTMassWeights(TFile *fileZMassPtWeights, TString ZMassPtWeightsFileName, TH2D* histZMassPtWeights, TString ZMassPtWeightsHistName){
 
   if (fileZMassPtWeights->IsZombie()) {
      std::cout << "File /src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
      exit(-1);
   }
   if (histZMassPtWeights==NULL) {
      std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file /src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName
                << std::endl;
      exit(-1);
   }
}

void  FillHistGenWeights( TTree *_inittree, TH1D* histWeightsH, bool isData){

   if (_inittree!=NULL) {
      Float_t genweight;
      Float_t genweight_;
      if (!isData)
         _inittree->SetBranchAddress("genweight",&genweight);
      Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
      std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
      for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
         _inittree->GetEntry(iEntry);
         if (isData)
            histWeightsH->Fill(0.,1.);
         else{
           if (genweight<0) genweight_ = -1;
           else             genweight_ = 1;
           histWeightsH->Fill(0.,genweight_);
         }
      }
   }
}

bool FillHistInputEvents(TFile *file_, TH1D*inputEventsH){
   
   bool hasInputEvents = true;
   TH1D * histoInputEvents = NULL;
   histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
   if (histoInputEvents==NULL) hasInputEvents = false;
   int NE = int(histoInputEvents->GetEntries());
   std::cout << "      number of input events    = " << NE << std::endl;
   for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
   return hasInputEvents;
}

Float_t GetEmbeddedWeight(std::vector<TLorentzVector> promptTausFirstCopy, RooWorkspace *correctionWS_embedded){
   
   Float_t embeddedWeight = 1.0;
   double gt1_pt  = promptTausFirstCopy[0].Pt();
   double gt1_eta = promptTausFirstCopy[0].Eta();
   double gt2_pt  = promptTausFirstCopy[1].Pt();
   double gt2_eta = promptTausFirstCopy[1].Eta();
   correctionWS_embedded->var("gt_pt")->setVal(gt1_pt);
   correctionWS_embedded->var("gt_eta")->setVal(gt1_eta);
   double id1_embed = correctionWS_embedded->function("m_sel_idEmb_ratio")->getVal();
   correctionWS_embedded->var("gt_pt")->setVal(gt2_pt);
   correctionWS_embedded->var("gt_eta")->setVal(gt2_eta);
   double id2_embed = correctionWS_embedded->function("m_sel_idEmb_ratio")->getVal();
   correctionWS_embedded->var("gt1_pt")->setVal(gt1_pt);
   correctionWS_embedded->var("gt1_eta")->setVal(gt1_eta);
   correctionWS_embedded->var("gt2_pt")->setVal(gt2_pt);
   correctionWS_embedded->var("gt2_eta")->setVal(gt2_eta);
   double trg_embed = correctionWS_embedded->function("m_sel_trg_ratio")->getVal();
   embeddedWeight = id1_embed * id2_embed * trg_embed;

   return embeddedWeight;
}

Float_t GetZptMassWeight(Float_t bosonMass, Float_t bosonPt, TH2D * histZMassPtWeights){
   
   Float_t zptmassweight = 1;
   if (bosonMass>50.0) {
      float bosonMassX = bosonMass;
      float bosonPtX = bosonPt;
      if (bosonMassX>histZMassPtWeights->GetXaxis()->GetXmax()) {
         bosonMassX = histZMassPtWeights->GetXaxis()->GetXmax() - 0.00001;
      }
      double minpt;
      if (histZMassPtWeights->GetYaxis()->GetXmin() > 1.) minpt = histZMassPtWeights->GetYaxis()->GetXmin();
      else minpt =1.;
      if (bosonPtX<minpt) bosonPtX = minpt+0.00001;
      if (bosonPtX>histZMassPtWeights->GetYaxis()->GetXmax()) bosonPtX = histZMassPtWeights->GetYaxis()->GetXmax()- 0.00001;

      zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
                                                        histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
   }
   return zptmassweight;
}


void SearchForBtagDiscriminant(AC1B &analysisTree, TString BTagDiscriminator1, TString BTagDiscriminator2, TString BTagDiscriminator3, unsigned int &nBTagDiscriminant1, unsigned int &nBTagDiscriminant2, unsigned int &nBTagDiscriminant3, TString era){

   for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
      TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
      if (discr == BTagDiscriminator1)
         nBTagDiscriminant1 = iBTag;
      if (era!="2016" && discr == BTagDiscriminator2)
         nBTagDiscriminant2 = iBTag;
      if (era!="2016" && discr == BTagDiscriminator3)
         nBTagDiscriminant3 = iBTag;
   }
}

void SelectMuonElePair(AC1B &analysisTree, vector<int> muons, vector<int> electrons, bool isMuonIsoR03, bool isElectronIsoR03, float dRleptonsCut, float ptMuonHighCut, float ptElectronHighCut, int &electronIndex, int &muonIndex, float &isoMuMin, float &isoEleMin, TString era){
         
   for (unsigned int im=0; im<muons.size(); ++im) {
      unsigned int mIndex  = muons.at(im);
      float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
      float photonIsoMu = analysisTree.muon_photonIso[mIndex];
      float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
      float puIsoMu = analysisTree.muon_puIso[mIndex];
      if (isMuonIsoR03) {
         neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
         photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
         chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
         puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
      }
      float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
      neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
      float absIsoMu = chargedHadIsoMu + neutralIsoMu;
      float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];
      
      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
         unsigned int eIndex = electrons.at(ie);
         float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
                           analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);
         
         if (dR<dRleptonsCut) continue;
         
         bool trigMatch =
            (analysisTree.muon_pt[mIndex]>ptMuonHighCut) ||
            (analysisTree.electron_pt[eIndex]>ptElectronHighCut);
         if (!trigMatch) continue;
         
         float absIsoEle; 
         float relIsoEle;
         if (era=="2016"){
            float neutralHadIsoEle = analysisTree.electron_neutralHadIso[eIndex];
            float photonIsoEle = analysisTree.electron_photonIso[eIndex];
            float chargedHadIsoEle = analysisTree.electron_chargedHadIso[eIndex];
            float puIsoEle = analysisTree.electron_puIso[eIndex];
            if (isElectronIsoR03) {
               neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
               photonIsoEle = analysisTree.electron_r03_sumPhotonEt[eIndex];
               chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
               puIsoEle = analysisTree.electron_r03_sumPUPt[eIndex];
            }
            float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
            neutralIsoEle = TMath::Max(float(0),neutralIsoEle);
            absIsoEle =  chargedHadIsoEle + neutralIsoEle;
            relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];
         }
         else {
            float rhoNeutral = analysisTree.rho;
            float  eA = getEffectiveArea( fabs(analysisTree.electron_superclusterEta[eIndex]) );
            absIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex] +
               TMath::Max(0.0f,analysisTree.electron_r03_sumNeutralHadronEt[eIndex]+analysisTree.electron_r03_sumPhotonEt[eIndex]-eA*rhoNeutral);
            relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];
         }
         if (int(mIndex)!=muonIndex) {
            if (relIsoMu==isoMuMin) {
               if (analysisTree.muon_pt[mIndex]>analysisTree.muon_pt[muonIndex]) {
                  isoMuMin  = relIsoMu;
                  muonIndex = int(mIndex);
                  isoEleMin = relIsoEle;
                  electronIndex = int(eIndex);
               }
            }
            else if (relIsoMu<isoMuMin) {
               isoMuMin  = relIsoMu;
               muonIndex = int(mIndex);
               isoEleMin = relIsoEle;
               electronIndex = int(eIndex);
            }
         }
         else {
            if (relIsoEle==isoEleMin) {
               if (analysisTree.electron_pt[eIndex]>analysisTree.electron_pt[electronIndex]) {
                  isoEleMin = relIsoEle;
                  electronIndex = int(eIndex);
               }
            }
            else if (relIsoEle<isoEleMin) {
               isoEleMin = relIsoEle;
               electronIndex = int(eIndex);
            }
         }
      }
   }
}


bool ElectronVeto(AC1B &analysisTree, int electronIndex, float ptVetoElectronCut, float etaVetoElectronCut, float dxyVetoElectronCut, float dzVetoElectronCut, TString era, bool applyVetoElectronId, bool isElectronIsoR03, float isoVetoElectronCut){
 
   bool foundExtraElectron = false;
   for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
      if (int(ie)==electronIndex) continue;
      if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
      if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
      if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
      if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue; 
      bool electronMvaId = true;
      if (era=="2016") {
         electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie]>0.5; 
      }
      else if (era == "2017"){
         electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[ie]>0.5;
      }
      else if (era == "2018"){
         electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[ie]>0.5;
      }
      if (!electronMvaId&&applyVetoElectronId) continue;
      if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
 
      if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
      float relIsoEle;
      float absIsoEle;
      if (era=="2016"){
         float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
         float photonIsoEle = analysisTree.electron_photonIso[ie];
         float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
         float puIsoEle = analysisTree.electron_puIso[ie];
               if (isElectronIsoR03) {                                               
                  neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
                  photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
                  chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
                  puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
               }
               float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
               neutralIsoEle = TMath::Max(float(0),neutralIsoEle);
               absIsoEle =  chargedHadIsoEle + neutralIsoEle;
               relIsoEle = absIsoEle/analysisTree.electron_pt[ie];
            }
      else {
         float rhoNeutral = analysisTree.rho;
         float  eA = getEffectiveArea( fabs(analysisTree.electron_superclusterEta[ie]) );
         absIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie] +
            TMath::Max(0.0f,analysisTree.electron_r03_sumNeutralHadronEt[ie]+analysisTree.electron_r03_sumPhotonEt[ie]-eA*rhoNeutral);
         relIsoEle = absIsoEle/analysisTree.electron_pt[ie];
      }
      if (relIsoEle>isoVetoElectronCut) continue;
 
      foundExtraElectron = true;
   }
   return foundExtraElectron;
}

bool MuonVeto(AC1B &analysisTree, int muonIndex, float ptVetoMuonCut, float etaVetoMuonCut, float dxyVetoMuonCut, float dzVetoMuonCut, bool applyICHEPMuonId, bool applyVetoMuonId, bool isMuonIsoR03, float isoVetoMuonCut){
   
   bool foundExtraMuon = false;
   for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
      if (int(im)==muonIndex) continue;
      if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
      if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
      if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
      if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
      bool muonId = analysisTree.muon_isMedium[im];
      if (applyICHEPMuonId) muonId = analysisTree.muon_isICHEP[im];
      if (!muonId&&applyVetoMuonId) continue;
      float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
      float photonIsoMu = analysisTree.muon_photonIso[im];
      float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
      float puIsoMu = analysisTree.muon_puIso[im];
      if (isMuonIsoR03) {
         neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
         photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
         chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
         puIsoMu = analysisTree.muon_r03_sumPUPt[im];
      }
      float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
      neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
      float absIsoMu = chargedHadIsoMu + neutralIsoMu;
      float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
      if (relIsoMu>isoVetoMuonCut) continue;
      foundExtraMuon = true;
   }
   
   return foundExtraMuon;
}
