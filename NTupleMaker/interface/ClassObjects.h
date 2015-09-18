#ifndef CLOBJECTS_H
#define CLOBJECTS_H

#include "NtupleTools3.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>

// shouldn't be using 'using namespace' in header files
//using namespace std;


// Extend the EasyChain class with getObjects functions
class ObjectChain: EasyChain {

};

// Object classes

// common object variable class ~ simple TLorentzVector -> inherit from TLV??
class ParticleObject: public TLorentzVector{
//private:
//    ParticleObject() {}

public:
// inherit baseclass constructors (C++11)
    using TLorentzVector::TLorentzVector;
    ParticleObject() {}
// conflict with default constructor?
    ParticleObject(Double_t pt, Double_t eta, Double_t phi, Double_t mass){
	SetPtEtaPhiM(pt,eta,phi,mass);
    }

};

/*
class ParticleObject{

private:
// constructor
    ParticleObject() {}

public:
    TLorentzVector vect;

    void SetVect(Double_t pt, Double_t eta, Double_t phi, Double_t mass){
	vect.SetPtEtaPhiM(pt,eta,phi,mass);
    }

// for fast access to kinematic variables
    Double_t Pt(){ return vect.Pt() }
    Double_t Eta(){ return vect.Eta() }
    Double_t Phi(){ return vect.Phi() }

};

// constuctor as TLV SetPtEtaPhiM
ParticleObject::ParticleObject(Double_t pt, Double_t eta, Double_t phi, Double_t mass){
    SetVect(pt,eta,phi,mass);
}

*/

// Lepton
class Lepton: public ParticleObject{
public:
    using ParticleObject::ParticleObject;

    Int_t pdgID;
    Double_t relIso03;
    Int_t tightID;
    Int_t charge;
};


// Generator Lepton
class GenLepton: public ParticleObject{
public:
    using ParticleObject::ParticleObject;

    int status;
};

// a bit more generic, generator partices
class GenParticle: public ParticleObject{
public:
    using ParticleObject::ParticleObject;

    int pdgid;
    int motherid;
    int grandmaid;
};

// Jet
class FatJet: public ParticleObject{
private:
// determine btag in constructor?
//    Jet() { if(btagCSV > 0.676) btag = true; }

public:
  Double_t prunedMass;
  Double_t trimmedMass;
  Double_t filteredMass;
  Double_t  tau1;
  Double_t tau2;
  Double_t tau3;
  Double_t topMass;
  Double_t minMass; 
  Double_t nSubJets;

  bool topTagged;


};


class Jet: public ParticleObject{
private:
public:
    bool btag;
    Double_t btagCSV;
};


// MET
class MET: public ParticleObject{
public:

//    Double_t met(){ return TLorentzVector::Pt() }
    Double_t met(){ return Pt(); }
};

class GetObjects{                                                                                              
 public:                                                                       
  void GetJets(EasyChain * tree);
  void GetFatJets(EasyChain * tree);
  void GetMET(EasyChain * tree);
  void GetGenMET(EasyChain * tree);

  void GetMETnoPU(EasyChain * tree);
  void GetLeptons(EasyChain * tree);
  void GetGenLeptons(EasyChain * tree);
  void GetGenParticles(EasyChain * tree);
  void GetGenLeptonsFromTau(EasyChain * tree);
  void GetGenTaus(EasyChain * tree);
  void GetEventCleaning(EasyChain * tree);
  void GetKinVariables();
  //void GetKinVariables(std::vector<Lepton> SelectedLep, std::vector<Jet> goodJet, TLorentzVector MET);

   std::vector<Jet> goodJet;
   std::vector<Jet> goodBJet;
   std::vector<FatJet> goodFatJet;
   std::vector<FatJet> goodTopTagJet;
   std::vector<FatJet> goodWTagJet;
   std::vector<FatJet> goodWmassTagJet;
  
   TLorentzVector MET;
   TLorentzVector genMET;
   TLorentzVector METnoPU;
  
   std::vector<Lepton> goodLep;
   std::vector<Lepton> goodEl;
   std::vector<Lepton> goodMu;
  
   std::vector<Lepton> SelectedLep;

   std::vector<Lepton> softLep;
   std::vector<Lepton> softEl;
   std::vector<Lepton> softMu;

   std::vector<Lepton> vetoLep;
   std::vector<Lepton> vetoEl;
   std::vector<Lepton> vetoMu;
  
   std::vector<Lepton> SoftvetoLep;
   std::vector<Lepton> SoftvetoEl;
   std::vector<Lepton> SoftvetoMu;

   std::vector<GenLepton> genLep;
   std::vector<GenLepton> genEl;
   std::vector<GenLepton> genMu;
   std::vector<GenLepton> genTau;

   std::vector<GenLepton> genLepFromTau;
   std::vector<GenLepton> genElFromTau;
   std::vector<GenLepton> genMuFromTau;

   std::vector<GenParticle> genPart;    


  // objects number can be aslo detemined as object.size()
  
   Int_t nLepGood;
   Int_t nMuGood;
   Int_t nElGood;
  
   Int_t nSoftLepGood;
   Int_t nSoftMuGood;
   Int_t nSoftElGood;

   Int_t nLepVeto;
   Int_t nElVeto;
   Int_t nMuVeto;
  
   Int_t nSoftLepVeto;
   Int_t nSoftElVeto;
   Int_t nSoftMuVeto;


   Int_t nJetGood;
   Int_t nFatJetGood;
   Int_t nTopTagJetGood;
   Int_t nWmassTagJetGood;
   Int_t nWTagJetGood;
   Int_t nBJetGood;
  
   Int_t nGenTau;
   Int_t nGenLep;
   Int_t nGenPart;
   Int_t nGenLepFromTau;
   
   
   bool Noise;
   bool goodvertex;

   Double_t HT40;
   Double_t ST;  
   Double_t DelPhiWLep;  
   Double_t DelPhiMetLep;  
   Double_t minDelPhibMet;  
   Double_t minDelPhiJMet;  
   Double_t minDelPhibW;  
   Double_t minDelPhibLep;  
   
   Double_t MTMetLep;
   Double_t MTbMet;
   Double_t MTbLep;
   Double_t MTbW;
   Double_t  DelRJMet0;
   Double_t  DelRJMet1;
   Double_t  DelRJMet2;
   Double_t  DelRJMet01;
   Double_t  minDelRJLep;
   Double_t  minDelRbL ;


};


#endif
