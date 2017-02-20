// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef LeptonScaleSys_h
#define LeptonScaleSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)

#include "TH2D.h"

#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics.h>

using namespace utils;

class LeptonScaleSys : public Systematics {
public:
  
  LeptonScaleSys(){};
  
  LeptonScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "lScale";
    
    this->Init(cenTree);
  };
  
  virtual ~LeptonScaleSys(){
    if(sf_up != 0)
      delete sf_up;
    if(sf_down != 0)
      delete sf_down;

    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      delete it->second;
  };
  
  virtual void Eval(utils::channel ch){
    this->Central();

    this->ScaleUp(ch);
    this->ScaleDown(ch);
  };

  virtual void Write(){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      it->second->Write();
  };
  
  void SetSvFitVisPtResolution(TFile* f){
    svFit_visPtResolution = f;
  };
  
protected:
  virtual void Init(Spring15Tree* c){
    cenTree = c;

    this->InitSF();

    this->InitTree("Up");
    this->InitTree("Down");
  };

  virtual void InitSF() = 0;
  virtual void ScaleUp(utils::channel ch) = 0;
  virtual void ScaleDown(utils::channel ch) = 0;

  virtual void InitTree(const char* shift){
    std::cout<<label+shift<<std::endl;
    outTree[shift] = cenTree->fChain->CloneTree(0);
    outTree[shift]->SetName(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetTitle(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetDirectory(cenTree->fChain->GetDirectory());
  };
  
  virtual void Central(){
    lep1 = TLorentzVector(); lep1.SetXYZM(cenTree->pt_1 * cos(cenTree->phi_1),
					  cenTree->pt_1 * sin(cenTree->phi_1),
					  cenTree->pt_1 * sinh(cenTree->eta_1),
					  cenTree->m_1);
    
    lep2 = TLorentzVector(); lep2.SetXYZM(cenTree->pt_2 * cos(cenTree->phi_2),
					  cenTree->pt_2 * sin(cenTree->phi_2),
					  cenTree->pt_2 * sinh(cenTree->eta_2),
					  cenTree->m_2);
  } 

  virtual void Fill(utils::channel ch, const char* shift){
    if (lep1.Pt() == lep1_scaled.Pt() && lep2.Pt() == lep2_scaled.Pt()){
      outTree[shift]->Fill();
      return;
    }

    // store cen values
    float pt_1_cen = cenTree->pt_1;
    float pt_2_cen = cenTree->pt_2;
    float met_cen = cenTree->met;

    float mt_1_cen = cenTree->mt_1;
    float pfmt_1_cen = cenTree->pfmt_1;    
    float mt_2_cen = cenTree->mt_2;
    float pfmt_2_cen = cenTree->pfmt_2;    

    float pt_tt_cen = cenTree->pt_tt;
    float pzetavis_cen = cenTree->pzetavis;
    float pzetamiss_cen = cenTree->pzetamiss;
    float pfpzetamiss_cen = cenTree->pfpzetamiss;    
    
    float m_vis_cen = cenTree->m_vis;

    float m_sv_cen = cenTree->m_sv;
    float pt_sv_cen = cenTree->pt_sv;
    float eta_sv_cen = cenTree->eta_sv;
    float phi_sv_cen = cenTree->phi_sv;
    float met_sv_cen = cenTree->met_sv;
    float mt_sv_cen = cenTree->mt_sv;

    // calc shifted values
    cenTree->pt_1 = lep1_scaled.Pt();
    cenTree->pt_2 = lep2_scaled.Pt();

    // central value of the met    
    TLorentzVector pfmetLV; pfmetLV.SetXYZT(cenTree->met * cos(cenTree->metphi),
					    cenTree->met * sin(cenTree->metphi),
					    0.,
					    TMath::Sqrt( cenTree->met * cos(cenTree->metphi) * cenTree->met * cos(cenTree->metphi) +
							 cenTree->met * sin(cenTree->metphi) * cenTree->met * sin(cenTree->metphi) ) ) ;

    TLorentzVector mvametLV; mvametLV.SetXYZT(cenTree->mvamet * cos(cenTree->mvametphi),
					      cenTree->mvamet * sin(cenTree->mvametphi),
					      0.,
					      TMath::Sqrt( cenTree->mvamet * cos(cenTree->mvametphi) * cenTree->mvamet * cos(cenTree->mvametphi) +
							   cenTree->mvamet * sin(cenTree->mvametphi) * cenTree->mvamet * sin(cenTree->mvametphi) ) ) ;    

    // propagate the lepton pt shift to the MET 
	pfmetLV.SetPx(pfmetLV.Px()- (lep1_scaled.Px()-lep1.Px()));
	pfmetLV.SetPy(pfmetLV.Py()- (lep1_scaled.Px()-lep1.Py()));

	mvametLV.SetPx(mvametLV.Px()- (lep1_scaled.Px()-lep1.Px()));
	mvametLV.SetPy(mvametLV.Py()- (lep1_scaled.Px()-lep1.Py()));

    // calc shifted values
	cenTree->met = pfmetLV.Pt();
	cenTree->mvamet = mvametLV.Pt();

    //cenTree->mt_1 = sqrt(2*lep1_scaled.Pt()*mvametLV.Pt()*(1.-cos(lep1_scaled.Phi()-mvametLV.Phi())));
    cenTree->pfmt_1 = sqrt(2*lep1_scaled.Pt()*pfmetLV.Pt()*(1.-cos(lep1_scaled.Phi()-pfmetLV.Phi())));
    //cenTree->mt_2 = sqrt(2*lep2_scaled.Pt()*mvametLV.Pt()*(1.-cos(lep2_scaled.Phi()-cenTree->mvametLV.Phi())));
    cenTree->pfmt_2 = sqrt(2*lep2_scaled.Pt()*pfmetLV.Pt()*(1.-cos(lep2_scaled.Phi()-pfmetLV.Phi())));
    cenTree->mt_1 = sqrt(2*lep1_scaled.Pt()*pfmetLV.Pt()*(1.-cos(lep1_scaled.Phi()-pfmetLV.Phi())));
    cenTree->mt_2 = sqrt(2*lep2_scaled.Pt()*pfmetLV.Pt()*(1.-cos(lep2_scaled.Phi()-pfmetLV.Phi())));
    
    TMatrixD covMET(2, 2);
    covMET[0][0] = cenTree->metcov00;
    covMET[1][0] = cenTree->metcov10;
    covMET[0][1] = cenTree->metcov01;
    covMET[1][1] = cenTree->metcov11;

    TLorentzVector dileptonLV = lep1_scaled + lep2_scaled;
    
    cenTree->pt_tt = (lep1_scaled+lep2_scaled+pfmetLV).Pt();
    cenTree->pzetavis = calc::pzetavis(lep1_scaled, lep2_scaled);
    //cenTree->pzetamiss = calc::pzetamiss( lep1_scaled, lep2_scaled, mvametLV);
    cenTree->pzetamiss = calc::pzetamiss( lep1_scaled, lep2_scaled, pfmetLV);
    cenTree->pfpzetamiss = calc::pzetamiss( lep1_scaled, lep2_scaled, pfmetLV);    
    
    cenTree->m_vis = dileptonLV.M();

    cenTree->m_sv = 0.;
    cenTree->pt_sv = 0.;
    cenTree->eta_sv = 0.;
    cenTree->phi_sv = 0.;
    cenTree->met_sv = 0.;
    cenTree->mt_sv = 0.;

    if(ch != UNKNOWN){
      std::shared_ptr<SVfitStandaloneAlgorithm> algo = calc::svFit(lep1_scaled, cenTree->tau_decay_mode_1, 
								   lep2_scaled, cenTree->tau_decay_mode_2, 
								   pfmetLV, covMET, 
								   ch, 
								   svFit_visPtResolution); 
      if (algo != 0){
	cenTree->m_sv = algo->mass();
	cenTree->pt_sv = algo->pt();
	cenTree->eta_sv = algo->eta();
	cenTree->phi_sv = algo->phi();      
	cenTree->met_sv = algo->fittedMET().Rho();
	cenTree->mt_sv = algo->transverseMass();
      }
    }
    outTree[shift]->Fill();


    // restore cen values
    cenTree->pt_1 = pt_1_cen;
    cenTree->pt_2 = pt_2_cen;
    cenTree->met = met_cen;
    cenTree->mt_1 = mt_1_cen;
    cenTree->pfmt_1 = pfmt_1_cen;    
    cenTree->mt_2 = mt_2_cen;
    cenTree->pfmt_2 = pfmt_2_cen;    

    cenTree->pt_tt = pt_tt_cen;
    cenTree->pzetavis = pzetavis_cen;
    cenTree->pzetamiss = pzetamiss_cen;
    cenTree->pfpzetamiss = pfpzetamiss_cen;    
    
    cenTree->m_vis = m_vis_cen;

    cenTree->m_sv = m_sv_cen;
    cenTree->pt_sv = pt_sv_cen;
    cenTree->eta_sv = eta_sv_cen;
    cenTree->phi_sv = phi_sv_cen;
    cenTree->met_sv = met_sv_cen;
    cenTree->mt_sv = mt_sv_cen;    
  }

  //std::map< std::string, Float_t >  mt_sv;
  
  TH2D* sf_up;
  TH2D* sf_down;  

  TLorentzVector lep1, lep2;
  TLorentzVector lep1_scaled, lep2_scaled;

  TFile* svFit_visPtResolution;

  std::map< std::string, TTree* >  outTree;
  //std::map< std::string, SpringTree* >  outTree;
};

class MuonScaleSys : public LeptonScaleSys { 
public:
  MuonScaleSys(){};
  
  MuonScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "CMS_scale_m_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~MuonScaleSys(){};
  
protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};

    const int eta_bins = 1;
    double eta_edges[eta_bins + 1] = {0., 2.5};

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);
  };  

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    
    if (ch == EMU)
      lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M());
    else if (ch == MUTAU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == MUMU){
      lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M());
      lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M());
    }
    
    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    
    if (ch == EMU)
      lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M());
    else if (ch == MUTAU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == MUMU){
      lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M());
      lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M());
    }
    
    this->Fill(ch, "Down");
  };
};


class ElectronScaleSys : public LeptonScaleSys { 
public:
  ElectronScaleSys(){};

  ElectronScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "eScale";
    
    this->Init(cenTree);
  };
  
  virtual ~ElectronScaleSys(){};

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 2;
    double eta_edges[eta_bins + 1] = {0., 1.479, 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.01);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), 0.025);   
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.01);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), 0.025);
  };
  
  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    
    if (ch == EMU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == ETAU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == EE){
      lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M());
      lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M());
    }

    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    
    if (ch == EMU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == MUTAU)
      lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M());
    else if (ch == MUMU){
      lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M());
      lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M());
    }
    
    this->Fill(ch, "Down");
  };
};

class ElectronEBScaleSys : public ElectronScaleSys {
public:
  ElectronEBScaleSys(){};
  
  ElectronEBScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "CMS_scale_eEB_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~ElectronEBScaleSys(){};

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 2;
    double eta_edges[eta_bins + 1] = {0., 1.5, 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.01);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), 0.);   
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.01);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), 0.);
  };
};

class ElectronEEScaleSys : public ElectronScaleSys {
public:
  ElectronEEScaleSys(){};
  
  ElectronEEScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "CMS_scale_eEE_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~ElectronEEScaleSys(){};

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 2;
    double eta_edges[eta_bins + 1] = {0., 1.5, 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), 0.025);   
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.00);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), 0.025);
  };
};


class TauScaleSys : public LeptonScaleSys { 
public:
  TauScaleSys(){};
  
  TauScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "CMS_shape_t_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~TauScaleSys(){

  };

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 1;
    double eta_edges[eta_bins + 1] = {0., 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.03);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.03);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == TAUTAU){
        lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))));
      
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      }
    }
    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == TAUTAU){
        lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))));
			  
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      }
    }
    this->Fill(ch, "Down");
  };
};


class LepTauFakeScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeScaleSys(){};
  
  LepTauFakeScaleSys(Spring15Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeScaleSys(){

  };

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 1;
    double eta_edges[eta_bins + 1] = {0., 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.01);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.01);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2<5){
      if (ch == ETAU && (cenTree->gen_match_2 == 1 || cenTree->gen_match_2==3) )
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU && (cenTree->gen_match_2 == 2 || cenTree->gen_match_2==4))
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
    }
    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2<5){
      if (ch == ETAU && (cenTree->gen_match_2 == 1 || cenTree->gen_match_2==3))
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU && (cenTree->gen_match_2 == 2 || cenTree->gen_match_2==4))
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
    }
    this->Fill(ch, "Down");
  };
};

#undef addvar

#endif //!endif LeptonScaleSys_h
