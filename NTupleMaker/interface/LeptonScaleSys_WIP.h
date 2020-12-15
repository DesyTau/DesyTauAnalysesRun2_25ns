// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef LeptonScaleSys_h
#define LeptonScaleSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)

#include "TH2D.h"

#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics.h>
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

using namespace utils;

class LeptonScaleSys : public Systematics {
public:
  
  LeptonScaleSys(){};
  
  LeptonScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "lScale";
    usePuppiMET = true;
    useSVFit = false;
    useFastMTT = false;
    errorBarrelUp = 0.0;
    errorBarrelDown = 0.0;
    errorEndcapUp = 0.0;
    errorEndcapDown = 0.0;
    barrelEdge = 1.5;
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

  virtual void Write(const char *name="", Int_t option=0){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      it->second->Write(name,option);
  };

  void SetLabel(TString l){
    label = l;
  };

  void SetSvFitVisPtResolution(TFile* f){
    svFit_visPtResolution = f;
  };

  void SetUseSVFit(bool isSV){
    useSVFit = isSV;	
  }

  void SetUseFastMTT(bool isFastMTT){
    useFastMTT = isFastMTT;	
  }

  void SetUsePuppiMET(bool isPuppiMET) {
    usePuppiMET = isPuppiMET;
  }
  
  void SetScale(float Central, float Error) {
    central = Central;
    error = Error;
    this->Init(cenTree);
  }  
  
  void SetBarrelEdge(float BarrelEdge) {
    barrelEdge = BarrelEdge;
  }  
  
  void SetScaleBarrelUp(float Central, float ErrorBarrelUp) {
    errorBarrelUp = ErrorBarrelUp;
    this->Init(cenTree);
  }  
  
  void SetScaleBarrelDown(float Central, float ErrorBarrelDown) {
    errorBarrelDown = ErrorBarrelDown;
    this->Init(cenTree);
  }  
  
  void SetScaleEndcapUp(float Central, float ErrorEndcapUp) {
    errorEndcapUp = ErrorEndcapUp;
    this->Init(cenTree);
  }  
  
  void SetScaleEndcapDown(float Central, float ErrorEndcapDown) {
    errorEndcapDown = ErrorEndcapDown;
    this->Init(cenTree);
  }  

protected:
  virtual void Init(Synch17Tree* c){
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
    float metphi_cen = cenTree->metphi;
    float puppimet_cen = cenTree->puppimet;
    float puppimetphi_cen = cenTree->puppimetphi;

    float mt_1_cen = cenTree->mt_1;
    float mt_2_cen = cenTree->mt_2;
    float puppimt_1_cen = cenTree->puppimt_1;
    float puppimt_2_cen = cenTree->puppimt_2;

    float pt_tt_cen = cenTree->pt_tt;
    float pzetavis_cen = cenTree->pzetavis;
    float pzetamiss_cen = cenTree->pzetamiss;
    float pzeta_cen = cenTree->pzeta;
    float mt_tot_cen = cenTree->mt_tot;

    float m_vis_cen = cenTree->m_vis;

    float m_sv_cen = cenTree->m_sv;
    float pt_sv_cen = cenTree->pt_sv;
    float eta_sv_cen = cenTree->eta_sv;
    float phi_sv_cen = cenTree->phi_sv;
    float met_sv_cen = cenTree->met_sv;
    float mt_sv_cen = cenTree->mt_sv;
    
    float mt_fast_cen = cenTree->mt_fast; 
    float m_fast_cen = cenTree->m_fast;
    float pt_fast_cen = cenTree->pt_fast;
    float eta_fast_cen = cenTree->eta_fast;
    float phi_fast_cen = cenTree->phi_fast; 

    // calc shifted values
    cenTree->pt_1 = lep1_scaled.Pt();
    cenTree->pt_2 = lep2_scaled.Pt();

    // central value of the met    
    TLorentzVector pfmetLV; pfmetLV.SetXYZT(cenTree->met * cos(cenTree->metphi),
					    cenTree->met * sin(cenTree->metphi),
					    0.,
					    cenTree->met);


    // central value of the puppi met
    TLorentzVector puppimetLV; puppimetLV.SetXYZT(cenTree->puppimet * cos(cenTree->puppimetphi),
						  cenTree->puppimet * sin(cenTree->puppimetphi),
						  0.,
						  cenTree->puppimet);

    // propagate the tau pt shift to the MET 
    pfmetLV.SetPx(pfmetLV.Px()- (lep2_scaled.Px()-lep2.Px()));
    pfmetLV.SetPy(pfmetLV.Py()- (lep2_scaled.Py()-lep2.Py()));

    pfmetLV.SetPx(pfmetLV.Px()- (lep1_scaled.Px()-lep1.Px()));
    pfmetLV.SetPy(pfmetLV.Py()- (lep1_scaled.Py()-lep1.Py()));
    
    // propagate the lepton pt shift to the puppiMET
    puppimetLV.SetPx(puppimetLV.Px()- (lep2_scaled.Px()-lep2.Px()));
    puppimetLV.SetPy(puppimetLV.Py()- (lep2_scaled.Py()-lep2.Py()));

    puppimetLV.SetPx(puppimetLV.Px()- (lep1_scaled.Px()-lep1.Px()));
    puppimetLV.SetPy(puppimetLV.Py()- (lep1_scaled.Py()-lep1.Py()));

    // calc shifted values
    cenTree->met = pfmetLV.Pt();
    cenTree->metphi = pfmetLV.Phi();

    cenTree->mt_1 = mT(lep1_scaled,pfmetLV);
    cenTree->mt_2 = mT(lep2_scaled,pfmetLV);

    // changing puppi MET
    cenTree->puppimet = puppimetLV.Pt();
    cenTree->puppimetphi = puppimetLV.Phi();

    cenTree->puppimt_1 = mT(lep1_scaled,puppimetLV);
    cenTree->puppimt_2 = mT(lep2_scaled,puppimetLV);
    
    TMatrixD covMET(2, 2);
    covMET[0][0] = cenTree->metcov00;
    covMET[1][0] = cenTree->metcov10;
    covMET[0][1] = cenTree->metcov01;
    covMET[1][1] = cenTree->metcov11;

    double measuredMETx =  cenTree->met * cos(cenTree->metphi);
    double measuredMETy =  cenTree->met * sin(cenTree->metphi);

    if (usePuppiMET) {
      measuredMETx = cenTree->puppimet * cos(cenTree->puppimetphi);
      measuredMETy = cenTree->puppimet * sin(cenTree->puppimetphi);
    }

    if (usePuppiMET) {
      covMET[0][0] = cenTree->puppimetcov00;
      covMET[1][0] = cenTree->puppimetcov10;
      covMET[0][1] = cenTree->puppimetcov01;
      covMET[1][1] = cenTree->puppimetcov11;
    }

    TLorentzVector dileptonLV = lep1_scaled + lep2_scaled;

    TLorentzVector metxLV = pfmetLV;
    if (usePuppiMET) 
      metxLV = puppimetLV;
    
    cenTree->pt_tt = (lep1_scaled+lep2_scaled+metxLV).Pt();
    
    cenTree->pzetavis = calc::pzetavis(lep1_scaled, lep2_scaled);
    cenTree->pzetamiss = calc::pzetamiss( lep1_scaled, lep2_scaled, metxLV);
    cenTree->pzeta = calc::pzeta(lep1_scaled, lep2_scaled, metxLV);
    
    cenTree->m_vis = dileptonLV.M();
    cenTree->mt_tot = calc::mTtot(lep1_scaled,lep2_scaled,metxLV);    

    // add flag for svfit
    if(ch != UNKNOWN && (useSVFit || useFastMTT)){

      std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
      classic_svFit::MeasuredTauLepton::kDecayType type_ = classic_svFit::MeasuredTauLepton::kUndefinedDecayType;

      if (ch == MUTAU)  type_ = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
      if (ch == ETAU)   type_ = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
      if (ch == TAUTAU) type_ = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
      if (ch == EMU)    type_ = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
      // define lepton four vectors                                                                                                             

      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type_,
								    lep1_scaled.Pt(),
								    lep1_scaled.Eta(),
								    lep1_scaled.Phi(),
								    lep1_scaled.M()));
      if (ch == EMU) {
	measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay,
								      lep2_scaled.Pt(),
								      lep2_scaled.Eta(),
								      lep2_scaled.Phi(),
								      lep2_scaled.M())); 
	
      }
      else {
	measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,
								      lep2_scaled.Pt(),
								      lep2_scaled.Eta(),
								      lep2_scaled.Phi(),
								      lep2_scaled.M(), 
  								      cenTree->tau_decay_mode_2));
      }
      if (useSVFit) {	  
				       
	int verbosity = 1;
	ClassicSVfit svFitAlgo(verbosity);
	double kappa = 4.; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
	
	svFitAlgo.addLogM_fixed(true, kappa);
	svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	bool isValidSolution = svFitAlgo.isValidSolution();
	
	cenTree->m_sv   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
	cenTree->pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPt();
	cenTree->eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getEta();
	cenTree->phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPhi();      
	//otree->met_sv = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getfittedMET().Rho();
	cenTree->mt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
      
      }
      if (useFastMTT) {

	// FasMTT
	LorentzVector tau1P4;
	LorentzVector tau2P4;
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	tau1P4 = aFastMTTAlgo.getTau1P4();
	tau2P4 = aFastMTTAlgo.getTau2P4();
    
	double dPhiTT = dPhiFrom2P( tau1P4.Px(), tau1P4.Py(), tau2P4.Px(), tau2P4.Py() );
	cenTree->mt_fast = TMath::Sqrt(2*tau1P4.Pt()*tau2P4.Pt()*(1 - TMath::Cos(dPhiTT)));
	cenTree->m_fast = ttP4.M();
	cenTree->pt_fast = ttP4.Pt();
	cenTree->eta_fast = ttP4.Eta();
	cenTree->phi_fast = ttP4.Phi();
      }
    }

    outTree[shift]->Fill();
/*
  std::cout << " SVFit in shift  - Inputs -  " << std::endl;
  std::cout << " Lep 1 : pt = " << lep1_scaled.Pt() << " eta : " << lep1_scaled.Eta() << " phi : " << lep1_scaled.Phi() << " M : " << lep1_scaled.M() <<   " DM : " << cenTree->tau_decay_mode_1 << std::endl;
  std::cout << " Lep 2 : pt = " << lep2_scaled.Pt() << " eta : " << lep2_scaled.Eta() << " phi : " << lep2_scaled.Phi() << " M : " << lep2_scaled.M() <<   " DM : " << cenTree->tau_decay_mode_2 << std::endl;
  std::cout << " MET : px = " <<  pfmetLV.Px() << " py = " <<  pfmetLV.Py() << std::endl;
  std::cout << " MET COV : 00 = " << cenTree->metcov00 << " 01 " << cenTree->metcov01  << " 10 " << cenTree->metcov10 << " 11 " << cenTree->metcov11 << std::endl;
  std::cout << " SVFit in shift  - Output -  " << std::endl;
  std::cout << " Mass sv = " << cenTree->m_sv << "pt sv = " << cenTree->pt_sv << std::endl;  
  std::cout << " SVFit in shift  - END -  " << std::endl;
*/
    // restore cen values
    cenTree->pt_1 = pt_1_cen;
    cenTree->pt_2 = pt_2_cen;
    cenTree->met = met_cen;
    cenTree->metphi = metphi_cen;
    cenTree->puppimet = puppimet_cen;
    cenTree->puppimetphi = puppimetphi_cen;
    
    cenTree->mt_1 = mt_1_cen;
    cenTree->mt_2 = mt_2_cen;

    cenTree->puppimt_1 = puppimt_1_cen;
    cenTree->puppimt_2 = puppimt_2_cen;

    cenTree->pt_tt = pt_tt_cen;
    cenTree->pzetavis = pzetavis_cen;
    cenTree->pzetamiss = pzetamiss_cen;
    
    cenTree->m_vis = m_vis_cen;

    cenTree->m_sv = m_sv_cen;
    cenTree->pt_sv = pt_sv_cen;
    cenTree->eta_sv = eta_sv_cen;
    cenTree->phi_sv = phi_sv_cen;
    cenTree->met_sv = met_sv_cen;
    cenTree->mt_sv = mt_sv_cen;    
	
    cenTree->mt_fast = mt_fast_cen; 
    cenTree->m_fast = m_fast_cen;
    cenTree->pt_fast = pt_fast_cen;
    cenTree->eta_fast = eta_fast_cen;
    cenTree->phi_fast = phi_fast_cen;

    cenTree->mt_tot = mt_tot_cen;

  }

  //std::map< std::string, Float_t >  mt_sv;
  
  TH2D* sf_up;
  TH2D* sf_down;  

  TLorentzVector lep1, lep2;
  TLorentzVector lep1_scaled, lep2_scaled;

  TFile* svFit_visPtResolution;
  bool useSVFit;
  bool usePuppiMET;
  bool useFastMTT;
  float central;
  float error;
  float errorBarrelUp;
  float errorBarrelDown;
  float errorEndcapUp;
  float errorEndcapDown;
  float barrelEdge;
  std::map< std::string, TTree* >  outTree;
  //std::map< std::string, SpringTree* >  outTree;
};

class MuonScaleSys : public LeptonScaleSys { 
public:
  MuonScaleSys(){};
  
  MuonScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_mu_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~MuonScaleSys(){};
  
protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};

    const int eta_bins = 3;
    double eta_edges[eta_bins + 1] = {0., 1.2, 2.1,  2.4};

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent(sf_up->FindBin(100.,0.9),0.004);
    sf_down->SetBinContent(sf_down->FindBin(100.,0.9),0.004);

    sf_up->SetBinContent(sf_up->FindBin(100.,1.5),0.009);
    sf_down->SetBinContent(sf_down->FindBin(100.,1.5),0.009);

    sf_up->SetBinContent(sf_up->FindBin(100.,2.2),0.017);
    sf_down->SetBinContent(sf_down->FindBin(100.,2.2),0.017);

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

  ElectronScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_e_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~ElectronScaleSys(){};

  void SetAC1B(const AC1B * tree){
    analysisTree = tree;
  };
  
  void SetElectronIndex(int index){
    electronIndex = index;
  };
  
  void SetIsEmbedded(bool isEmbeddedFlag){
    isEmbedded = isEmbeddedFlag;
  };

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
  
    // this is used only for Embedded samples: ±0.50% (barrel) and  ±1.25%(endcap) for all the years (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016Legacy#Lepton_energy_scale_corrections)
    // otherwise read the up/down variations from BigNTuple tree
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), 0.005);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), 0.0125);   
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), 0.005);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), 0.0125);
  };
  
  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;
    
    if (isEmbedded) {
      
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

    } else {
      if (ch == EMU)
      lep1_scaled.SetXYZM(analysisTree->electron_px_energyscale_up[electronIndex],
        analysisTree->electron_py_energyscale_up[electronIndex],
        analysisTree->electron_pz_energyscale_up[electronIndex],
        lep1.M());
      else if (ch == ETAU)
        lep1_scaled.SetXYZM(analysisTree->electron_px_energyscale_up[electronIndex],
          analysisTree->electron_py_energyscale_up[electronIndex],
          analysisTree->electron_pz_energyscale_up[electronIndex],
          lep1.M());
    }
    // 
    // float pt1 = lep1.Pt();
    // float absEta1 = fabs(lep1.Eta());
    // float pt2 = lep2.Pt();
    // float absEta2 = fabs(lep2.Eta());

    // if (ch == EMU)
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    // else if (ch == ETAU)
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    // else if (ch == EE){
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    //   lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
		// 	  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
		// 	  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
		// 	  lep2.M());
    // }

    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    if (isEmbedded) {
      // for Embedded read variations from SF hist
      float pt1 = lep1.Pt();
      float absEta1 = fabs(lep1.Eta());
      float pt2 = lep2.Pt();
      float absEta2 = fabs(lep2.Eta());
      
      if (ch == EMU)
        lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.M());
      else if (ch == ETAU)
        lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
  			  lep1.M());
    } else {
      // otherwie reading scaled momentum values from BigNTuple tree 
      if (ch == EMU)
      lep1_scaled.SetXYZM(analysisTree->electron_px_energyscale_down[electronIndex],
        analysisTree->electron_py_energyscale_down[electronIndex],
        analysisTree->electron_pz_energyscale_down[electronIndex],
        lep1.M());
      else if (ch == ETAU)
        lep1_scaled.SetXYZM(analysisTree->electron_px_energyscale_down[electronIndex],
          analysisTree->electron_py_energyscale_down[electronIndex],
          analysisTree->electron_pz_energyscale_down[electronIndex],
          lep1.M());
    }

    // float pt1 = lep1.Pt();
    // float absEta1 = fabs(lep1.Eta());
    // float pt2 = lep2.Pt();
    // float absEta2 = fabs(lep2.Eta());
    // 
    // if (ch == EMU)
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    // else if (ch == MUTAU)
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    // else if (ch == MUMU){
    //   lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
		// 	  lep1.M());
    //   lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
		// 	  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
		// 	  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
		// 	  lep2.M());
    // }
    
    this->Fill(ch, "Down");
  };
  
  const AC1B * analysisTree;
  int electronIndex;
  bool isEmbedded;
};

class ElectronEBScaleSys : public ElectronScaleSys {
public:
  ElectronEBScaleSys(){};
  
  ElectronEBScaleSys(Synch17Tree* c){
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
  
  ElectronEEScaleSys(Synch17Tree* c){
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

// don't split shapes by decay mode.
class TauScaleSys : public LeptonScaleSys { 
public:
  TauScaleSys(){};
  
  TauScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_shape_t_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~TauScaleSys(){};

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

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), error);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), error);
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

// split by decay mode. One prong shape, for tau-decay_mode_2 ==0
class TauOneProngScaleSys : public LeptonScaleSys { 
public:
  TauOneProngScaleSys(){};
  
  TauOneProngScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_t_1prong_13TeV";
    
    this->Init(cenTree);
  };

  virtual ~TauOneProngScaleSys(){};

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

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), error);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), error);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2 ==0){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  );
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor up : " << 1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2)) << endl;
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.M()  );
      }
      else if (ch == TAUTAU){
        lep1_scaled.SetXYZM(lep1.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt1, absEta1))),
			  lep1.M()  );
      
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  );
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
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2==0){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M() );
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor dwon : " << 1. + sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2)) << endl;
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.M() );
      }
      else if (ch == TAUTAU){
        lep1_scaled.SetXYZM(lep1.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt1, absEta1))),
			  lep1.M() );
			  
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M() );
      }
    }
    this->Fill(ch, "Down");
  };
};

// split by decay mode. OneProngOnePi0 shape, for tau_decay_mode_2 ==1
class TauOneProngOnePi0ScaleSys : public LeptonScaleSys { 
public:
  TauOneProngOnePi0ScaleSys(){};
  
  TauOneProngOnePi0ScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_t_1prong1pizero_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~TauOneProngOnePi0ScaleSys(){};

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

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), error);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), error);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2 ==1){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor up : " << 1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2)) << endl;
	lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      }
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
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2==1){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor down : " << 1. + sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2)) << endl;
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      }
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


// split by decay mode. Three prong shape, for tau_decay_mode_2 ==10
class TauThreeProngScaleSys : public LeptonScaleSys { 
public:
  TauThreeProngScaleSys(){};
  
  TauThreeProngScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_t_3prong_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~TauThreeProngScaleSys(){};

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

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), error);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), error);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2 ==10){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor up : " << 1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2)) << endl;
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      }
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
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2==10){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//        cout << "scale factor down : " << 1. + sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2)) << endl;
	lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      }
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

// split by decay mode. Three prong shape, for tau_decay_mode_2 ==11
class TauThreeProngOnePi0ScaleSys : public LeptonScaleSys { 
public:
  TauThreeProngOnePi0ScaleSys(){};
  
  TauThreeProngOnePi0ScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_scale_t_3prong1pizero_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~TauThreeProngOnePi0ScaleSys(){};

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

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), error);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), error);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2 ==11){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//	cout << "scale factor up : " << 1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2)) << endl;
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      }
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
    if (cenTree->gen_match_2==5 && cenTree->tau_decay_mode_2==11){
      if (ch == ETAU)
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU) {
	//	cout << "decay mode : " << cenTree->tau_decay_mode_2 << endl;
	//        cout << "scale factor down : " << 1. + sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2)) << endl;
	lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			    lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      }
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

// shift tau pt and propagate to MEt and svfit - for all fake taus. 
class LepTauFakeScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeScaleSys(){};
  
  LepTauFakeScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeScaleSys(){};

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

// shift tau pt and propagate to MEt and svfit - for 1prong fake taus
//!!! it is assumed that the choice of mt or et channels is done and the corresponding values are passed *outside* of the class
//!!! the same applies also for barrel/endcap choice
class LepTauFakeOneProngScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeOneProngScaleSys(){};
  
  LepTauFakeOneProngScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_1prong_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeOneProngScaleSys(){};

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 2;
    double eta_edges[eta_bins + 1] = {0., barrelEdge, 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), errorBarrelUp);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), errorBarrelDown);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), errorEndcapUp);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), errorEndcapDown);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->tau_decay_mode_2==0)
      lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                    		  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                    		  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                    		  lep2.M()  );
    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->tau_decay_mode_2==0)
      lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.M() );
    this->Fill(ch, "Down");
  };
};


// shift tau pt and propagate to MEt and svfit - for 1prong1pi0 fake taus
//!!! it is assumed that the choice of mt or et channels is done and the corresponding values are passed *outside* of the class
//!!! the same applies also for barrel/endcap choice
class LepTauFakeOneProngOnePi0ScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeOneProngOnePi0ScaleSys(){};
  
  LepTauFakeOneProngOnePi0ScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_1prong1pi_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeOneProngOnePi0ScaleSys(){};

protected:
  virtual void InitSF(){
    const int pt_bins = 1;
    double pt_edges[pt_bins + 1] = {0., 14001.};
    double pt_central[pt_bins] = {};
    for(int ibin = 0; ibin < pt_bins; ibin++)
      pt_central[ibin] = (pt_edges[ibin+1] + pt_edges[ibin])/2.; 

    const int eta_bins = 2;
    double eta_edges[eta_bins + 1] = {0., barrelEdge, 2.5};
    double eta_central[eta_bins] = {};
    for(int ibin = 0; ibin < eta_bins; ibin++)
      eta_central[ibin] = (eta_edges[ibin+1] + eta_edges[ibin])/2.; 

    sf_up = new TH2D(label+"_sf_up", label+"_sf_up", pt_bins, pt_edges, eta_bins, eta_edges);
    sf_down = new TH2D(label+"_sf_down", label+"_sf_down", pt_bins, pt_edges, eta_bins, eta_edges);

    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[0] ), errorBarrelUp);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[0] ), errorBarrelDown);
    sf_up->SetBinContent( sf_up->FindBin( pt_central[0], eta_central[1] ), errorEndcapUp);
    sf_down->SetBinContent( sf_down->FindBin( pt_central[0], eta_central[1] ), errorEndcapDown);
  };

  virtual void ScaleUp(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->tau_decay_mode_2==1)
      lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                  			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                  			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
                  			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
    this->Fill(ch, "Up");
  };
  
  virtual void ScaleDown(utils::channel ch){
    lep1_scaled = lep1;
    lep2_scaled = lep2;

    float pt1 = lep1.Pt();
    float absEta1 = fabs(lep1.Eta());
    float pt2 = lep2.Pt();
    float absEta2 = fabs(lep2.Eta());
    if (cenTree->tau_decay_mode_2==1)
      lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
                  			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
    this->Fill(ch, "Down");
  };
};


// shift tau pt and propagate to MEt and svfit - for 3prong fake taus. Shift all lep->tau fakes. 
// for mu-tau, 1.5% shift. For e-tau, 3% shift.
class LepTauFakeThreeProngScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeThreeProngScaleSys(){};
  
  LepTauFakeThreeProngScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_3prong_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeThreeProngScaleSys(){};

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
    if (cenTree->gen_match_2<5 && cenTree->tau_decay_mode_2==10){
      if (ch == ETAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
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
    if (cenTree->gen_match_2<5 && cenTree->tau_decay_mode_2==10){
      if (ch == ETAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
    }
    this->Fill(ch, "Down");
  };
};


// shift tau pt and propagate to MEt and svfit - for 3prong fake taus. Shift all lep->tau fakes. 
// for mu-tau, 1.5% shift. For e-tau, 3% shift.
class LepTauFakeThreeProngOnePi0ScaleSys : public LeptonScaleSys { 
public:
  LepTauFakeThreeProngOnePi0ScaleSys(){};
  
  LepTauFakeThreeProngOnePi0ScaleSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_htt_ZLShape_3prong1pizero_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~LepTauFakeThreeProngOnePi0ScaleSys(){};

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
    if (cenTree->gen_match_2<5 && cenTree->tau_decay_mode_2==11){
      if (ch == ETAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
      else if (ch == MUTAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. + 0.5*sf_up->GetBinContent( sf_up->FindBin(pt2, absEta2))));
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
    if (cenTree->gen_match_2<5 && cenTree->tau_decay_mode_2==11){
      if (ch == ETAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
      else if (ch == MUTAU )
        lep2_scaled.SetXYZM(lep2.Px() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Py() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.Pz() * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))),
			  lep2.M()  * (1. - 0.5*sf_down->GetBinContent( sf_down->FindBin(pt2, absEta2))));
    }
    this->Fill(ch, "Down");
  };
};





#undef addvar

#endif //!endif LeptonScaleSys_h
