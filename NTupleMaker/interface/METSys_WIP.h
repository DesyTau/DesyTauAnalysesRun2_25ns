// MET Systematics evaluator
// Author: Alexei Raspereza  (rasp@mail.desy.de)

#ifndef METSys_h
#define METSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)

#include <TLorentzVector.h>
#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h>
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include <TString.h>


class METSys : public Systematics {
 public:
  
  METSys() {};
  METSys(Synch17Tree* c, TString name){
    cenTree = c;
    label = "CMS_met" + name + "_13TeV";
    this->Init(cenTree);
  };

  virtual ~METSys() {
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      delete it->second;
  };


  virtual void Write(const char *name="", Int_t option=0){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      it->second->Write(name,option);
  };
  
  void SetSvFitVisPtResolution(TFile* f){
    svFit_visPtResolution = f;
  };

  void SetUseSVFit(bool isSV){
    useSVFit = isSV;	
  };

  void SetAC1B(const AC1B * tree){
    analysisTree = tree;
  };

 protected:

  struct Observables {
    float met;
    float metphi;
    float mt_1;
    float mt_2;
    float metx;
    float mety;
    float pt_tt;
    float mt_tot;
    float m_sv;
    float pt_sv;
    float eta_sv;
    float phi_sv;
    float met_sv;
    float mt_sv;
    float pzetamiss;
    float pzeta;
    
  } obs;


  virtual void Init(Synch17Tree* c){
    cenTree = c;
    this->InitTree("Up");
    this->InitTree("Down");
  };

  virtual void InitTree(const char* shift){
    std::cout<<label+shift<<std::endl;
    outTree[shift] = cenTree->fChain->CloneTree(0);
    outTree[shift]->SetName(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetTitle(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetDirectory(cenTree->fChain->GetDirectory());
  };

  TFile* svFit_visPtResolution;
  bool useSVFit;
  const AC1B * analysisTree;

  virtual void PropagateUncertainty() {
    obs.metx = obs.met*cos(obs.metphi);
    obs.mety = obs.met*sin(obs.metphi);
    TLorentzVector metLV; metLV.SetXYZT(obs.metx,obs.mety,0,obs.met);

    TLorentzVector lep1LV; lep1LV.SetXYZM(cenTree->pt_1 * cos(cenTree->phi_1),
					  cenTree->pt_1 * sin(cenTree->phi_1),
					  cenTree->pt_1 * sinh(cenTree->eta_1),
					  cenTree->m_1);

    TLorentzVector lep2LV; lep2LV.SetXYZM(cenTree->pt_2 * cos(cenTree->phi_2),
					  cenTree->pt_2 * sin(cenTree->phi_2),
					  cenTree->pt_2 * sinh(cenTree->eta_2),
					  cenTree->m_2);

    obs.mt_1 = mT(lep1LV,metLV);
    obs.mt_2 = mT(lep2LV,metLV);
    obs.pt_tt = (lep1LV+lep2LV+metLV).Pt();
    obs.mt_tot  = calc::mTtot(lep1LV,lep2LV,metLV);

    obs.pzetamiss = calc::pzetamiss(lep1LV,lep2LV,metLV);
    obs.pzeta = calc::pzeta(lep1LV,lep2LV,metLV);
    
    // initialize
    obs.m_sv = cenTree->m_sv;
    obs.pt_sv = cenTree->pt_sv;
    obs.eta_sv = cenTree->eta_sv;
    obs.phi_sv = cenTree->phi_sv;
    obs.met_sv = cenTree->met_sv;
    obs.mt_sv = cenTree->mt_sv;

  };

  virtual void FillPuppiMET(const char* shift) {

    if (cenTree->puppimet==obs.met&&cenTree->puppimetphi==obs.metphi) {
      outTree[shift]->Fill();
      return;
    }

    float puppimt_1_cen = cenTree->puppimt_1;
    float puppimt_2_cen = cenTree->puppimt_2;
    float puppimet_cen = cenTree->puppimet;
    float puppimetphi_cen = cenTree->puppimetphi;
    float pt_tt_cen = cenTree->pt_tt;
    float mt_tot_cen = cenTree->mt_tot;
    float pzeta_cen = cenTree->pzeta;
    float pzetamiss_cen = cenTree->pzetamiss;
    float m_sv_cen = cenTree->m_sv;
    float pt_sv_cen = cenTree->pt_sv;
    float eta_sv_cen = cenTree->eta_sv;
    float phi_sv_cen = cenTree->phi_sv;
    float met_sv_cen = cenTree->met_sv;
    float mt_sv_cen  = cenTree->mt_sv;

    PropagateUncertainty();

    cenTree->puppimt_1 = obs.mt_1;
    cenTree->puppimt_2 = obs.mt_2;
    cenTree->puppimet  = obs.met;
    cenTree->puppimetphi = obs.metphi;
    cenTree->pt_tt = obs.pt_tt;
    cenTree->mt_tot = obs.mt_tot;
    cenTree->pzeta = obs.pzeta;
    cenTree->pzetamiss = obs.pzetamiss;

    if (useSVFit) {
      cenTree->m_sv = obs.m_sv;
      cenTree->pt_sv = obs.pt_sv;
      cenTree->eta_sv = obs.eta_sv;
      cenTree->phi_sv = obs.phi_sv;
      cenTree->met_sv = obs.met_sv;
      cenTree->mt_sv  = obs.mt_sv;
    }

    // Filling
    outTree[shift]->Fill();
    // restore central values

    cenTree->puppimt_1 = puppimt_1_cen;
    cenTree->puppimt_2 = puppimt_2_cen;
    cenTree->puppimet = puppimet_cen;
    cenTree->puppimetphi = puppimetphi_cen;
    cenTree->pt_tt = pt_tt_cen;
    cenTree->mt_tot = mt_tot_cen;
    cenTree->pzeta = pzeta_cen;
    cenTree->pzetamiss = pzetamiss_cen;
    cenTree->m_sv = m_sv_cen;
    cenTree->pt_sv = pt_sv_cen;
    cenTree->eta_sv = eta_sv_cen;
    cenTree->phi_sv = phi_sv_cen;
    cenTree->met_sv = met_sv_cen;
    cenTree->mt_sv = mt_sv_cen;

  };

  virtual void FillPFMET(const char* shift) {

    if (cenTree->met==obs.met&&cenTree->metphi==obs.metphi) {
      outTree[shift]->Fill();
      return;
    }

    // save central values
    float mt_1_cen = cenTree->mt_1;
    float mt_2_cen = cenTree->mt_2;
    float met_cen = cenTree->met;
    float metphi_cen = cenTree->metphi;
    float pt_tt_cen = cenTree->pt_tt;
    float mt_tot_cen = cenTree->mt_tot;
    float m_sv_cen = cenTree->m_sv;
    float pt_sv_cen = cenTree->pt_sv;
    float eta_sv_cen = cenTree->eta_sv;
    float phi_sv_cen = cenTree->phi_sv;
    float met_sv_cen = cenTree->met_sv;
    float mt_sv_cen  = cenTree->mt_sv;
    float pzeta_cen = cenTree->pzeta;
    float pzetamiss_cen = cenTree->pzetamiss;
    
    PropagateUncertainty();

    cenTree->mt_1 = obs.mt_1;
    cenTree->mt_2 = obs.mt_2;
    cenTree->met  = obs.met;
    cenTree->metphi = obs.metphi;
    cenTree->pt_tt = obs.pt_tt;
    cenTree->mt_tot = obs.mt_tot;
    cenTree->pzeta = obs.pzeta;
    cenTree->pzetamiss = obs.pzetamiss;

    if (useSVFit) {
      cenTree->m_sv = obs.m_sv;
      cenTree->pt_sv = obs.pt_sv;
      cenTree->eta_sv = obs.eta_sv;
      cenTree->phi_sv = obs.phi_sv;
      cenTree->met_sv = obs.met_sv;
      cenTree->mt_sv  = obs.mt_sv;
    }

    // Filling
    outTree[shift]->Fill();
    // restore central values
    cenTree->mt_1 = mt_1_cen;
    cenTree->mt_2 = mt_2_cen;
    cenTree->met = met_cen;
    cenTree->metphi = metphi_cen;
    cenTree->pt_tt = pt_tt_cen;
    cenTree->mt_tot = mt_tot_cen;
    cenTree->m_sv = m_sv_cen;
    cenTree->pt_sv = pt_sv_cen;
    cenTree->eta_sv = eta_sv_cen;
    cenTree->phi_sv = phi_sv_cen;
    cenTree->met_sv = met_sv_cen;
    cenTree->mt_sv  = mt_sv_cen;
    cenTree->pzeta = pzeta_cen;
    cenTree->pzetamiss = pzetamiss_cen;

  };

  std::map< std::string, TTree* >  outTree;

};


#undef addvar

#endif
