// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef JetEnergyScaleSys_h
#define JetEnergyScaleSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics.h>

using namespace utils;

class JetEnergyScaleSys : public Systematics {
public:
  
  JetEnergyScaleSys(){};
  
  JetEnergyScaleSys(Synch17Tree* c, TString name){
    cenTree = c;
    label = "CMS_scale_j_13TeV";    
	this->SetUncertaintyName(name);
    this->Init(cenTree);
  };
  
  virtual ~JetEnergyScaleSys(){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      delete it->second;
  };
  
  virtual void Eval(utils::channel ch = utils::UNKNOWN){
    this->Central();
    this->ScaleUp();
    this->ScaleDown();
  };

  virtual void Write(){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      it->second->Write();
  };

void SetAC1B(const AC1B * tree){ 
  analysisTree = tree;
}; 

void SetConfig(Config * config){
  cfg = config;
};

void SetBtagScaling(const btag_scaling_inputs * _InputsBtagScaling){
	inputs_btag_scaling = _InputsBtagScaling;
};

void SetUncertaintyName(TString name){
	uncertainty_name = name;
	label = "CMS_scale_j_"+uncertainty_name+"13TeV";
}

void SetJESUncertainties(JESUncertainties * jec){
	jecUncertainties = jec;
}

protected:

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
  
  virtual void Central(){
  };

  virtual void ScaleUp(){
    met_cen = cenTree->met;
    metphi_cen = cenTree->metphi;
    puppimet_cen = cenTree->puppimet;
    puppimetphi_cen = cenTree->puppimetphi;
    mt_tot_cen = cenTree->mt_tot;
    pt_tt_cen = cenTree->pt_tt;
    mt_1_cen = cenTree->mt_1;
    mt_2_cen = cenTree->mt_2;
    puppimt_1_cen = cenTree->puppimt_1;
    puppimt_2_cen = cenTree->puppimt_2;    
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Up", jecUncertainties);
    this->Fill("Up");
  };
  
  virtual void ScaleDown(){    
    met_cen = cenTree->met;
    metphi_cen = cenTree->metphi;
    puppimet_cen = cenTree->puppimet;
    puppimetphi_cen = cenTree->puppimetphi;
    mt_tot_cen = cenTree->mt_tot;
    pt_tt_cen = cenTree->pt_tt;
    mt_1_cen = cenTree->mt_1;
    mt_2_cen = cenTree->mt_2;
    puppimt_1_cen = cenTree->puppimt_1;
    puppimt_2_cen = cenTree->puppimt_2;
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Down", jecUncertainties);
    this->Fill("Down");
  };
  

  virtual void Fill(const char* shift){
    outTree[shift]->Fill();
    cenTree->met = met_cen;
    cenTree->metphi = metphi_cen;
    cenTree->puppimet = puppimet_cen;
    cenTree->puppimetphi = puppimetphi_cen;
    cenTree->mt_1 = mt_1_cen;
    cenTree->mt_2 = mt_2_cen;
    cenTree->puppimt_1 = puppimt_1_cen;
    cenTree->puppimt_2 = puppimt_2_cen;
    cenTree->mt_tot = mt_tot_cen;
    cenTree->pt_tt = pt_tt_cen;
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, "central");
  }

  //const AC1B * analysisTree;
 const AC1B * analysisTree; //Merijn 2019 5 7: addapted to have consistent behaviour w.r.t. jets.h 2019 8 2 adapted back to const..
  Config * cfg;
  const btag_scaling_inputs * inputs_btag_scaling;
  std::map< std::string, TTree* >  outTree;
  TString uncertainty_name;
  JESUncertainties * jecUncertainties;

  float met_cen; 
  float metphi_cen; 
  float puppimet_cen; 
  float puppimetphi_cen; 
  float mt_tot_cen; 
  float pt_tt_cen; 
  float mt_1_cen; 
  float mt_2_cen; 
  float puppimt_1_cen; 
  float puppimt_2_cen; 
  
};

#undef addvar

#endif //!endif TopPtWeightSys_h
