// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef JetEnergyScaleSys_h
#define JetEnergyScaleSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/leptau_jets.h>
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

//void SetAC1B(const AC1B * tree){ //Merijn 2019 5 7: needed to have consistent behaviour w.r.t. jets.h
void SetAC1B(AC1B * tree){ 
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
	jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Up", jecUncertainties);
    this->Fill("Up");
  };
  
  virtual void ScaleDown(){
	jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Down", jecUncertainties);
    this->Fill("Down");
  };
  
  virtual void Fill(const char* shift){
    outTree[shift]->Fill();
	jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, "central");
  }

  //const AC1B * analysisTree;
  AC1B * analysisTree; //Merijn 2019 5 7: needed to have consistent behaviour w.r.t. jets.h
  Config * cfg;
  const btag_scaling_inputs * inputs_btag_scaling;
  std::map< std::string, TTree* >  outTree;
  TString uncertainty_name;
  JESUncertainties * jecUncertainties;
  
};

#undef addvar

#endif //!endif TopPtWeightSys_h
