// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef BtagSys_h
#define BtagSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h>

using namespace utils;

class BtagSys : public Systematics {
public:
  
  BtagSys(){};
  
  BtagSys(Synch17Tree* c, TString name){
    cenTree = c;
    label = "CMS_eff_b_13TeV";    
    this->SetUncertaintyName(name);
    this->Init(cenTree);
  };
  
  virtual ~BtagSys(){
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
};

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
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Up");
    this->Fill("Up");
  };
  
  virtual void ScaleDown(){    
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Down");
    this->Fill("Down");
  };
  
  virtual void Fill(const char* shift){
    outTree[shift]->Fill();
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, "central");
  }

  //const AC1B * analysisTree;
 const AC1B * analysisTree; //Merijn 2019 5 7: addapted to have consistent behaviour w.r.t. jets.h 2019 8 2 adapted back to const..
  Config * cfg;
  const btag_scaling_inputs * inputs_btag_scaling;
  std::map< std::string, TTree* >  outTree;
  TString uncertainty_name;
  
};

#undef addvar

#endif //!endif TopPtWeightSys_h
