// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef ZPtWeightSys_h
#define ZPtWeightSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)

#include "TH2D.h"

#include <DesyTauAnalyses/NTupleMaker/interface/functions.h>
#include <DesyTauAnalyses/NTupleMaker/interface/Systematics.h>

using namespace utils;

class ZPtWeightSys : public Systematics {
public:
  
  ZPtWeightSys(){};
  
  ZPtWeightSys(Synch17Tree* c){
    cenTree = c;
    label = "CMS_shape_dyShape_13TeV";
    
    this->Init(cenTree);
  };
  
  virtual ~ZPtWeightSys(){
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
    w_cen = cenTree->zptweight;
  };

  virtual void ScaleUp(){
    w_scaled = cenTree->zptweight * cenTree->zptweight;
    this->Fill("Up");
  };
  
  virtual void ScaleDown(){
    w_scaled = 1.;
    this->Fill("Down");
  };
  
  virtual void Fill(const char* shift){
    cenTree->zptweight = w_scaled;
    outTree[shift]->Fill();

    cenTree->zptweight = w_cen;
  }

  float w_scaled, w_cen;

  std::map< std::string, TTree* >  outTree;
};

#undef addvar

#endif //!endif ZPtWeightSys_h
