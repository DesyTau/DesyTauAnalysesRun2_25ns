#ifndef PFMETSys_h
#define PFMETSys_h

#include "DesyTauAnalyses/NTupleMaker/interface/METSys_WIP.h"

class PFMETSys : public METSys {

 public:
  PFMETSys(){};
  PFMETSys(Synch17Tree* c, TString name) {
    cenTree = c;
    label = "CMS_met_"+name+"_13TeV";
    this->Init(cenTree);
  };
  virtual void Eval(utils::channel ch = utils::UNKNOWN) {
    this->Central();
    this->ScaleUp();
    this->ScaleDown();
  };
  
  virtual ~PFMETSys(){};

 protected:

  virtual void Central() {};
  virtual void ScaleDown() {
    bool uncertaintyFound = false;
    if (label.Contains("UnclusteredEn")) {
      obs.metx = analysisTree->pfmetcorr_ex_UnclusteredEnDown;
      obs.mety = analysisTree->pfmetcorr_ey_UnclusteredEnDown;
      uncertaintyFound = true;
    }
    else if (label.Contains("JetRes")) {
      obs.metx = analysisTree->pfmetcorr_ex_JetResDown;
      obs.mety = analysisTree->pfmetcorr_ey_JetResDown;
      uncertaintyFound = true;
   }    
    else if (label.Contains("JetEn")) {
      obs.metx = analysisTree->pfmetcorr_ex_JetEnDown;
      obs.mety = analysisTree->pfmetcorr_ey_JetEnDown;
      uncertaintyFound = true;
   }    
    else {
      std::cout << "Systematic uncertainty " << label << std::endl;
      obs.met = cenTree->met;
      obs.metphi = cenTree->metphi;
    }
    if (uncertaintyFound) {
      obs.metphi = atan2(obs.mety,obs.metx);
      obs.met = sqrt(obs.metx*obs.metx+obs.mety*obs.mety);
      //      std::cout << label << " (met,metphi) : central = (" 
      //		<< cenTree->met << ","
      //		<< cenTree->metphi << ")" 
      //		<< " down = (" << obs.met << ","
      //                << obs.metphi << ")" << std::endl;
    }
    FillPFMET("Down");
  };

  virtual void ScaleUp() {
    bool uncertaintyFound = false;
    if (label.Contains("UnclusteredEn")) {
      obs.metx = analysisTree->pfmetcorr_ex_UnclusteredEnUp;
      obs.mety = analysisTree->pfmetcorr_ey_UnclusteredEnUp;
      uncertaintyFound = true;
    }
    else if (label.Contains("JetRes")) {
      obs.metx = analysisTree->pfmetcorr_ex_JetResUp;
      obs.mety = analysisTree->pfmetcorr_ey_JetResUp;
      uncertaintyFound = true;
   }    
    else if (label.Contains("JetEn")) {
      obs.metx = analysisTree->pfmetcorr_ex_JetEnUp;
      obs.mety = analysisTree->pfmetcorr_ey_JetEnUp;
      uncertaintyFound = true;
   }    
    else {
      std::cout << "Systematic uncertainty " << label << std::endl;
      obs.met = cenTree->met;
      obs.metphi = cenTree->metphi;
    }
    if (uncertaintyFound) {
      obs.metphi = atan2(obs.mety,obs.metx);
      obs.met = sqrt(obs.metx*obs.metx+obs.mety*obs.mety);
      obs.metphi = atan2(obs.mety,obs.metx);
      obs.met = sqrt(obs.metx*obs.metx+obs.mety*obs.mety);
      //      std::cout << label << " (met,metphi) : central = (" 
      //		<< cenTree->met << ","
      //		<< cenTree->metphi << ")" 
      //		<< " up = (" << obs.met << ","
      //                << obs.metphi << ")" << std::endl;
    }
    FillPFMET("Up");
  };

};

#endif
