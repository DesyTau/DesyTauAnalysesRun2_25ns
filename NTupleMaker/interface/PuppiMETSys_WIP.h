#ifndef PuppiMETSys_h
#define PuppiMETSys_h

#include "DesyTauAnalyses/NTupleMaker/interface/METSys_WIP.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"

using namespace kit;

class PuppiMETSys : public METSys {

 public:
  PuppiMETSys(){};
  PuppiMETSys(Synch17Tree* c, TString name) {
    cenTree = c;
    label = "CMS_met_"+name+"_13TeV";
    this->Init(cenTree);
  };
  void SetMEtSys(kit::MEtSys * sys) {
    metSys = sys;
  };
  virtual void Eval(utils::channel ch = utils::UNKNOWN) {
    this->Central();
    this->ScaleUp();
    this->ScaleDown();
  };
  
  virtual ~PuppiMETSys(){};

 protected:

  kit::MEtSys * metSys;

  virtual void Central() {};
  virtual void ScaleDown() {
    bool uncertaintyFound = false;
    if (label.Contains("UnclusteredEn")) {
      obs.metx = analysisTree->puppimet_ex_UnclusteredEnDown;
      obs.mety = analysisTree->puppimet_ey_UnclusteredEnDown;
      uncertaintyFound = true;
    }
    else if (label.Contains("JetRes")) {
      obs.metx = analysisTree->puppimet_ex_JetResDown;
      obs.mety = analysisTree->puppimet_ey_JetResDown;
      uncertaintyFound = true;
    }    
    else if (label.Contains("JetEn")) {
      obs.metx = analysisTree->puppimet_ex_JetEnDown;
      obs.mety = analysisTree->puppimet_ey_JetEnDown;
      uncertaintyFound = true;
    }    
    else if (label.Contains("boson")) {
      int systype = 0;
      if (label.Contains("response"))
	systype = MEtSys::SysType::Response;
      else if (label.Contains("resolution"))
	systype = MEtSys::SysType::Resolution;
      TLorentzVector genV = genTools::genV(*analysisTree);
      TLorentzVector genL = genTools::genL(*analysisTree);
      float puppimet_ex = cenTree->puppimet*TMath::Cos(cenTree->puppimetphi);
      float puppimet_ey = cenTree->puppimet*TMath::Sin(cenTree->puppimetphi);
      metSys->ApplyMEtSys(puppimet_ex,
			  puppimet_ey,
			  genV.Px(),
			  genV.Py(),
			  genL.Px(),
			  genL.Py(),
			  cenTree->njets,
			  systype,
			  MEtSys::SysShift::Down,
			  obs.metx,
			  obs.mety);
      uncertaintyFound = true;
    }
    else {
      std::cout << "Systematic uncertainty " << label << std::endl;
      obs.met = cenTree->puppimet;
      obs.metphi = cenTree->puppimetphi;
    }
    if (uncertaintyFound) {
      obs.metphi = atan2(obs.mety,obs.metx);
      obs.met = sqrt(obs.metx*obs.metx+obs.mety*obs.mety);
      //      std::cout << label << " (met,metphi) : central = (" 
      //      		<< cenTree->puppimet << ","
      //      		<< cenTree->puppimetphi << ")" 
      //      		<< " down = (" << obs.met << ","
      //		<< obs.metphi << ")" << std::endl;
    }
    FillPuppiMET("Down");
  };
  
  virtual void ScaleUp() {
    bool uncertaintyFound = false;
    if (label.Contains("UnclusteredEn")) {
      obs.metx = analysisTree->puppimet_ex_UnclusteredEnUp;
      obs.mety = analysisTree->puppimet_ey_UnclusteredEnUp;
      uncertaintyFound = true;
    }
    else if (label.Contains("JetRes")) {
      obs.metx = analysisTree->puppimet_ex_JetResUp;
      obs.mety = analysisTree->puppimet_ey_JetResUp;
      uncertaintyFound = true;
   }    
    else if (label.Contains("JetEn")) {
      obs.metx = analysisTree->puppimet_ex_JetEnUp;
      obs.mety = analysisTree->puppimet_ey_JetEnUp;
      uncertaintyFound = true;
   }    
    else if (label.Contains("boson")) {
      int systype = 0;
      if (label.Contains("response"))
	systype = MEtSys::SysType::Response;
      else if (label.Contains("resolution"))
	systype = MEtSys::SysType::Resolution;
      TLorentzVector genV = genTools::genV(*analysisTree);
      TLorentzVector genL = genTools::genL(*analysisTree);
      float puppimet_ex = cenTree->puppimet*TMath::Cos(cenTree->puppimetphi);
      float puppimet_ey = cenTree->puppimet*TMath::Sin(cenTree->puppimetphi);
      metSys->ApplyMEtSys(puppimet_ex,
			  puppimet_ey,
			  genV.Px(),
			  genV.Py(),
			  genL.Px(),
			  genL.Py(),
			  cenTree->njets,
			  systype,
			  MEtSys::SysShift::Up,
			  obs.metx,
			  obs.mety);
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
      //      		<< cenTree->puppimet << ","
      //      		<< cenTree->puppimetphi << ")" 
      //      		<< " up = (" << obs.met << ","
      //		<< obs.metphi << ")" << std::endl;
    }
    FillPuppiMET("Up");
  };

};

#endif
