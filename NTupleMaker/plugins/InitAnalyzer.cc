// -*- C++ -*-
//
// Package:    Demo/InitAnalyzer
// Class:      InitAnalyzer
// 
/**\class InitAnalyzer InitAnalyzer.cc Demo/InitAnalyzer/plugins/InitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexis Kalogeropoulos
//         Created:  Mon, 28 Sep 2015 11:33:21 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

//
// class declaration
//

using namespace reco;
using namespace std;
class InitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit InitAnalyzer(const edm::ParameterSet&);
      ~InitAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  TTree* tree0;
  TH1D*  nEvents;
  TH1D*  nPU_TrueNumInteractions;
  TH1D*  nPU_NumInteractionsBX0;
  TH1D*  nPU_NumInteractionsBXm1;
  TH1D*  nPU_NumInteractionsBXp1;

  bool cdata;
  bool cFastSim;
  bool cgen;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct > GenToken_;
  Float_t genweight;
  Float_t fill_weight;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
InitAnalyzer::InitAnalyzer(const edm::ParameterSet& iConfig) :
  cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
  cFastSim(iConfig.getUntrackedParameter<bool>("IsFastSim", false)),
  cgen(iConfig.getUntrackedParameter<bool>("GenParticles", false)),
  PUInfoToken_( consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"))),
  GenToken_( consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
	
}

InitAnalyzer::~InitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
InitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   fill_weight = +1;
   if(!cdata && genweight<0) fill_weight = -1;
   nEvents->Fill(0.,fill_weight);

  if (!cdata) {
     edm::Handle<vector<PileupSummaryInfo> > PUInfo;
     iEvent.getByToken( PUInfoToken_, PUInfo);
      if(PUInfo.isValid())
	{
	  for(vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	    {
	      int BX = PVI->getBunchCrossing();
	      if(BX == -1)
		{ 
		 nPU_NumInteractionsBXm1->Fill(PVI->getPU_NumInteractions());
		}
	      else if(BX == 0)
		{ 
		  nPU_NumInteractionsBX0->Fill(PVI->getPU_NumInteractions());	        
		  nPU_TrueNumInteractions->Fill(PVI->getTrueNumInteractions());
		}
	      else if(BX == 1)
		{ 
		 nPU_NumInteractionsBXp1->Fill(PVI->getPU_NumInteractions());
		}
	      
	    }
	}
  }

  genweight = 1.;
  if(cgen && !cdata)
    {
      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByToken( GenToken_, HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  //  	  cout << "Event weight from HEPMC : " << genweight << endl;
	}
    }
  tree0->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
InitAnalyzer::beginJob()
{
  edm::Service<TFileService> FS;
  tree0 = FS->make<TTree>("AC1B", "AC1B", 1);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);
  if (!cdata){
    nPU_TrueNumInteractions = FS->make<TH1D>("nPU_TrueNumInteractions","nPU_TrueNumInteractions",100, 0., 100);
    nPU_NumInteractionsBX0 = FS->make<TH1D>("nPU_NumInteractionsBX0","nPU_NumInteractionsBX0",100, 0., 100);
    nPU_NumInteractionsBXm1 = FS->make<TH1D>("nPU_NumInteractionsBXm1","nPU_NumInteractionsBXm1",100, 0., 100);
    nPU_NumInteractionsBXp1 = FS->make<TH1D>("nPU_NumInteractionsBXp1","nPU_NumInteractionsBXp1",100, 0., 100);
  }
  if (cgen && !cdata) {
    tree0->Branch("genweight", &genweight, "genweight/F");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void 
InitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
InitAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
InitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
InitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
InitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InitAnalyzer);
