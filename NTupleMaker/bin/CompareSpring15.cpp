/////////////////////////////////////////////////////////////
// Analysis Macro for Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
//
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
/////////////////////////////////////////////////////////////

#include <iostream>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TMath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Spring15Tree.h"

bool areEqual(double x1, double x2){
  if (x1*x2 < 0.) return false;
  if (fabs(log10(fabs(x1))-log10(fabs(x2))) > 1.e-5) return false;
  return true;
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

struct sort_second {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};

void load_entries( Spring15Tree& t, std::vector< std::pair< int, int> >& entry){
  while(t.GetEntry(-1) > 0)
    entry.push_back(std::make_pair(t.LoadedEntryId(), t.evt));

  std::sort(entry.begin(), entry.end(), sort_second());
}


// main

int main(int argc, char * argv[]) {

  // first argument - ref file
  // second argument - test file
  
  using namespace std;

  if(argc < 2) return -1;  

  TChain* tref = new TChain("TauCheck");
  TChain* ttest = new TChain("TauCheck");
  
  if(argc < 2) return -1;

  std::cout<<"ref: "<< argv[1]<<std::endl;
  std::cout<<"test: "<< argv[2]<<std::endl;  
  
  TFile* fref = new TFile(argv[1], "read");
  if(!fref->FindObjectAny("TauCheck")) return -2;
  fref->Close();
  
  tref->Add(argv[1]);
  if(!tref) return -3;
  
  TFile* ftest = new TFile(argv[2], "read");
  if(!ftest->FindObjectAny("TauCheck")) return -2;
  ftest->Close();
  
  ttest->Add(argv[2]);
  if(!ttest) return -3;
  
  Spring15Tree* ref = new Spring15Tree(tref);
  std::vector<std::pair<int, int> > ref_entry;  
  load_entries(*ref, ref_entry);
  
  Spring15Tree* test = new Spring15Tree(ttest);
  std::vector<std::pair<int, int> > test_entry;  
  load_entries(*test, test_entry);

  if(test_entry.size() != ref_entry.size())
    std::cout<<"Ntry mismatch: ref = "<<ref_entry.size()<<"; test = "<<test_entry.size()<<";"<<std::endl;

  UInt_t iref = 0;
  UInt_t itest = 0;

  int ref_id = 0;
  int test_id = 0;

  std::vector<std::pair<int, int> > ref_missing; 
  std::vector<std::pair<int, int> > test_missing; 

  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_events;
  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_leptons;
  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_met;
  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_mvamet;
  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_jets;
  
  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > good_events;  

  bool check_mvamet = true;
  
  Float_t ref_min_pt_1 = 9999999999.;
  Float_t ref_min_pt_2 = 9999999999.;

  Float_t test_min_pt_1 = 9999999999.;
  Float_t test_min_pt_2 = 9999999999.;  

  for( ; iref < ref_entry.size(); iref++){
    ref->GetEntry(ref_entry.at(iref).first);

    if(ref->pt_1 < ref_min_pt_1)
      ref_min_pt_1 = ref->pt_1;

    if(ref->pt_2 < ref_min_pt_2)
      ref_min_pt_2 = ref->pt_2;	
  }
  
  for( ; itest < test_entry.size(); itest++){
    test->GetEntry(test_entry.at(itest).first);

    if(test->pt_1 < test_min_pt_1)
      test_min_pt_1 = test->pt_1;

    if(test->pt_2 < test_min_pt_2)
      test_min_pt_2 = test->pt_2;	
  }

  std::cout<<"ref: "<<ref_min_pt_1<<" "<<ref_min_pt_2<<std::endl;
  std::cout<<"test: "<<test_min_pt_1<<" "<<test_min_pt_2<<std::endl;  

  itest = 0;
  iref = 0;
  
  for( ; iref < ref_entry.size(); iref++){
    ref_id = ref_entry.at(iref).second;

    for(; itest < test_entry.size(); itest++){
      test_id = test_entry.at(itest).second;

      if(test_id < ref_id)
	ref_missing.push_back(test_entry.at(itest));
      else if(test_id > ref_id){
	test_missing.push_back(ref_entry.at(iref));
	break;
      }
      else
	break;
    }

    if (test_id != ref_id)
      continue;
    
    // lets compare!
    bool isGood = true;
    
    ref->GetEntry(ref_entry.at(iref).first);
    test->GetEntry(test_entry.at(itest).first);
    
    // compare leptons
    if( !areEqual(ref->pt_1, test->pt_1) ||
	!areEqual(ref->iso_1, test->iso_1) ||
	!areEqual(ref->pt_2, test->pt_2) ||
	!areEqual(ref->iso_2, test->iso_2)){

      isGood = false;
      
      std::cout<<"Event "<<ref_entry.at(iref).second<<":"<<std::endl;
      std::cout<<"    lep1: pt="<<ref->pt_1<<" "<<test->pt_1<<" eta="<<ref->eta_1<<" "<<test->eta_1
	       <<" phi="<<ref->phi_1<<" "<<test->phi_1<<" iso="<<ref->iso_1<<" "<<test->iso_1
	       <<" q="<<ref->q_1<<" "<<test->q_1<<std::endl;

      std::cout<<"    lep2: pt="<<ref->pt_2<<" "<<test->pt_2<<" eta="<<ref->eta_2<<" "<<test->eta_2
	       <<" phi="<<ref->phi_2<<" "<<test->phi_2<<" iso="<<ref->iso_2<<" "<<test->iso_2
	       <<" q="<<ref->q_2<<" "<<test->q_2<<std::endl;


      //std::cout<<ref->againstElectronVLooseMVA5_2<<" "<<ref->againstMuonTight3_2<<std::endl;
      //std::cout<<test->againstElectronVLooseMVA5_2<<" "<<test->againstMuonTight3_2<<std::endl;
      
      std::cout<<"    deltaR(lep1, lep2)="<<deltaR(ref->eta_1, ref->phi_1, ref->eta_2, ref->phi_2)
	       <<" "<<deltaR(test->eta_1, test->phi_1, test->eta_2, test->phi_2)<<std::endl;

      bad_leptons.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));
    }

    
    // compare pfmet
    if( !areEqual(ref->met, test->met) ||
	!areEqual(ref->metphi, test->metphi) ||
	!areEqual(ref->metcov00, test->metcov00) ||
	!areEqual(ref->metcov01, test->metcov01) ||
	!areEqual(ref->metcov10, test->metcov10) ||
	!areEqual(ref->metcov11, test->metcov11)) {
      isGood = false;
      
      std::cout<<"Event "<<ref_entry.at(iref).second<<":"<<std::endl;
      std::cout<<"    met: |met|="<<ref->met<<" "<<test->met<<" phi="<<ref->metphi<<" "<<test->metphi
	       <<" cov00="<<ref->metcov00<<" "<<test->metcov00<<" cov01="<<ref->metcov01<<" "<<test->metcov01
	       <<" cov10="<<ref->metcov10<<" "<<test->metcov10<<" cov11="<<ref->metcov11<<" "<<test->metcov11
	       <<std::endl;
      
      bad_met.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));
    }

    // compare mva met
    if( check_mvamet &&
	(!areEqual(ref->mvamet, test->mvamet) ||
	 !areEqual(ref->mvametphi, test->mvametphi) ||
	 !areEqual(ref->mvacov00, test->mvacov00) ||
	 !areEqual(ref->mvacov01, test->mvacov01) ||
	 !areEqual(ref->mvacov10, test->mvacov10) ||
	 !areEqual(ref->mvacov11, test->mvacov11))) {
      isGood = false;
      
      std::cout<<"Event "<<ref_entry.at(iref).second<<":"<<std::endl;
      std::cout<<"    mvamet: |met|="<<ref->mvamet<<" "<<test->mvamet<<" phi="<<ref->mvametphi<<" "<<test->mvametphi
	       <<" cov00="<<ref->mvacov00<<" "<<test->mvacov00<<" cov01="<<ref->mvacov01<<" "<<test->mvacov01
	       <<" cov10="<<ref->mvacov10<<" "<<test->mvacov10<<" cov11="<<ref->mvacov11<<" "<<test->mvacov11
	       <<std::endl;
      
      bad_mvamet.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));
    }

    if(isGood)
      good_events.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));
    else
      bad_events.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));
    
    // compare jets
    /*if (!areEqual(ref->njets, test->njets) ||
	!areEqual(ref->njetspt20, test->njetspt20) ||
	!areEqual(ref->nbtag, test->nbtag) ||
	!areEqual(ref->jpt_1, test->jpt_1) ||
	!areEqual(ref->jpt_2, test->jpt_2)) {*/
    /*if (ref->njets > test->njets){
      std::cout<<"Event "<<ref_entry.at(iref).second<<":"<<std::endl;
      std::cout<<"    njets = "<<ref->njets<<" "<<test->njets<<" njetspt20="<<ref->njetspt20<<" "<<test->njetspt20
	       <<" nbtag="<<ref->nbtag<<" "<<test->nbtag<<std::endl;
      std::cout<<"    jpt_1 = "<<ref->jpt_1<<" "<<test->jpt_1<<" eta_1="<<ref->jeta_1<<" "<<test->jeta_1<<std::endl;
      std::cout<<"    jpt_2 = "<<ref->jpt_2<<" "<<test->jpt_2<<" eta_2="<<ref->jeta_2<<" "<<test->jeta_2<<std::endl;

      
      std::cout<<"    deltaR(jet1, lep1)="<<deltaR(ref->jeta_1, ref->jphi_1, ref->eta_1, ref->phi_1)
	       <<" "<<deltaR(test->jeta_1, test->jphi_1, test->eta_1, test->phi_1)<<std::endl;

      std::cout<<"    deltaR(jet1, lep2)="<<deltaR(ref->jeta_1, ref->jphi_1, ref->eta_2, ref->phi_2)
	       <<" "<<deltaR(test->jeta_1, test->jphi_1, test->eta_2, test->phi_2)<<std::endl;
      
      std::cout<<"    deltaR(jet2, lep1)="<<deltaR(ref->jeta_2, ref->jphi_2, ref->eta_1, ref->phi_1)
	       <<" "<<deltaR(test->jeta_2, test->jphi_2, test->eta_1, test->phi_1)<<std::endl;

      std::cout<<"    deltaR(jet2, lep2)="<<deltaR(ref->jeta_2, ref->jphi_2, ref->eta_2, ref->phi_2)
	       <<" "<<deltaR(test->jeta_2, test->jphi_2, test->eta_2, test->phi_2)<<std::endl;
      
      bad_jets.push_back(std::make_pair( ref_entry.at(iref), test_entry.at(itest)));      
    }




    
    
    /*
    // compare jet1
    if( !areEqual(ref->jpt_1, test->jpt_1))
      std::cout<<"    jpt_1 = "<<ref->jpt_1<<" "<<test->jpt_1<<std::endl;

    // compare jet2
    if( !areEqual(ref->jpt_2, test->jpt_2))
      std::cout<<"    jpt_2 = "<<ref->jpt_2<<" "<<test->jpt_2<<std::endl;    
    */
    
    itest++;
  }
  std::cout<<std::endl;

  std::cout<<test_entry.size()-ref_missing.size()-bad_events.size()<<" events with perfect match."<<std::endl;
  std::cout<<ref_entry.size()-test_missing.size()-bad_events.size()<<" events with perfect match."<<std::endl;  
  std::cout<<std::endl;
  
  std::cout<<test_missing.size()<<" events missing in test tree:"<<std::endl;
  for( UInt_t i = 0; i < test_missing.size(); i++)
    std::cout<<test_missing.at(i).second<<", ";
  std::cout<<std::endl<<std::endl;

  std::cout<<ref_missing.size()<<" events missing in ref tree:"<<std::endl;
  for( UInt_t i = 0; i < ref_missing.size(); i++)
    std::cout<<ref_missing.at(i).second<<", ";
  std::cout<<std::endl<<std::endl;

  std::cout<<bad_events.size()<<" events show differences:"<<std::endl;
  for( UInt_t i = 0; i < bad_events.size(); i++)
    std::cout<<bad_events.at(i).second.second<<", ";
  std::cout<<std::endl<<std::endl;
    
  return 0;
}
