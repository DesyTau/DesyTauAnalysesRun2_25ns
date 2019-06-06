#ifndef NTupleMakerLepTauFunctions_h
#define NTupleMakerLepTauFunctions_h

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"

struct btag_scaling_inputs{
  BTagCalibrationReader reader_B;
  BTagCalibrationReader reader_C;
  BTagCalibrationReader reader_Light;
  TH2F *tagEff_B;
  TH2F *tagEff_C;
  TH2F *tagEff_Light;
  TRandom3 *rand;
};

namespace jets{

JESUncertainties * dummyJEC = 0;

float get_jetPt(const AC1B *analysisTree, int jetIndex, TString JESname, TString direction, JESUncertainties * jecUncertainties){ // direction can be Up or Down

	float jetPt = -9999;
	float shift = 0.;
	// get the relative shift
	if (JESname == "central")  { jetPt = analysisTree->pfjet_pt[jetIndex]; return jetPt;}
	else if (JESname== "JES")  shift = analysisTree->pfjet_jecUncertainty[jetIndex];
	else if (std::find(jecUncertainties->getUncertNames().begin(), jecUncertainties->getUncertNames().end(), JESname) != jecUncertainties->getUncertNames().end()) 
		{
         shift = jecUncertainties->getUncertainty(std::string(JESname), analysisTree->pfjet_pt[jetIndex],analysisTree->pfjet_eta[jetIndex]);
		 // cout for debugging
		 /*std::cout << "------------------------------------------------" << endl ;			    
         std::cout << "in get jet pt -- " << " name: " << JESname ;			
 	     std::cout << " pt of jet #" << jetIndex << " :: before correction = " <<  analysisTree->pfjet_pt[jetIndex] << std::endl ;
		 std::cout << " shift  = " << shift << std::endl; */
		}
	// calculate shifted pt 
	if (direction == "Up") jetPt = (analysisTree->pfjet_pt[jetIndex])*(1+shift);
	else if (direction == "Down") jetPt = (analysisTree->pfjet_pt[jetIndex])*(1-shift);
	//std::cout << " shifted by " << shift << " in direction" << direction << std::endl;
    //std::cout << "------------------------------------------------" << endl ;			    
	return jetPt;
};


float get_jetE(const AC1B *analysisTree, int jetIndex, TString JESname, TString direction, JESUncertainties * jecUncertainties){
	float jetE = -9999;
	float shift = 0.;
	// get the relative shift
	if (JESname == "central")  { jetE = analysisTree->pfjet_e[jetIndex]; return jetE;}
	else if (JESname== "JES")   shift = analysisTree->pfjet_jecUncertainty[jetIndex];
	else if (std::find(jecUncertainties->getUncertNames().begin(), jecUncertainties->getUncertNames().end(), JESname) != jecUncertainties->getUncertNames().end()) 
		shift = jecUncertainties->getUncertainty(std::string(JESname), analysisTree->pfjet_pt[jetIndex],analysisTree->pfjet_eta[jetIndex]);
	// calculate shifted energy
	if (direction == "Up")  jetE =  (analysisTree->pfjet_e[jetIndex])*(1+shift);
	if (direction== "Down") jetE =  (analysisTree->pfjet_e[jetIndex])*(1-shift);
	return jetE;
};


void counting_jets(const AC1B *analysisTree, Synch17Tree *otree, const Config *cfg, const btag_scaling_inputs *inputs_btag_scaling, TString JESname = "central", TString direction = "None",  JESUncertainties * jecUncertainties = dummyJEC){

  //std::cout << "JESname: " << JESname << "| direction " << direction << std::endl;

  vector<unsigned int> jets; jets.clear();
  vector<unsigned int> jetspt20; jetspt20.clear();
  vector<unsigned int> bjets; bjets.clear();

  int indexLeadingJet = -1;
  float ptLeadingJet = -1;

  int indexSubLeadingJet = -1;
  float ptSubLeadingJet = -1;

  int indexLeadingBJet = -1;
  float ptLeadingBJet = -1;

  int indexSubLeadingBJet = -1;
  float ptSubLeadingBJet = -1;

  TH2F* histo_tageff_ = 0;

  for (unsigned int jet=0; jet<analysisTree->pfjet_count; ++jet) {

    float jetEta    = analysisTree->pfjet_eta[jet];
    float absJetEta = fabs(analysisTree->pfjet_eta[jet]);
    if (absJetEta>=cfg->get<float>("JetEtaCut")) continue;

    //float jetPt = analysisTree->pfjet_pt[jet];
    float jetPt = get_jetPt(analysisTree, jet, JESname, direction, jecUncertainties);
    if (jetPt<=cfg->get<float>("JetPtLowCut")) continue;

    float dR1 = deltaR(analysisTree->pfjet_eta[jet],analysisTree->pfjet_phi[jet],otree->eta_1,otree->phi_1);
    if (dR1<=cfg->get<float>("dRJetLeptonCut")) continue;

    float dR2 = deltaR(analysisTree->pfjet_eta[jet],analysisTree->pfjet_phi[jet],otree->eta_2,otree->phi_2);
    if (dR2<=cfg->get<float>("dRJetLeptonCut")) continue;

    // jetId
	// ALSO PROPAGATED HERE?
    float energy = analysisTree->pfjet_e[jet];
    energy *= analysisTree->pfjet_energycorr[jet];
    float chf = analysisTree->pfjet_chargedhadronicenergy[jet]/energy;
    float nhf = analysisTree->pfjet_neutralhadronicenergy[jet]/energy;
    float phf = analysisTree->pfjet_neutralemenergy[jet]/energy;
    float elf = analysisTree->pfjet_chargedemenergy[jet]/energy;
    float muf = analysisTree->pfjet_muonenergy[jet]/energy;
    float chm = analysisTree->pfjet_chargedmulti[jet];
    float nm  = analysisTree->pfjet_neutralmulti[jet];
    float npr = analysisTree->pfjet_chargedmulti[jet] + analysisTree->pfjet_neutralmulti[jet];
    
    bool isPFJetId = false;
    if (absJetEta<=2.7)
      isPFJetId = (nhf < 0.99 && phf < 0.99 && npr > 1) && (absJetEta>2.4 || (chf>0 && chm > 0 && elf < 0.99));
    else if (absJetEta<=3.0)
      isPFJetId = (phf < 0.9 && nm > 2);
    else
      isPFJetId = phf < 0.9 && nm > 10;
    
    if (!isPFJetId) continue;

    jetspt20.push_back(jet);

    if (absJetEta < cfg->get<float>("bJetEtaCut")) { // jet within b-tagging acceptance
      
      bool tagged = ( analysisTree->pfjet_btag[jet][0]>cfg->get<float>("btagCut") );

      if(!cfg->get<bool>("isData")  && cfg->get<bool>("ApplyBTagScaling")) {

	int flavor = abs(analysisTree->pfjet_flavour[jet]);
	double jet_scalefactor = 1;
	double JetPtForBTag    = jetPt;
	double tageff          = 1;
	 // std::cout << "before inputs_btag_scaling" << std::endl;
	if (flavor==5) {
	  jet_scalefactor = inputs_btag_scaling->reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, jetEta, JetPtForBTag);

	  histo_tageff_=inputs_btag_scaling->tagEff_B;
	}
	else if (flavor==4) {
	  jet_scalefactor = inputs_btag_scaling->reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, jetEta, JetPtForBTag);
	  histo_tageff_=inputs_btag_scaling->tagEff_C;

	}
	else {
	  jet_scalefactor = inputs_btag_scaling->reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetEta, JetPtForBTag);
	  histo_tageff_=inputs_btag_scaling->tagEff_Light;
	}
	
	if(JetPtForBTag > histo_tageff_->GetXaxis()->GetBinLowEdge(histo_tageff_->GetNbinsX()+1)){
	  tageff = histo_tageff_->GetBinContent(histo_tageff_->GetNbinsX(),histo_tageff_->GetYaxis()->FindBin(absJetEta));
	}
	else{
	  tageff = histo_tageff_->GetBinContent(histo_tageff_->GetXaxis()->FindBin(JetPtForBTag), histo_tageff_->GetYaxis()->FindBin(absJetEta));
	}

	if (tageff<1e-5)      tageff = 1e-5;
	if (tageff>0.99999)   tageff = 0.99999;
	inputs_btag_scaling->rand->SetSeed((int)((jetEta+5)*100000));
	double rannum = inputs_btag_scaling->rand->Rndm();

	if (jet_scalefactor<1 && tagged)  { // downgrade - demote
	  if (rannum<1-jet_scalefactor)  tagged = false;
	}
	if (jet_scalefactor>1 && !tagged) { // upgrade - promote
	  double fraction = (1.0-jet_scalefactor)/(1.0-1.0/tageff);
	  if (rannum<fraction) tagged = true;
	}
      }
      if (tagged) {

	bjets.push_back(jet);

	if (indexLeadingBJet>=0) {
	  if (jetPt<ptLeadingBJet && jetPt>ptSubLeadingBJet) {
	    indexSubLeadingBJet = jet;
	    ptSubLeadingBJet = jetPt;
	  }
	}
	if (jetPt>ptLeadingBJet) {
	  ptLeadingBJet = jetPt;
	  indexLeadingBJet = jet;
	}
      }
    } //if (absJetEta < cfg->get<float>("bJetEtaCut"))
	  //std::cout << "after inputs_btag_scaling" << std::endl;
    if (indexLeadingJet>=0) {
      if (jetPt<ptLeadingJet && jetPt>ptSubLeadingJet) {
        indexSubLeadingJet = jet;
        ptSubLeadingJet = jetPt;
      }
    }

    if (jetPt>ptLeadingJet) {
      indexLeadingJet = jet;
      ptLeadingJet = jetPt;
    }

    if (jetPt<cfg->get<float>("JetPtHighCut")) continue;
    jets.push_back(jet);
  }

  otree->njets = jets.size();
  otree->njetspt20 = jetspt20.size();
  otree->nbtag = bjets.size();
  
  otree->bpt_1   = -9999;
  otree->beta_1  = -9999;
  otree->bphi_1  = -9999;
  otree->bcsv_1  = -9999;
  
  if (indexLeadingBJet>=0) {
    otree->bpt_1   = get_jetPt(analysisTree, indexLeadingBJet, JESname, direction, jecUncertainties);
    otree->beta_1  = analysisTree->pfjet_eta[indexLeadingBJet];
    otree->bphi_1  = analysisTree->pfjet_phi[indexLeadingBJet];
  
    otree->bcsv_1  = analysisTree->pfjet_btag[indexLeadingBJet][0];
  }

  otree->bpt_2   = -9999;
  otree->beta_2  = -9999;
  otree->bphi_2  = -9999;
  otree->bcsv_2  = -9999;
  
  if (indexSubLeadingBJet>=0) {
    otree->bpt_2   = get_jetPt(analysisTree, indexSubLeadingBJet, JESname, direction, jecUncertainties);
    otree->beta_2  = analysisTree->pfjet_eta[indexSubLeadingBJet];
    otree->bphi_2  = analysisTree->pfjet_phi[indexSubLeadingBJet];
    otree->bcsv_2  = analysisTree->pfjet_btag[indexSubLeadingBJet][0];
  }

  otree->jpt_1 = -9999;
  otree->jeta_1 = -9999;
  otree->jphi_1 = -9999;
  
  if ( indexLeadingJet>=0 && indexSubLeadingJet>=0 && indexLeadingJet==indexSubLeadingJet )
cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;

  if (indexLeadingJet>=0) {
	//std::cout << "pt of jet #" << indexLeadingJet << " :: before correction = " <<  analysisTree->pfjet_pt[indexLeadingJet] ;
	//std::cout << " Leading jet --------------------- " << std::endl;
    otree->jpt_1 = get_jetPt(analysisTree, indexLeadingJet, JESname, direction, jecUncertainties); //analysisTree->pfjet_pt[indexLeadingJet];
	//std::cout << " after correction = " << otree->jpt_1 << std::endl;
	//std::cout << "  --------------------- EOF Leading jet " << std::endl;
    otree->jeta_1 = analysisTree->pfjet_eta[indexLeadingJet];
    otree->jphi_1 = analysisTree->pfjet_phi[indexLeadingJet];
  }

  otree->jpt_2 = -9999;
  otree->jeta_2 = -9999;
  otree->jphi_2 = -9999;
  
  if (indexSubLeadingJet>=0) {
    otree->jpt_2 = get_jetPt(analysisTree, indexSubLeadingJet, JESname, direction, jecUncertainties);
    otree->jeta_2 = analysisTree->pfjet_eta[indexSubLeadingJet];
    otree->jphi_2 = analysisTree->pfjet_phi[indexSubLeadingJet];
  }

  otree->mjj =  -9999;
  otree->jdeta =  -9999;
  otree->dijetpt =  -9999;
  otree->jdphi =  -9999;
  otree->njetingap = -1;

  if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
    otree->njetingap = 0;
    TLorentzVector jet1; 
    TLorentzVector jet2; 

	jet1.SetPtEtaPhiE(otree->jpt_1, otree->jeta_1, otree->jphi_1, get_jetE(analysisTree, indexLeadingJet, JESname, direction, jecUncertainties));
	jet2.SetPtEtaPhiE(otree->jpt_2, otree->jeta_2, otree->jphi_2, get_jetE(analysisTree, indexSubLeadingJet, JESname, direction, jecUncertainties));

    otree->mjj = (jet1+jet2).M();
    otree->dijetpt = (jet1+jet2).Pt();//Merijn: added, needed for DNN inputs in mt
    otree->jdeta = abs(analysisTree->pfjet_eta[indexLeadingJet]-
          analysisTree->pfjet_eta[indexSubLeadingJet]);
    otree->jdphi = dPhiFrom2P(jet1.Px(),jet1.Py(),
			      jet2.Px(),jet2.Py());
   
    float etamax = analysisTree->pfjet_eta[indexLeadingJet];
    float etamin = analysisTree->pfjet_eta[indexSubLeadingJet];
    if (etamax<etamin) {
      float tmp = etamax;
      etamax = etamin;
      etamin = tmp;
    }
    for (unsigned int jet=0; jet<jets.size(); ++jet) {
      int index = jets.at(jet);
      float etaX = analysisTree->pfjet_eta[index];
      if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax) 
        otree->njetingap++;
    }


  }

}




} // end of leptau namespace 

#endif
