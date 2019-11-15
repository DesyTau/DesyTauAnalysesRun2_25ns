#ifndef NTupleMakerLepTauFunctions_h
#define NTupleMakerLepTauFunctions_h

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"

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

float get_jetPt(const AC1B *analysisTree, int jetIndex, TString JESname, TString direction, JESUncertainties * jecUncertainties){ // direction can be Up or Down
	float jetPt = -9999;
	float shift = 0.;
  
	// get the relative shift
	if (JESname == "central"){
    jetPt = analysisTree->pfjet_pt[jetIndex];
    return jetPt;
  }
	else if (JESname== "JES")  shift = analysisTree->pfjet_jecUncertainty[jetIndex];
	else if (std::find(jecUncertainties->getUncertNames().begin(), jecUncertainties->getUncertNames().end(), JESname) != jecUncertainties->getUncertNames().end()) 
      shift = jecUncertainties->getUncertainty(std::string(JESname), analysisTree->pfjet_pt[jetIndex],analysisTree->pfjet_eta[jetIndex]);
	
  // calculate shifted pt 
	if (direction == "Up") jetPt = (analysisTree->pfjet_pt[jetIndex])*(1+shift);
	else if (direction == "Down") jetPt = (analysisTree->pfjet_pt[jetIndex])*(1-shift);
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


void counting_jets(const AC1B *analysisTree, Synch17Tree *otree, const Config *cfg, const btag_scaling_inputs *inputs_btag_scaling, 
                   TString JESname = "central", TString direction = "None",  JESUncertainties * jecUncertainties = 0){

  float MaxBJetPt = 1000.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;

  vector<unsigned int> jets; jets.clear();
  vector<unsigned int> jetspt20; jetspt20.clear();
  vector<unsigned int> bjets; bjets.clear();
  vector<unsigned int> bjetsRaw; bjetsRaw.clear();
 
  int indexLeadingJet = -1;
  float ptLeadingJet = -1;

  int indexSubLeadingJet = -1;
  float ptSubLeadingJet = -1;

  int indexLeadingBJet = -1;
  float ptLeadingBJet = -1;

  int indexSubLeadingBJet = -1;
  float ptSubLeadingBJet = -1;

  TH2F* histo_tageff_ = 0;

  bool isData = cfg->get<bool>("isData");
  bool ApplyBTagScaling = cfg->get<bool>("ApplyBTagScaling");

  bool is2016 = false;
  bool is2017 = false;
  bool is2018 = false;
  int era = cfg->get<int>("era");
  if(era == 2016) is2016 = true;
  else if(era == 2017) is2017 = true;
  else if(era == 2018) is2018 = true;
  else{cout<<"no era found in cfg file, exiting"<<endl; exit(0);}
  
  const float JetEtaCut = cfg->get<float>("JetEtaCut");
  const float JetPtLowCut = cfg->get<float>("JetPtLowCut");
  const float JetPtHighCut = cfg->get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg->get<float>("dRJetLeptonCut");
  const float bJetEtaCut = cfg->get<float>("bJetEtaCut");
  const float btagCut = cfg->get<float>("btagCut");
  TString BTagDiscriminator1(cfg->get<string>("BTagDiscriminator1"));
  TString BTagDiscriminator2(cfg->get<string>("BTagDiscriminator2"));
  TString BTagDiscriminator3(cfg->get<string>("BTagDiscriminator3"));
    
  int nBTagDiscriminant1 = -1;
  int nBTagDiscriminant2 = -1;
  int nBTagDiscriminant3 = -1;
  
  for (unsigned int iBTag = 0; iBTag < analysisTree->run_btagdiscriminators->size(); ++iBTag) {
    TString discr(analysisTree->run_btagdiscriminators->at(iBTag));
            
    if (discr == BTagDiscriminator1)
      nBTagDiscriminant1 = iBTag;
    if (!is2016 && discr == BTagDiscriminator2)
      nBTagDiscriminant2 = iBTag;
    if (!is2016 && discr == BTagDiscriminator3)
      nBTagDiscriminant3 = iBTag;
  }
    
  if (nBTagDiscriminant1 == -1) {
    cout << "couldn\'t find "<< BTagDiscriminator1 << " in run_btagdiscriminators, exiting" <<endl;
    exit(-1);
  }
  if (nBTagDiscriminant2 == -1) {
    cout << "couldn\'t find "<< BTagDiscriminator2 << " in run_btagdiscriminators, exiting" <<endl;
    exit(-1);
  }
  if (nBTagDiscriminant3 == -1) {
    cout << "couldn\'t find "<< BTagDiscriminator3 << " in run_btagdiscriminators, exiting" <<endl;
    exit(-1);
  }

  for (unsigned int jet = 0; jet < analysisTree->pfjet_count; ++jet) {
    
    float jetEta    = analysisTree->pfjet_eta[jet];
    float absJetEta = fabs(analysisTree->pfjet_eta[jet]);
    if (absJetEta >= JetEtaCut) continue;
  
    float jetPt = get_jetPt(analysisTree, jet, JESname, direction, jecUncertainties);
    if (jetPt <= JetPtLowCut) continue;

    float dR1 = deltaR(analysisTree->pfjet_eta[jet], analysisTree->pfjet_phi[jet], otree->eta_1, otree->phi_1);
    if (dR1 <= dRJetLeptonCut) continue;

    float dR2 = deltaR(analysisTree->pfjet_eta[jet], analysisTree->pfjet_phi[jet], otree->eta_2, otree->phi_2);
    if (dR2 <= dRJetLeptonCut) continue;

    // skip prefiring region for 2017:
    if(is2017 && jetPt < 50 && absJetEta > 2.65 && absJetEta < 3.139) continue; 

    // see definition in Jets.h
    bool isPFJetId = tightJetID((*analysisTree), int(jet), era);
    if (!isPFJetId) continue;
     
    jetspt20.push_back(jet);

    if (absJetEta < bJetEtaCut) { // jet within b-tagging acceptance
                
      // check if meets working point cut <=> tagged      
      bool tagged = false;
      if (is2016)
      	tagged = analysisTree->pfjet_btag[jet][nBTagDiscriminant1] > btagCut; 
      else
      	// tagged = (analysisTree->pfjet_btag[jet][nBTagDiscriminant1] + analysisTree->pfjet_btag[jet][nBTagDiscriminant2]) > btagCut;
      	tagged = (analysisTree->pfjet_btag[jet][nBTagDiscriminant1] + analysisTree->pfjet_btag[jet][nBTagDiscriminant2] + analysisTree->pfjet_btag[jet][nBTagDiscriminant3]) > btagCut;
      bool taggedRaw = tagged;
      
      if(!isData && ApplyBTagScaling) {
      	int flavor = abs(analysisTree->pfjet_flavour[jet]);
      	double jet_scalefactor = 1;
      	double JetPtForBTag    = jetPt;
      	double tageff          = 1;

      	if (JetPtForBTag > MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
      	if (JetPtForBTag < MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
  	
        // getting SFs and efficiencies
      	if (flavor == 5) {
      	  jet_scalefactor = inputs_btag_scaling->reader_B.eval_auto_bounds("central", BTagEntry::FLAV_B, jetEta, JetPtForBTag);
      	  histo_tageff_= inputs_btag_scaling->tagEff_B;
      	}
      	else if (flavor == 4) {
      	  jet_scalefactor = inputs_btag_scaling->reader_C.eval_auto_bounds("central", BTagEntry::FLAV_C, jetEta, JetPtForBTag);
      	  histo_tageff_= inputs_btag_scaling->tagEff_C;
      	}
      	else {
      	  jet_scalefactor = inputs_btag_scaling->reader_Light.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, jetEta, JetPtForBTag);
      	  histo_tageff_= inputs_btag_scaling->tagEff_Light;
      	}
      	tageff = histo_tageff_->Interpolate(JetPtForBTag, absJetEta);
      	if (tageff < 1e-5)      tageff = 1e-5;
      	if (tageff > 0.99999)   tageff = 0.99999;
        
        // random seed
      	inputs_btag_scaling->rand->SetSeed((int)((jetEta+5)*100000));
      	double rannum = inputs_btag_scaling->rand->Rndm();

        // promote-demote method
      	if (jet_scalefactor < 1 && tagged)  { // downgrade - demote
      	  if (rannum < 1 - jet_scalefactor)  
            tagged = false;
      	}
      	if (jet_scalefactor > 1 && !tagged) { // upgrade - promote
      	  double fraction = (1.0 - jet_scalefactor)/(1.0 - 1.0 / tageff);
      	  if (rannum < fraction) tagged = true;
      	}
      }

      if (taggedRaw) bjetsRaw.push_back(jet); 
      if (tagged) {
      	bjets.push_back(jet);
        
      	if (indexLeadingBJet >= 0) {
      	  if (jetPt < ptLeadingBJet && jetPt > ptSubLeadingBJet) {
      	    indexSubLeadingBJet = jet;
      	    ptSubLeadingBJet = jetPt;
      	  }
      	}
      	if (jetPt > ptLeadingBJet) {
          indexSubLeadingBJet = indexLeadingBJet;
          ptSubLeadingBJet = ptLeadingBJet;
          indexLeadingBJet = jet;
          ptLeadingBJet = jetPt;
      	}
      }
    } // jet within b-tagging acceptance 
    
    if (indexLeadingJet >= 0) {
      if (jetPt < ptLeadingJet && jetPt > ptSubLeadingJet) {
        indexSubLeadingJet = jet;
        ptSubLeadingJet = jetPt;
      }
    }

    if (jetPt > ptLeadingJet) {
      indexSubLeadingJet = indexLeadingJet;
      ptSubLeadingJet = ptLeadingJet;
      indexLeadingJet = jet;
      ptLeadingJet = jetPt;
    }

    if (jetPt < JetPtHighCut) continue;
    jets.push_back(jet);
  } // jets loop

  otree->njets = jets.size();
  otree->njetspt20 = jetspt20.size();
  otree->nbtag = bjets.size();

  // leading b-jet variables
  if (indexLeadingBJet >= 0) {
    otree->bpt_1   = get_jetPt(analysisTree, indexLeadingBJet, JESname, direction, jecUncertainties);
    otree->beta_1  = analysisTree->pfjet_eta[indexLeadingBJet];
    otree->bphi_1  = analysisTree->pfjet_phi[indexLeadingBJet];
    otree->bcsv_1  = analysisTree->pfjet_btag[indexLeadingBJet][1] + analysisTree->pfjet_btag[indexLeadingBJet][2];
  }
  else {
    otree->bpt_1   = -10;
    otree->beta_1  = -10;
    otree->bphi_1  = -10;
    otree->bcsv_1  = -10;
  }

  // subleading b-jet variables
  if (indexSubLeadingBJet >= 0) {
    otree->bpt_2   = get_jetPt(analysisTree, indexSubLeadingBJet, JESname, direction, jecUncertainties);
    otree->beta_2  = analysisTree->pfjet_eta[indexSubLeadingBJet];
    otree->bphi_2  = analysisTree->pfjet_phi[indexSubLeadingBJet];
    otree->bcsv_2  = analysisTree->pfjet_btag[indexSubLeadingBJet][1] + analysisTree->pfjet_btag[indexSubLeadingBJet][2]; // 3 - pfDeepFlavourJetTags:probb
  }
  else {
    otree->bpt_2   = -10;
    otree->beta_2  = -10;
    otree->bphi_2  = -10;
    otree->bcsv_2  = -10;
  }
  
  // leading jet variables
  if ( indexLeadingJet >= 0 && indexSubLeadingJet >= 0 && indexLeadingJet == indexSubLeadingJet )
    cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
  if (indexLeadingJet >= 0) {
    otree->jpt_1 = get_jetPt(analysisTree, indexLeadingJet, JESname, direction, jecUncertainties);
    otree->jeta_1 = analysisTree->pfjet_eta[indexLeadingJet];
    otree->jphi_1 = analysisTree->pfjet_phi[indexLeadingJet];
    otree->jcsv_1 = analysisTree->pfjet_btag[indexLeadingJet][1] + analysisTree->pfjet_btag[indexLeadingJet][2];
  }
  else {
    otree->jpt_1 = -10;
    otree->jeta_1 = -10;
    otree->jphi_1 = -10;
    otree->jcsv_1 = -10;
  }
  
  // subleading jet variables
  if (indexSubLeadingJet >= 0) {
    otree->jpt_2 = get_jetPt(analysisTree, indexSubLeadingJet, JESname, direction, jecUncertainties);
    otree->jeta_2 = analysisTree->pfjet_eta[indexSubLeadingJet];
    otree->jphi_2 = analysisTree->pfjet_phi[indexSubLeadingJet];
    otree->jcsv_2 = analysisTree->pfjet_btag[indexSubLeadingJet][1] + analysisTree->pfjet_btag[indexSubLeadingJet][2];
  }
  else {
    otree->jpt_2 = -10;
    otree->jeta_2 = -10;
    otree->jphi_2 = -10;
    otree->jcsv_2 = -10;
  }

  // dijet variables
  if (indexLeadingJet >= 0 && indexSubLeadingJet >= 0) {
    otree->njetingap = 0;
    TLorentzVector jet1; 
    TLorentzVector jet2; 

  	jet1.SetPtEtaPhiE(otree->jpt_1, otree->jeta_1, otree->jphi_1, get_jetE(analysisTree, indexLeadingJet, JESname, direction, jecUncertainties));
  	jet2.SetPtEtaPhiE(otree->jpt_2, otree->jeta_2, otree->jphi_2, get_jetE(analysisTree, indexSubLeadingJet, JESname, direction, jecUncertainties));

    otree->mjj = (jet1 + jet2).M();
    otree->dijetpt = (jet1 + jet2).Pt();
    otree->dijeteta = (jet1 + jet2).Eta();
    otree->dijetphi = (jet1 + jet2).Phi();
    otree->jdeta = abs(analysisTree->pfjet_eta[indexLeadingJet] - analysisTree->pfjet_eta[indexSubLeadingJet]);
    otree->jdphi = dPhiFrom2P(jet1.Px(), jet1.Py(), jet2.Px(),jet2.Py());
   
    float etamax = analysisTree->pfjet_eta[indexLeadingJet];
    float etamin = analysisTree->pfjet_eta[indexSubLeadingJet];
    if (etamax < etamin) {
      float tmp = etamax;
      etamax = etamin;
      etamin = tmp;
    }
    for (unsigned int jet = 0; jet < jets.size(); ++jet) {
      int index = jets.at(jet);
      float etaX = analysisTree->pfjet_eta[index];
      if ((index != indexLeadingJet) && (index != indexSubLeadingJet) && (etaX > etamin) && (etaX < etamax)) 
        otree->njetingap++;
    }
  }
  else {
    otree->mjj =  -10;
    otree->jdeta =  -10;
    otree->dijetpt =  -10;
    otree->jdphi =  -10;
    otree->njetingap = -10;
  }
 } // counting_jets
} // jets namespace 

#endif
