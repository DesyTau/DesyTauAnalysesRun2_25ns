#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTBase.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TSystem.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

//#include "RecoilCorrector.hh"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
//#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
//#include "TauAnalysis/CandidateTools/interface/neuralMtautauAuxFunctions.h"
//#include "BtagSF.hh" 
//#include "readJSONFile.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/VectorUtil.h"

#include <vector>
#include <utility>
#include <map>
#include <algorithm>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#define SAVE   true
#define MINPt1 20.0 
#define MINPt2 20.0
//#define PtVETO 20.0
#define PtVETO 30.0
//MB#define MAXEta  5.0 
#define MAXEta  4.7 //MB
#define MINJetID 0.5
#define DEBUG false

using namespace ROOT::Math;
using namespace std;

struct DiTauInfo 
{ 
  DiTauInfo(){}; 
  int diTauCharge_; 
  double sumPt_; 
  double sumIso_;
  int index1_;
  int index2_;
}; 

struct SortDiTauPairs 
{ 
  bool operator() (const DiTauInfo t1, const DiTauInfo t2) 
  { 
    // 1st criterion: OS 
    //if ( t1.diTauCharge_ < t2.diTauCharge_ ) return true; 
    //if ( t1.diTauCharge_ > t2.diTauCharge_ ) return false; 
    if(t1.sumIso_ < t2.sumIso_ ) return true;
    if(t1.sumIso_ > t2.sumIso_ ) return false;
    // 2nd criterion: sumPt of diTau pair 
    return (t1.sumPt_ > t2.sumPt_);  
  } 
}; 

typedef map< float , int > MAPDITAU_etaL2;
typedef map< float , MAPDITAU_etaL2 > MAPDITAU_etaL1;
typedef map< float , MAPDITAU_etaL1 > MAPDITAU_ptL2;
typedef map< float , MAPDITAU_ptL2 > MAPDITAU_ptL1;
typedef map< int , MAPDITAU_ptL1 > MAPDITAU_event;
typedef map< int , MAPDITAU_event > MAPDITAU_lumi;
typedef map< int , MAPDITAU_lumi > MAPDITAU_run;

typedef std::vector<std::string> vstring;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

//edm::LumiReWeighting *LumiWeights_ = new edm::LumiReWeighting("/nfs/dust/cms/user/anayak/CMS/MyHTTAnalysis/data/MC_Summer12_PU_S10-600bins.root",
//"/nfs/dust/cms/user/anayak/CMS/MyHTTAnalysis/data/Data_Pileup_2012_ReRecoPixel-600bins.root","pileup","pileup");

//RecoilCorrector *RecoilCorrector_;
enum MotherNames{HIGGS=1, WBOSON, ZBOSON, TAU};
enum MvaMetChannel{EMU=1, ETAU, MUTAU, TAUTAU, MUMU, EE, UNKNOWN};
enum MCEVENTTYPE{HADHAD, LEPHAD, LEPLEP, LL, LNU, LEPNU, TAUNU, OTHERS};

enum BVariation{kNo = 0, kDown = 1, kUp = 2};
//BtagSF* btsf = new BtagSF(12345);

//Reading json file and filter lumi section
std::map<int, std::vector<std::pair<int, int> > >
readJSONFile(const std::string& inFileName)
{
  std::ifstream inFile(inFileName.c_str(), std::ios::in);
  
  std::string line;
  while(!inFile.eof())
    {
      std::string buffer;
      inFile >> buffer;
      line += buffer;
    }
  
  // define map with result
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  
  // loop on JSON file
  for(std::string::const_iterator it = line.begin(); it < line.end(); ++it)
    {
      // find run number
      if( (*(it) == '"') && (*(it+7) == '"') )   
        {
	  std::string run(it+1, it+7);
          //std::cout << "found run " << run << std::endl;
	  
          // find lumi sections
	  std::vector<std::pair<int, int> > lumisections;
          for(std::string::const_iterator it2 = it+10; it2 < line.end(); ++it2)
            {
              if( (*(it2) == ']') && (*(it2-1) == ']') ) break;
              if( *(it2) != '[' ) continue;
              
	      std::string::const_iterator it_beg = it2;
	      std::string::const_iterator it_mid;
	      std::string::const_iterator it_end;
	      
              for(std::string::const_iterator it3 = it_beg; it3 < line.end(); ++it3)
                {
                  if( *(it3) == ',' ) it_mid = it3;
                  if( *(it3) == ']' )
                    {
                      it_end = it3;
                      break;
                    }
                }
	      std::string lumi_beg(it_beg+1, it_mid);
	      std::string lumi_end(it_mid+1, it_end);
              //std::cout << "[" << lumi_beg;
              //std::cout << ",";
              //std::cout << lumi_end << "]" << std::endl;
              
	      std::pair<int, int> tempLS(atoi(lumi_beg.c_str()), atoi(lumi_end.c_str()));
              lumisections.push_back(tempLS);
              
              it2 = it_end;
            }
          
          jsonMap[atoi(run.c_str())] = lumisections;
        } // find run number
      
    } // loop on JSON file
  
  return jsonMap;
}

bool AcceptEventByRunAndLumiSection(const int& runId, const int& lumiId,
                                    std::map<int, std::vector<std::pair<int, int> > >& jsonMap)
{
  // select by runId
  if( jsonMap.find(runId) == jsonMap.end() ) return false;
  
  // select by lumiId
  std::vector<std::pair<int, int> > lumisections = jsonMap[runId];
  
  int skipEvent = true;
  for(unsigned int i = 0; i < lumisections.size(); ++i)
    if( (lumiId >= lumisections.at(i).first) &&
        (lumiId <= lumisections.at(i).second) )
      skipEvent = false;
  
  if( skipEvent == true ) return false;
  
  return true;
}

//compute DeltaR between two vectors
double deltaR(LV v1, LV v2) {

  double deta = v1.Eta() - v2.Eta();
  //double dphi = v1.Phi() - v2.Phi();
  double dphi = TVector2::Phi_mpi_pi(v1.Phi() - v2.Phi());
  return TMath::Sqrt( TMath::Power(deta,2) + TMath::Power(dphi,2) );

}

//Weights for N-Jet MC samples
float reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  
  //with new high stat sample
  if(nJets==0)      return 0.492871535;
  else if(nJets==1) return 0.100267473;
  else if(nJets==2) return 0.031238278;
  else if(nJets==3) return 0.019961315;
  else if(nJets>=4) return 0.018980202;
  else return 1 ;

}

float reweightHEPNUPDYJets(int hepNUP) {

  int nJets = hepNUP-5;

  if(nJets==0)      return 0.11906618;
  else if(nJets==1) return 0.022478688;
  else if(nJets==2) return 0.009092586;
  else if(nJets==3) return 0.005278795;
  else if(nJets>=4) return 0.004118808;
  else return 1 ;

}
/*
float pileupWeight( float intimepileup){
  float weight_ = LumiWeights_->weight(intimepileup); 
  return weight_;
}
*/
int getJetIDMVALoose(double pt, double eta, double rawMVA)
{
  float eta_bin[] = {0,2.5,2.75,3.0,5.0};
  float Pt010_Loose[]    = {-0.95,-0.96,-0.94,-0.95};
  float Pt1020_Loose[]   = {-0.95,-0.96,-0.94,-0.95};
  float Pt2030_Loose[]   = {-0.63,-0.60,-0.55,-0.45};
  float Pt3050_Loose[]   = {-0.63,-0.60,-0.55,-0.45};

  int passId = 0;
  for(int i = 0; i < 4; i++){
    if(TMath::Abs(eta) >= eta_bin[i] && TMath::Abs(eta) < eta_bin[i+1]){
      if(pt < 10){
        if(rawMVA > Pt010_Loose[i])passId = 1;
      }
      else if(pt >= 10 && pt < 20){
        if(rawMVA > Pt1020_Loose[i])passId = 1;
      }
      else if(pt >= 20 && pt < 30){ 
        if(rawMVA > Pt2030_Loose[i])passId = 1; 
      }
      else if(pt >= 30){  
        if(rawMVA > Pt3050_Loose[i])passId = 1;  
      }
    }
  }
  return passId;
}
//Remove duplication of events in ntuples
bool checkEventIsDuplicated(MAPDITAU_run &mapDiTau, int run, int lumi, int event, float ptL1, float ptL2, float etaL1, float etaL2)
{
  
  MAPDITAU_run::const_iterator iter_run = mapDiTau.find(run) ;
  
  if( iter_run != mapDiTau.end() ) {
    MAPDITAU_lumi::const_iterator iter_lumi = iter_run->second.find(lumi) ;
    
    if( iter_lumi != iter_run->second.end() ) {
      MAPDITAU_event::const_iterator iter_event = iter_lumi->second.find(event) ;
      
      if( iter_event != iter_lumi->second.end() ) {
	MAPDITAU_ptL1::const_iterator iter_ptL1 = iter_event->second.find(ptL1) ;
	
	if( iter_ptL1 != iter_event->second.end() ) {
	  MAPDITAU_ptL2::const_iterator iter_ptL2 = iter_ptL1->second.find(ptL2) ;
	    
	  if( iter_ptL2 != iter_ptL1->second.end() ) {
	    MAPDITAU_etaL1::const_iterator iter_etaL1 = iter_ptL2->second.find(etaL1) ;
	        
	    if( iter_etaL1 != iter_ptL2->second.end() ) {
	      MAPDITAU_etaL2::const_iterator iter_etaL2 = iter_etaL1->second.find(etaL2) ;

	      if( iter_etaL2 != iter_etaL1->second.end() ) {
		return true;
	      }
	    }
	  }
	}
      }
    }
  }
  mapDiTau[run][lumi][event][ptL1][ptL2][etaL1][etaL2] = 1;

  return false;
}
void splitstring(string input, vector<string>& output)
{
  UInt_t posstart = 0;
  for(UInt_t i = 0 ; i < input.size() ; i++)
    {
      if(input[i] == ' ')
	{
	  output.push_back(input.substr(posstart, i-posstart));
	  posstart=i+1;
	}
    }
}

string combinestring(const vector<string>& input)
{
  string output;
  for(UInt_t i = 0 ; i < input.size() ; i++)
    {
      output += input[i] + string(" ");
    }
  return(output);
}
Int_t GetTauDiscriminator(string disname, string disnamelist, ULong64_t dishps) 
{
  vector<string> taudiscriminators;
  splitstring(disnamelist, taudiscriminators);
  Int_t pos = -1;
  for(UInt_t i = 0 ; i < taudiscriminators.size() ; i++){
    if(disname == taudiscriminators[i]){
      pos = i;
      break;
    }
  }
  if(pos == -1) return(-1);
  if((dishps & ((ULong64_t)1<<(ULong64_t)pos)) != 0) return(1);
  return(0);
}
Float_t GetBTagDiscriminator(string disname, std::vector<string> disnamelist, float *bDiscrValues)
{
  Int_t pos = -1;
  for(UInt_t i = 0 ; i < disnamelist.size() ; i++){
    if(disname == disnamelist[i]){
      pos = i;
      break;
    }
  }
  if(pos == -1) return(-99.);
  else return bDiscrValues[pos];
}
/* //old code
Int_t GetTriggerMatch(string pathName_, string triglist, ULong64_t trigResult_)
{
  vector<string> trigPathNames;
  splitstring(triglist, trigPathNames);
  Int_t index = -1;
  for(UInt_t i = 0 ; i < trigPathNames.size() ; i++)
    {
      std::string indexedTrigger = trigPathNames[i];
      size_t pos = pathName_.find(":");
      string pathName_name = pathName_;
      string pathName_filter = "";
      if(pos != string::npos){
	pathName_name = pathName_.substr(0, pos);
	pathName_filter = pathName_.substr(pos+1);
      }
      if(indexedTrigger.find(pathName_name) != string::npos){
	if(pathName_filter.size() > 0){
	  if(indexedTrigger.find(pathName_filter) != string::npos)
	    index = i;
	}
	else index = i;
      }
    }

  if(index == -1) return(0);
  bool result = (trigResult_ & 1<<index) != 0;
  return int(result);
}
*/
bool IsHLTMatched(string filterName_, std::vector<string> hltfilters, LV trigCandP4_, bool trigobject_filters[50], LV LeptonP4_){
  bool trigMatch_ = false;
  int filter_index = -1;
  for(size_t i = 0; i < hltfilters.size(); i++){
    if(hltfilters[i].find(filterName_) != std::string::npos) {
      filter_index = i;
      break;
    }
  }
  //cout<<filterName_<<" index "<<filter_index<<endl;
  if(filter_index >= 0){
    if(trigobject_filters[filter_index]){
      if(deltaR(trigCandP4_, LeptonP4_) < 0.5) trigMatch_ = true;
    }
  }

  return trigMatch_;
}

LV TauEnergyCorrector(LV tauP4_, LV genTauP4_, int decaymode_)
{
  double shiftP = 1.;
  double shiftMass = 1.;
  if(genTauP4_.pt() > 8.0 && ROOT::Math::VectorUtil::DeltaR(tauP4_, genTauP4_) < 0.5){
    
    if(decaymode_ == 1){ //1prong 0pi0
      shiftP = 1.01; //New correction for Winter2013
      shiftMass = 1.0;
    }
    else if(decaymode_ == 2){ //1prong 1pi0
      shiftP = 1.01; //New correction for Winter2013
      shiftMass = 1.01;
    }
    else if(decaymode_ == 3){ //3prong
      shiftP = 1.01; //New correction for Winter2013
      shiftMass = 1.01;
    }
  }

  double pxS = tauP4_.px()*shiftP;
  double pyS = tauP4_.py()*shiftP;
  double pzS = tauP4_.pz()*shiftP;
  double massS = tauP4_.mass()*shiftMass;
  double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
  LV p4S_( pxS, pyS, pzS, enS );

  return p4S_;
}

LV TauEnergyRescaler(LV tauP4_, int decaymode_)
{
  double shift = 0.03;
  double scale = 1+shift;
  if(decaymode_ == 1){ //1prong 0pi0
    scale = sqrt( tauP4_.energy()*(1+shift)*tauP4_.energy()*(1+shift) - tauP4_.mass()*tauP4_.mass() )/tauP4_.P();
  }

  LV p4S_(tauP4_.px()*scale , tauP4_.py()*scale, tauP4_.pz()*scale, tauP4_.energy()*(1+shift) ); 

  return p4S_;
}

LV RescaledMET(LV pfMet_, std::vector<LV> leptons_, double shift_)
{
  float dPx = 0; float dPy=0;
  for(size_t i = 0; i < leptons_.size(); i++){
    dPx += leptons_[i].px() * shift_;
    dPy += leptons_[i].py() * shift_;
  }
  
  double scaledMETPx = pfMet_.px() - dPx;
  double scaledMETPy = pfMet_.py() - dPy;
  
  LV scaledMET(scaledMETPx, scaledMETPy, 0, sqrt(scaledMETPx*scaledMETPx + scaledMETPy*scaledMETPy));

  return scaledMET;
}
/*
bool GetTriggerResult(std::map<std::string, int> hltriggerresults_, std::string pathName_)
{
  cout<<" required path "<<pathName_<<endl;
  bool result_ = false;
  for(std::map<std::string, int>::const_iterator it = hltriggerresults_.begin(); it != hltriggerresults_.end(); it++){
    std::string thisPathName_ = it->first;
    cout<<" stored path "<<thisPathName_<<endl;
    if(thisPathName_.find(pathName_) != std::string::npos){
      if(it->second > 0.5)result_ = true;
    }
  }

  return result_;
}
*/
bool GetTriggerResult(std::vector<std::string> hltriggerresults_, std::string pathName_)
{
  //cout<<" required path "<<pathName_<<endl;
  bool result_ = false;
  for(std::vector<std::string>::const_iterator it = hltriggerresults_.begin(); it != hltriggerresults_.end(); it++){
    std::string thisPathName_ = (*it);
    //cout<<" stored path "<<thisPathName_<<endl;
    if(thisPathName_.find(pathName_) != std::string::npos){
      result_ = true;
    }
  }

  return result_;
}

void fillTrees_TauTauStream(TChain* currentTree,
			   TTree* outTree,
			   double nEventsRead =0.,
			   string analysis_ = "", 
			   string sample_ = "",
			   float xsec_ = 0., 
			   float skimEff_ = 0., 
			   int iJson_=-1,
			   int iDiv = 0,
			   int nDiv = 1
			   )
{


  // normalization Lumi
  Float_t Lumi=1000;
  
  //Event info
  UInt_t event_nr, event_run, event_luminosityblock;
  Float_t genweight;

  //PV
  unsigned int primvertex_count;
  float primvertex_x, primvertex_y, primvertex_z, primvertex_ndof;
  //PU
  Float_t numtruepileupinteractions, rhoNeutral;

  //Muons
  UInt_t muon_count; 
  float muon_px[40], muon_py[40], muon_pz[40];
  float muon_dxy[40], muon_dz[40];
  float //muon_chi2[40], muon_ndof[40], // muon_innertrack_dxy[40], muon_innertrack_dz[40], 
    muon_r03_sumChargedHadronPt[40], muon_r03_sumNeutralHadronEt[40],
    muon_r03_sumPhotonEt[40],  muon_r03_sumPUPt[40];
  
  //UInt_t muon_nMuonStations[40], muon_nMuonHits[40], 
  //muon_nTrackerHits[40], muon_nPixelHits[40];
  Float_t muon_charge[40];
  bool muon_isGlobal[40], muon_isTracker[40], muon_isTight[40],
    muon_isLoose[40], muon_isMedium[40];
  
  //Electron
  UInt_t electron_count;
  float electron_px[40], electron_py[40], electron_pz[40],
    electron_superclusterEta[40],
    //electron_trackchi2[40], electron_trackndof[40],
    electron_dxy[40], electron_dz[40], electron_charge[40],
    //electron_deltaetasuperclustertrack[40], electron_deltaphisuperclustertrack[40], 
    //electron_full5x5_sigmaietaieta[40], electron_ehcaloverecal[40],
    //electron_ooemoop[40],
    electron_r03_sumChargedHadronPt[40], electron_r03_sumNeutralHadronEt[40], 
    electron_r03_sumPhotonEt[40], electron_r03_sumPUPt[40];
    //electron_mva_id_nontrigPhys14[40];
    
  bool electron_pass_conversion[40], electron_mva_wp90_nontrig_Spring15_v1[40];
  UChar_t electron_nmissinginnerhits[40];
  
  //Taus
  UInt_t tau_count;
  float tau_e[40], tau_px[40], tau_py[40], tau_pz[40], tau_mass[40],
    tau_dxy[40], tau_dz[40], tau_vertexx[40], tau_vertexy[40], tau_vertexz[40]; 
  unsigned int tau_signalChargedHadrCands_size[40], tau_signalGammaCands_size[40];
  Float_t tau_charge[40]; int tau_decayMode[40];
  Float_t tau_leadchargedhadrcand_dz[40], tau_leadchargedhadrcand_dxy[40];
  //unsigned long int tau_dishps[40] ;
  //Char_t run_taudiscriminators[10000];
  bool tau_L1trigger_match[40];
  //string tau_genTaudecayMode[40];
  float tau_decayModeFinding[40], tau_decayModeFindingNewDMs[40], 
    tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[40], tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[40],
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[40], tau_byTightCombinedIsolationDeltaBetaCorr3Hits[40],
    tau_againstMuonLoose3[40], tau_againstMuonTight3[40],
    tau_againstElectronVLooseMVA5[40], tau_againstElectronLooseMVA5[40], tau_againstElectronMediumMVA5[40],
    tau_againstElectronTightMVA5[40];
  float tau_genjet_px[40], tau_genjet_py[40], tau_genjet_pz[40], tau_genjet_e[40];
  
  //Jets
  UInt_t pfjet_count;
  float pfjet_e[100], pfjet_px[100], pfjet_py[100], pfjet_pz[100],
    pfjet_neutralhadronicenergy[100], pfjet_chargedhadronicenergy[100],
    pfjet_neutralemenergy[100], pfjet_chargedemenergy[100],
    pfjet_btag[100][10], pfjet_energycorr[100];
  UInt_t pfjet_neutralmulti[100], pfjet_chargedmulti[100], pfjet_chargedhadronmulti[100];
  float pfjet_pu_jet_full_mva[100];
  int pfjet_flavour[100];
  std::vector<std::string> *run_btagdiscriminators = new std::vector<std::string> ();
  
  //MET
  float pfmet_ex, pfmet_ey, pfmet_sigxx, pfmet_sigxy, pfmet_sigyx, pfmet_sigyy;
  UInt_t mvamet_count;
  float mvamet_ex[200], mvamet_ey[200], mvamet_sigxx[200], mvamet_sigxy[200], mvamet_sigyx[200], mvamet_sigyy[200];
  UChar_t mvamet_channel[200];
  UInt_t mvamet_lep1[200], mvamet_lep2[200];

  //GetParticles
  unsigned int genparticles_count;
  float genparticles_e[1000], genparticles_px[1000], genparticles_py[1000],genparticles_pz[1000],
    genparticles_vx[1000], genparticles_vy[1000], genparticles_vz[1000]; 
  Int_t genparticles_pdgid[1000], genparticles_status[1000];
  UChar_t genparticles_mother[1000];
  Int_t genparticles_isPrompt[1000], genparticles_isDirectPromptTauDecayProduct[1000];

  UInt_t gentau_count;
  float gentau_e[20], gentau_px[20], gentau_py[20], gentau_pz[20],
    gentau_visible_e[20], gentau_visible_px[20], gentau_visible_py[20],
    gentau_visible_pz[20];
  int  gentau_decayMode[20], gentau_isPrompt[20]; 
  UChar_t gentau_mother[20];

  std::vector<std::string> *hltriggerresultsV = new std::vector<std::string> ();
  std::vector<std::string> *run_hltfilters = new std::vector<std::string>();
  UInt_t trigobject_count;
  float trigobject_px[100], trigobject_py[100], trigobject_pz[100];
  bool trigobject_filters[100][50];

  currentTree->SetBranchStatus("*"        ,0);
  currentTree->SetBranchStatus("event_nr"        ,1);
  currentTree->SetBranchStatus("event_run"        ,1);
  currentTree->SetBranchStatus("event_luminosityblock"        ,1);
  currentTree->SetBranchStatus("genweight"        ,1);
  currentTree->SetBranchStatus("primvertex_count"      ,1);
  currentTree->SetBranchStatus("primvertex_x"        ,1);
  currentTree->SetBranchStatus("primvertex_y"        ,1);
  currentTree->SetBranchStatus("primvertex_z"        ,1);
  currentTree->SetBranchStatus("primvertex_ndof"        ,1);
  currentTree->SetBranchStatus("numtruepileupinteractions",  1);
  currentTree->SetBranchStatus("rho",    1);
  currentTree->SetBranchStatus("muon_count"        ,1);
  currentTree->SetBranchStatus("muon_px"        ,1);
  currentTree->SetBranchStatus("muon_py"        ,1);
  currentTree->SetBranchStatus("muon_pz"        ,1);
  currentTree->SetBranchStatus("muon_dxy"       ,1);
  currentTree->SetBranchStatus("muon_dz"        ,1);
  //currentTree->SetBranchStatus("muon_chi2"        ,1);
  //currentTree->SetBranchStatus("muon_ndof"        ,1);
  //currentTree->SetBranchStatus("muon_innertrack_dxy"        ,1);
  //currentTree->SetBranchStatus("muon_innertrack_dz"        ,1);
  currentTree->SetBranchStatus("muon_r03_sumChargedHadronPt"        ,1);
  currentTree->SetBranchStatus("muon_r03_sumNeutralHadronEt"        ,1);
  currentTree->SetBranchStatus("muon_r03_sumPhotonEt"        ,1);
  currentTree->SetBranchStatus("muon_r03_sumPUPt"        ,1);
  currentTree->SetBranchStatus("muon_charge"        ,1);
  //currentTree->SetBranchStatus("muon_nMuonHits"        ,1);
  //currentTree->SetBranchStatus("muon_nPixelHits"        ,1);
  //currentTree->SetBranchStatus("muon_nTrackerHits"        ,1);
  //currentTree->SetBranchStatus("muon_nMuonStations"        ,1);
  currentTree->SetBranchStatus("muon_isGlobal", 1);
  currentTree->SetBranchStatus("muon_isTracker", 1);
  currentTree->SetBranchStatus("muon_isTight", 1);
  currentTree->SetBranchStatus("muon_isLoose", 1);
  currentTree->SetBranchStatus("muon_isMedium", 1);
  currentTree->SetBranchStatus("electron_count"        ,1);
  currentTree->SetBranchStatus("electron_px"        ,1);
  currentTree->SetBranchStatus("electron_py"        ,1);
  currentTree->SetBranchStatus("electron_pz"        ,1);
  currentTree->SetBranchStatus("electron_superclusterEta"   ,1);
  //currentTree->SetBranchStatus("electron_trackchi2"        ,1);
  //currentTree->SetBranchStatus("electron_trackndof"        ,1);
  currentTree->SetBranchStatus("electron_dxy"        ,1);
  currentTree->SetBranchStatus("electron_dz"        ,1);
  currentTree->SetBranchStatus("electron_charge"        ,1);
  //currentTree->SetBranchStatus("electron_deltaetasuperclustertrack", 1);
  //currentTree->SetBranchStatus("electron_deltaphisuperclustertrack", 1);
  //currentTree->SetBranchStatus("electron_full5x5_sigmaietaieta", 1);
  //currentTree->SetBranchStatus("electron_ehcaloverecal", 1);
  //currentTree->SetBranchStatus("electron_ooemoop", 1);
  //currentTree->SetBranchStatus("electron_r03_sumChargedHadronPt", 1);
  currentTree->SetBranchStatus("electron_r03_sumNeutralHadronEt", 1);
  currentTree->SetBranchStatus("electron_r03_sumPhotonEt", 1);
  currentTree->SetBranchStatus("electron_r03_sumPUPt", 1);
  //currentTree->SetBranchStatus("electron_mva_id_nontrigPhys14", 1);
  currentTree->SetBranchStatus("electron_mva_wp90_nontrig_Spring15_v1", 1);
  currentTree->SetBranchStatus("electron_pass_conversion", 1);
  currentTree->SetBranchStatus("electron_nmissinginnerhits", 1);
  currentTree->SetBranchStatus("tau_count"        ,1);
  currentTree->SetBranchStatus("tau_e"        ,1);
  currentTree->SetBranchStatus("tau_px"        ,1);
  currentTree->SetBranchStatus("tau_py"        ,1);
  currentTree->SetBranchStatus("tau_pz"        ,1);
  currentTree->SetBranchStatus("tau_mass"       ,1);
  currentTree->SetBranchStatus("tau_dxy"        ,1);
  currentTree->SetBranchStatus("tau_dz"        ,1);
  currentTree->SetBranchStatus("tau_vertexx"        ,1);
  currentTree->SetBranchStatus("tau_vertexy"        ,1);
  currentTree->SetBranchStatus("tau_vertexz"        ,1);
  currentTree->SetBranchStatus("tau_leadchargedhadrcand_dz"        ,1);
  currentTree->SetBranchStatus("tau_leadchargedhadrcand_dxy"        ,1);
  currentTree->SetBranchStatus("tau_signalChargedHadrCands_size"        ,1);
  currentTree->SetBranchStatus("tau_signalGammaCands_size"        ,1);
  currentTree->SetBranchStatus("tau_charge"        ,1);
  currentTree->SetBranchStatus("tau_decayMode"        ,1);
  currentTree->SetBranchStatus("tau_L1trigger_match"        ,1);
  currentTree->SetBranchStatus("tau_decayModeFinding",    1);
  currentTree->SetBranchStatus("tau_decayModeFindingNewDMs",      1);
  currentTree->SetBranchStatus("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits"      ,1);
  currentTree->SetBranchStatus("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits"      ,1);
  currentTree->SetBranchStatus("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",   1);
  currentTree->SetBranchStatus("tau_byTightCombinedIsolationDeltaBetaCorr3Hits",   1);
  currentTree->SetBranchStatus("tau_againstMuonLoose3",   1);
  currentTree->SetBranchStatus("tau_againstMuonTight3",   1);
  currentTree->SetBranchStatus("tau_againstElectronVLooseMVA5",   1);
  currentTree->SetBranchStatus("tau_againstElectronLooseMVA5",   1);
  currentTree->SetBranchStatus("tau_againstElectronMediumMVA5",   1);
  currentTree->SetBranchStatus("tau_againstElectronTightMVA5",   1);
  currentTree->SetBranchStatus("tau_genjet_e"        ,1);
  currentTree->SetBranchStatus("tau_genjet_px"        ,1);
  currentTree->SetBranchStatus("tau_genjet_py"        ,1);
  currentTree->SetBranchStatus("tau_genjet_pz"        ,1);							    
  currentTree->SetBranchStatus("pfjet_count"        ,1);
  currentTree->SetBranchStatus("pfjet_e"        ,1);
  currentTree->SetBranchStatus("pfjet_px"        ,1);
  currentTree->SetBranchStatus("pfjet_py"        ,1);
  currentTree->SetBranchStatus("pfjet_pz"        ,1);
  currentTree->SetBranchStatus("pfjet_neutralhadronicenergy"        ,1);
  currentTree->SetBranchStatus("pfjet_chargedhadronicenergy"        ,1);
  currentTree->SetBranchStatus("pfjet_neutralemenergy"        ,1);
  currentTree->SetBranchStatus("pfjet_chargedemenergy"        ,1);
  currentTree->SetBranchStatus("pfjet_chargedmulti"        ,1);
  currentTree->SetBranchStatus("pfjet_neutralmulti"        ,1);
  currentTree->SetBranchStatus("pfjet_chargedhadronmulti", 1);
  currentTree->SetBranchStatus("pfjet_btag"        ,1);
  currentTree->SetBranchStatus("pfjet_energycorr"     ,1);
  currentTree->SetBranchStatus("pfjet_flavour"        ,1);
  currentTree->SetBranchStatus("pfjet_pu_jet_full_mva"        ,1);
  currentTree->SetBranchStatus("run_btagdiscriminators"        ,1);
  currentTree->SetBranchStatus("pfmet_ex",            1);
  currentTree->SetBranchStatus("pfmet_ey",            1);
  currentTree->SetBranchStatus("pfmet_sigxx",        1);
  currentTree->SetBranchStatus("pfmet_sigxy",        1);
  currentTree->SetBranchStatus("pfmet_sigyx",        1);
  currentTree->SetBranchStatus("pfmet_sigyy",        1);
  currentTree->SetBranchStatus("mvamet_count",       1);
  currentTree->SetBranchStatus("mvamet_ex",        1);
  currentTree->SetBranchStatus("mvamet_ey",        1);
  currentTree->SetBranchStatus("mvamet_sigxx",        1);
  currentTree->SetBranchStatus("mvamet_sigxy",        1);
  currentTree->SetBranchStatus("mvamet_sigyx",        1);
  currentTree->SetBranchStatus("mvamet_sigyy",        1);
  currentTree->SetBranchStatus("mvamet_channel",        1);
  currentTree->SetBranchStatus("mvamet_lep1",        1);
  currentTree->SetBranchStatus("mvamet_lep2",        1);
  currentTree->SetBranchStatus("genparticles_count"      ,1);
  currentTree->SetBranchStatus("genparticles_e"      ,1);
  currentTree->SetBranchStatus("genparticles_px"      ,1);
  currentTree->SetBranchStatus("genparticles_py"      ,1);
  currentTree->SetBranchStatus("genparticles_pz"      ,1);
  currentTree->SetBranchStatus("genparticles_vx"      ,1);
  currentTree->SetBranchStatus("genparticles_vy"      ,1);
  currentTree->SetBranchStatus("genparticles_vz"      ,1);
  currentTree->SetBranchStatus("genparticles_pdgid"      ,1);
  currentTree->SetBranchStatus("genparticles_status"      ,1);
  currentTree->SetBranchStatus("genparticles_mother"      ,1);
  currentTree->SetBranchStatus("genparticles_isPrompt"   ,1);
  currentTree->SetBranchStatus("genparticles_isDirectPromptTauDecayProduct" ,1);
  currentTree->SetBranchStatus("gentau_count", 1);
  currentTree->SetBranchStatus("gentau_e", 1);
  currentTree->SetBranchStatus("gentau_px", 1);
  currentTree->SetBranchStatus("gentau_py", 1);
  currentTree->SetBranchStatus("gentau_pz", 1);
  currentTree->SetBranchStatus("gentau_visible_e", 1);
  currentTree->SetBranchStatus("gentau_visible_px", 1);
  currentTree->SetBranchStatus("gentau_visible_py", 1);
  currentTree->SetBranchStatus("gentau_visible_pz", 1);
  currentTree->SetBranchStatus("gentau_decayMode", 1);
  currentTree->SetBranchStatus("gentau_mother", 1);
  currentTree->SetBranchStatus("gentau_isPrompt", 1);
  currentTree->SetBranchStatus("run_hltfilters"      ,1);
  currentTree->SetBranchStatus("trigobject_count",  1);
  currentTree->SetBranchStatus("trigobject_px",  1);
  currentTree->SetBranchStatus("trigobject_py",    1);
  currentTree->SetBranchStatus("trigobject_pz",    1);
  currentTree->SetBranchStatus("trigobject_filters",    1);
  currentTree->SetBranchStatus("hltriggerresultsV"      ,1);
  
  //Set branch address
  currentTree->SetBranchAddress("event_nr"        ,&event_nr);
  currentTree->SetBranchAddress("event_run"        ,&event_run);
  currentTree->SetBranchAddress("event_luminosityblock"        ,&event_luminosityblock);
  currentTree->SetBranchAddress("genweight"        ,&genweight);
  currentTree->SetBranchAddress("primvertex_count"      ,&primvertex_count);
  currentTree->SetBranchAddress("primvertex_x"        ,&primvertex_x);
  currentTree->SetBranchAddress("primvertex_y"        ,&primvertex_y);
  currentTree->SetBranchAddress("primvertex_z"        ,&primvertex_z);
  currentTree->SetBranchAddress("primvertex_ndof"        ,&primvertex_ndof);
  currentTree->SetBranchAddress("numtruepileupinteractions",  &numtruepileupinteractions);
  currentTree->SetBranchAddress("rho",    &rhoNeutral);
  currentTree->SetBranchAddress("muon_count"        ,&muon_count);
  currentTree->SetBranchAddress("muon_px"        ,muon_px);
  currentTree->SetBranchAddress("muon_py"        ,muon_py);
  currentTree->SetBranchAddress("muon_pz"        ,muon_pz);
  currentTree->SetBranchAddress("muon_dxy"       ,muon_dxy);
  currentTree->SetBranchAddress("muon_dz"        ,muon_dz);
  //currentTree->SetBranchAddress("muon_chi2"        ,muon_chi2);
  //currentTree->SetBranchAddress("muon_ndof"        ,muon_ndof);
  //currentTree->SetBranchAddress("muon_innertrack_dxy"        ,muon_innertrack_dxy);
  //currentTree->SetBranchAddress("muon_innertrack_dz"        ,muon_innertrack_dz);
  currentTree->SetBranchAddress("muon_r03_sumChargedHadronPt"        ,muon_r03_sumChargedHadronPt);
  currentTree->SetBranchAddress("muon_r03_sumNeutralHadronEt"        ,muon_r03_sumNeutralHadronEt);
  currentTree->SetBranchAddress("muon_r03_sumPhotonEt"        ,muon_r03_sumPhotonEt);
  currentTree->SetBranchAddress("muon_r03_sumPUPt"        ,muon_r03_sumPUPt);
  currentTree->SetBranchAddress("muon_charge"        ,muon_charge);
  //currentTree->SetBranchAddress("muon_nTrackerHits"        ,muon_nTrackerHits);
  //currentTree->SetBranchAddress("muon_nMuonHits"        ,muon_nMuonHits);
  //currentTree->SetBranchAddress("muon_nPixelHits"        ,muon_nPixelHits);
  //currentTree->SetBranchAddress("muon_nMuonStations"        ,muon_nMuonStations);
  currentTree->SetBranchAddress("muon_isGlobal", muon_isGlobal);
  currentTree->SetBranchAddress("muon_isTracker", muon_isTracker);
  currentTree->SetBranchAddress("muon_isTight", muon_isTight);
  currentTree->SetBranchAddress("muon_isLoose", muon_isLoose);
  currentTree->SetBranchAddress("muon_isMedium", muon_isMedium);
  currentTree->SetBranchAddress("electron_count"        ,&electron_count);
  currentTree->SetBranchAddress("electron_px"        ,electron_px);
  currentTree->SetBranchAddress("electron_py"        ,electron_py);
  currentTree->SetBranchAddress("electron_pz"        ,electron_pz);
  currentTree->SetBranchAddress("electron_superclusterEta"   ,electron_superclusterEta);
  //currentTree->SetBranchAddress("electron_trackchi2"        ,electron_trackchi2);
  //currentTree->SetBranchAddress("electron_trackndof"        ,electron_trackndof);
  currentTree->SetBranchAddress("electron_dxy"        ,electron_dxy);
  currentTree->SetBranchAddress("electron_dz"        ,electron_dz);
  currentTree->SetBranchAddress("electron_charge"        ,electron_charge);
  //currentTree->SetBranchAddress("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack);
  //currentTree->SetBranchAddress("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack);
  //currentTree->SetBranchAddress("electron_full5x5_sigmaietaieta", electron_full5x5_sigmaietaieta);
  //currentTree->SetBranchAddress("electron_ehcaloverecal", electron_ehcaloverecal);
  //currentTree->SetBranchAddress("electron_ooemoop", electron_ooemoop);
  currentTree->SetBranchAddress("electron_r03_sumChargedHadronPt", electron_r03_sumChargedHadronPt);
  currentTree->SetBranchAddress("electron_r03_sumNeutralHadronEt", electron_r03_sumNeutralHadronEt);
  currentTree->SetBranchAddress("electron_r03_sumPhotonEt", electron_r03_sumPhotonEt);
  currentTree->SetBranchAddress("electron_r03_sumPUPt", electron_r03_sumPUPt);
  //currentTree->SetBranchAddress("electron_mva_id_nontrigPhys14", electron_mva_id_nontrigPhys14);
  currentTree->SetBranchAddress("electron_mva_wp90_nontrig_Spring15_v1", electron_mva_wp90_nontrig_Spring15_v1);
  currentTree->SetBranchAddress("electron_pass_conversion", electron_pass_conversion);
  currentTree->SetBranchAddress("electron_nmissinginnerhits", electron_nmissinginnerhits);
  currentTree->SetBranchAddress("tau_count"        ,&tau_count);
  currentTree->SetBranchAddress("tau_e"        ,tau_e);
  currentTree->SetBranchAddress("tau_px"        ,tau_px);
  currentTree->SetBranchAddress("tau_py"        ,tau_py);
  currentTree->SetBranchAddress("tau_pz"        ,tau_pz);
  currentTree->SetBranchAddress("tau_mass"      ,tau_mass);
  currentTree->SetBranchAddress("tau_dxy"        ,tau_dxy);
  currentTree->SetBranchAddress("tau_dz"        ,tau_dz);
  currentTree->SetBranchAddress("tau_vertexx"        ,tau_vertexx);
  currentTree->SetBranchAddress("tau_vertexy"        ,tau_vertexy);
  currentTree->SetBranchAddress("tau_vertexz"        ,tau_vertexz);
  currentTree->SetBranchAddress("tau_leadchargedhadrcand_dz"        ,tau_leadchargedhadrcand_dz);
  currentTree->SetBranchAddress("tau_leadchargedhadrcand_dxy"        ,tau_leadchargedhadrcand_dxy);
  currentTree->SetBranchAddress("tau_signalChargedHadrCands_size"        ,tau_signalChargedHadrCands_size);
  currentTree->SetBranchAddress("tau_signalGammaCands_size"        ,tau_signalGammaCands_size);
  currentTree->SetBranchAddress("tau_charge"        ,tau_charge);
  currentTree->SetBranchAddress("tau_decayMode"        ,tau_decayMode);
  currentTree->SetBranchAddress("tau_L1trigger_match"        ,tau_L1trigger_match);
  currentTree->SetBranchAddress("tau_decayModeFinding",    tau_decayModeFinding);
  currentTree->SetBranchAddress("tau_decayModeFindingNewDMs",      tau_decayModeFindingNewDMs);
  currentTree->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits"      ,tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  currentTree->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits"      ,tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
  currentTree->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",   tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
  currentTree->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits",   tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
  currentTree->SetBranchAddress("tau_againstMuonLoose3",   tau_againstMuonLoose3);
  currentTree->SetBranchAddress("tau_againstMuonTight3",   tau_againstMuonTight3);
  currentTree->SetBranchAddress("tau_againstElectronVLooseMVA5",   tau_againstElectronVLooseMVA5);
  currentTree->SetBranchAddress("tau_againstElectronLooseMVA5",   tau_againstElectronLooseMVA5);
  currentTree->SetBranchAddress("tau_againstElectronMediumMVA5",   tau_againstElectronMediumMVA5);
  currentTree->SetBranchAddress("tau_againstElectronTightMVA5",   tau_againstElectronTightMVA5);
  currentTree->SetBranchAddress("tau_genjet_e"        ,tau_genjet_e);
  currentTree->SetBranchAddress("tau_genjet_px"        ,tau_genjet_px);
  currentTree->SetBranchAddress("tau_genjet_py"        ,tau_genjet_py);
  currentTree->SetBranchAddress("tau_genjet_pz"        ,tau_genjet_pz);
  currentTree->SetBranchAddress("pfjet_count"        ,&pfjet_count);
  currentTree->SetBranchAddress("pfjet_e"        ,pfjet_e);
  currentTree->SetBranchAddress("pfjet_px"        ,pfjet_px);
  currentTree->SetBranchAddress("pfjet_py"        ,pfjet_py);
  currentTree->SetBranchAddress("pfjet_pz"        ,pfjet_pz);
  currentTree->SetBranchAddress("pfjet_neutralhadronicenergy"        ,pfjet_neutralhadronicenergy);
  currentTree->SetBranchAddress("pfjet_neutralemenergy"        ,pfjet_neutralemenergy);
  currentTree->SetBranchAddress("pfjet_chargedhadronmulti", pfjet_chargedhadronmulti);
  currentTree->SetBranchAddress("pfjet_chargedhadronicenergy"        ,pfjet_chargedhadronicenergy);
  currentTree->SetBranchAddress("pfjet_chargedemenergy"        ,pfjet_chargedemenergy);
  currentTree->SetBranchAddress("pfjet_chargedmulti"        ,pfjet_chargedmulti);
  currentTree->SetBranchAddress("pfjet_neutralmulti"        ,pfjet_neutralmulti);
  currentTree->SetBranchAddress("pfjet_btag"        ,pfjet_btag);
  currentTree->SetBranchAddress("pfjet_energycorr"     ,pfjet_energycorr);
  currentTree->SetBranchAddress("pfjet_pu_jet_full_mva"        ,pfjet_pu_jet_full_mva);
  currentTree->SetBranchAddress("pfjet_flavour"        ,pfjet_flavour);
  currentTree->SetBranchAddress("run_btagdiscriminators"        ,&run_btagdiscriminators);
  currentTree->SetBranchAddress("pfmet_ex",          &pfmet_ex);
  currentTree->SetBranchAddress("pfmet_ey",          &pfmet_ey);
  currentTree->SetBranchAddress("pfmet_sigxx",        &pfmet_sigxx);
  currentTree->SetBranchAddress("pfmet_sigxy",        &pfmet_sigxy);
  currentTree->SetBranchAddress("pfmet_sigyx",        &pfmet_sigyx);
  currentTree->SetBranchAddress("pfmet_sigyy",        &pfmet_sigyy);
  currentTree->SetBranchAddress("mvamet_count",     &mvamet_count);
  currentTree->SetBranchAddress("mvamet_ex",      mvamet_ex);
  currentTree->SetBranchAddress("mvamet_ey",      mvamet_ey);
  currentTree->SetBranchAddress("mvamet_sigxx",        mvamet_sigxx);
  currentTree->SetBranchAddress("mvamet_sigxy",        mvamet_sigxy);
  currentTree->SetBranchAddress("mvamet_sigyx",        mvamet_sigyx);
  currentTree->SetBranchAddress("mvamet_sigyy",        mvamet_sigyy);
  currentTree->SetBranchAddress("mvamet_channel",        mvamet_channel);
  currentTree->SetBranchAddress("mvamet_lep1",        mvamet_lep1);
  currentTree->SetBranchAddress("mvamet_lep2",        mvamet_lep2);
  currentTree->SetBranchAddress("genparticles_count"      ,&genparticles_count);
  currentTree->SetBranchAddress("genparticles_e"      ,genparticles_e);
  currentTree->SetBranchAddress("genparticles_px"      ,genparticles_px);
  currentTree->SetBranchAddress("genparticles_py"      ,genparticles_py);
  currentTree->SetBranchAddress("genparticles_pz"      ,genparticles_pz);
  currentTree->SetBranchAddress("genparticles_vx"      ,genparticles_vx);
  currentTree->SetBranchAddress("genparticles_vy"      ,genparticles_vy);
  currentTree->SetBranchAddress("genparticles_vz"      ,genparticles_vz);
  currentTree->SetBranchAddress("genparticles_pdgid"      ,genparticles_pdgid);
  currentTree->SetBranchAddress("genparticles_status"      ,genparticles_status);
  currentTree->SetBranchAddress("genparticles_mother"      ,genparticles_mother);
  currentTree->SetBranchAddress("genparticles_isPrompt"    ,genparticles_isPrompt);
  currentTree->SetBranchAddress("genparticles_isDirectPromptTauDecayProduct" ,genparticles_isDirectPromptTauDecayProduct);
  currentTree->SetBranchAddress("gentau_count", &gentau_count);
  currentTree->SetBranchAddress("gentau_e", gentau_e);
  currentTree->SetBranchAddress("gentau_px", gentau_px);
  currentTree->SetBranchAddress("gentau_py", gentau_py);
  currentTree->SetBranchAddress("gentau_pz", gentau_pz);
  currentTree->SetBranchAddress("gentau_visible_e", gentau_visible_e);
  currentTree->SetBranchAddress("gentau_visible_px", gentau_visible_px);
  currentTree->SetBranchAddress("gentau_visible_py", gentau_visible_py);
  currentTree->SetBranchAddress("gentau_visible_pz", gentau_visible_pz);
  currentTree->SetBranchAddress("gentau_decayMode", gentau_decayMode);
  currentTree->SetBranchAddress("gentau_mother", gentau_mother);
  currentTree->SetBranchAddress("gentau_isPrompt", gentau_isPrompt);
  currentTree->SetBranchAddress("run_hltfilters"      ,&run_hltfilters);
  currentTree->SetBranchAddress("trigobject_count",  &trigobject_count);
  currentTree->SetBranchAddress("trigobject_px",  trigobject_px);
  currentTree->SetBranchAddress("trigobject_py",    trigobject_py);
  currentTree->SetBranchAddress("trigobject_pz",    trigobject_pz);
  currentTree->SetBranchAddress("trigobject_filters",    trigobject_filters);
  currentTree->SetBranchAddress("hltriggerresultsV"      ,&hltriggerresultsV);

  //OutTree
  // kinematical variables of first 2 jets  
  float ptj1,ptj2,etaj1,etaj2,Detajj,Mjj,Dphijj,phij1,phij2;
  float pumvaj1, pumvaj2, csvj1, csvj2, ptrawj1,ptrawj2;
  float ptB1, etaB1, phiB1;
  float diJetPt, diJetPhi;
  int nJets30, nJets20, nJets20BTagged, nJets20BTaggedBUp, nJets20BTaggedBDown, nJets20BTaggedLUp, nJets20BTaggedLDown, nJets20BTaggedLoose;
  // kinematical variables of the veto jet
  float ptVeto, etaVeto, phiVeto;
  int nVetoJets;

  // diTau related variables
  float diTauNSVfitMass_,diTauNSVfitMassErrUp_,diTauNSVfitMassErrDown_;
  //float diTauNSVfitMassZLUp_,diTauNSVfitMassZLDown_;
  float diTauNSVfitPt_,diTauNSVfitPtErrUp_,diTauNSVfitPtErrDown_;
  float diTauNSVfitEta_, diTauNSVfitPhi_;
  float diTauVisMass,diTauVisPt,diTauVisEta,diTauVisPhi;
  float diTauRecoPt,diTauRecoPhi;
  //float diTauVisMassZLUp, diTauVisMassZLDown;

  // taus/MET related variables
  float ptL1,ptL2,etaL1,etaL2,phiL1,phiL2,dRL1L2,dEtaL1L2,dPhiL1L2,dPhiL1J1,dPhiL1J2,dPhiL2J1,dPhiL2J2;
  float dzL1, dzL2, dxyL1, dxyL2; 
  float diTauCharge_, chargeL1_;
  float MEt, MEtPhi, MtLeg1_, MtLeg2_;
  float MEtMVA, MEtMVAPhi, MtLeg1MVA_, MtLeg2MVA_;
  float MEtCov00,MEtCov01,MEtCov10,MEtCov11;
  float MEtMVACov00,MEtMVACov01,MEtMVACov10,MEtMVACov11;

  //TauID
  int decayModeL1_,decayModeFindingL1_,decayModeFindingNewDML1_,decayModeFindingOldDML1_;
  int AntiEDeadEcalL1_,tightestAntiECutWPL1_,tightestAntiEMVA5WPL1_,AntiEMVA5categoryL1_;
  int tightestAntiMuWPL1_,tightestAntiMu2WPL1_,tightestAntiMu3WPL1_,tightestAntiMuMVAWPL1_;
  float AntiEMVA5rawL1_,AntiMuMVArawL1_;
  int tightestHPSDBWPL1_,tightestHPSDB3HWPL1_,tightestHPSMVA3newDMwLTWPL1_,tightestHPSMVA3newDMwoLTWPL1_,tightestHPSMVA3oldDMwLTWPL1_,tightestHPSMVA3oldDMwoLTWPL1_;
  float hpsDB3HL1_,hpsMVA3newDMwLTL1_,hpsMVA3newDMwoLTL1_,hpsMVA3oldDMwLTL1_,hpsMVA3oldDMwoLTL1_;

  int decayModeL2_,decayModeFindingL2_,decayModeFindingNewDML2_,decayModeFindingOldDML2_;
  int AntiEDeadEcalL2_,tightestAntiECutWPL2_,tightestAntiEMVA5WPL2_,AntiEMVA5categoryL2_;
  int tightestAntiMuWPL2_,tightestAntiMu2WPL2_,tightestAntiMu3WPL2_,tightestAntiMuMVAWPL2_;
  float AntiEMVA5rawL2_,AntiMuMVArawL2_;
  int tightestHPSDBWPL2_,tightestHPSDB3HWPL2_,tightestHPSMVA3newDMwLTWPL2_,tightestHPSMVA3newDMwoLTWPL2_,tightestHPSMVA3oldDMwLTWPL2_,tightestHPSMVA3oldDMwoLTWPL2_;
  float hpsDB3HL2_,hpsMVA3newDMwLTL2_,hpsMVA3newDMwoLTL2_,hpsMVA3oldDMwLTL2_,hpsMVA3oldDMwoLTL2_;

  //tau related variables
  float visibleTauMassL1, visibleTauMassL2;
  float genVMass, genVPt;
  float visGenTauMassL1, genTauPtL1, genTauEtaL1;
  int genDecayModeL1_;
  float visGenTauMassL2, genTauPtL2, genTauEtaL2;
  int genDecayModeL2_;
  int genMatchL1_, genMatchL2_;
  
  // event-related variables
  int nVetoMuon_, nVetoElectron_;
  float numPV_, npu_, rho_;

  // event id
  ULong64_t event_,run_,lumi_;
  int index_, pairIndex; //, counter;

  //trigger
  int HLTx, HLTmatchL1, HLTmatchL2;
  
  //event weights
  float evtweight, mcweight, puweight; 
  //clasify DY events for splitting
  bool isZtt_, isZl_, isZj_;

  outTree->Branch("ptj1",  &ptj1,"ptj1/F");
  outTree->Branch("ptj2",  &ptj2,"ptj2/F");
  outTree->Branch("etaj1", &etaj1,"etaj1/F");
  outTree->Branch("etaj2", &etaj2,"etaj2/F");
  outTree->Branch("phij1", &phij1,"phij1/F");
  outTree->Branch("phij2", &phij2,"phij2/F");
  outTree->Branch("pumvaj1",  &pumvaj1,"pumvaj1/F");
  outTree->Branch("pumvaj2",  &pumvaj2,"pumvaj2/F");
  outTree->Branch("csvj1",  &csvj1,"csvj1/F");
  outTree->Branch("csvj2",  &csvj2,"csvj2/F");
  outTree->Branch("ptrawj1",  &ptrawj1,"ptrawj1/F");
  outTree->Branch("ptrawj2",  &ptrawj2,"ptrawj2/F");
  outTree->Branch("Detajj", &Detajj,"Detajj/F");
  outTree->Branch("Dphijj", &Dphijj,"Dphijj/F");
  outTree->Branch("Mjj",  &Mjj,"Mjj/F");
  outTree->Branch("diJetPt",  &diJetPt , "diJetPt/F");
  outTree->Branch("diJetPhi", &diJetPhi , "diJetPhi/F");
  outTree->Branch("ptB1",  &ptB1, "ptB1/F");
  outTree->Branch("etaB1", &etaB1,"etaB1/F");
  outTree->Branch("phiB1", &phiB1,"phiB1/F");
  outTree->Branch("nJets30",       &nJets30,  "nJets30/I");  
  outTree->Branch("nJets20",       &nJets20,  "nJets20/I");
  outTree->Branch("nJets20BTagged",&nJets20BTagged,  "nJets20BTagged/I");  
  outTree->Branch("nJets20BTaggedLoose",&nJets20BTaggedLoose,  "nJets20BTaggedLoose/I");  
  outTree->Branch("nJets20BTaggedBUp",&nJets20BTaggedBUp,  "nJets20BTaggedBUp/I");
  outTree->Branch("nJets20BTaggedBDown",&nJets20BTaggedBDown,  "nJets20BTaggedBDown/I");
  outTree->Branch("nJets20BTaggedLUp",&nJets20BTaggedLUp,  "nJets20BTaggedLUp/I");
  outTree->Branch("nJets20BTaggedLDown",&nJets20BTaggedLDown,  "nJets20BTaggedLDown/I");
  outTree->Branch("ptVeto",  &ptVeto, "ptVeto/F");
  outTree->Branch("phiVeto", &phiVeto,"phiVeto/F");
  outTree->Branch("etaVeto", &etaVeto,"etaVeto/F");
  outTree->Branch("nVetoJets", &nVetoJets,"nVetoJets/I");
  
  outTree->Branch("diTauNSVfitMass",       &diTauNSVfitMass_,       "diTauNSVfitMass/F");
  outTree->Branch("diTauNSVfitMassErrUp",  &diTauNSVfitMassErrUp_,  "diTauNSVfitMassErrUp/F");
  outTree->Branch("diTauNSVfitMassErrDown",&diTauNSVfitMassErrDown_,"diTauNSVfitMassErrDown/F");
  outTree->Branch("diTauNSVfitPt",       &diTauNSVfitPt_,       "diTauNSVfitPt/F");
  outTree->Branch("diTauNSVfitPtErrUp",  &diTauNSVfitPtErrUp_,  "diTauNSVfitPtErrUp/F");
  outTree->Branch("diTauNSVfitPtErrDown",&diTauNSVfitPtErrDown_,"diTauNSVfitPtErrDown/F");
  outTree->Branch("diTauNSVfitEta",       &diTauNSVfitEta_,       "diTauNSVfitEta/F");
  outTree->Branch("diTauNSVfitPhi",       &diTauNSVfitPhi_,       "diTauNSVfitPhi/F");
  outTree->Branch("diTauRecoPt",   &diTauRecoPt,  "diTauRecoPt/F");
  outTree->Branch("diTauRecoPhi",  &diTauRecoPhi,  "diTauRecoPhi/F");
  outTree->Branch("diTauVisMass",&diTauVisMass,"diTauVisMass/F");
  outTree->Branch("diTauVisPt",  &diTauVisPt,"diTauVisPt/F");
  outTree->Branch("diTauVisEta", &diTauVisEta,"diTauVisEta/F");
  outTree->Branch("diTauVisPhi", &diTauVisPhi,"diTauVisPhi/F");

  outTree->Branch("etaL1",   &etaL1,"etaL1/F");
  outTree->Branch("etaL2",   &etaL2,"etaL2/F");
  outTree->Branch("ptL1",    &ptL1,"ptL1/F");
  outTree->Branch("ptL2",    &ptL2,"ptL2/F");
  outTree->Branch("phiL1",   &phiL1,"phiL1/F");
  outTree->Branch("phiL2",   &phiL2,"phiL2/F");
  outTree->Branch("dRL1L2"  ,&dRL1L2,"dRL1L2/F");
  outTree->Branch("dEtaL1L2",&dEtaL1L2,"dEtaL1L2/F");
  outTree->Branch("dPhiL1L2",&dPhiL1L2,"dPhiL1L2/F");
  outTree->Branch("dPhiL1J1",&dPhiL1J1,"dPhiL1J1/F");
  outTree->Branch("dPhiL1J2",&dPhiL1J2,"dPhiL1J2/F");
  outTree->Branch("dPhiL2J1",&dPhiL2J1,"dPhiL2J1/F");
  outTree->Branch("dPhiL2J2",&dPhiL2J2,"dPhiL2J2/F");
  outTree->Branch("dzL1", &dzL1, "dzL1/F");
  outTree->Branch("dzL2", &dzL2, "dzL2/F");
  outTree->Branch("dxyL1", &dxyL1, "dxyL1/F");
  outTree->Branch("dxyL2", &dxyL2, "dxyL2/F");
  outTree->Branch("visibleTauMassL1",          &visibleTauMassL1,"visibleTauMassL1/F");
  outTree->Branch("visGenTauMassL1",           &visGenTauMassL1, "visGenTauMassL1/F");
  outTree->Branch("genTauPtL1",                &genTauPtL1, "genTauPtL1/F");
  outTree->Branch("genTauEtaL1",               &genTauEtaL1, "genTauEtaL1/F");
  outTree->Branch("genDecayModeL1",            &genDecayModeL1_, "genDecayModeL1/I");
  outTree->Branch("genVMass",                  &genVMass,     "genVMass/F");
  outTree->Branch("genVPt",                    &genVPt,     "genVPt/F");
  outTree->Branch("visibleTauMassL2",          &visibleTauMassL2,"visibleTauMassL2/F");
  outTree->Branch("visGenTauMassL2",           &visGenTauMassL2, "visGenTauMassL2/F");
  outTree->Branch("genTauPtL2",                &genTauPtL2, "genTauPtL2/F");
  outTree->Branch("genTauEtaL2",               &genTauEtaL2, "genTauEtaL2/F");
  outTree->Branch("genDecayModeL2",            &genDecayModeL2_, "genDecayModeL2/I");
  outTree->Branch("genMatchL1", &genMatchL1_, "genMatchL1/I");
  outTree->Branch("genMatchL2", &genMatchL2_, "genMatchL2/I");
  
  outTree->Branch("diTauCharge", &diTauCharge_,"diTauCharge/F");
  outTree->Branch("chargeL1", &chargeL1_,"chargeL1/F");

  outTree->Branch("MtLeg1MVA",   &MtLeg1MVA_,"MtLeg1MVA/F");
  outTree->Branch("MtLeg2MVA",   &MtLeg2MVA_,"MtLeg2MVA/F");
  outTree->Branch("MEtMVA",      &MEtMVA,     "MEtMVA/F");
  outTree->Branch("MEtMVAPhi",   &MEtMVAPhi,  "MEtMVAPhi/F");
  outTree->Branch("MEtMVACov00", &MEtMVACov00, "MEtMVACov00/F");
  outTree->Branch("MEtMVACov01", &MEtMVACov01, "MEtMVACov01/F");
  outTree->Branch("MEtMVACov10", &MEtMVACov10, "MEtMVACov10/F");
  outTree->Branch("MEtMVACov11", &MEtMVACov11, "MEtMVACov11/F");
  outTree->Branch("MEtCov00", &MEtCov00, "MEtCov00/F");
  outTree->Branch("MEtCov01", &MEtCov01, "MEtCov01/F");
  outTree->Branch("MEtCov10", &MEtCov10, "MEtCov10/F");
  outTree->Branch("MEtCov11", &MEtCov11, "MEtCov11/F");
  outTree->Branch("MEt",      &MEt,     "MEt/F");
  outTree->Branch("MEtPhi",   &MEtPhi,  "MEtPhi/F");
  outTree->Branch("MtLeg1",   &MtLeg1_,"MtLeg1/F");
  outTree->Branch("MtLeg2",   &MtLeg2_,"MtLeg2/F");

  //TauID
  outTree->Branch("decayModeFindingL1",&decayModeFindingL1_,"decayModeFindingL1/I");
  outTree->Branch("decayModeFindingNewDML1",&decayModeFindingNewDML1_,"decayModeFindingNewDML1/I");
  outTree->Branch("decayModeFindingOldDML1",&decayModeFindingOldDML1_,"decayModeFindingOldDML1/I");
  outTree->Branch("AntiEDeadEcalL1",&AntiEDeadEcalL1_,"AntiEDeadEcalL1/I");
  outTree->Branch("tightestAntiECutWPL1",&tightestAntiECutWPL1_,"tightestAntiECutWPL1/I");
  outTree->Branch("tightestAntiEMVA5WPL1",&tightestAntiEMVA5WPL1_,"tightestAntiEMVA5WPL1/I");
  outTree->Branch("AntiEMVA5categoryL1",&AntiEMVA5categoryL1_,"AntiEMVA5categoryL1/I");
  outTree->Branch("AntiEMVA5rawL1",&AntiEMVA5rawL1_,"AntiEMVA5rawL1/F");
  outTree->Branch("tightestAntiMuWPL1",&tightestAntiMuWPL1_,"tightestAntiMuWPL1/I");
  outTree->Branch("tightestAntiMu2WPL1",&tightestAntiMu2WPL1_,"tightestAntiMu2WPL1/I");
  outTree->Branch("tightestAntiMu3WPL1",&tightestAntiMu3WPL1_,"tightestAntiMu3WPL1/I");
  outTree->Branch("tightestAntiMuMVAWPL1",&tightestAntiMuMVAWPL1_,"tightestAntiMuMVAWPL1/I");
  outTree->Branch("AntiMuMVArawL1",&AntiMuMVArawL1_,"AntiMuMVArawL1/F");
  outTree->Branch("tightestHPSDBWPL1",&tightestHPSDBWPL1_,"tightestHPSDBWPL1/I");
  outTree->Branch("tightestHPSDB3HWPL1",&tightestHPSDB3HWPL1_,"tightestHPSDB3HWPL1/I");
  outTree->Branch("hpsDB3HL1",&hpsDB3HL1_,"hpsDB3HL1/F");
  outTree->Branch("tightestHPSMVA3newDMwLTWPL1",&tightestHPSMVA3newDMwLTWPL1_,"tightestHPSMVA3newDMwLTWPL1/I");
  outTree->Branch("hpsMVA3newDMwLTL1",&hpsMVA3newDMwLTL1_,"hpsMVA3newDMwLTL1/F");
  outTree->Branch("tightestHPSMVA3newDMwoLTWPL1",&tightestHPSMVA3newDMwoLTWPL1_,"tightestHPSMVA3newDMwoLTWPL1/I");
  outTree->Branch("hpsMVA3newDMwoLTL1",&hpsMVA3newDMwoLTL1_,"hpsMVA3newDMwoLTL1/F");
  outTree->Branch("tightestHPSMVA3oldDMwLTWPL1",&tightestHPSMVA3oldDMwLTWPL1_,"tightestHPSMVA3oldDMwLTWPL1/I");
  outTree->Branch("hpsMVA3oldDMwLTL1",&hpsMVA3oldDMwLTL1_,"hpsMVA3oldDMwLTL1/F");
  outTree->Branch("tightestHPSMVA3oldDMwoLTWPL1",&tightestHPSMVA3oldDMwoLTWPL1_,"tightestHPSMVA3oldDMwoLTWPL1/I");
  outTree->Branch("hpsMVA3oldDMwoLTL1",&hpsMVA3oldDMwoLTL1_,"hpsMVA3oldDMwoLTL1/F");
  outTree->Branch("decayModeL1",          &decayModeL1_,"decayModeL1/I");
  outTree->Branch("decayModeFindingL2",&decayModeFindingL2_,"decayModeFindingL2/I");
  outTree->Branch("decayModeFindingNewDML2",&decayModeFindingNewDML2_,"decayModeFindingNewDML2/I");
  outTree->Branch("decayModeFindingOldDML2",&decayModeFindingOldDML2_,"decayModeFindingOldDML2/I");
  outTree->Branch("AntiEDeadEcalL2",&AntiEDeadEcalL2_,"AntiEDeadEcalL2/I");
  outTree->Branch("tightestAntiECutWPL2",&tightestAntiECutWPL2_,"tightestAntiECutWPL2/I");
  outTree->Branch("tightestAntiEMVA5WPL2",&tightestAntiEMVA5WPL2_,"tightestAntiEMVA5WPL2/I");
  outTree->Branch("AntiEMVA5categoryL2",&AntiEMVA5categoryL2_,"AntiEMVA5categoryL2/I");
  outTree->Branch("AntiEMVA5rawL2",&AntiEMVA5rawL2_,"AntiEMVA5rawL2/F");
  outTree->Branch("tightestAntiMuWPL2",&tightestAntiMuWPL2_,"tightestAntiMuWPL2/I");
  outTree->Branch("tightestAntiMu2WPL2",&tightestAntiMu2WPL2_,"tightestAntiMu2WPL2/I");
  outTree->Branch("tightestAntiMu3WPL2",&tightestAntiMu3WPL2_,"tightestAntiMu3WPL2/I");
  outTree->Branch("tightestAntiMuMVAWPL2",&tightestAntiMuMVAWPL2_,"tightestAntiMuMVAWPL2/I");
  outTree->Branch("AntiMuMVArawL2",&AntiMuMVArawL2_,"AntiMuMVArawL2/F");
  outTree->Branch("tightestHPSDBWPL2",&tightestHPSDBWPL2_,"tightestHPSDBWPL2/I");
  outTree->Branch("tightestHPSDB3HWPL2",&tightestHPSDB3HWPL2_,"tightestHPSDB3HWPL2/I");
  outTree->Branch("hpsDB3HL2",&hpsDB3HL2_,"hpsDB3HL2/F");
  outTree->Branch("tightestHPSMVA3newDMwLTWPL2",&tightestHPSMVA3newDMwLTWPL2_,"tightestHPSMVA3newDMwLTWPL2/I");
  outTree->Branch("hpsMVA3newDMwLTL2",&hpsMVA3newDMwLTL2_,"hpsMVA3newDMwLTL2/F");
  outTree->Branch("tightestHPSMVA3newDMwoLTWPL2",&tightestHPSMVA3newDMwoLTWPL2_,"tightestHPSMVA3newDMwoLTWPL2/I");
  outTree->Branch("hpsMVA3newDMwoLTL2",&hpsMVA3newDMwoLTL2_,"hpsMVA3newDMwoLTL2/F");
  outTree->Branch("tightestHPSMVA3oldDMwLTWPL2",&tightestHPSMVA3oldDMwLTWPL2_,"tightestHPSMVA3oldDMwLTWPL2/I");
  outTree->Branch("hpsMVA3oldDMwLTL2",&hpsMVA3oldDMwLTL2_,"hpsMVA3oldDMwLTL2/F");
  outTree->Branch("tightestHPSMVA3oldDMwoLTWPL2",&tightestHPSMVA3oldDMwoLTWPL2_,"tightestHPSMVA3oldDMwoLTWPL2/I");
  outTree->Branch("hpsMVA3oldDMwoLTL2",&hpsMVA3oldDMwoLTL2_,"hpsMVA3oldDMwoLTL2/F");
  outTree->Branch("decayModeL2",          &decayModeL2_,"decayModeL2/I");

  outTree->Branch("nVetoMuon",          &nVetoMuon_,    "nVetoMuon/I");
  outTree->Branch("nVetoElectron",      &nVetoElectron_,    "nVetoElectron/I");
  outTree->Branch("numPV",              &numPV_,        "numPV/F");
  outTree->Branch("npu",       &npu_,      "npu/F");
  outTree->Branch("rho",       &rho_,      "rho/F");
  
  outTree->Branch("event",&event_,"event/l");
  outTree->Branch("run",  &run_,  "run/l");
  outTree->Branch("lumi", &lumi_, "lumi/l");
  outTree->Branch("index", &index_, "index/I");
  outTree->Branch("pairIndex", &pairIndex, "pairIndex/I");
  outTree->Branch("HLTx", &HLTx, "HLTx/I");
  outTree->Branch("HLTmatchL1", &HLTmatchL1, "HLTmatchL1/I");
  outTree->Branch("HLTmatchL2", &HLTmatchL2, "HLTmatchL2/I");
  
  outTree->Branch("weight", &evtweight, "weight/F");
  outTree->Branch("mcweight", &mcweight,"mcweight/F");
  outTree->Branch("puweight", &puweight, "puweight/F");
  outTree->Branch("isZtt", &isZtt_, "isZtt/O");
  outTree->Branch("isZl", &isZl_, "isZl/O");
  outTree->Branch("isZj", &isZj_, "isZj/O");
  
  // define JSON selector //
  cout << "JSON selection" << endl;
  int nJson=1;
  string jsonFile[nJson];
  string dirJson = "/nfs/dust/cms/user/anayak/CMS/OnSLC6/CMSSW_7414_htt/src/DesyTauAnalyses/NTupleMaker/test/"; 
  jsonFile[0] = dirJson+"/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"; // promptReco 
  map<int, vector<pair<int, int> > > jsonMap[nJson] ;  
  for(int iJ=0 ; iJ<nJson ; iJ++)
    jsonMap[iJ] = readJSONFile(jsonFile[iJ]);
  bool isGoodRun=false;

  //isData
  TString sample(sample_.c_str());
  bool isData = sample.Contains("2015");

  int nEntries    = currentTree->GetEntries() ;
  float crossSection = xsec_;
  //float scaleFactor = (crossSection != 0) ? Lumi / (  float(nEventsRead)/(crossSection*skimEff_) )  : 1.0;
  //first loop over the whole events once to get the total sum of gen weights for normalization
  float totalGenWeight_ = 0;
  if(!isData){
    for(int n = 0 ; n < nEntries ; n++) {
      currentTree->GetEntry(n);
      
      totalGenWeight_ += genweight;
    }
  }
  float scaleFactor = (crossSection != 0) ? (Lumi*crossSection) / float(totalGenWeight_)  : 1.0;

  int nProc,n1,n2 = nEntries;
  //cout<<"nDiv = "<<nDiv<<endl;
  //cout<<"nEntries = "<<nEntries<<endl;
  nProc = nEntries/ nDiv ;
  n1 = iDiv * nProc ;
  if( iDiv < (nDiv-1) )
    n2 = (iDiv+1) * nProc ;
  else if( iDiv == nDiv-1 )
    n2 = nEntries;
  
  //cout<<"n1 = "<<n1<<endl;
  //cout<<"n2 = "<<n2<<endl;

  //TString sample(sample_.c_str());
  cout << "Processing sample " << sample << endl;
  cout<< "nEventsRead = " << nEventsRead << endl;
  cout<< "nEntries    = " << nEntries << endl;
  cout<< "crossSection " << crossSection << " pb ==> scaleFactor " << scaleFactor << endl;

  //bool dyFinalState=false;

  ///////////////////////
  // LOOP OVER ENTRIES //
  ///////////////////////

  for(int n = n1 ; n < n2 ; n++) {
  //for(int n = n1 ; n < 1000 ; n++) { 
    //     cout<<"n/n2 = "<<n<<"/"<<n2<<endl;
    
    currentTree->GetEntry(n);
    if(n%1000==0) cout << n <<"/"<<(n2-n1)<< endl;
    //     if(n%1000==0) cout << n <<"/"<<nEntries<< endl;

    // APPLY JSON SELECTION //
    isGoodRun=true;

    if(iJson_>=0)
      isGoodRun = AcceptEventByRunAndLumiSection(event_run, event_luminosityblock, jsonMap[iJson_]);
    
    if(!isGoodRun) continue;
    if(event_nr == 30706)std::cout<<"pass none: event "<<event_nr<<std::endl;
    //cut on PV
    if(primvertex_count <= 0)continue;
    if(event_nr == 30706)std::cout<<"pass PV1: event "<<event_nr<<std::endl;
    //if(primvertex_ndof < 4)continue;
    //if(primvertex_z < -24 || primvertex_z > 24) continue;
    if(event_nr == 30706)std::cout<<"pass PV: event "<<event_nr<<std::endl;
    evtweight = scaleFactor;
    mcweight = genweight;
    puweight = 1.0; //no Pileup weight for now

    run_ = event_run;
    event_ = event_nr;
    lumi_ = event_luminosityblock;

    numPV_ = primvertex_count;
    npu_ = numtruepileupinteractions;
    rho_ = rhoNeutral;
    
    // Get Gen boson and daughters to clasify the event
    LV genVP4_(0, 0, 0, 0); unsigned int genVType_ = 100;
    std::vector<LV> genLeptonP4_; genLeptonP4_.clear();
    std::vector<LV> genTauP4_; genTauP4_.clear();
    std::vector<LV> genTauVisP4_; genTauVisP4_.clear();
    std::vector<int>genLeptonPdgId_; genLeptonPdgId_.clear();
    int nLepton_ = 0, TauLepDecay_ = 0, TauHadDecay_ = 0;

    /* //not needed to define this way
    if(!isData){ //for mc only

      if( (sample_.find("WJets")!=string::npos && sample_.find("WWJets")==string::npos ) ||
	  sample_.find("DYJets")!=string::npos || sample_.find("HToTauTau")!=string::npos ||
	  sample_.find("SUSYGGH")!=string::npos || sample_.find("SUSYBBH")!=string::npos ||
          sample_.find("GGFH")!=string::npos || sample_.find("VBFH")!=string::npos
	  )
	{
	  for(unsigned int ig = 0; ig < genparticles_count; ig++){

	    int GenPdgId = genparticles_pdgid[ig];
	    int GenStatus = genparticles_status[ig];
	    if(GenStatus == 1 && (fabs(GenPdgId) == 11 || fabs(GenPdgId) == 13)) {
	      
	      genLeptonP4_.push_back(LV(genparticles_px[ig], genparticles_py[ig], genparticles_pz[ig], genparticles_e[ig]));
	      genLeptonPdgId_.push_back(GenPdgId);
	      nLepton_++;
	      if(genparticles_mother[ig] == WBOSON || genparticles_mother[ig] == ZBOSON || genparticles_mother[ig] == HIGGS)
		genVType_ = genparticles_mother[ig];
	    }
	  }

	  for(unsigned int ig = 0; ig < gentau_count; ig++){
	    genTauP4_.push_back(LV(gentau_px[ig], gentau_py[ig], gentau_pz[ig], gentau_e[ig]));
	    genTauVisP4_.push_back(LV(gentau_visible_px[ig], gentau_visible_py[ig], gentau_visible_pz[ig], gentau_visible_e[ig]));
	    if(gentau_decayMode[ig] == 8 || gentau_decayMode[ig] == 9)
	      TauLepDecay_++;
	    else TauHadDecay_++;
	    if(gentau_mother[ig] == WBOSON || gentau_mother[ig] == ZBOSON || gentau_mother[ig] == HIGGS)
	      genVType_ = gentau_mother[ig];
	  }

	}
    }//end of MCInfo

    //Classify event, HADHAD, LEPHAD, LEPLEP, LL, LNU, LEPNU, TAUNU, OTHERS
    int EventType_ = OTHERS;
    if(nLepton_ >= 2 && TauLepDecay_ == 0 && TauHadDecay_<=1) 
      EventType_ = LL;
    else if (nLepton_ >= 1 && TauLepDecay_ == 0 && TauHadDecay_<=1)
      EventType_ = LNU;
    else if(TauLepDecay_ >= 2 && TauLepDecay_ == nLepton_ && TauHadDecay_<=1)
      EventType_ = LEPLEP;
    else if(TauLepDecay_ == 1 && TauLepDecay_ == nLepton_ && TauHadDecay_ == 0)
      EventType_ = LEPNU;
    else if(TauLepDecay_ == 1 && TauHadDecay_ >= 1 )
      EventType_ = LEPHAD;
    else if(TauLepDecay_ == 0 && TauHadDecay_ >= 2 )
      EventType_ = HADHAD;
    else if(TauLepDecay_ == 0 && nLepton_ == 0 && TauHadDecay_ == 1)
      EventType_ = TAUNU;
    */

    //require HLT
    string hltPath_(""); string hltFilter_("");
    if(!isData){
      hltPath_ = "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v";
      hltFilter_ = "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg";
    }
    else{
      if(sample.Contains("2015D")){
	hltPath_ = "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v";
	hltFilter_ = "hltDoublePFTau35TrackPt1MediumIsolationDz02Reg";
      }
      else{
	hltPath_ = "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v";
	hltFilter_ = "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg";
      }
    }
    if(GetTriggerResult((*hltriggerresultsV), hltPath_) < 0.5) continue;
    if(event_nr== 30706)std::cout<<"pass tigger: event "<<event_nr<<std::endl;
    std::vector<DiTauInfo>sortDiTauInfos; sortDiTauInfos.clear();
    //Loop over taus
    for(unsigned int it = 0; it < tau_count; it++){ //tauL1
      
      LV tauLeg1_(tau_px[it], tau_py[it], tau_pz[it], tau_e[it]);
      if(tauLeg1_.pt() < 45 || TMath::Abs(tauLeg1_.eta()) > 2.1) continue;
      //if(TMath::Abs(tau_vertexz[it] - primvertex_z) > 0.2) continue;
      //if(tau_vertexz[it] != primvertex_z) continue;
      if(TMath::Abs(tau_leadchargedhadrcand_dz[it]) >= 0.2) continue;

      if(tau_decayModeFindingNewDMs[it] < 0.5) continue;
      //if(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it] > 1.0) continue;
      //if(tau_againstElectronVLooseMVA5[it] < 0.5) continue;
      //if(tau_againstMuonLoose3[it] < 0.5) continue;
      if(TMath::Abs(tau_charge[it]) != 1) continue;
      if(event_nr== 30706)std::cout<<" pass offline sel. Leg1: event "<<event_nr<<std::endl;
      //HLT match
      bool HLTmatchLeg1_ = false;
      bool matchLeg1Level1_ = false; bool matchLeg1Level2_ = false; bool matchLeg1Level3_ = false;
      for(unsigned int it = 0; it < trigobject_count; it++){
        LV trigCandP4_(trigobject_px[it], trigobject_py[it], trigobject_pz[it], sqrt(trigobject_px[it]*trigobject_px[it] + trigobject_py[it]*trigobject_py[it] + trigobject_pz[it]*trigobject_pz[it]));

        //if(IsHLTMatched("hltL1sDoubleTauJet36erORDoubleTauJet68er", (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg1_) )matchLeg1Level1_ = true;
        //if(IsHLTMatched("hltDoubleL2IsoTau35eta2p1", (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg1_) )matchLeg1Level2_ = true;
        if(IsHLTMatched(hltFilter_, (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg1_) )matchLeg1Level3_ = true;
      }
      HLTmatchLeg1_ = matchLeg1Level3_; //(matchLeg1Level1_ && matchLeg1Level2_ && matchLeg1Level3_);
      if(!HLTmatchLeg1_) continue;
      if(event_nr== 30706)std::cout<<" found HLTLeg1 : event "<<event_nr<<std::endl;
      for(unsigned int jt = it+1; jt < tau_count; jt++){ //tauL2   
	    
	LV tauLeg2_(tau_px[jt], tau_py[jt], tau_pz[jt], tau_e[jt]);
	if(tauLeg2_.pt() < 45 || TMath::Abs(tauLeg2_.eta()) > 2.1) continue;
	if(TMath::Abs(tau_leadchargedhadrcand_dz[jt]) >= 0.2) continue;
	//if(TMath::Abs(tau_vertexz[jt] - primvertex_z) > 0.2) continue;
	//if(tau_vertexz[jt] != primvertex_z) continue;

	if(tau_decayModeFindingNewDMs[jt] < 0.5) continue;
	//if(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[jt] > 1.0) continue;
	//if(tau_againstElectronVLooseMVA5[jt] < 0.5) continue;
	//if(tau_againstMuonLoose3[jt] < 0.5) continue;
	if(TMath::Abs(tau_charge[jt]) != 1) continue;
	if(event_nr== 30706)std::cout<<" pass offline sel. Leg2: event "<<event_nr<<std::endl;
	//HLT match
	bool HLTmatchLeg2_ = false;
	bool matchLeg2Level1_ = false; bool matchLeg2Level2_ = false; bool matchLeg2Level3_ = false;
	for(unsigned int it = 0; it < trigobject_count; it++){
	  LV trigCandP4_(trigobject_px[it], trigobject_py[it], trigobject_pz[it], sqrt(trigobject_px[it]*trigobject_px[it] + trigobject_py[it]*trigobject_py[it] + trigobject_pz[it]*trigobject_pz[it]));

	  //if(IsHLTMatched("hltL1sDoubleTauJet36erORDoubleTauJet68er", (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg2_) )matchLeg2Level1_ = true;
	  //if(IsHLTMatched("hltDoubleL2IsoTau35eta2p1", (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg2_) )matchLeg2Level2_ = true;
	  if(IsHLTMatched(hltFilter_, (*run_hltfilters), trigCandP4_, trigobject_filters[it], tauLeg2_) )matchLeg2Level3_ = true;
	}
	HLTmatchLeg2_ = matchLeg2Level3_; //(matchLeg2Level1_ && matchLeg2Level2_ && matchLeg2Level3_);
	if(!HLTmatchLeg2_) continue;
	if(event_nr== 30706)std::cout<<" found HLTLeg2 : event "<<event_nr<<std::endl;
	if(ROOT::Math::VectorUtil::DeltaR(tauLeg1_, tauLeg2_) < 0.5) continue;
	if(event_nr== 30706)std::cout<<" pass DeltaR filter : event "<<event_nr<<std::endl;
	float sumPt = tauLeg1_.pt() + tauLeg2_.pt();
	float sumIso = tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it] + tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[jt];
	int pairCharge = tau_charge[it]*tau_charge[jt];
	DiTauInfo sortDiTauInfo; 
	sortDiTauInfo.index1_ = it;
	sortDiTauInfo.index2_ = jt;
	sortDiTauInfo.sumPt_ = sumPt; 
	sortDiTauInfo.sumIso_ = sumIso;
	sortDiTauInfo.diTauCharge_ = pairCharge; 
	sortDiTauInfos.push_back(sortDiTauInfo); 
      }
    }
  

    //sort diTaus, 1st according to sumIso, then sumPt  
    std::sort(sortDiTauInfos.begin(), sortDiTauInfos.end(), SortDiTauPairs()); 

    int diTauCounter = -1;
    for(std::vector<DiTauInfo>::iterator iter = sortDiTauInfos.begin();  iter != sortDiTauInfos.end() ; iter++){ 
      if(diTauCounter >= 0) continue;
      //diTauCounter++;
      //pairIndex   = diTauCounter;
      if(event_nr== 30706)std::cout<<"pass selection: event "<<event_nr<<std::endl;
      int tau1 = iter->index1_;
      int tau2 = iter->index2_;

      LV temp_Leg1_(tau_px[tau1], tau_py[tau1], tau_pz[tau1], tau_e[tau1]);
      LV temp_Leg2_(tau_px[tau2], tau_py[tau2], tau_pz[tau2], tau_e[tau2]);
      LV Leg1P4_, Leg2P4_;
      if(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau1] < tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau2]){
      //if(temp_Leg1_.pt() > temp_Leg2_.pt()){
	Leg1P4_ = temp_Leg1_; 
	Leg2P4_ = temp_Leg2_;
      }
      else if(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau1] > tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau2]){
	Leg1P4_ = temp_Leg2_; tau1 = iter->index2_;
        Leg2P4_ = temp_Leg1_; tau2 = iter->index1_;
      }
      else if(temp_Leg1_.pt() > temp_Leg2_.pt()){
	Leg1P4_ = temp_Leg1_;
        Leg2P4_ = temp_Leg2_;
      }
      else {
	Leg1P4_ = temp_Leg2_; tau1 = iter->index2_;
	Leg2P4_ = temp_Leg1_; tau2 = iter->index1_;
      }
	  

      // FIND matched gen lepton
      LV Leg1GenP4_(0, 0, 0, 0);
      LV Leg2GenP4_(0, 0, 0, 0);
      //if(EventType_ == HTTHADHAD || EventType_ == ZTTHADHAD || EventType_ == WTAUHAD){
      for(size_t ij = 0; ij < genTauVisP4_.size(); ij++){
	if(ROOT::Math::VectorUtil::DeltaR(genTauVisP4_[ij], Leg1P4_) < 0.3 )
	  Leg1GenP4_ = genTauVisP4_[ij];
	else if(ROOT::Math::VectorUtil::DeltaR(genTauVisP4_[ij], Leg2P4_) < 0.3 )
	  Leg2GenP4_ = genTauVisP4_[ij];
      }
      //}
      
      ///////Check Gen Lepton matching
      LV Leg1GenLepP4_(0, 0, 0, 0);
      LV Leg2GenLepP4_(0, 0, 0, 0);
      for(size_t il = 0; il < genLeptonP4_.size(); il++){
	if(ROOT::Math::VectorUtil::DeltaR(genLeptonP4_[il], Leg1P4_) < 0.3 )
	  Leg1GenLepP4_ = genLeptonP4_[il];
	else if(ROOT::Math::VectorUtil::DeltaR(genLeptonP4_[il], Leg2P4_) < 0.3 )
	  Leg2GenLepP4_ = genLeptonP4_[il];
      }

      //////Select and Fill Jets
      ptj1 = -999; ptj2 = -999.; etaj1 = -999.; etaj2 = -999.; Detajj = -999.; Mjj = -999.; Dphijj = -999.; phij1 = -999.; phij2 = -999.;
      ptB1 = -999.; etaB1 = -999.; phiB1 = -999.;
      diJetPt = -999.; diJetPhi = -999.;
      nJets30 = -99; nJets20 = -99; nJets20BTagged = -99; 
      nJets20BTaggedBUp = -99; nJets20BTaggedBDown = -99; nJets20BTaggedLUp = -99; nJets20BTaggedLDown = -99; nJets20BTaggedLoose = -99;
      ptVeto = -999.; etaVeto = -999.; phiVeto = -999.;
      nVetoJets  = 0;


      nJets30 = 0; nJets20 = 0;
      nJets20BTagged = 0; nJets20BTaggedLoose = 0;
      std::vector<LV> pfJetsP4_; pfJetsP4_.clear();
      std::vector<int> pfJetsIndex_; pfJetsIndex_.clear();
      for(unsigned int ijet = 0; ijet < pfjet_count; ijet++){
	LV JetP4_(pfjet_px[ijet], pfjet_py[ijet], pfjet_pz[ijet], pfjet_e[ijet]);

	//Exclude leptons/taus
	if(ROOT::Math::VectorUtil::DeltaR(Leg1P4_, JetP4_) < 0.5 ||
	   ROOT::Math::VectorUtil::DeltaR(Leg2P4_, JetP4_) < 0.5) continue;

	//apply PF JetID
	bool passLooseJetID = false;
	float totalPFEnergy = (JetP4_*pfjet_energycorr[ijet]).energy(); 
	float NHF = pfjet_neutralhadronicenergy[ijet]/totalPFEnergy;
	float NEMF = pfjet_neutralemenergy[ijet]/totalPFEnergy;
	float NumConstituent = pfjet_chargedmulti[ijet]+pfjet_neutralmulti[ijet];
	float CHF = pfjet_chargedhadronicenergy[ijet]/totalPFEnergy;
	float CHM = pfjet_chargedmulti[ijet];
	float CEMF = pfjet_chargedemenergy[ijet]/totalPFEnergy;
	float NumNeutral = pfjet_neutralmulti[ijet];

	if((TMath::Abs(JetP4_.eta())>3.0 && 
	    NEMF < 0.90 && NumNeutral > 10) || 
	   (TMath::Abs(JetP4_.eta())<=3.0 && 
	    NHF<0.99 && NEMF<0.99 && NumConstituent > 1 && 
	    ((TMath::Abs(JetP4_.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || 
	      TMath::Abs(JetP4_.eta())>2.4)
	    )
	   ) passLooseJetID = true;
	if(!passLooseJetID) continue;
	//Apply PU JetID
	//if(getJetIDMVALoose(JetP4_.pt(), JetP4_.eta(), pfjet_pu_jet_full_mva[ijet]) < 0.5) continue;
	if(JetP4_.pt() < 20 || TMath::Abs(JetP4_.eta()) > 4.7) continue;
	nJets20++;
	if(JetP4_.pt() > 30)nJets30++;
	pfJetsP4_.push_back(JetP4_);
	pfJetsIndex_.push_back(ijet);
	//count b-tagged jets
	if(TMath::Abs(JetP4_.eta()) < 2.4){
	  //int jetFlavour = pfjet_flavour[ijet]; //get it from tree
	  //bool isBtag = btsf->isbtagged(JetP4_.Pt(), JetP4_.Eta(), pfjet_btag[ijet][2], jetFlavour, isData ,kNo, kNo, true); //use CSV Medium WP
	  float BDiscr = GetBTagDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags", (*run_btagdiscriminators), pfjet_btag[ijet]);
	  bool isBtag = (BDiscr > 0.89); //pfCombinedInclusiveSecondaryVertexV2BJetTags
	  //bool isBtag = (pfjet_btag[ijet][8] > 0.814); //pfCombinedInclusiveSecondaryVertexV2BJetTags

	  if(isBtag){
	    nJets20BTagged++;
	    if(nJets20BTagged < 2){
	      ptB1  = JetP4_.Pt();
	      phiB1 = JetP4_.Phi();
	      etaB1 = JetP4_.Eta();
	    }
	  }

	  if(pfjet_btag[ijet][8] > 0.244){ //Loose WP
	    nJets20BTaggedLoose++;
	  }
	}
      }
      
      //compute VBF info
      if(pfJetsP4_.size() >= 1){
	ptj1 = pfJetsP4_[0].pt();
	etaj1 = pfJetsP4_[0].eta();
	phij1 = pfJetsP4_[0].phi();
	pumvaj1 = pfjet_pu_jet_full_mva[pfJetsIndex_[0]];
	csvj1 = pfjet_btag[pfJetsIndex_[0]][6];
	ptrawj1 = pfjet_energycorr[pfJetsIndex_[0]]*ptj1;
	if(pfJetsP4_.size() >= 2){
	  ptj2 = pfJetsP4_[1].pt();
	  etaj2 = pfJetsP4_[1].eta();
	  phij2 = pfJetsP4_[1].phi();
	  pumvaj2 = pfjet_pu_jet_full_mva[pfJetsIndex_[1]];
	  csvj2 = pfjet_btag[pfJetsIndex_[1]][6];
	  ptrawj2 = pfjet_energycorr[pfJetsIndex_[1]]*ptj2;

	  Detajj = TMath::Abs(pfJetsP4_[0].eta() - pfJetsP4_[1].eta());
	  Mjj = (pfJetsP4_[0] + pfJetsP4_[1]).mass();
	  Dphijj = ROOT::Math::VectorUtil::DeltaPhi(pfJetsP4_[0], pfJetsP4_[1]);
	  diJetPt = (pfJetsP4_[0] + pfJetsP4_[1]).pt();
	  diJetPhi = (pfJetsP4_[0] + pfJetsP4_[1]).phi();
	  
	  //check CJV
	  if(pfJetsP4_.size() > 2 ){
	    ptVeto = pfJetsP4_[2].pt();
	    etaVeto = pfJetsP4_[2].eta();
	    phiVeto = pfJetsP4_[2].phi();

	    nVetoJets = 0;
	    for(size_t ijet = 2; ijet < pfJetsP4_.size(); ijet++){
	      if(pfJetsP4_[ijet].pt() > 30 && (pfJetsP4_[ijet].eta() - pfJetsP4_[0].eta())*(pfJetsP4_[ijet].eta() - pfJetsP4_[1].eta()) <= 0)
		nVetoJets++;
	    }
	  }
	}
      }
      

      //Fill Leg infos
      decayModeL1_ = 0;
      if(tau_signalChargedHadrCands_size[tau1] == 1 && tau_signalGammaCands_size[tau1]<=0)decayModeL1_ = 1; 
      else if(tau_signalChargedHadrCands_size[tau1] == 1 && tau_signalGammaCands_size[tau1]>0)decayModeL1_ = 2;
      else if(tau_signalChargedHadrCands_size[tau1] == 3) decayModeL1_ = 3;
      decayModeFindingL1_ = tau_decayModeFinding[tau1];
      decayModeFindingNewDML1_ = tau_decayModeFindingNewDMs[tau1];
      //decayModeFindingOldDML1_ = GetTauDiscriminator("decayModeFindingOldDMs", run_taudiscriminators, tau_dishps[tau1]);
      tightestAntiEMVA5WPL1_ = 0;
      if(tau_againstElectronVLooseMVA5[tau1] > 0.5)tightestAntiEMVA5WPL1_ = 1;
      if(tau_againstElectronLooseMVA5[tau1] > 0.5)tightestAntiEMVA5WPL1_ = 2;
      if(tau_againstElectronMediumMVA5[tau1] > 0.5)tightestAntiEMVA5WPL1_ = 3;
      if(tau_againstElectronTightMVA5[tau1] > 0.5)tightestAntiEMVA5WPL1_ = 4;
      tightestAntiMu3WPL1_ = 0;
      if(tau_againstMuonLoose3[tau1] > 0.5)tightestAntiMu3WPL1_ = 1;
      if(tau_againstMuonTight3[tau1] > 0.5)tightestAntiMu3WPL1_ = 2;
      tightestHPSDB3HWPL1_ = 0;
      if(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tau1] > 0.5) tightestHPSDB3HWPL1_ = 1;
      if(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau1] > 0.5) tightestHPSDB3HWPL1_ = 2;
      if(tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tau1] > 0.5) tightestHPSDB3HWPL1_ = 3;
      tightestHPSMVA3oldDMwLTWPL1_ = 0;
      hpsDB3HL1_ = tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau1];

      decayModeL2_ = 0;
      if(tau_signalChargedHadrCands_size[tau2] == 1 && tau_signalGammaCands_size[tau2]<=0)decayModeL2_ = 1;
      else if(tau_signalChargedHadrCands_size[tau2] == 1 && tau_signalGammaCands_size[tau2]>0)decayModeL2_ = 2;
      else if(tau_signalChargedHadrCands_size[tau2] == 3) decayModeL2_ = 3;
      decayModeFindingL2_ = tau_decayModeFinding[tau2];
      decayModeFindingNewDML2_ = tau_decayModeFindingNewDMs[tau2];
      tightestAntiEMVA5WPL2_ = 0;
      if(tau_againstElectronVLooseMVA5[tau2] > 0.5)tightestAntiEMVA5WPL2_ = 1;
      if(tau_againstElectronLooseMVA5[tau2] > 0.5)tightestAntiEMVA5WPL2_ = 2;
      if(tau_againstElectronMediumMVA5[tau2] > 0.5)tightestAntiEMVA5WPL2_ = 3;
      if(tau_againstElectronTightMVA5[tau2] > 0.5)tightestAntiEMVA5WPL2_ = 4;
      tightestAntiMu3WPL2_ = 0;
      if(tau_againstMuonLoose3[tau2] > 0.5)tightestAntiMu3WPL2_ = 1;
      if(tau_againstMuonTight3[tau2] > 0.5)tightestAntiMu3WPL2_ = 2;
      tightestHPSDB3HWPL2_ = 0;
      if(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tau2] > 0.5) tightestHPSDB3HWPL2_ = 1;
      if(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau2] > 0.5) tightestHPSDB3HWPL2_ = 2;
      if(tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tau2] > 0.5) tightestHPSDB3HWPL2_ = 3;
      tightestHPSMVA3oldDMwLTWPL2_ = 0;
      hpsDB3HL2_ = tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau2];

      //Apply Tau ES Correction, switch off for now
      //LV Leg1GenJetP4_(tau_genjet_px[tau1], tau_genjet_py[tau1], tau_genjet_pz[tau1], tau_genjet_e[tau1]);
      //LV Leg2GenJetP4_(tau_genjet_px[tau2], tau_genjet_py[tau2], tau_genjet_pz[tau2], tau_genjet_e[tau2]);
      //Leg1P4_ = TauEnergyCorrector(Leg1P4_, Leg1GenJetP4_, decayModeL1_);
      //Leg2P4_ = TauEnergyCorrector(Leg2P4_, Leg2GenJetP4_, decayModeL2_);
      
      ptL1 = Leg1P4_.pt(); etaL1 = Leg1P4_.eta(); phiL1 = Leg1P4_.phi();
      ptL2 = Leg2P4_.pt(); etaL2 = Leg2P4_.eta(); phiL2 = Leg2P4_.phi();
      visibleTauMassL1 = tau_mass[tau1]; visibleTauMassL2 = tau_mass[tau2];
      diTauCharge_ = iter->diTauCharge_;
      chargeL1_ = tau_charge[tau1];
      dxyL1 = tau_leadchargedhadrcand_dxy[tau1]; dxyL2 = tau_leadchargedhadrcand_dxy[tau2];
      dzL1 = tau_leadchargedhadrcand_dz[tau1]; dzL2 = tau_leadchargedhadrcand_dz[tau2];

      //Apply MET Recoil Correction
      //For Z+Jets, W+jets and Higgs
      int mvamet_index = -1;
      for(UInt_t imet = 0; imet < mvamet_count; imet++){
      	//std::cout<<"mva met channel "<<mvamet_channel[imet]<<std::endl;
      	if(mvamet_channel[imet] == TAUTAU){
	  if((int(mvamet_lep1[imet]) == tau1 && int(mvamet_lep2[imet]) == tau2) ||
	     (int(mvamet_lep1[imet]) == tau2 && int(mvamet_lep2[imet]) == tau1))
	    mvamet_index = imet;
	}
      }
      if(mvamet_index < 0) { std::cout<<"mvamet_index < 0"<<std::endl; break;}

      LV mvaMetP4_(mvamet_ex[mvamet_index], mvamet_ey[mvamet_index], 0, 
		   sqrt(mvamet_ex[mvamet_index]*mvamet_ex[mvamet_index] + mvamet_ey[mvamet_index]*mvamet_ey[mvamet_index]));
      double newMvaMetPt_ = mvaMetP4_.Pt();
      double newMvaMetPhi_ = mvaMetP4_.Phi();
      //No recoil correction for now

      TLorentzVector newMvaMetP4_;
      newMvaMetP4_.SetPtEtaPhiM(newMvaMetPt_, 0, newMvaMetPhi_, 0);
      double scaledMetPx_ = newMvaMetP4_.Px(); double scaledMetPy_ = newMvaMetP4_.Py();
      LV corMvaMetP4_(scaledMetPx_, scaledMetPy_, 0, sqrt(scaledMetPx_*scaledMetPx_ + scaledMetPy_*scaledMetPy_));
      
      float scalarSumPtLeg1MVA  = Leg1P4_.pt() + corMvaMetP4_.pt();
      float vectorSumPtLeg1MVA  = (Leg1P4_ +  corMvaMetP4_).pt() ;
      MtLeg1MVA_  = TMath::Sqrt( scalarSumPtLeg1MVA*scalarSumPtLeg1MVA - vectorSumPtLeg1MVA*vectorSumPtLeg1MVA ) ;

      float scalarSumPtLeg2MVA  = Leg2P4_.pt() + corMvaMetP4_.pt();
      float vectorSumPtLeg2MVA  = (Leg2P4_ +  corMvaMetP4_).pt() ;
      MtLeg2MVA_  = TMath::Sqrt( scalarSumPtLeg2MVA*scalarSumPtLeg2MVA - vectorSumPtLeg2MVA*vectorSumPtLeg2MVA ) ;

      MEtMVA = newMvaMetPt_;
      MEtMVAPhi = newMvaMetPhi_;
      MEtMVACov00 = mvamet_sigxx[mvamet_index];
      MEtMVACov01 = mvamet_sigxy[mvamet_index];
      MEtMVACov10 = mvamet_sigyx[mvamet_index];
      MEtMVACov11 = mvamet_sigyy[mvamet_index];
      //also store pf met
      LV pfMetP4_(pfmet_ex, pfmet_ey, 0, sqrt(pfmet_ex*pfmet_ex + pfmet_ey*pfmet_ey));
      MEt = pfMetP4_.Pt();
      MEtPhi = pfMetP4_.Phi();
      float scalarSumPtLeg1  = Leg1P4_.pt() + pfMetP4_.pt();
      float vectorSumPtLeg1  = (Leg1P4_ +  pfMetP4_).pt() ;
      MtLeg1_  = TMath::Sqrt( scalarSumPtLeg1*scalarSumPtLeg1 - vectorSumPtLeg1*vectorSumPtLeg1 ) ;
      
      float scalarSumPtLeg2  = Leg2P4_.pt() + pfMetP4_.pt();
      float vectorSumPtLeg2  = (Leg2P4_ +  pfMetP4_).pt() ;
      MtLeg2_  = TMath::Sqrt( scalarSumPtLeg2*scalarSumPtLeg2 - vectorSumPtLeg2*vectorSumPtLeg2 ) ;

      MEtCov00 = pfmet_sigxx;
      MEtCov01 = pfmet_sigxy;
      MEtCov10 = pfmet_sigyx;
      MEtCov11 = pfmet_sigyy;

      /* no SVfit for now
      /////Perform SVFit to compute di-tau mass
      std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, Leg1P4_.pt(), Leg1P4_.eta(), Leg1P4_.phi(), Leg1P4_.mass(), tau_decayMode[tau1]));
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, Leg2P4_.pt(), Leg2P4_.eta(), Leg2P4_.phi(), Leg2P4_.mass(), tau_decayMode[tau2]));
      //NSVfitStandalone::Vector measuredMET( corMvaMetP4_.Px(), corMvaMetP4_.Py(), 0); 
      TMatrixD covMET(2,2); 
      covMET[0][0] = mvamet_sigxx[mvamet_index];  
      covMET[0][1] = mvamet_sigxy[mvamet_index];
      covMET[1][0] = mvamet_sigyx[mvamet_index];  
      covMET[1][1] = mvamet_sigyy[mvamet_index];  

      SVfitStandaloneAlgorithm algo(measuredTauLeptons, corMvaMetP4_.Px(), corMvaMetP4_.Py(), covMET, 0);
      algo.addLogM(false); 

      edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
      TH1::AddDirectory(false);  
      TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
      algo.shiftVisPt(true, inputFile_visPtResolution);
      algo.integrateMarkovChain();

      diTauNSVfitMass_ = algo.getMass(); 
      diTauNSVfitMassErrUp_ = (diTauNSVfitMass_ + algo.massUncert());
      diTauNSVfitMassErrDown_ = (diTauNSVfitMass_ - algo.massUncert());
      diTauNSVfitPt_ = algo.pt();  
      diTauNSVfitPtErrUp_ = (diTauNSVfitPt_ + algo.ptUncert());
      diTauNSVfitPtErrDown_ = (diTauNSVfitPt_ - algo.ptUncert());
      diTauNSVfitEta_ = algo.eta();
      diTauNSVfitPhi_ = algo.phi();

      delete inputFile_visPtResolution;
      */

      //DiTau Visible Mass
      LV diTauVisP4_ = Leg1P4_+Leg2P4_;
      diTauVisMass = diTauVisP4_.mass();
      diTauVisPt = diTauVisP4_.pt();
      diTauVisEta = diTauVisP4_.eta();
      diTauVisPhi = diTauVisP4_.phi();
      diTauRecoPt = (diTauVisP4_ + corMvaMetP4_).pt();
      diTauRecoPhi = (diTauVisP4_ + corMvaMetP4_).phi();
      
      
      //compute 3rd lepton veto
      nVetoElectron_ = 0;
      for(unsigned int ie = 0; ie < electron_count; ie++){
	LV eleP4_(electron_px[ie], electron_py[ie], electron_pz[ie], sqrt(electron_px[ie]*electron_px[ie] + electron_py[ie]*electron_py[ie] + electron_pz[ie]*electron_pz[ie]));
	float scEta = electron_superclusterEta[ie];
	if(eleP4_.pt() < 10 || TMath::Abs(eleP4_.eta()) > 2.5) continue;
	if(TMath::Abs(electron_dxy[ie]) > 0.045 || TMath::Abs(electron_dz[ie]) > 0.2) continue;
	
	
	//Isolation
	float eleIso_ = (electron_r03_sumChargedHadronPt[ie] + std::max(electron_r03_sumPhotonEt[ie]+electron_r03_sumNeutralHadronEt[ie] - 0.5*electron_r03_sumPUPt[ie], 0.0))/eleP4_.pt(); 
	//MVA Id, 90% Eff WP
	/*bool pass_eid_ = ((TMath::Abs(scEta)< 0.8 && electron_mva_id_nontrigPhys14[ie] > 0.933 ) || 
			  (TMath::Abs(scEta)> 0.8 && TMath::Abs(scEta)< 1.479 && electron_mva_id_nontrigPhys14[ie] > 0.825 ) ||
			  (TMath::Abs(scEta)> 1.479 && electron_mva_id_nontrigPhys14[ie] > 0.337 )
			  );
	*/
	bool pass_eid_ = (electron_mva_wp90_nontrig_Spring15_v1[ie] > 0);

	bool pass_conversion_ = (electron_pass_conversion[ie] == true && electron_nmissinginnerhits[ie] <= 1);

	if(pass_eid_ && pass_conversion_ && eleIso_ < 0.3)nVetoElectron_++;
      } 

      nVetoMuon_ = 0;
      for(unsigned int im = 0; im < muon_count; im++){
	LV muP4_(muon_px[im], muon_py[im], muon_pz[im], sqrt(muon_px[im]*muon_px[im] + muon_py[im]*muon_py[im] + muon_pz[im]*muon_pz[im]));
	if(muP4_.pt() < 10 || TMath::Abs(muP4_.eta()) > 2.4) continue;
	if(TMath::Abs(muon_dxy[im]) > 0.045 || TMath::Abs(muon_dz[im]) > 0.2) continue;

	bool pass_muid_ = muon_isMedium[im]; 
	float muIso_ = (muon_r03_sumChargedHadronPt[im] + std::max(muon_r03_sumPhotonEt[im]+muon_r03_sumNeutralHadronEt[im] - 0.5*muon_r03_sumPUPt[im], 0.0)) / muP4_.pt();

	if(pass_muid_ && muIso_ < 0.3)nVetoMuon_++;
      }

      //Get Trigger Informations
      HLTx =  GetTriggerResult((*hltriggerresultsV), "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v");
      HLTmatchL1 = 0; HLTmatchL2 = 0;
      bool matchLeg1Level1_ = false; bool matchLeg1Level2_ = false; bool matchLeg1Level3_ = false;
      bool matchLeg2Level1_ = false; bool matchLeg2Level2_ = false; bool matchLeg2Level3_ = false;
      for(unsigned int it = 0; it < trigobject_count; it++){
        LV trigCandP4_(trigobject_px[it], trigobject_py[it], trigobject_pz[it], sqrt(trigobject_px[it]*trigobject_px[it] + trigobject_py[it]*trigobject_py[it] + trigobject_pz[it]*trigobject_pz[it]));
	if(IsHLTMatched(hltFilter_, (*run_hltfilters), trigCandP4_, trigobject_filters[it], Leg1P4_) )matchLeg1Level3_ = true;
	if(IsHLTMatched(hltFilter_, (*run_hltfilters), trigCandP4_, trigobject_filters[it], Leg2P4_) )matchLeg2Level3_ = true;

      }
      HLTmatchL1 = matchLeg1Level3_; 
      HLTmatchL2 = matchLeg2Level3_; 
      
      //Define MC DY event type (used to split DY events)
      isZtt_=false; isZl_=false; isZj_=false;
      genMatchL1_ = 99; genMatchL2_ = 99;
      if(sample_.find("DYJets") != std::string::npos){
	//std::cout<<"sample "<<sample_<<" type "<<EventType_<<std::endl;
	int Leg1GenMatch_ = 6; int Leg2GenMatch_ = 6;
	for(unsigned int ig = 0; ig < genparticles_count; ig++){
	  int GenPdgId = genparticles_pdgid[ig];
	  int GenStatus = genparticles_status[ig];
	  //check electon
	  if(fabs(GenPdgId) == 11 && (genparticles_isPrompt[ig] > 0 || genparticles_isDirectPromptTauDecayProduct[ig] > 0)) {
	    LV genElecP4_(genparticles_px[ig], genparticles_py[ig], genparticles_pz[ig], genparticles_e[ig]);
	    if(genElecP4_.pt() > 8){
	      if(deltaR(Leg1P4_, genElecP4_) < 0.2){
		if(genparticles_isPrompt[ig] > 0)Leg1GenMatch_ = 1;
		else if(genparticles_isDirectPromptTauDecayProduct[ig] > 0)Leg1GenMatch_ = 3;
	      }
	      if(deltaR(Leg2P4_, genElecP4_) < 0.2){
                if(genparticles_isPrompt[ig] > 0)Leg2GenMatch_ = 1;
		else if(genparticles_isDirectPromptTauDecayProduct[ig] > 0)Leg2GenMatch_ = 3;
              }
	    }
	  }
	  //check muon
	  if(fabs(GenPdgId) == 13 && GenStatus == 1 && (genparticles_isPrompt[ig] > 0 || genparticles_isDirectPromptTauDecayProduct[ig] > 0)) {
            LV genMuonP4_(genparticles_px[ig], genparticles_py[ig], genparticles_pz[ig], genparticles_e[ig]);
            if(genMuonP4_.pt() > 8){
              if(deltaR(Leg1P4_, genMuonP4_) < 0.2){
                if(genparticles_isPrompt[ig] > 0)Leg1GenMatch_ = 2;
		else if(genparticles_isDirectPromptTauDecayProduct[ig] > 0)Leg1GenMatch_ = 4;
              }
	      if(deltaR(Leg2P4_, genMuonP4_) < 0.2){
		if(genparticles_isPrompt[ig] > 0)Leg2GenMatch_ = 2;
                else if(genparticles_isDirectPromptTauDecayProduct[ig] > 0)Leg2GenMatch_ = 4;
              }
            }
          }
	}
	//check tau
	for(unsigned int ig = 0; ig < gentau_count; ig++){                                         
	  if(gentau_isPrompt[ig] <= 0) continue;
	  LV genTauVisP4_(gentau_visible_px[ig], gentau_visible_py[ig], gentau_visible_pz[ig], gentau_visible_e[ig]); 
	  if(genTauVisP4_.pt() > 15){
	    if(deltaR(Leg1P4_, genTauVisP4_) < 0.2) Leg1GenMatch_ = 5;
	    if(deltaR(Leg2P4_, genTauVisP4_) < 0.2) Leg2GenMatch_ = 5;
	  }
	}

	if(Leg1GenMatch_ == 5 && Leg2GenMatch_ == 5) isZtt_=true;
	else if(Leg1GenMatch_ < 5 || Leg2GenMatch_ < 5) isZl_= true;
	else if(Leg1GenMatch_ == 6 || Leg2GenMatch_ == 6) isZj_ = true;

	genMatchL1_ = Leg1GenMatch_; genMatchL2_ = Leg2GenMatch_;
      }//end of gen matching

      pairIndex = -1;
      if(HLTx && HLTmatchL1 && HLTmatchL2){
	diTauCounter++; 
	pairIndex   = diTauCounter;
      }
      
      outTree->Fill(); //fill tree for each pair

    }//end of di-tau pair

  }

  delete hltriggerresultsV; delete run_hltfilters;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[])
{
  
  std::cout << "preAnalyzerTauTau_Summer15" << std::endl;
  gROOT->SetBatch(true);
  
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  

  //--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  
  //--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("PreAnalyzerTauTau") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgPreAnalyzerTauTau = cfg.getParameter<edm::ParameterSet>("preAnalyzerTauTau");

  std::string sample = cfgPreAnalyzerTauTau.getParameter<std::string>("sample");
  std::string analysis = cfgPreAnalyzerTauTau.getParameter<std::string>("analysis");
  double xSection = cfgPreAnalyzerTauTau.getParameter<double>("xSection");
  double skimEff = cfgPreAnalyzerTauTau.getParameter<double>("skimEff");
  int iJson = cfgPreAnalyzerTauTau.getParameter<int>("iJson");
  int iDiv = cfgPreAnalyzerTauTau.getParameter<int>("iDiv");
  int nDiv = cfgPreAnalyzerTauTau.getParameter<int>("nDiv");

  fwlite::InputSource inputFiles(cfg); 
  //int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  string analysisFileName = analysis;
  if( !(analysis.find("Up")!=string::npos || analysis.find("Down")!=string::npos) &&  analysis.find("Raw")==string::npos)
    analysisFileName = "Nominal";
  if( !(analysis.find("Up")!=string::npos || analysis.find("Down")!=string::npos) &&  analysis.find("Raw")!=string::npos)
    analysisFileName = "RawNominal";

  cout << "Now skimming analysis " << analysis << endl;
  if(analysis=="nominal") analysis="";

  TTree* outTree = fs.make<TTree>(TString(("outTree"+analysis).c_str()),"tree jets pT-ord");

  double nEventsRead = 0;
  /*
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end(); ++inputFileName ) {
    //--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("PreAnalyzerTauTau") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    TString histoName("allEventsFilter/totalEvents");
    TH1D* histo =(TH1D*)inputFile->Get(histoName);
    //     cout<<"histo "<<histo->GetEntries();
    if(histo)nEventsRead += histo->GetBinContent(1) ;
    else throw cms::Exception("PreAnalyzerTauTau") 
    	   << "Failed to read histogram "<<histoName<<" from inputFile = " << (*inputFileName) << " !!\n";
    //--- close input file
    delete inputFile;
  }
  cout<< "nEventsRead " << nEventsRead << endl;
  */
  /*string anlyzerName = analysis;
  if( analysis.find("Jet")!=string::npos && analysis.find("Raw")==string::npos)
    anlyzerName = "";
  if( analysis.find("Jet")!=string::npos && analysis.find("Raw")!=string::npos)
    anlyzerName = "Raw";
  */
  TString treeName("makeroottree/AC1B");
  TChain* currentTree = new TChain (treeName);
  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {
    currentTree->Add(inputFileName->data());
  }
  nEventsRead = currentTree->GetEntries();
  cout<< "nEventsRead " << nEventsRead << endl;
  cout<<"iDiv = "<<iDiv<<endl;
  cout<<"nDiv = "<<nDiv<<endl;


  fillTrees_TauTauStream(currentTree,outTree,nEventsRead,analysis,sample,xSection,skimEff,iJson,iDiv,nDiv);

  //delete outTree;
  //delete currentTree;
  return 0;
}
