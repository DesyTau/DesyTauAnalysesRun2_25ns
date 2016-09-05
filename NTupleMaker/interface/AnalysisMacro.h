
//#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/Basic_Mt2_332_Calculator.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/Basic_MPairProd_Calculator.h"

#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/lester_mt2_bisect.h"
#include "DesyTauAnalyses/NTupleMaker/interface/mTBound.h"
#include "TTree.h"
#include <algorithm>
using namespace std;


const  int CutN=21;



unsigned int RunMin = 9999999;
unsigned int RunMax = 0;
     
int tau_index=-1;
int tau_loose=-1;
int tau_tight=-1;
int mu_index=-1;
int el_index=-1;




   Int_t          primvert_count;
   Float_t         primvert_x;
   Float_t         primvert_y;
   Float_t         primvert_z;
   Float_t 	   SusyMother;
   Float_t 	   SusyLSP;


   Float_t 	   CFCounter_[30];
   Int_t	   muon_index;
   Int_t	   electron_index;
   Int_t	   taus_index;
   Int_t           mu_count;
   Int_t	   nbtag;
   Int_t	   njets;
   Int_t	   npv;
   Int_t 	   npu;
   Float_t         mu_px[10];   //[mu_count]
   Float_t         mu_py[10];   //[mu_count]
   Float_t         mu_pz[10];   //[mu_count]
   Float_t         mu_pt[10];   //[mu_count]
   Float_t         mu_eta[10];   //[mu_count]
   Float_t         mu_phi[10];   //[mu_count]
   Float_t         mu_charge[10];   //[mu_count]
   Float_t         mu_miniISO[10];   //[mu_count]
   Float_t         mu_dxy[10];   //[mu_count]
   Float_t         mu_dz[10];   //[mu_count]
   Float_t         mu_relIso[10];   //[mu_count]
 
   Float_t     mu_neutralHadIso[10]; 
   Float_t     mu_photonIso[10]; 
   Float_t     mu_chargedHadIso[10]; 
   Float_t     mu_puIso[10]; 
   Float_t     mu_neutralIso[10];
   Float_t     mu_absIsoMu[10]; 
   Float_t     mu_relIsoMu[10]; 

   Float_t     el_neutralHadIso[10]; 
   Float_t     el_photonIso[10]; 
   Float_t     el_chargedHadIso[10]; 
   Float_t     el_puIso[10]; 
   Float_t     el_neutralIso[10];
   Float_t     el_absIsoEl[10]; 
   Float_t     el_relIsoEl[10]; 







   Int_t           jet_count;
   Int_t           jets_cleaned[30];
   Float_t 	   jet_jecUn[30];
   Float_t         jet_e[30];   //[jet_count]
   Float_t         jet_px[30];   //[jet_count]
   Float_t         jet_py[30];   //[jet_count]
   Float_t         jet_pz[30];   //[jet_count]
   Float_t         jet_pt[30];   //[jet_count]
   Float_t         jet_eta[30];   //[jet_count]
   Float_t         jet_phi[30];   //[jet_count]
   Int_t           jet_flavour[30];   //[jet_count]
   Float_t         jet_btag[30];   //[jet_count]
   Int_t	   jet_isLoose[30];
   string	   datasetName;
   string	   regionName;
   Int_t           el_count;
   Float_t         el_px[10];   //[el_count]
   Float_t         el_py[10];   //[el_count]
   Float_t         el_pz[10];   //[el_count]
   Float_t         el_pt[10];   //[el_count]
   Float_t         el_eta[10];   //[el_count]
   Float_t         el_phi[10];   //[el_count]
   Float_t         el_miniISO[10];   //[el_count]
   Float_t         el_dxy[10];   //[el_count]
   Float_t         el_dz[10];   //[el_count]
   Float_t         el_charge[10];   //[el_count]
   Float_t         el_relIso[10];   //[el_count]


   Int_t           ta_count;
   Float_t         ta_px[30];   //[ta_count]
   Float_t         ta_py[30];   //[ta_count]
   Float_t         ta_pz[30];   //[ta_count]
   Float_t         ta_mass[30];   //[ta_count]
   Float_t         ta_eta[30];   //[ta_count]
   Float_t         ta_phi[30];   //[ta_count]
   Float_t         ta_pt[30];   //[ta_count]
   Float_t         ta_dxy[30];   //[ta_count]
   Float_t         ta_dz[30];   //[ta_count]
   Float_t         ta_charge[30];   //[ta_count]
   Float_t         ta_IsoFlag;   //[ta_count]
   Float_t         ta_relIso[10];   //[ta_count]
   Float_t         ta_puCorrPtSum[30];   //[ta_count]
   Float_t         ta_chargedIsoPtSum[30];   //[ta_count]
   Float_t         ta_neutralIsoPtSum[30];   //[ta_count]



   Float_t         met_ex;
   Float_t         met_ey;
   Float_t         met_ez;
   Float_t         met_ex_JetEnUp;
   Float_t         met_ey_JetEnUp;
   Float_t         met_ex_JetEnDown;
   Float_t         met_ey_JetEnDown;
   Float_t         met_ex_UnclusteredEnUp;
   Float_t         met_ey_UnclusteredEnUp;
   Float_t         met_ex_UnclusteredEnDown;
   Float_t         met_ey_UnclusteredEnDown;
   Float_t         met_pt;
   Float_t         met_phi;

   Float_t         gen_weight;
   Float_t 	   pu_weight;
   Float_t 	   LSF_weight;
   Float_t 	   TFR_weight;
   Float_t 	   top_weight;
   Float_t 	   all_weight;
   Float_t 	   trig_weight;
   Float_t 	   xsecs;
   Float_t 	   event_sign;
   Float_t 	   event_secondLeptonVeto;
   Float_t 	   met_flag;
   Float_t 	   event_thirdLeptonVeto;
   Float_t 	   event_leptonDrTrigger;
   Float_t	   genTauMatched;
   Float_t	   genLeptonMatched;
   Float_t	   qcdweight;
   Float_t	   qcdweightup;
   Float_t	   qcdweightdown;
   Int_t 	   npartons;



//bool isData = false;


double ChiMass=0;
double mIntermediate = tauMass;
double sumpT = 0;
double XSec=-1;
double xs,fact,fact2;
  
/*
  int nPtBins = 8;
  float ptBins[9] = {10,13,16,20,25,30,40,60,1000};

  int nPtBinsTrig = 16;
  float ptBinsTrig[17] = {10,
			  13,
			  16,
			  19,
			  22,
			  25,
			  28,
			  31,
			  34,
			  37,
			  40,
			  45,
			  50,
			  60,
			  70,
			  100,
			  1000};  
  
  int nEtaBins = 3;
  float etaBins[4] = {0,0.9,1.2,2.4}; 
  
  TString PtBins[8] = {"Pt10to13",
		       "Pt13to16",
		       "Pt16to20",
		       "Pt20to25",
		       "Pt25to30",
		       "Pt30to40",
		       "Pt40to60",
		       "PtGt60"};
  
  TString PtBinsTrig[16] = {"Pt10to13",
			    "Pt13to16",
			    "Pt16to19",
			    "Pt19to22",
			    "Pt22to25",
			    "Pt25to28",
			    "Pt28to31",
			    "Pt31to34",
			    "Pt34to37",
			    "Pt37to40",
			    "Pt40to45",
			    "Pt45to50",
			    "Pt50to60",
			    "Pt60to70",
			    "Pt70to100",
			    "PtGt100"};

  TString EtaBins[3] = {"EtaLt0p9",
			"Eta0p9to1p2",
			"EtaGt1p2"};
*/
float topPtWeight(float pt1,
		  float pt2) {

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter
 //   float a =  0.159; //l+jets
 //   float b =  -0.00141;
  
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
    
  if (pt1>400) w1 = 1.;
  if (pt2>400) w2 = 1.;
//cout<<" w1  "<<w1<<"  "<<w2<<endl;
  return TMath::Sqrt(w1*w2);

}

double TauFakeRateOld(float pt,float eta){

float SF = 1;

if (  fabs(eta) < 0.9 ) 
	{
		if (pt>20 && pt<30) SF = 1.12569;
		if (pt>30 && pt<50) SF = 1.17716;
		if (pt>50 && pt<60) SF = 1.1107;
		if (pt>60 )	    SF = 0.929284;
	}
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 ) 
	{

		if (pt>20 && pt<30) SF = 1.12176;
		if (pt>30 && pt<50) SF = 1.14511;
		if (pt>50 && pt<60) SF = 1.13235;
		if (pt>60 )	    SF = 0.803793;
	}

if (  fabs(eta) > 1.2 && fabs(eta) < 2.4 ) 
	{

		if (pt>20 && pt<30) SF = 1.19201;
		if (pt>30 && pt<50) SF = 1.36684;
		if (pt>50 && pt<60) SF = 0.822642;
		if (pt>60 )	    SF = 0.935916;
	}

return SF;


}



//string CutList[10];
vector<string> var;
vector < string > vec;
double var_[1000];

vector<string> CutList;

TH1D * histRuns = new TH1D("histRuns","",6000,24000,30000);

TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
TH1D * histTopPt = new TH1D("histTopPt","",1,-0.5,0.5);

TH1D * hWeights [CutN];

TH1D * hmu_dxy [CutN];
TH1D * hmu_dz [CutN];



TH1D * hel_dxy [CutN];
TH1D * hel_dz [CutN];

TH1D * htau_dxy [CutN];
TH1D * htau_dz [CutN];

TH1D *hMeffMuon[CutN];
TH1D *hMeffEl[CutN];
TH1D *hMeffTau[CutN];
TH1D *hHTOsqrMET[CutN];
TH1D *hPtOHT[CutN];
TH1D *hMeffMuonOsqrMET[CutN];
TH1D *hMeffElOsqrMET[CutN];
TH1D *hMeffTauOsqrMET[CutN];

TH1D *hHT[CutN];
TH1D *hHText[CutN];
TH1D *hRht[CutN];
TH1D *hHT2[CutN];
TH1D *hHT3[CutN];
TH1D *hHT4[CutN];
//TH1D *hST[CutN];
//TH1D *h0JetpT[CutN];
TH1D *hnJet[CutN];
TH1D *hnBJet[CutN];

TH1D *hCentrality[CutN];

TH1D *hPtJ0[CutN];
TH1D *hPtJ1[CutN];
TH1D *hPtJ2[CutN];
TH1D *hPtJ3[CutN];

TH1D *hInvMassMuTau[CutN];
TH1D *hInvMassMuEl[CutN];
TH1D *hInvMassTauTau[CutN];
TH1D *hInvMassElEl[CutN];
TH1D *hInvMassElTau[CutN];
TH1D *hInvMassMuMu[CutN];


TH1D *hnEl[CutN];
TH1D *hElpt[CutN];
TH1D *hEleta[CutN];
TH1D *hel_relISO[CutN];
TH1D *hel_relISOL[CutN];
TH1D *hel_miniISO[CutN];
TH1D *hel_miniISOL[CutN];

TH1D *hnLep[CutN];
TH1D *hLeppt[CutN];
TH1D *hLepeta[CutN];

TH1D *hnMu[CutN];
TH1D *hMupt[CutN];
TH1D *hMueta[CutN];

TH1D *hmu_relISO[CutN];
TH1D *hmu_relISOL[CutN];
TH1D *hmu_miniISO[CutN];
TH1D *hmu_miniISOL[CutN];

TH1D *htau_ISO[CutN];
TH1D *htau_ISOL[CutN];



TH1D *hnTau[CutN];
TH1D *hTaupt[CutN];
TH1D *hTaueta[CutN];


TH1D *hMET[CutN];
//TH1D *hnOver[CutN];
TH1D *hdPhiMETLep[CutN];
TH1D *hdPhiMETMu[CutN];
TH1D *hdPhiMETEl[CutN];
TH1D *hdPhiMETTau[CutN];
TH1D *hdPhiJMET[CutN];
TH1D *hdPhiJ0MET[CutN];
TH1D *hdPhiJ1MET[CutN];
TH1D *hdPhiJ2MET[CutN];
TH1D *hdPhiJ3MET[CutN];
TH1D *hdPhiMuMET[CutN];
TH1D *hdPhiElMET[CutN];
TH1D *hdPhiTauMET[CutN];

TH1D *hMT[CutN];
TH1D *hMTel[CutN];
TH1D *hMTmu[CutN];
TH1D *hMTtau[CutN];
TH1D *hDZeta[CutN];

TH1D *hMTtautau[CutN];


/////////// mutau
TH1D *hMTmutau[CutN];
TH1D *hMt2lestermutau[CutN];
TH1D *hMt2mutau[CutN];
TH1D *hMCTmutau[CutN];
TH1D *hMCTxmutau[CutN];
TH1D *hMCTymutau[CutN];
TH1D *hMCTbmutau[CutN];
TH1D *hTBoundmutau[CutN];
TH1D *hdR_mutau[CutN];

TH1D *hMTmuel[CutN];
TH1D *hMt2lestermuel[CutN];
TH1D *hMt2muel[CutN];
TH1D *hMCTmuel[CutN];
TH1D *hMCTxmuel[CutN];
TH1D *hMCTymuel[CutN];
TH1D *hMCTbmuel[CutN];
TH1D *hTBoundmuel[CutN];
TH1D *hdR_muel[CutN];

TH1D *hMTeltau[CutN];
TH1D *hMt2lestereltau[CutN];
TH1D *hMt2eltau[CutN];
TH1D *hMCTeltau[CutN];
TH1D *hMCTxeltau[CutN];
TH1D *hMCTyeltau[CutN];
TH1D *hMCTbeltau[CutN];
TH1D *hdR_eltau[CutN];
TH1D *hTBoundeltau[CutN];


TH1D *hdR_tautau[CutN];


TH1D *hnpv[CutN];
TH1D *hnpu[CutN];
TH1D *hnrho[CutN];


TH2D *hmet_MT[CutN];
TH2D *hmet_MTel[CutN];
TH2D *hmet_MTmu[CutN];
TH2D *hmet_MTtau[CutN];

TH2D *hdR_dPhi[CutN];
TH2D *hdRmt_dPhi[CutN];
TH2D *hdRet_dPhi[CutN];
TH2D *hdRme_dPhi[CutN];
TH2D *hdRtt_dPhi[CutN];

TH2D *hmet_dPhi[CutN];
TH2D *hmet_dPhiel[CutN];
TH2D *hmet_dPhimu[CutN];
TH2D *hmet_dPhitau[CutN];

TH2D *hMT_dPhi[CutN];
TH2D *hMT_dPhiel[CutN];
TH2D *hMT_dPhimu[CutN];
TH2D *hMT_dPhitau[CutN];


TH2D *hIso_sign[CutN];
  
TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutN,1,CutN+1);
TH1D *CutFlowUnW= new TH1D("CutFlowUnW","Cut Flow",CutN,1,CutN+1);
TH1D *CutFlowUnWNorm= new TH1D("CutFlowUnWNorm","Cut Flow",CutN,1,CutN+1);

TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);


TH1D * muonPtAllH = new TH1D("muonPtAllH","",10,0,200);
TH1D * electronPtAllH = new TH1D("electronPtAllH","",10,0,200);
TH1D * tauPtAllH = new TH1D("tauPtAllH","",10,0,200);

// histograms (dilepton selection)
TH1D * electronPtH  = new TH1D("electronPtH","",10,0,200);
TH1D * electronEtaH = new TH1D("electronEtaH","",50,-2.5,2.5); 
TH1D * muonPtH  = new TH1D("muonPtH","",10,0,200);
TH1D * muonEtaH = new TH1D("muonEtaH","",50,-2.5,2.5); 
TH1D * tauEtaAllH = new TH1D("tauEtaAllH","",50,-2.5,2.5); 
 
TH1D * dileptonMassH = new TH1D("dileptonMassH","",10,0,200);
TH1D * dileptonPtH = new TH1D("dileptonPtH","",10,0,200);
TH1D * dileptonEtaH = new TH1D("dileptonEtaH","",100,-5,5);
TH1D * dileptondRH = new TH1D("dileptondRH","",60,0,6);
TH1D * ETmissH = new TH1D("ETmissH","",10,0,200);
TH1D * MtH = new TH1D("MtH_2l","",10,0,200);
TH1D * DZetaH = new TH1D("DZetaH","",60,-400,200);

// histograms (dilepton selection + DZeta cut DZeta)
TH1D * electronPtSelH  = new TH1D("electronPtSelH","",10,0,200);
TH1D * electronEtaSelH = new TH1D("electronEtaSelH","",50,-2.5,2.5); 
TH1D * muonPtSelH  = new TH1D("muonPtSelH","",10,0,200);
TH1D * muonEtaSelH = new TH1D("muonEtaSelH","",50,-2.5,2.5); 

TH1D * dileptonMassSelH = new TH1D("dileptonMassSelH","",10,0,200);
TH1D * dileptonPtSelH = new TH1D("dileptonPtSelH","",10,0,200);
TH1D * dileptonEtaSelH = new TH1D("dileptonEtaSelH","",100,-5,5);
TH1D * dileptondRSelH = new TH1D("dileptondRSelH","",60,0,6);
TH1D * ETmissSelH = new TH1D("ETmissSelH","",10,0,200);
TH1D * MtSelH = new TH1D("MtSelH_2l","",10,0,200);
TH1D * DZetaSelH = new TH1D("DZetaSelH","",60,-400,200);

TLorentzVector ElV, MuV, TauV, JetsV, METV;

vector<TLorentzVector> AllJets_Lepton_noMet;
vector<TLorentzVector> JetsMV;
vector<TLorentzVector>  ElMV;
vector<TLorentzVector>  MuMV;
vector<TLorentzVector>  TauMV;
vector<TLorentzVector>  LeptMV;


std::vector<pair<string,float> > variables_;
std::vector<pair<string,float> > variablesMC_;

TTree *T;

void WriteHists(int CutNer, TFile *in, TString dir){

  in->cd(dir);

  for(int cj = 0; cj < CutNer; cj++)
    {

      hHTOsqrMET[cj]->Write();
    }
}


//void SetupTree(string & fout_,string &Sel){
void SetupTree(){
/*
TFile * filetree = new TFile(fout_.c_str(),"recreate");
	
filetree->mkdir(Sel.c_str());
filetree->cd(Sel.c_str());
*/
T  = new TTree("T","T");

  T->Branch("met_ex", &met_ex, "met_ex/F");
  T->Branch("met_ey", &met_ey, "met_ey/F");
  T->Branch("met_ez", &met_ez, "met_ez/F");

  T->Branch("met_ex_JetEnUp", &met_ex_JetEnUp, "met_ex_JetEnUp/F");
  T->Branch("met_ey_JetEnUp", &met_ey_JetEnUp, "met_ey_JetEnUp/F");

  T->Branch("met_ex_JetEnDown", &met_ex_JetEnDown, "met_ex_JetEnDown/F");
  T->Branch("met_ey_JetEnDown", &met_ey_JetEnDown, "met_ey_JetEnDown/F");

  T->Branch("met_ex_UnclusteredEnDown", &met_ex_UnclusteredEnDown, "met_ex_UnclusteredEnDown/F");
  T->Branch("met_ey_UnclusteredEnDown", &met_ey_UnclusteredEnDown, "met_ey_UnclusteredEnDown/F");

  T->Branch("met_ex_UnclusteredEnUp", &met_ex_UnclusteredEnUp, "met_ex_UnclusteredEnUp/F");
  T->Branch("met_ey_UnclusteredEnUp", &met_ey_UnclusteredEnUp, "met_ey_UnclusteredEnUp/F");
  T->Branch("met_pt", &met_pt, "met_pt/F");
  T->Branch("met_phi", &met_phi, "met_phi/F");
 
  T->Branch("gen_weight", &gen_weight, "gen_weight/F");
  T->Branch("pu_weight", &pu_weight, "pu_weight/F");
  T->Branch("LSF_weight", &LSF_weight, "LSF_weight/F");
  T->Branch("TFR_weight", &TFR_weight, "TFR_weight/F");
  T->Branch("top_weight", &top_weight, "top_weight/F");
  T->Branch("all_weight", &all_weight, "all_weight/F");
  T->Branch("trig_weight", &trig_weight, "trig_weight/F");

  T->Branch("xsecs", &xsecs, "xsecs/F");
  T->Branch("event_sign", &event_sign, "event_sign/F");
  T->Branch("met_flag", &met_flag, "met_flag/F");
  T->Branch("event_secondLeptonVeto", &event_secondLeptonVeto, "event_secondLeptonVeto/F");
  T->Branch("event_thirdLeptonVeto", &event_thirdLeptonVeto, "event_thirdLeptonVeto/F");
  T->Branch("event_leptonDrTrigger", &event_leptonDrTrigger, "event_leptonDrTrigger/F");

  T->Branch("muon_index", &muon_index, "muon_index/I");
  T->Branch("electron_index", &electron_index, "electron_index/I");
  T->Branch("taus_index", &taus_index, "taus_index/I");

  T->Branch("primvert_count", &primvert_count, "primvert_count/I");
  T->Branch("primvert_x", &primvert_x, "primvert_x/F");
  T->Branch("primvert_y", &primvert_x, "primvert_y/F");
  T->Branch("primvert_z", &primvert_x, "primvert_z/F");

  T->Branch("mu_count", &mu_count, "mu_count/I");
  T->Branch("mu_px", mu_px, "mu_px[10]/F");
  T->Branch("mu_py", mu_py, "mu_py[10]/F");
  T->Branch("mu_pz", mu_pz, "mu_pz[10]/F");
  T->Branch("mu_pt", mu_pt, "mu_pt[10]/F");
  T->Branch("mu_eta", mu_eta, "mu_eta[10]/F");
  T->Branch("mu_phi", mu_phi, "mu_phi[10]/F");
  T->Branch("mu_charge", mu_charge, "mu_charge[10]/F");
  T->Branch("mu_miniISO", mu_miniISO, "mu_miniISO[10]/F");
  T->Branch("mu_dxy", mu_dxy, "mu_dxy[10]/F");
  T->Branch("mu_dz", mu_dz, "mu_dz[10]/F");
  T->Branch("mu_relIso", mu_relIso, "mu_relIso[10]/F");
 
  T->Branch("mu_neutralHadIso", mu_neutralHadIso, "mu_neutralHadIso[10]/F");
  T->Branch("mu_photonIso", mu_photonIso, "mu_photonIso[10]/F");
  T->Branch("mu_chargedHadIso", mu_chargedHadIso, "mu_chargedHadIso[10]/F");
  T->Branch("mu_puIso", mu_puIso, "mu_puIso[10]/F");
  T->Branch("mu_neutralIso", mu_neutralIso, "mu_neutralIso[10]/F");
  T->Branch("mu_absIsoMu", mu_absIsoMu, "mu_absIsoMu[10]/F");
  T->Branch("mu_relIsoMu", mu_relIsoMu, "mu_relIsoMu[10]/F");

  T->Branch("el_neutralHadIso", el_neutralHadIso, "el_neutralHadIso[10]/F");
  T->Branch("el_photonIso", el_photonIso, "el_photonIso[10]/F");
  T->Branch("el_chargedHadIso", el_chargedHadIso, "el_chargedHadIso[10]/F");
  T->Branch("el_puIso", el_puIso, "el_puIso[10]/F");
  T->Branch("el_neutralIso", el_neutralIso, "el_neutralIso[10]/F");
  T->Branch("el_absIsoEl", el_absIsoEl, "el_absIsoEl[10]/F");
  T->Branch("el_relIsoEl", el_relIsoEl, "el_relIsoEl[10]/F");

  T->Branch("jet_count", &jet_count, "jet_count/I");
  T->Branch("njets", &njets, "njets/I");
  T->Branch("npv", &npv, "npv/I");
  T->Branch("npu", &npu, "npu/I");
  T->Branch("jets_cleaned", &jets_cleaned, "jets_cleaned[30]/I");
  T->Branch("jet_jecUn", jet_jecUn, "jet_jecUn[30]/F");
  T->Branch("jet_e", jet_e, "jet_e[30]/F");
  T->Branch("jet_px", jet_px, "jet_px[30]/F");
  T->Branch("jet_py", jet_py, "jet_py[30]/F");
  T->Branch("jet_pz", jet_pz, "jet_pz[30]/F");
  T->Branch("jet_pt", jet_pt, "jet_pt[30]/F");
  T->Branch("jet_eta", jet_eta, "jet_eta[30]/F");
  T->Branch("jet_phi", jet_phi, "jet_phi[30]/F");
  T->Branch("jet_flavour", jet_flavour, "jet_flavour[30]/F");
  T->Branch("jet_btag", jet_btag, "jet_btag[30]/F");
  T->Branch("jet_isLoose", jet_isLoose, "jet_isLoose[30]/I");
  
  T->Branch("CFCounter_", CFCounter_, "CFCounter_[30]/F");


  T->Branch("el_count", &el_count, "el_count/I");
  T->Branch("el_px", el_px, "el_px[10]/F");
  T->Branch("el_py", el_py, "el_py[10]/F");
  T->Branch("el_pz", el_pz, "el_pz[10]/F");
  T->Branch("el_pt", el_pt, "el_pt[10]/F");
  T->Branch("el_eta", el_eta, "el_eta[10]/F");
  T->Branch("el_phi", el_phi, "el_phi[10]/F");
  T->Branch("el_miniISO", el_miniISO, "el_miniISO[10]/F");
  T->Branch("el_dxy", el_dxy, "el_dxy[10]/F");
  T->Branch("el_dz", el_dz, "el_dz[10]/F");
  T->Branch("el_charge", el_charge, "el_charge[10]/F");
  T->Branch("el_relIso", el_relIso, "el_relIso[10]/F");


  T->Branch("ta_count", &ta_count, "ta_count/I");
  T->Branch("ta_px", ta_px, "ta_px[10]/F");
  T->Branch("ta_py", ta_py, "ta_py[10]/F");
  T->Branch("ta_pz", ta_pz, "ta_pz[10]/F");
  T->Branch("ta_mass", ta_mass, "ta_mass[10]/F");
  T->Branch("ta_eta", ta_eta, "ta_eta[10]/F");
  T->Branch("ta_phi", ta_phi, "ta_phi[10]/F");
  T->Branch("ta_pt", ta_pt, "ta_pt[10]/F");
  T->Branch("ta_dxy", ta_dxy, "ta_dxy[10]/F");
  T->Branch("ta_dz", ta_dz, "ta_dz[10]/F");
  T->Branch("ta_charge", ta_charge, "ta_charge[10]/F");
  T->Branch("ta_relIso", ta_relIso, "ta_relIso[10]/F");
  T->Branch("ta_IsoFlag", &ta_IsoFlag, "ta_IsoFlag/F");
  T->Branch("ta_chargedIsoPtSum", &ta_chargedIsoPtSum, "ta_chargedIsoPtSum/F");
  T->Branch("ta_neutralIsoPtSum", &ta_neutralIsoPtSum, "ta_neutralIsoPtSum/F");
  T->Branch("ta_puCorrPtSum", &ta_puCorrPtSum, "ta_puCorrPtSum/F");
  T->Branch("qcdweight", &qcdweight, "qcdweight/F");
  T->Branch("qcdweightup", &qcdweightup, "qcdweightup/F");
  T->Branch("qcdweightdown", &qcdweightdown, "qcdweightdown/F");
  T->Branch("nbtag", &nbtag, "nbtag/I");

  T->Branch("datasetName", &datasetName);
  T->Branch("regionName", &regionName);
  T->Branch("genTauMatched", &genTauMatched);
  T->Branch("genLeptonMatched", &genLeptonMatched);
  T->Branch("npartons",&npartons,"npartons/I");
  T->Branch("SusyMother",&SusyMother,"SusyMother/F");
  T->Branch("SusyLSP",&SusyLSP,"SusyLSP/F");



      char arg[100];
     for (unsigned int i = 0; i < vec.size (); i++)
    {
      var_[i] = -8888.;
      //string name = vec[i].c_str()+"_"+ds;
      sprintf (arg, "%s/F", vec[i].c_str ());
      T->Branch (vec[i].c_str (), &var_[i], arg);
      //T->Branch("pt_tt", &pt_tt, "pt_tt/F");
      //cout << " creating the TTree. " << vec[i].c_str()<< "  "<<i<<endl;
    }
     //T->Fill();
     //T->Write();
     
}

/*
void SetupHistsCut(int CutNer){


for(int cj = 0; cj < CutNer; cj++)
    {
      CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      CutFlowUnW->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
    }
}
*/
//string CutList[CutN];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNer){


for(int cj = 0; cj < CutNer; cj++)
    {
      CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      CutFlowUnW->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
    }
}
/*
 
 
      TString cutName=CutList[cj];
      TString nCut;
      nCut.Form("%i",cj);
      ///generic variables
      
      //T->Branch("hHT_"+nCut,"TH1D", &hHT[cj], 32000,0);
      
      hmu_dxy[cj] = new TH1D ("hmu_dxy_"+nCut,"hmu_dxy "+cutName,100,-5.,5.);
      
      hel_dxy[cj] = new TH1D ("hel_dxy_"+nCut,"hel_dxy "+cutName,100,-5.,5.);
      htau_dxy[cj] = new TH1D ("htau_dxy_"+nCut,"htau_dxy "+cutName,100,-5.,5.);
      
      htau_dz[cj] = new TH1D ("htau_dz_"+nCut,"htau_dz "+cutName,100,-50.,50.);
      hel_dz[cj] = new TH1D ("hel_dz_"+nCut,"hel_dz "+cutName,100,-50.,50.);
      hmu_dz[cj] = new TH1D ("hmu_dz_"+nCut,"hmu_dz "+cutName,100,-50.,50.);

      hHTOsqrMET[cj] = new TH1D ("hHTOsqrMET_"+nCut,"hHTOsqrMET "+cutName,60,0.0,1600.0);
      hHTOsqrMET[cj]->Sumw2();
      hPtOHT[cj] = new TH1D ("hPtOHT_"+nCut,"hPtOHT "+cutName,10,0.0,10.0);
      hPtOHT[cj]->Sumw2();
      
      hMeffMuonOsqrMET[cj] = new TH1D ("hMeffMuonOsqrMET_"+nCut,"hMeffMuonOsqrMET "+cutName,60,0.0,1600.0);
      hMeffMuonOsqrMET[cj]->Sumw2();
      hMeffElOsqrMET[cj] = new TH1D ("hMeffElOsqrMET_"+nCut,"hMeffElOsqrMET "+cutName,60,0.0,1600.0);
      hMeffElOsqrMET[cj]->Sumw2();
      hMeffTauOsqrMET[cj] = new TH1D ("hMeffTauOsqrMET_"+nCut,"hMeffTauOsqrMET "+cutName,60,0.0,1600.0);
      hMeffTauOsqrMET[cj]->Sumw2();

      hMeffMuon[cj] = new TH1D ("hMeffMuon_"+nCut,"hMeffMuon "+cutName,60,0.0,1600.0);
      hMeffMuon[cj]->Sumw2();
      hMeffEl[cj] = new TH1D ("hMeffEl_"+nCut,"hMeffEl "+cutName,60,0.0,1600.0); 
      hMeffEl[cj]->Sumw2();
      hMeffTau[cj] = new TH1D ("hMeffTau_"+nCut,"hMeffTau "+cutName,60,0.0,1600.0); 
      hMeffTau[cj]->Sumw2();

      hCentrality[cj]  = new TH1D ("hCentrality_"+nCut,"hCentrality "+cutName,20,0.,2.);
      hCentrality[cj]->Sumw2();

      hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,60,0.0,1600.0);
      hHT[cj]->Sumw2();
 
      hRht[cj] = new TH1D ("Rht_"+nCut,"Rht "+cutName,10,0.0,2.0);
      hRht[cj]->Sumw2();
    
      hHText[cj] = new TH1D ("HText_"+nCut,"HText "+cutName,60,0.0,1600.0);
      hHText[cj]->Sumw2();
      hHT2[cj] = new TH1D ("HT2_"+nCut,"HT2 "+cutName,60,0.0,1600.0);
      hHT2[cj]->Sumw2();
      hHT3[cj] = new TH1D ("HT3_"+nCut,"HT3 "+cutName,60,0.0,1600.0);
      hHT3[cj]->Sumw2();
      hHT4[cj] = new TH1D ("HT4_"+nCut,"HT4 "+cutName,60,0.0,1600.0);
      hHT4[cj]->Sumw2();
          
      hPtJ0[cj] = new TH1D ("hPtJ0_"+nCut,"hPtJ0 "+cutName,60,0.0,1600.0);
      hPtJ0[cj]->Sumw2();
      hPtJ1[cj] = new TH1D ("hPtJ1_"+nCut,"hPtJ1 "+cutName,60,0.0,1600.0);
      hPtJ1[cj]->Sumw2();
      hPtJ2[cj] = new TH1D ("hPtJ2_"+nCut,"hPtJ2 "+cutName,60,0.0,1600.0);
      hPtJ2[cj]->Sumw2();
      hPtJ3[cj] = new TH1D ("hPtJ3_"+nCut,"hPtJ3 "+cutName,60,0.0,1600.0);
      hPtJ3[cj]->Sumw2();
     

      //h0JetpT[cj] = new TH1D ("0JetpT_"+nCut,"0JetpT "+cutName,60,0.0,1600.0);
      //h0JetpT[cj]->Sumw2();
      hnJet[cj] = new TH1D ("nJet_"+nCut,"nJet "+cutName,25,-0.5,24.5);
      hnJet[cj]->Sumw2();
      hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,10,-0.5,9.5);
      hnBJet[cj]->Sumw2();

      hWeights[cj] = new TH1D ("hWeights_"+nCut,"hWeights "+cutName,10,-1,9);
      hWeights[cj]->Sumw2();
	
      hDZeta[cj] = new TH1D("hDZeta_"+nCut,"hDZeta"+cutName,60,-400,200);
      hDZeta[cj]->Sumw2();
      hInvMassMuTau[cj] = new TH1D ("hInvMassMuTau_"+nCut,"hInvMassMuTau "+cutName,25,0,500);
      hInvMassMuEl[cj] = new TH1D ("hInvMassMuEl_"+nCut,"hInvMassMuel "+cutName,25,0,500);
      hInvMassMuMu[cj] = new TH1D ("hInvMassMuMu_"+nCut,"hInvMassMuMu "+cutName,25,0,500);
      hInvMassElTau[cj] = new TH1D ("hInvMassElTau_"+nCut,"hInvMassElTau "+cutName,25,0,500);
      hInvMassElEl[cj] = new TH1D ("hInvMassElEl_"+nCut,"hInvMassElEl "+cutName,25,0,500);
      hInvMassTauTau[cj] = new TH1D ("hInvMassTauTau_"+nCut,"hInvMassTauTau "+cutName,25,0,500);
        
      //Leptons
      //
      //
      hnLep[cj] = new TH1D ("nLep_"+nCut,"nLep "+cutName,10,-0.5,9.5);
      hnLep[cj]->Sumw2();
      hLeppt[cj] = new TH1D ("LeppT_"+nCut,"Lep pT "+cutName,25,0,500);
      hLeppt[cj]->Sumw2();
      hLepeta[cj] = new TH1D ("Lepeta_"+nCut,"Lep eta "+cutName,40,-4,4);
      hLepeta[cj]->Sumw2();
      //hST[cj] = new TH1D ("ST_"+nCut,"ST "+cutName,60,0.0,1600.0);
      //hST[cj]->Sumw2();
        
      //Muons
      //
      //
      hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,-0.5,9.5);
      hnMu[cj]->Sumw2();
      hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,25,0,500);
      hMupt[cj]->Sumw2();
      hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,40,-4,4);
      hMueta[cj]->Sumw2();
        
      //Taus
      //
      //
      hnTau[cj] = new TH1D ("nTau_"+nCut,"nTau "+cutName,10,-0.5,9.5);
      hnTau[cj]->Sumw2();
      hTaupt[cj] = new TH1D ("TaupT_"+nCut,"Tau pT "+cutName,25,0,500);
      hTaupt[cj]->Sumw2();
      hTaueta[cj] = new TH1D ("Taueta_"+nCut,"Tau eta "+cutName,40,-4,4);
      hTaueta[cj]->Sumw2();
	
      //hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
      //Electrons
      //
      //
      hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,-0.5,9.5);
      hnEl[cj]->Sumw2();
      hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,25,0,500);
      hElpt[cj]->Sumw2();
      hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,40,-4,4);
      hEleta[cj]->Sumw2();
       
       
      hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,40.0,0.0,800.0);
      hMET[cj]->Sumw2();
	

      //dPhi
      //
      //
      hdPhiMETLep[cj] = new TH1D("dPhiMETLep_"+nCut,"dPhiMETLep "+cutName,64,0.0,3.2);
      hdPhiMETLep[cj]->Sumw2();
    
      hdPhiMETMu[cj] = new TH1D("dPhiMETMu_"+nCut,"dPhiMETMu "+cutName,64,0.0,3.2);
      hdPhiMETMu[cj]->Sumw2();
      hdPhiMETEl[cj] = new TH1D("dPhiMETEl_"+nCut,"dPhiMETEl "+cutName,64,0.0,3.2);
      hdPhiMETEl[cj]->Sumw2();
      hdPhiMETTau[cj] = new TH1D("dPhiMETTau_"+nCut,"dPhiMETTau "+cutName,64,0.0,3.2);
      hdPhiMETTau[cj]->Sumw2();
    
      hdPhiJMET[cj] = new TH1D("dPhiJMET_"+nCut,"dPhiJMET "+cutName,64,0.0,3.2);
      hdPhiJMET[cj]->Sumw2();
      hdPhiJ0MET[cj] = new TH1D("dPhiJ0MET_"+nCut,"dPhiJ0MET "+cutName,64,0.0,3.2);
      hdPhiJ0MET[cj]->Sumw2();
      hdPhiJ1MET[cj] = new TH1D("dPhiJ1MET_"+nCut,"dPhiJ1MET "+cutName,64,0.0,3.2);
      hdPhiJ1MET[cj]->Sumw2();
      hdPhiJ2MET[cj] = new TH1D("dPhiJ2MET_"+nCut,"dPhiJ2MET "+cutName,64,0.0,3.2);
      hdPhiJ2MET[cj]->Sumw2();
      hdPhiJ3MET[cj] = new TH1D("dPhiJ3MET_"+nCut,"dPhiJ3MET "+cutName,64,0.0,3.2);
      hdPhiJ3MET[cj]->Sumw2();

      hdPhiMuMET[cj] = new TH1D("dPhiMuMET_"+nCut,"dPhiMuMET "+cutName,64,0.0,3.2);
      hdPhiMuMET[cj]->Sumw2();

      hdPhiElMET[cj] = new TH1D("dPhiElMET_"+nCut,"dPhiElMET "+cutName,64,0.0,3.2);
      hdPhiElMET[cj]->Sumw2();
   
      hdPhiTauMET[cj] = new TH1D("dPhiTauMET_"+nCut,"dPhiTauMET "+cutName,64,0.0,3.2);
      hdPhiTauMET[cj]->Sumw2();
   
      //MT
      //
      //
      hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,10,0,200);
      hMT[cj]->Sumw2();
      hMTel[cj] = new TH1D ("MTel_"+nCut,"MTel "+cutName,10,0,200);
      hMTel[cj]->Sumw2();
      hMTmu[cj] = new TH1D ("MTmu_"+nCut,"MTmu "+cutName,10,0,200);
      hMTmu[cj]->Sumw2();

      hMTtau[cj] = new TH1D ("MTtau_"+nCut,"MTtau "+cutName,10,0,200);
      hMTtau[cj]->Sumw2();

      hMTmutau[cj] = new TH1D ("MTmutau_"+nCut,"MTmutau "+cutName,10,0,200);
      hMTmutau[cj]->Sumw2();
      hMt2mutau[cj] = new TH1D ("Mt2mutau_"+nCut,"Mt2mutau "+cutName,50,0,1000);
      hMt2mutau[cj]->Sumw2();
      hMt2lestermutau[cj] = new TH1D ("Mt2lestermutau_"+nCut,"Mt2lestermutau "+cutName,30,0,600);
      hMt2lestermutau[cj]->Sumw2();
      hMCTmutau[cj] = new TH1D ("MCTmutau_"+nCut,"MCTmutau "+cutName,25,0,500);
      hMCTmutau[cj]->Sumw2();
      hMCTxmutau[cj] = new TH1D ("MCTxmutau_"+nCut,"MCTxmutau "+cutName,25,0,500);
      hMCTxmutau[cj]->Sumw2();
      hMCTymutau[cj] = new TH1D ("MCTymutau_"+nCut,"MCTymutau "+cutName,25,0,500);
      hMCTymutau[cj]->Sumw2();
      hMCTbmutau[cj] = new TH1D ("MCTbmutau_"+nCut,"MCTbmutau "+cutName,25,0,500);
      hMCTbmutau[cj]->Sumw2();
      hTBoundmutau[cj] = new TH1D ("TBoundmutau_"+nCut,"TBoundmutau "+cutName,100,0,2000);
      hTBoundmutau[cj] ->Sumw2();
      hdR_mutau[cj]= new TH1D ("dR_mutau_"+nCut,"dR_mutau "+cutName,60,0,6);;
      hdR_mutau[cj]->Sumw2();


      hMTeltau[cj] = new TH1D ("MTeltau_"+nCut,"MTeltau "+cutName,10,0,200);
      hMTeltau[cj]->Sumw2();
      hMt2eltau[cj] = new TH1D ("Mt2eltau_"+nCut,"Mt2eltau "+cutName,50,0,1000);
      hMt2eltau[cj]->Sumw2();
      hMt2lestereltau[cj] = new TH1D ("Mt2lestereltau_"+nCut,"Mt2lestereltau "+cutName,30,0,600);
      hMt2lestereltau[cj]->Sumw2();
      hMCTeltau[cj] = new TH1D ("MCTeltau_"+nCut,"MCTeltau "+cutName,25,0,500);
      hMCTeltau[cj]->Sumw2();
      hMCTxeltau[cj] = new TH1D ("MCTxeltau_"+nCut,"MCTxeltau "+cutName,25,0,500);
      hMCTxeltau[cj]->Sumw2();
      hMCTyeltau[cj] = new TH1D ("MCTyeltau_"+nCut,"MCTyeltau "+cutName,25,0,500);
      hMCTyeltau[cj]->Sumw2();
      hMCTbeltau[cj] = new TH1D ("MCTbeltau_"+nCut,"MCTbeltau "+cutName,25,0,500);
      hMCTbeltau[cj]->Sumw2();
      hTBoundeltau[cj] = new TH1D ("TBoundeltau_"+nCut,"TBoundeltau "+cutName,100,0,2000);
      hTBoundeltau[cj] ->Sumw2();
      hdR_eltau[cj]= new TH1D ("dR_eltau_"+nCut,"dR_eltau "+cutName,60,0,6);;
      hdR_eltau[cj]->Sumw2();


      hMTmuel[cj] = new TH1D ("MTmuel_"+nCut,"MTmuel "+cutName,10,0,200);
      hMTmuel[cj]->Sumw2();
      hMt2muel[cj] = new TH1D ("Mt2muel_"+nCut,"Mt2muel "+cutName,50,0,1000);
      hMt2muel[cj]->Sumw2();
      hMt2lestermuel[cj] = new TH1D ("Mt2lestermuel_"+nCut,"Mt2lestermuel "+cutName,30,0,600);
      hMt2lestermuel[cj]->Sumw2();
      hMCTmuel[cj] = new TH1D ("MCTmuel_"+nCut,"MCTmuel "+cutName,25,0,500);
      hMCTmuel[cj]->Sumw2();
      hMCTxmuel[cj] = new TH1D ("MCTxmuel_"+nCut,"MCTxmuel "+cutName,25,0,500);
      hMCTxmuel[cj]->Sumw2();
      hMCTymuel[cj] = new TH1D ("MCTymuel_"+nCut,"MCTymuel "+cutName,25,0,500);
      hMCTymuel[cj]->Sumw2();
      hMCTbmuel[cj] = new TH1D ("MCTbmuel_"+nCut,"MCTbmuel "+cutName,25,0,500);
      hMCTbmuel[cj]->Sumw2();
      hTBoundmuel[cj] = new TH1D ("TBoundmuel_"+nCut,"TBoundmuel "+cutName,100,0,2000);
      hTBoundmuel[cj] ->Sumw2();
      hdR_muel[cj]= new TH1D ("dR_muel_"+nCut,"dR_muel "+cutName,60,0,6);;
      hdR_muel[cj]->Sumw2();




      htau_ISOL[cj]= new TH1D ("tauISOL_"+nCut,"tauISOL "+cutName,50,0,5);;
      htau_ISO[cj]= new TH1D ("tauISO_"+nCut,"tauISO "+cutName,50,0,5);;
      htau_ISO[cj] ->Sumw2(); 
      htau_ISOL[cj] ->Sumw2(); 

      hel_miniISO[cj]= new TH1D ("elminiISO_"+nCut,"elminiISO "+cutName,50,0,5);;
      hel_miniISO[cj]->Sumw2();
      hel_miniISOL[cj]= new TH1D ("elminiISOL_"+nCut,"elminiISOL "+cutName,50,0,5);;
      hel_miniISOL[cj]->Sumw2();

      hel_relISO[cj]= new TH1D ("elrelISO_"+nCut,"elrelISO "+cutName,50,0,5);;
      hel_relISO[cj]->Sumw2();
      hel_relISOL[cj]= new TH1D ("elrelISOL_"+nCut,"elrelISOL "+cutName,50,0,5);;
      hel_relISOL[cj]->Sumw2();
        
        
      hmu_miniISO[cj]= new TH1D ("muminiISO_"+nCut,"muminiISO "+cutName,50,0,5);;
      hmu_miniISO[cj]->Sumw2();
      hmu_miniISOL[cj]= new TH1D ("muminiISOL_"+nCut,"muminiISOL "+cutName,50,0,5);;
      hmu_miniISOL[cj]->Sumw2();

      hmu_relISO[cj]= new TH1D ("murelISO_"+nCut,"murelISO "+cutName,50,0,5);;
      hmu_relISO[cj]->Sumw2();
      hmu_relISOL[cj]= new TH1D ("murelISOL_"+nCut,"murelISOL "+cutName,50,0,5);;
      hmu_relISOL[cj]->Sumw2();
 
        

      hdR_tautau[cj]= new TH1D ("dR_tautau_"+nCut,"dR_tautau "+cutName,60,0,6);;
      hdR_tautau[cj]->Sumw2();
	
      

      hnpv[cj]= new TH1D ("npv_"+nCut,"npv "+cutName,100,-0.5,99.5);;
      hnpv[cj]->Sumw2();
      hnpu[cj]= new TH1D ("npu_"+nCut,"npu "+cutName,100,-0.5,99.5);;
      hnpu[cj]->Sumw2();
      hnrho[cj]= new TH1D ("nrho_"+nCut,"nrho "+cutName,100,-0.5,99.5);;
      hnrho[cj]->Sumw2();
 
      
      hdRmt_dPhi[cj] = new TH2D ("dRm_dPhi_"+nCut,"dRm_dPhi "+cutName,60,0.0,6.0,64,0.0,3.2);
      hdRmt_dPhi[cj]->Sumw2();
      hdRet_dPhi[cj] = new TH2D ("dRe_dPhi_"+nCut,"dRe_dPhi "+cutName,60.0,0.0,6.0,64,0.0,3.2);
      hdRet_dPhi[cj]->Sumw2();
      hdRme_dPhi[cj] = new TH2D ("dRme_dPhi_"+nCut,"dRme_dPhi "+cutName,60.0,0.0,6.0,64,0.0,3.2);
      hdRme_dPhi[cj]->Sumw2();
      hdRtt_dPhi[cj] = new TH2D ("dRtt_dPhi_"+nCut,"dRtt_dPhi "+cutName,60.0,0.0,6.0,64,0.0,3.2);
      hdRtt_dPhi[cj]->Sumw2();
 
      hmet_dPhi[cj] = new TH2D ("met_dPhi_"+nCut,"met_dPhi "+cutName,40.0,0.0,800.0,64,0.0,3.2);
      hmet_dPhi[cj]->Sumw2();
      hmet_MT[cj] = new TH2D ("met_MT_"+nCut,"met_MT "+cutName,40.0,0.0,800.0,10,0,200);
      hmet_MT[cj]->Sumw2();

      hmet_dPhiel[cj] = new TH2D ("met_dPhiel_"+nCut,"met_dPhiel "+cutName,40.0,0.0,800.0,64,0.0,3.2);
      hmet_dPhiel[cj]->Sumw2();
      hmet_MTel[cj] = new TH2D ("met_MTel_"+nCut,"met_MTel "+cutName,40.0,0.0,800.0,10,0,200);
      hmet_MTel[cj]->Sumw2();
 

      hmet_dPhimu[cj] = new TH2D ("met_dPhimu_"+nCut,"met_dPhimu "+cutName,40.0,0.0,800.0,64,0.0,3.2);
      hmet_dPhimu[cj]->Sumw2();
      hmet_MTmu[cj] = new TH2D ("met_MTmu_"+nCut,"met_MTmu "+cutName,40.0,0.0,800.0,10,0,200);
      hmet_MTmu[cj]->Sumw2();

      hmet_dPhitau[cj] = new TH2D ("met_dPhitau_"+nCut,"met_dPhitau "+cutName,40.0,0.0,800.0,64,0.0,3.2);
      hmet_dPhitau[cj]->Sumw2();
      hmet_MTtau[cj] = new TH2D ("met_MTtau_"+nCut,"met_MTtau "+cutName,40.0,0.0,800.0,10,0,200);
      hmet_MTtau[cj]->Sumw2();


 
      hMT_dPhi[cj]= new TH2D ("MT_dPhi_"+nCut,"MT_dPhi "+cutName,10,0,200,64,0.0,3.2);
      hMT_dPhi[cj]->Sumw2();
      hMT_dPhiel[cj]= new TH2D ("MTel_dPhi_"+nCut,"MTel_dPhi "+cutName,10,0,200,64,0.0,3.2);
      hMT_dPhiel[cj]->Sumw2();
      hMT_dPhimu[cj]= new TH2D ("MTmu_dPhi_"+nCut,"MTmu_dPhi "+cutName,10,0,200,64,0.0,3.2);
      hMT_dPhimu[cj]->Sumw2();
      hMT_dPhitau[cj]= new TH2D ("MTtau_dPhi_"+nCut,"MTtau_dPhi "+cutName,10,0,200,64,0.0,3.2);
      hMT_dPhitau[cj]->Sumw2();
      hIso_sign[cj]= new TH2D ("Iso_sign_"+nCut,"Iso_sign "+cutName,50,0,5,2,0,2);
      hIso_sign[cj]->GetYaxis()->SetBinLabel(1,"SS");
      hIso_sign[cj]->GetYaxis()->SetBinLabel(2,"OS");
      hIso_sign[cj]->Sumw2();






    }

}
*/
void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV,vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_){}


void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel){}

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex){};


double Centrality (vector < TLorentzVector > AllJets_Lepton_noMet_){
    // Centrality
    double Centrality = 0, Centrality_N = 0, Centrality_D = 0;
    for (unsigned int i = 0; i < AllJets_Lepton_noMet_.size (); i++)
    {
    Centrality_N += AllJets_Lepton_noMet_[i].Pt ();
    Centrality_D += AllJets_Lepton_noMet_[i].P ();
    }
    double Centr_=-1;

    Centrality_D > 0 ? Centr_ = Centrality_N / Centrality_D : Centr_=-1;

    return Centr_ ;
  }


void WriteTree() 

  {
	  T->Print();
	  //T->Write();
	  T->AutoSave();
  }


void FillTree() {

	  T->Fill();
}

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, double Chimass, double mintermediate,AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex){
  }

