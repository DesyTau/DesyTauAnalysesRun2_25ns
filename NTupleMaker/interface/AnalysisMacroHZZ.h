
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
int tau_index1=-1;
int tau_index2=-1;
int tau_loose=-1;
int tau_tight=-1;
int tau_loose2=-1;
int tau_tight2=-1;
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
   Int_t	   muon_index_1;
   Int_t	   muon_index_2;
   Int_t	   electron_index;
   Int_t	   taus_index;
   Int_t	   taus_index2;
   Int_t           mu_count;
   Int_t	   nbtag;
   Int_t	   njets;
   Int_t	   npv;
   Int_t 	   npu;
   int 	   	   genTauDecayMode;
   int 	   	   event_lumi;
   int 	   	   event_run;
   int 	   	   genTauDecayMode1;
   int 	   	   genTauDecayMode2;
   int 	   	   genTauDecayModeAll;
   Float_t         mu_px[20];   //[mu_count]
   Float_t         mu_py[20];   //[mu_count]
   Float_t         mu_pz[20];   //[mu_count]
   Float_t         mu_pt[20];   //[mu_count]
   Float_t         mu_eta[20];   //[mu_count]
   Float_t         mu_phi[20];   //[mu_count]
   Float_t         mu_charge[20];   //[mu_count]
   Float_t         mu_miniISO[20];   //[mu_count]
   Float_t         mu_dxy[20];   //[mu_count]
   Float_t         mu_dxyerr[20];   //[mu_count]
   Float_t         mu_dz[20];   //[mu_count]
   Float_t         mu_dzerr[20];   //[mu_count]
   Float_t         mu_relIso[20];   //[mu_count]
 
   Float_t     mu_neutralHadIso[20]; 
   Float_t     mu_photonIso[20]; 
   Float_t     mu_chargedHadIso[20]; 
   Float_t     mu_puIso[20]; 
   Float_t     mu_neutralIso[20];
   Float_t     mu_absIsoMu[20]; 
   Float_t     mu_relIsoMu[20]; 

   Float_t     el_neutralHadIso[20]; 
   Float_t     el_photonIso[20]; 
   Float_t     el_chargedHadIso[20]; 
   Float_t     el_puIso[20]; 
   Float_t     el_neutralIso[20];
   Float_t     el_absIsoEl[20]; 
   Float_t     el_relIsoEl[20]; 

   Float_t         wScale0;
   Float_t         wScale1;
   Float_t         wScale2;
   Float_t         wScale3;
   Float_t         wScale4;
   Float_t         wScale5;
   Float_t         wScale6;
   Float_t         wScale7;
   Float_t         wScale8;
   Float_t         wPDFUp;
   Float_t         wPDFDown;
   Float_t         wPDFdev;






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
   Float_t         el_px[20];   //[el_count]
   Float_t         el_py[20];   //[el_count]
   Float_t         el_pz[20];   //[el_count]
   Float_t         el_pt[20];   //[el_count]
   Float_t         el_eta[20];   //[el_count]
   Float_t         el_phi[20];   //[el_count]
   Float_t         el_miniISO[20];   //[el_count]
   Float_t         el_dxy[20];   //[el_count]
   Float_t         el_dxyerr[20];   //[el_count]
   Float_t         el_dz[20];   //[el_count]
   Float_t         el_dzerr[20];   //[el_count]
   Float_t         el_charge[20];   //[el_count]
   Float_t         el_relIso[20];   //[el_count]
   Float_t         el_isMVA[20];   //[el_count]
   Float_t         el_isnotrig_MVA80[20];   //[el_count]
   Float_t         el_isnotrig_MVA90[20];   //[el_count]
   Float_t         el_istrig_MVA80[20];   //[el_count]
   Float_t         el_istrig_MVA90[20];   //[el_count]


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
   Float_t         ta_relIso[20];   //[ta_count]
   Float_t         ta_puCorrPtSum[30];   //[ta_count]
   Float_t         ta_chargedIsoPtSum[30];   //[ta_count]
   Float_t         ta_neutralIsoPtSum[30];   //[ta_count]

   Float_t 	   ta_IsoFlagVTight[5];
   Float_t 	   ta_IsoFlagTight[5];
   Float_t 	   ta_IsoFlagLoose[5];
   Float_t 	   ta_IsoFlagVLoose[5];
   Float_t 	   ta_IsoFlagMedium[5];
   Float_t 	   isElTau;
   Float_t 	   isMuTau;
   Float_t 	   isTauTau;

   Float_t 	   ta_isLoose;
   Float_t 	   ta_isTight;


   Float_t         genmet;
   Float_t         genmet_Ex;
   Float_t         genmet_Ey;
   Float_t         genmetphi;
   Float_t         met_scaleUp;
   Float_t         metphi_scaleUp;
   Float_t         met_scaleDown;
   Float_t         metphi_scaleDown;
   Float_t         met_resoUp;
   Float_t         met_resoDown;
   Float_t         metphi_resoUp;
   Float_t         metphi_resoDown;
   Float_t         NuPx;
   Float_t         NuPy;
   Float_t         NuPz;
   Float_t         NuPt;
   Float_t         NuPhi;
   Float_t         PtSystem;

   Float_t         met_ex;
   Float_t         met_ey;
   Float_t         met_ex_recoil;
   Float_t         met_ey_recoil;
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
   Float_t 	   LSF_weight_mu;
   Float_t 	   LSF_weight_el;
   Float_t 	   LSF_weight_1;
   Float_t 	   LSF_weight_2;
   Float_t 	   TFR_weight;
   Float_t 	   top_weight;
   Float_t 	   all_weight;
   Float_t 	   trig_weight;
   Float_t 	   trig_weight_1;
   Float_t 	   trig_weight_2;
   Float_t 	   zptmassweight;
   Float_t 	   xsecs;
   Float_t 	   event_sign;
   Float_t 	   event_secondLeptonVeto;
   Float_t 	   eleMVA;
   Float_t 	   met_flag;
   Float_t 	   event_thirdLeptonVeto;
   Float_t 	   event_leptonDrTrigger;
   Float_t	   genTauMatched;
   Float_t	   genTauMatched2;
   Float_t         genLeptonMatchedPrompEl;
   Float_t         genLeptonMatchedPrompMu;
   Float_t         genElMatchedToTauDecay;
   Float_t         genMuMatchedToTauDecay;
   Float_t         genTauMatchedToTauDecay;
   Float_t	   genLeptonMatchedPromptEl;
   Float_t	   genLeptonMatchedPromptMu;
   Float_t	   genLeptonMatchedPromptTau;
   Float_t	   genElMatchedHadrDecay;
   Float_t	   genMuMatchedHadrDecay;
   Float_t	   genTauMatchedHadrDecay;
   Float_t	   genLeptonMatchedGluon;
   Float_t	   genLeptonMatchedHFQ;
   Float_t	   genLeptonMatchedLFQ;

   Float_t         matchedTauToPromptEl;
   Float_t         matchedTauToPromptMu;
   Float_t         matchedTauToPromptTau;
   Float_t         matchedTauToTauDecEl;
   Float_t         matchedTauToTauDecMu;
   Float_t         matchedTauToTauDecTau;

   Float_t	   matchedTauToElHadronDec;
   Float_t	   matchedTauToMuHadronDec;
   Float_t	   matchedTauToTauHadronDec;
   Float_t	   matchedTauToGluon;
   Float_t	   matchedTauToNothing;
   Float_t	   matchedTauToHFQ;
   Float_t	   matchedTauToLFQ;
   Float_t	   matchedTauToUpQ;
   Float_t	   matchedTauToDownQ;
   Float_t	   matchedTauToStrangeQ;
   Float_t	   matchedTauToCharmQ;
   Float_t	   matchedTauToBottomQ;
   Float_t	   matchedTauToBottom;

   Float_t	   isDYTT;
   Float_t	   isDYLL;
   Float_t	   isDYEE;
   Float_t	   isDYMM;
   Float_t	   isDYNuNu;



   Float_t	   genLeptonPromptElMatched;
   Float_t	   genLeptonPromptMuMatched;
   Float_t	   genLeptonMatched;
   Float_t	   genTauDecayedElMatched;
   Float_t	   genTauDecayedMuMatched;
   Float_t	   genLeptonMatchedEl;
   Float_t	   genLeptonMatchedMu;
   Float_t	   genLeptonMatched2;
   Float_t	   qcdweight;
   Float_t	   qcdweightup;
   Float_t	   qcdweightdown;
   Int_t 	   npartons;


   Float_t         met_ex_JetEnUp_recoil;
   Float_t         met_ey_JetEnUp_recoil;
   Float_t         met_ex_JetEnDown_recoil;
   Float_t         met_ey_JetEnDown_recoil;
   Float_t         met_ex_UnclusteredEnUp_recoil;
   Float_t         met_ey_UnclusteredEnUp_recoil;
   Float_t         met_ex_UnclusteredEnDown_recoil;
   Float_t         met_ey_UnclusteredEnDown_recoil;





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

  float a = 0.0615;    // Run1 a parameter
  float b = -0.0005;  // Run1 b parameter
 //   float a =  0.159; //l+jets
 //   float b =  -0.00141;
  
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
    
  if (pt1>400) w1 = 1.;
  if (pt2>400) w2 = 1.;
//cout<<" w1  "<<w1<<"  "<<w2<<endl;
  return TMath::Sqrt(w1*w2);

}


//string CutList[20];
vector<string> var;
vector < string > vec;
double var_[2000];

vector<string> CutList;

TH1D * histRuns = new TH1D("histRuns","",6000,24000,30000);

int nBinsPt = 8;
double binsPt[9] ={0,50,100,150,200,300,400,600,1000};
TH1D * histPt = new TH1D("histPt","histPt",nBinsPt,binsPt);

TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
TH1D * histWeightsScalesUp = new TH1D("histWeightsScalesUp","",1,-0.5,0.5);
TH1D * histWeightsScalesDown = new TH1D("histWeightsScalesDown","",1,-0.5,0.5);
TH1D * histWeightsPDFUp = new TH1D("histWeightsPDFUp","",1,-0.5,0.5);
TH1D * histWeightsPDFDown = new TH1D("histWeightsPDFDown","",1,-0.5,0.5);
TH1D * histTopPt = new TH1D("histTopPt","",1,-0.5,0.5);
TH1D * histTopPtSq = new TH1D("histTopPtSq","",1,-0.5,0.5);

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

TLorentzVector ElV, MuV, TauV, JetsV, METV, LeptV1, LeptV2;

vector<TLorentzVector> AllJets_Lepton_noMet;
vector<TLorentzVector> JetsMV;
vector<TLorentzVector>  ElMV;
vector<TLorentzVector>  MuMV;
vector<TLorentzVector>  TauMV;
vector<TLorentzVector>  LeptMV;
vector<TLorentzVector>  LeptMV1;
vector<TLorentzVector>  LeptMV2;


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
/*
  T->Branch("genmet", &genmet, "genmet/F");
  T->Branch("genmetphi", &genmetphi, "genmetphi/F");


   Float_t         genmet;
   Float_t         genmetphi;
   Float_t         met_scaleUp;
   Float_t         metphi_scaleUp;
   Float_t         met_scaleDown;
   Float_t         metphi_scaleDown;
   Float_t         met_resoUp;
   Float_t         met_resoDown;
   Float_t         metphi_resoUp;
   Float_t         metphi_resoDown;
*/
  T->Branch("genmet", &genmet, "genmet/F");
  T->Branch("genmet_Ex", &genmet_Ex, "genmet_Ex/F");
  T->Branch("genmet_Ey", &genmet_Ey, "genmet_Ey/F");
  T->Branch("genmetphi", &genmetphi, "genmetphi/F");
  T->Branch("PtSystem", &PtSystem, "PtSystem/F");
  
  T->Branch("wScale0", &wScale0, "wScale0/F");
  T->Branch("wScale1", &wScale1, "wScale1/F");
  T->Branch("wScale2", &wScale2, "wScale2/F");
  T->Branch("wScale3", &wScale3, "wScale3/F");
  T->Branch("wScale4", &wScale4, "wScale4/F");
  T->Branch("wScale5", &wScale5, "wScale5/F");
  T->Branch("wScale6", &wScale6, "wScale6/F");
  T->Branch("wScale7", &wScale7, "wScale7/F");
  T->Branch("wScale8", &wScale8, "wScale8/F");
  T->Branch("wPDFUp", &wPDFUp, "wPDFUp/F");
  T->Branch("wPDFDown", &wPDFDown, "wPDFDown/F");

  T->Branch("NuPx", &NuPx, "NuPx/F");
  T->Branch("NuPy", &NuPy, "NuPy/F");
  T->Branch("NuPz", &NuPz, "NuPz/F");
  T->Branch("NuPt", &NuPt, "NuPt/F");
  T->Branch("NuPhi", &NuPhi, "NuPhi/F");

  T->Branch("met_scaleUp", &met_scaleUp, "met_scaleUp/F");
  T->Branch("met_scaleDown", &met_scaleDown, "met_scaleDown/F");
  T->Branch("metphi_scaleUp", &metphi_scaleUp, "metphi_scaleUp/F");
  T->Branch("metphi_scaleDown", &metphi_scaleDown, "metphi_scaleDown/F");
  T->Branch("met_resoUp", &met_resoUp, "met_resoUp/F");
  T->Branch("met_resoDown", &met_resoDown, "met_resoDown/F");
  T->Branch("metphi_resoUp", &metphi_resoUp, "metphi_resoUp/F");
  T->Branch("metphi_resoDown", &metphi_resoDown, "metphi_resoDown/F");

  T->Branch("met_ex", &met_ex, "met_ex/F");
  T->Branch("met_ey", &met_ey, "met_ey/F");
  T->Branch("met_ex_recoil", &met_ex_recoil, "met_ex_recoil/F");
  T->Branch("met_ey_recoil", &met_ey_recoil, "met_ey_recoil/F");
  T->Branch("met_ez", &met_ez, "met_ez/F");

  T->Branch("met_ex_JetEnUp", &met_ex_JetEnUp, "met_ex_JetEnUp/F");
  T->Branch("met_ey_JetEnUp", &met_ey_JetEnUp, "met_ey_JetEnUp/F");

  T->Branch("met_ex_JetEnDown", &met_ex_JetEnDown, "met_ex_JetEnDown/F");
  T->Branch("met_ey_JetEnDown", &met_ey_JetEnDown, "met_ey_JetEnDown/F");

  T->Branch("met_ex_UnclusteredEnDown", &met_ex_UnclusteredEnDown, "met_ex_UnclusteredEnDown/F");
  T->Branch("met_ey_UnclusteredEnDown", &met_ey_UnclusteredEnDown, "met_ey_UnclusteredEnDown/F");

  T->Branch("met_ex_UnclusteredEnUp", &met_ex_UnclusteredEnUp, "met_ex_UnclusteredEnUp/F");
  T->Branch("met_ey_UnclusteredEnUp", &met_ey_UnclusteredEnUp, "met_ey_UnclusteredEnUp/F");

  T->Branch("met_ex_JetEnUp_recoil", &met_ex_JetEnUp_recoil, "met_ex_JetEnUp/F_recoil");
  T->Branch("met_ey_JetEnUp_recoil", &met_ey_JetEnUp_recoil, "met_ey_JetEnUp/F_recoil");

  T->Branch("met_ex_JetEnDown_recoil", &met_ex_JetEnDown_recoil, "met_ex_JetEnDown/F_recoil");
  T->Branch("met_ey_JetEnDown_recoil", &met_ey_JetEnDown_recoil, "met_ey_JetEnDown/F_recoil");

  T->Branch("met_ex_UnclusteredEnDown_recoil", &met_ex_UnclusteredEnDown_recoil, "met_ex_UnclusteredEnDown/F_recoil");
  T->Branch("met_ey_UnclusteredEnDown_recoil", &met_ey_UnclusteredEnDown_recoil, "met_ey_UnclusteredEnDown/F_recoil");

  T->Branch("met_ex_UnclusteredEnUp_recoil", &met_ex_UnclusteredEnUp_recoil, "met_ex_UnclusteredEnUp/F_recoil");
  T->Branch("met_ey_UnclusteredEnUp_recoil", &met_ey_UnclusteredEnUp_recoil, "met_ey_UnclusteredEnUp/F_recoil");


  T->Branch("met_pt", &met_pt, "met_pt/F");
  T->Branch("met_phi", &met_phi, "met_phi/F");
 
  T->Branch("gen_weight", &gen_weight, "gen_weight/F");
  T->Branch("pu_weight", &pu_weight, "pu_weight/F");
  T->Branch("LSF_weight", &LSF_weight, "LSF_weight/F");
  T->Branch("LSF_weight_mu", &LSF_weight_mu, "LSF_weight_mu/F");
  T->Branch("LSF_weight_el", &LSF_weight_el, "LSF_weight_el/F");
  T->Branch("LSF_weight_1", &LSF_weight_1, "LSF_weight_1/F");
  T->Branch("LSF_weight_2", &LSF_weight_2, "LSF_weight_2/F");
  T->Branch("TFR_weight", &TFR_weight, "TFR_weight/F");
  T->Branch("top_weight", &top_weight, "top_weight/F");
  T->Branch("all_weight", &all_weight, "all_weight/F");
  T->Branch("trig_weight", &trig_weight, "trig_weight/F");
  T->Branch("trig_weight_1", &trig_weight_1, "trig_weight_1/F");
  T->Branch("trig_weight_2", &trig_weight_2, "trig_weight_2/F");
  T->Branch("zptmassweight", &zptmassweight, "zptmassweight/F");

  T->Branch("xsecs", &xsecs, "xsecs/F");
  T->Branch("event_sign", &event_sign, "event_sign/F");
  T->Branch("met_flag", &met_flag, "met_flag/F");
  T->Branch("event_secondLeptonVeto", &event_secondLeptonVeto, "event_secondLeptonVeto/F");
  T->Branch("event_thirdLeptonVeto", &event_thirdLeptonVeto, "event_thirdLeptonVeto/F");
  T->Branch("event_leptonDrTrigger", &event_leptonDrTrigger, "event_leptonDrTrigger/F");
  T->Branch("eleMVA", &eleMVA, "eleMVA/F");

  T->Branch("muon_index", &muon_index, "muon_index/I");
  T->Branch("muon_index_1", &muon_index_1, "muon_index_1/I");
  T->Branch("muon_index_2", &muon_index_2, "muon_index_2/I");
  T->Branch("electron_index", &electron_index, "electron_index/I");
  T->Branch("taus_index", &taus_index, "taus_index/I");
  T->Branch("taus_index2", &taus_index2, "taus_index2/I");

  T->Branch("primvert_count", &primvert_count, "primvert_count/I");
  T->Branch("primvert_x", &primvert_x, "primvert_x/F");
  T->Branch("primvert_y", &primvert_x, "primvert_y/F");
  T->Branch("primvert_z", &primvert_x, "primvert_z/F");

  T->Branch("mu_count", &mu_count, "mu_count/I");
  T->Branch("mu_px", mu_px, "mu_px[20]/F");
  T->Branch("mu_py", mu_py, "mu_py[20]/F");
  T->Branch("mu_pz", mu_pz, "mu_pz[20]/F");
  T->Branch("mu_pt", mu_pt, "mu_pt[20]/F");
  T->Branch("mu_eta", mu_eta, "mu_eta[20]/F");
  T->Branch("mu_phi", mu_phi, "mu_phi[20]/F");
  T->Branch("mu_charge", mu_charge, "mu_charge[20]/F");
  T->Branch("mu_miniISO", mu_miniISO, "mu_miniISO[20]/F");
  T->Branch("mu_dxy", mu_dxy, "mu_dxy[20]/F");
  T->Branch("mu_dz", mu_dz, "mu_dz[20]/F");
  T->Branch("mu_dxyerr", mu_dxyerr, "mu_dxyerr[20]/F");
  T->Branch("mu_dzerr", mu_dzerr, "mu_dzerr[20]/F");
  T->Branch("mu_relIso", mu_relIso, "mu_relIso[20]/F");
 
  T->Branch("mu_neutralHadIso", mu_neutralHadIso, "mu_neutralHadIso[20]/F");
  T->Branch("mu_photonIso", mu_photonIso, "mu_photonIso[20]/F");
  T->Branch("mu_chargedHadIso", mu_chargedHadIso, "mu_chargedHadIso[20]/F");
  T->Branch("mu_puIso", mu_puIso, "mu_puIso[20]/F");
  T->Branch("mu_neutralIso", mu_neutralIso, "mu_neutralIso[20]/F");
  T->Branch("mu_absIsoMu", mu_absIsoMu, "mu_absIsoMu[20]/F");
  T->Branch("mu_relIsoMu", mu_relIsoMu, "mu_relIsoMu[20]/F");

  T->Branch("el_neutralHadIso", el_neutralHadIso, "el_neutralHadIso[20]/F");
  T->Branch("el_photonIso", el_photonIso, "el_photonIso[20]/F");
  T->Branch("el_chargedHadIso", el_chargedHadIso, "el_chargedHadIso[20]/F");
  T->Branch("el_puIso", el_puIso, "el_puIso[20]/F");
  T->Branch("el_neutralIso", el_neutralIso, "el_neutralIso[20]/F");
  T->Branch("el_absIsoEl", el_absIsoEl, "el_absIsoEl[20]/F");
  T->Branch("el_relIsoEl", el_relIsoEl, "el_relIsoEl[20]/F");

  T->Branch("jet_count", &jet_count, "jet_count/I");
  T->Branch("njets", &njets, "njets/I");
  T->Branch("genTauDecayMode1", &genTauDecayMode1, "genTauDecayMode1/I");
  T->Branch("genTauDecayMode2", &genTauDecayMode2, "genTauDecayMode2/I");
  T->Branch("genTauDecayMode", &genTauDecayMode, "genTauDecayMode/I");
  T->Branch("genTauDecayModeAll", &genTauDecayModeAll, "genTauDecayModeAll/I");
  T->Branch("event_run", &event_run, "event_run/I");
  T->Branch("event_lumi", &event_lumi, "event_lumi/I");
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
  T->Branch("el_px", el_px, "el_px[20]/F");
  T->Branch("el_py", el_py, "el_py[20]/F");
  T->Branch("el_pz", el_pz, "el_pz[20]/F");
  T->Branch("el_pt", el_pt, "el_pt[20]/F");
  T->Branch("el_eta", el_eta, "el_eta[20]/F");
  T->Branch("el_phi", el_phi, "el_phi[20]/F");
  T->Branch("el_miniISO", el_miniISO, "el_miniISO[20]/F");
  T->Branch("el_dxy", el_dxy, "el_dxy[20]/F");
  T->Branch("el_dz", el_dz, "el_dz[20]/F");
  T->Branch("el_dxyerr", el_dxyerr, "el_dxyerr[20]/F");
  T->Branch("el_dzerr", el_dzerr, "el_dzerr[20]/F");
  T->Branch("el_charge", el_charge, "el_charge[20]/F");
  T->Branch("el_relIso", el_relIso, "el_relIso[20]/F");
  T->Branch("el_isMVA", el_isMVA, "el_isMVA[20]/F");
  T->Branch("el_isnotrig_MVA80", el_isnotrig_MVA80, "el_isnotrig_MVA80[20]/F");
  T->Branch("el_isnotrig_MVA90", el_isnotrig_MVA90, "el_isnotrig_MVA90[20]/F");
  T->Branch("el_istrig_MVA80", el_istrig_MVA80, "el_istrig_MVA80[20]/F");
  T->Branch("el_istrig_MVA90", el_istrig_MVA90, "el_istrig_MVA90[20]/F");


  T->Branch("ta_count", &ta_count, "ta_count/I");
  T->Branch("ta_px", ta_px, "ta_px[20]/F");
  T->Branch("ta_py", ta_py, "ta_py[20]/F");
  T->Branch("ta_pz", ta_pz, "ta_pz[20]/F");
  T->Branch("ta_mass", ta_mass, "ta_mass[20]/F");
  T->Branch("ta_eta", ta_eta, "ta_eta[20]/F");
  T->Branch("ta_phi", ta_phi, "ta_phi[20]/F");
  T->Branch("ta_pt", ta_pt, "ta_pt[20]/F");
  T->Branch("ta_dxy", ta_dxy, "ta_dxy[20]/F");
  T->Branch("ta_dz", ta_dz, "ta_dz[20]/F");
  T->Branch("ta_charge", ta_charge, "ta_charge[20]/F");
  T->Branch("ta_relIso", ta_relIso, "ta_relIso[20]/F");
  T->Branch("ta_IsoFlag", &ta_IsoFlag, "ta_IsoFlag/F");
  T->Branch("ta_chargedIsoPtSum", &ta_chargedIsoPtSum, "ta_chargedIsoPtSum/F");
  T->Branch("ta_neutralIsoPtSum", &ta_neutralIsoPtSum, "ta_neutralIsoPtSum/F");
  T->Branch("ta_puCorrPtSum", &ta_puCorrPtSum, "ta_puCorrPtSum/F");

  T->Branch("isMuTau", &isMuTau, "isMuTau/F");
  T->Branch("isElTau", &isElTau, "isElTau/F");
  T->Branch("isTauTau", &isTauTau, "isTauTau/F");

  T->Branch("ta_IsoFlagVTight", ta_IsoFlagVTight, "ta_IsoFlagVTight[5]/F");
  T->Branch("ta_IsoFlagTight", ta_IsoFlagTight, "ta_IsoFlagTight[5]/F");
  T->Branch("ta_IsoFlagLoose", ta_IsoFlagLoose, "ta_IsoFlagLoose[5]/F");
  T->Branch("ta_IsoFlagVLoose", ta_IsoFlagVLoose, "ta_IsoFlagVLoose[5]/F");
  T->Branch("ta_IsoFlagMedium", ta_IsoFlagMedium, "ta_IsoFlagMedium[5]/F");

  T->Branch("ta_isLoose", &ta_isLoose, "ta_isLoose/F");
  T->Branch("ta_isTight", &ta_isTight, "ta_isTigh/F");



  T->Branch("qcdweight", &qcdweight, "qcdweight/F");
  T->Branch("qcdweightup", &qcdweightup, "qcdweightup/F");
  T->Branch("qcdweightdown", &qcdweightdown, "qcdweightdown/F");
  T->Branch("nbtag", &nbtag, "nbtag/I");

  T->Branch("datasetName", &datasetName);
  T->Branch("regionName", &regionName);
  T->Branch("genTauMatched", &genTauMatched);
  T->Branch("genTauMatched2", &genTauMatched2);

  T->Branch("genLeptonMatchedPrompEl", &genLeptonMatchedPrompEl);
  T->Branch("genLeptonMatchedPrompMu", &genLeptonMatchedPrompMu);
  T->Branch("genElMatchedToTauDecay", &genElMatchedToTauDecay);
  T->Branch("genMuMatchedToTauDecay", &genMuMatchedToTauDecay);
  T->Branch("genTauMatchedToTauDecay", &genTauMatchedToTauDecay);
  T->Branch("matchedTauToPromptEl", &matchedTauToPromptEl);
  T->Branch("matchedTauToPromptMu", &matchedTauToPromptMu);
  T->Branch("matchedTauToPromptTau", &matchedTauToPromptTau);
  T->Branch("matchedTauToTauDecEl", &matchedTauToTauDecEl);
  T->Branch("matchedTauToTauDecMu", &matchedTauToTauDecMu);
  T->Branch("matchedTauToTauDecTau", &matchedTauToTauDecTau);
  T->Branch("matchedTauToElHadronDec", &matchedTauToElHadronDec);
  T->Branch("matchedTauToMuHadronDec", &matchedTauToMuHadronDec);
  T->Branch("matchedTauToTauHadronDec", &matchedTauToTauHadronDec);
  T->Branch("matchedTauToGluon", &matchedTauToGluon);
  T->Branch("matchedTauToHFQ", &matchedTauToHFQ);
  T->Branch("matchedTauToLFQ", &matchedTauToLFQ);
  T->Branch("genLeptonMatchedPromptEl", &genLeptonMatchedPromptEl);
  T->Branch("genLeptonMatchedPromptMu", &genLeptonMatchedPromptMu);
  T->Branch("genLeptonMatchedPromptTau", &genLeptonMatchedPromptTau);
  T->Branch("genElMatchedHadrDecay", &genElMatchedHadrDecay);
  T->Branch("genMuMatchedHadrDecay", &genMuMatchedHadrDecay);
  T->Branch("genTauMatchedHadrDecay", &genTauMatchedHadrDecay);
  T->Branch("genLeptonMatchedGluon", &genLeptonMatchedGluon);
  T->Branch("genLeptonMatchedLFQ", &genLeptonMatchedLFQ);
  T->Branch("genLeptonMatchedHFQ", &genLeptonMatchedHFQ);
  
  
  T->Branch("isDYTT", &isDYTT);
  T->Branch("isDYLL", &isDYLL);
  T->Branch("isDYEE", &isDYEE);
  T->Branch("isDYMM", &isDYMM);
  T->Branch("isDYNuNu", &isDYNuNu);



  T->Branch("genLeptonPromptElMatched", &genLeptonPromptElMatched);
  T->Branch("genLeptonPromptMuMatched", &genLeptonPromptMuMatched);
  T->Branch("genLeptonMatched", &genLeptonMatched);
  T->Branch("genLeptonMatchedEl", &genLeptonMatchedEl);
  T->Branch("genLeptonMatchedMu", &genLeptonMatchedMu);
  T->Branch("genLeptonMatched2", &genLeptonMatched2);
  T->Branch("genTauDecayedMuMatched", &genTauDecayedMuMatched);
  T->Branch("genTauDecayedElMatched", &genTauDecayedElMatched);
  T->Branch("npartons",&npartons,"npartons/I");
  T->Branch("SusyMother",&SusyMother,"SusyMother/F");
  T->Branch("SusyLSP",&SusyLSP,"SusyLSP/F");



      char arg[200];
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

void SetupHists(int CutNer){


for(int cj = 0; cj < CutNer; cj++)
    {
      CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      CutFlowUnW->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
    }
}

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV,vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_){}


void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel){}

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex){};

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex, int category_){};


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

