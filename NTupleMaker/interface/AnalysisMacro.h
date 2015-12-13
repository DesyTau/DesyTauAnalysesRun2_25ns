
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/Basic_Mt2_332_Calculator.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/Basic_MPairProd_Calculator.h"
#include "DesyTauAnalyses/NTupleMaker/interface/lester_mt2_bisect.h"
#include "DesyTauAnalyses/NTupleMaker/interface/mTBound.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
using namespace std;


const  int CutN=18;


unsigned int RunMin = 9999999;
unsigned int RunMax = 0;
     
unsigned int tau_index=-1;
unsigned int mu_index=-1;
unsigned int el_index=-1;

//bool isData = false;

double ChiMass=0;
double mIntermediate = tauMass;

double XSec=-1;
double xs,fact,fact2;
  

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

//string CutList[10];
vector<string> CutList;


TH1D * histRuns = new TH1D("histRuns","",6000,24000,30000);

TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);

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

TH1D *hMTmutau[CutN];
TH1D *hMTmuel[CutN];
TH1D *hMTeltau[CutN];
TH1D *hMTtautau[CutN];

TH1D *hMt2mutau[CutN];
TH1D *hTBoundmutau[CutN];

TH1D *hdR_mutau[CutN];
TH1D *hdR_eltau[CutN];
TH1D *hdR_tautau[CutN];
TH1D *hdR_muel[CutN];


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
  
TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutN,0.5,CutN+0.5);
TH1D *CutFlowUnW= new TH1D("CutFlowUnW","Cut Flow",CutN,0.5,CutN+0.5);
TH1D *CutFlowUnWNorm= new TH1D("CutFlowUnWNorm","Cut Flow",CutN,0.5,CutN+0.5);

TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);


TH1D * muonPtAllH = new TH1D("muonPtAllH","",40,0,200);
TH1D * electronPtAllH = new TH1D("electronPtAllH","",40,0,200);
TH1D * tauPtAllH = new TH1D("tauPtAllH","",40,0,200);

// histograms (dilepton selection)
TH1D * electronPtH  = new TH1D("electronPtH","",40,0,200);
TH1D * electronEtaH = new TH1D("electronEtaH","",50,-2.5,2.5); 
TH1D * muonPtH  = new TH1D("muonPtH","",40,0,200);
TH1D * muonEtaH = new TH1D("muonEtaH","",50,-2.5,2.5); 
TH1D * tauEtaAllH = new TH1D("tauEtaAllH","",50,-2.5,2.5); 
 
TH1D * dileptonMassH = new TH1D("dileptonMassH","",40,0,200);
TH1D * dileptonPtH = new TH1D("dileptonPtH","",40,0,200);
TH1D * dileptonEtaH = new TH1D("dileptonEtaH","",100,-5,5);
TH1D * dileptondRH = new TH1D("dileptondRH","",60,0,6);
TH1D * ETmissH = new TH1D("ETmissH","",40,0,200);
TH1D * MtH = new TH1D("MtH_2l","",40,0,200);
TH1D * DZetaH = new TH1D("DZetaH","",60,-400,200);

// histograms (dilepton selection + DZeta cut DZeta)
TH1D * electronPtSelH  = new TH1D("electronPtSelH","",40,0,200);
TH1D * electronEtaSelH = new TH1D("electronEtaSelH","",50,-2.5,2.5); 
TH1D * muonPtSelH  = new TH1D("muonPtSelH","",40,0,200);
TH1D * muonEtaSelH = new TH1D("muonEtaSelH","",50,-2.5,2.5); 

TH1D * dileptonMassSelH = new TH1D("dileptonMassSelH","",40,0,200);
TH1D * dileptonPtSelH = new TH1D("dileptonPtSelH","",40,0,200);
TH1D * dileptonEtaSelH = new TH1D("dileptonEtaSelH","",100,-5,5);
TH1D * dileptondRSelH = new TH1D("dileptondRSelH","",60,0,6);
TH1D * ETmissSelH = new TH1D("ETmissSelH","",40,0,200);
TH1D * MtSelH = new TH1D("MtSelH_2l","",40,0,200);
TH1D * DZetaSelH = new TH1D("DZetaSelH","",60,-400,200);

TLorentzVector ElV, MuV, TauV, JetsV, METV;

vector<TLorentzVector> AllJets_Lepton_noMet;
vector<TLorentzVector> JetsMV;
vector<TLorentzVector>  ElMV;
vector<TLorentzVector>  MuMV;
vector<TLorentzVector>  TauMV;
vector<TLorentzVector>  LeptMV;







void WriteHists(int CutNer, TFile *in, TString dir){

  in->cd(dir);

  for(int cj = 0; cj < CutNer; cj++)
    {

      hHTOsqrMET[cj]->Write();
    }
}

//string CutList[CutN];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNer){
  for(int cj = 0; cj < CutNer; cj++)
    {
      CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      CutFlowUnW->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      TString cutName=CutList[cj];
      TString nCut;
      nCut.Form("%i",cj);
      ///generic variables
      //
      //Jets
      hmu_dxy[cj] = new TH1D ("hmu_dxy_"+nCut,"hmu_dxy "+cutName,100,-5.,5.);
      hel_dxy[cj] = new TH1D ("hel_dxy_"+nCut,"hel_dxy "+cutName,100,-5.,5.);
      htau_dxy[cj] = new TH1D ("htau_dxy_"+nCut,"htau_dxy "+cutName,100,-5.,5.);
      
      htau_dz[cj] = new TH1D ("htau_dz_"+nCut,"htau_dz "+cutName,200,-50.,50.);
      hel_dz[cj] = new TH1D ("hel_dz_"+nCut,"hel_dz "+cutName,200,-50.,50.);
      hmu_dz[cj] = new TH1D ("hmu_dz_"+nCut,"hmu_dz "+cutName,200,-50.,50.);

      hHTOsqrMET[cj] = new TH1D ("hHTOsqrMET_"+nCut,"hHTOsqrMET "+cutName,400,0.0,4000.0);
      hHTOsqrMET[cj]->Sumw2();
      
      hMeffMuonOsqrMET[cj] = new TH1D ("hMeffMuonOsqrMET_"+nCut,"hMeffMuonOsqrMET "+cutName,400,0.0,4000.0);
      hMeffMuonOsqrMET[cj]->Sumw2();
      hMeffElOsqrMET[cj] = new TH1D ("hMeffElOsqrMET_"+nCut,"hMeffElOsqrMET "+cutName,400,0.0,4000.0);
      hMeffElOsqrMET[cj]->Sumw2();
      hMeffTauOsqrMET[cj] = new TH1D ("hMeffTauOsqrMET_"+nCut,"hMeffTauOsqrMET "+cutName,400,0.0,4000.0);
      hMeffTauOsqrMET[cj]->Sumw2();

      hMeffMuon[cj] = new TH1D ("hMeffMuon_"+nCut,"hMeffMuon "+cutName,400,0.0,4000.0);
      hMeffMuon[cj]->Sumw2();
      hMeffEl[cj] = new TH1D ("hMeffEl_"+nCut,"hMeffEl "+cutName,400,0.0,4000.0); 
      hMeffEl[cj]->Sumw2();
      hMeffTau[cj] = new TH1D ("hMeffTau_"+nCut,"hMeffTau "+cutName,400,0.0,4000.0); 
      hMeffTau[cj]->Sumw2();

      hCentrality[cj]  = new TH1D ("hCentrality_"+nCut,"hCentrality "+cutName,20,0.,2.);
      hCentrality[cj]->Sumw2();

      hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,400,0.0,4000.0);
      hHT[cj]->Sumw2();
 
      hRht[cj] = new TH1D ("Rht_"+nCut,"Rht "+cutName,20,0.0,2.0);
      hRht[cj]->Sumw2();
    
      hHText[cj] = new TH1D ("HText_"+nCut,"HText "+cutName,400,0.0,4000.0);
      hHText[cj]->Sumw2();
      hHT2[cj] = new TH1D ("HT2_"+nCut,"HT2 "+cutName,400,0.0,4000.0);
      hHT2[cj]->Sumw2();
      hHT3[cj] = new TH1D ("HT3_"+nCut,"HT3 "+cutName,400,0.0,4000.0);
      hHT3[cj]->Sumw2();
      hHT4[cj] = new TH1D ("HT4_"+nCut,"HT4 "+cutName,400,0.0,4000.0);
      hHT4[cj]->Sumw2();
          
      hPtJ0[cj] = new TH1D ("hPtJ0_"+nCut,"hPtJ0 "+cutName,200,0.0,2000.0);
      hPtJ0[cj]->Sumw2();
      hPtJ1[cj] = new TH1D ("hPtJ1_"+nCut,"hPtJ1 "+cutName,200,0.0,2000.0);
      hPtJ1[cj]->Sumw2();
      hPtJ2[cj] = new TH1D ("hPtJ2_"+nCut,"hPtJ2 "+cutName,200,0.0,2000.0);
      hPtJ2[cj]->Sumw2();
      hPtJ3[cj] = new TH1D ("hPtJ3_"+nCut,"hPtJ3 "+cutName,200,0.0,2000.0);
      hPtJ3[cj]->Sumw2();

      //h0JetpT[cj] = new TH1D ("0JetpT_"+nCut,"0JetpT "+cutName,200,0.0,2000.0);
      //h0JetpT[cj]->Sumw2();
      hnJet[cj] = new TH1D ("nJet_"+nCut,"nJet "+cutName,25,-0.5,24.5);
      hnJet[cj]->Sumw2();
      hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,10,-0.5,9.5);
      hnBJet[cj]->Sumw2();

      hWeights[cj] = new TH1D ("hWeights_"+nCut,"hWeights "+cutName,10,-1,9);
      hWeights[cj]->Sumw2();
	
      hInvMassMuTau[cj] = new TH1D ("hInvMassMuTau_"+nCut,"hInvMassMuTau "+cutName,80,0,160);
      hInvMassMuMu[cj] = new TH1D ("hInvMassMuMu_"+nCut,"hInvMassMuMu "+cutName,80,0,160);
      hInvMassElTau[cj] = new TH1D ("hInvMassElTau_"+nCut,"hInvMassElTau "+cutName,80,0,160);
      hInvMassElEl[cj] = new TH1D ("hInvMassElEl_"+nCut,"hInvMassElEl "+cutName,80,0,160);
      hInvMassTauTau[cj] = new TH1D ("hInvMassTauTau_"+nCut,"hInvMassTauTau "+cutName,80,0,160);
        
      //Leptons
      //
      //
      hnLep[cj] = new TH1D ("nLep_"+nCut,"nLep "+cutName,10,-0.5,9.5);
      hnLep[cj]->Sumw2();
      hLeppt[cj] = new TH1D ("LeppT_"+nCut,"Lep pT "+cutName,100,0,1000);
      hLeppt[cj]->Sumw2();
      hLepeta[cj] = new TH1D ("Lepeta_"+nCut,"Lep eta "+cutName,100,-4,4);
      hLepeta[cj]->Sumw2();
      //hST[cj] = new TH1D ("ST_"+nCut,"ST "+cutName,400,0.0,4000.0);
      //hST[cj]->Sumw2();
        
      //Muons
      //
      //
      hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,-0.5,9.5);
      hnMu[cj]->Sumw2();
      hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,100,0,1000);
      hMupt[cj]->Sumw2();
      hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,100,-4,4);
      hMueta[cj]->Sumw2();
        
      //Taus
      //
      //
      hnTau[cj] = new TH1D ("nTau_"+nCut,"nTau "+cutName,10,-0.5,9.5);
      hnTau[cj]->Sumw2();
      hTaupt[cj] = new TH1D ("TaupT_"+nCut,"Tau pT "+cutName,100,0,1000);
      hTaupt[cj]->Sumw2();
      hTaueta[cj] = new TH1D ("Taueta_"+nCut,"Tau eta "+cutName,100,-4,4);
      hTaueta[cj]->Sumw2();
	
      //hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
      //Electrons
      //
      //
      hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,-0.5,9.5);
      hnEl[cj]->Sumw2();
      hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,100,0,1000);
      hElpt[cj]->Sumw2();
      hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,100,-4,4);
      hEleta[cj]->Sumw2();
       
       
      hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,200.0,0.0,2000.0);
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
      hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,40,0,200);
      hMT[cj]->Sumw2();
      hMTel[cj] = new TH1D ("MTel_"+nCut,"MTel "+cutName,40,0,200);
      hMTel[cj]->Sumw2();
      hMTmu[cj] = new TH1D ("MTmu_"+nCut,"MTmu "+cutName,40,0,200);
      hMTmu[cj]->Sumw2();

      hMTtau[cj] = new TH1D ("MTtau_"+nCut,"MTtau "+cutName,40,0,200);
      hMTtau[cj]->Sumw2();

      hMTmutau[cj] = new TH1D ("MTmutau_"+nCut,"MTmutau "+cutName,40,0,200);
      hMTmutau[cj]->Sumw2();
     
      hMTmuel[cj] = new TH1D ("MTmuel_"+nCut,"MTmuel "+cutName,40,0,200);
      hMTmuel[cj]->Sumw2();
       
      hMTeltau[cj] = new TH1D ("MTeltau_"+nCut,"MTmeltau "+cutName,40,0,200);
      hMTeltau[cj]->Sumw2();
      
      hMTtautau[cj] = new TH1D ("MTtautau_"+nCut,"MTtautau "+cutName,40,0,200);
      hMTtautau[cj]->Sumw2();

      hMt2mutau[cj] = new TH1D ("Mt2mutau_"+nCut,"Mt2mutau "+cutName,200,0,1000);
      hMt2mutau[cj]->Sumw2();
      hTBoundmutau[cj] = new TH1D ("TBoundmutau_"+nCut,"TBoundmutau "+cutName,400,0,2000);
      hTBoundmutau[cj] ->Sumw2();

      /*       hDZeta[CutN]= new TH1D ("DZeta_"+nCut,"DZeta "+cutName,300,-400,200);
	       hDZeta[cj]->Sumw2();
      */
      htau_ISOL[cj]= new TH1D ("tauISOL_"+nCut,"tauISOL "+cutName,100,0,5);;
      htau_ISO[cj]= new TH1D ("tauISO_"+nCut,"tauISO "+cutName,100,0,5);;
      htau_ISO[cj] ->Sumw2(); 
      htau_ISOL[cj] ->Sumw2(); 

      hel_miniISO[cj]= new TH1D ("elminiISO_"+nCut,"elminiISO "+cutName,100,0,5);;
      hel_miniISO[cj]->Sumw2();
      hel_miniISOL[cj]= new TH1D ("elminiISOL_"+nCut,"elminiISOL "+cutName,100,0,5);;
      hel_miniISOL[cj]->Sumw2();

      hel_relISO[cj]= new TH1D ("elrelISO_"+nCut,"elrelISO "+cutName,100,0,5);;
      hel_relISO[cj]->Sumw2();
      hel_relISOL[cj]= new TH1D ("elrelISOL_"+nCut,"elrelISOL "+cutName,100,0,5);;
      hel_relISOL[cj]->Sumw2();
        
        
      hmu_miniISO[cj]= new TH1D ("muminiISO_"+nCut,"muminiISO "+cutName,100,0,5);;
      hmu_miniISO[cj]->Sumw2();
      hmu_miniISOL[cj]= new TH1D ("muminiISOL_"+nCut,"muminiISOL "+cutName,100,0,5);;
      hmu_miniISOL[cj]->Sumw2();

      hmu_relISO[cj]= new TH1D ("murelISO_"+nCut,"murelISO "+cutName,100,0,5);;
      hmu_relISO[cj]->Sumw2();
      hmu_relISOL[cj]= new TH1D ("murelISOL_"+nCut,"murelISOL "+cutName,100,0,5);;
      hmu_relISOL[cj]->Sumw2();
 
      hdR_eltau[cj]= new TH1D ("dR_eltau_"+nCut,"dR_eltau "+cutName,60,0,6);;
      hdR_eltau[cj]->Sumw2();
        
      hdR_mutau[cj]= new TH1D ("dR_mutau_"+nCut,"dR_mutau "+cutName,60,0,6);;
      hdR_mutau[cj]->Sumw2();

      hdR_tautau[cj]= new TH1D ("dR_tautau_"+nCut,"dR_tautau "+cutName,60,0,6);;
      hdR_tautau[cj]->Sumw2();
	
      hdR_muel[cj]= new TH1D ("dR_muel_"+nCut,"dR_muel "+cutName,60,0,6);;
      hdR_muel[cj]->Sumw2();

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
 
      hmet_dPhi[cj] = new TH2D ("met_dPhi_"+nCut,"met_dPhi "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
      hmet_dPhi[cj]->Sumw2();
      hmet_MT[cj] = new TH2D ("met_MT_"+nCut,"met_MT "+cutName,200.0,0.0,2000.0,40,0,200);
      hmet_MT[cj]->Sumw2();

      hmet_dPhiel[cj] = new TH2D ("met_dPhiel_"+nCut,"met_dPhiel "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
      hmet_dPhiel[cj]->Sumw2();
      hmet_MTel[cj] = new TH2D ("met_MTel_"+nCut,"met_MTel "+cutName,200.0,0.0,2000.0,40,0,200);
      hmet_MTel[cj]->Sumw2();
 

      hmet_dPhimu[cj] = new TH2D ("met_dPhimu_"+nCut,"met_dPhimu "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
      hmet_dPhimu[cj]->Sumw2();
      hmet_MTmu[cj] = new TH2D ("met_MTmu_"+nCut,"met_MTmu "+cutName,200.0,0.0,2000.0,40,0,200);
      hmet_MTmu[cj]->Sumw2();

      hmet_dPhitau[cj] = new TH2D ("met_dPhitau_"+nCut,"met_dPhitau "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
      hmet_dPhitau[cj]->Sumw2();
      hmet_MTtau[cj] = new TH2D ("met_MTtau_"+nCut,"met_MTtau "+cutName,200.0,0.0,2000.0,40,0,200);
      hmet_MTtau[cj]->Sumw2();


 
      hMT_dPhi[cj]= new TH2D ("MT_dPhi_"+nCut,"MT_dPhi "+cutName,40,0,200,64,0.0,3.2);
      hMT_dPhi[cj]->Sumw2();
      hMT_dPhiel[cj]= new TH2D ("MTel_dPhi_"+nCut,"MTel_dPhi "+cutName,40,0,200,64,0.0,3.2);
      hMT_dPhiel[cj]->Sumw2();
      hMT_dPhimu[cj]= new TH2D ("MTmu_dPhi_"+nCut,"MTmu_dPhi "+cutName,40,0,200,64,0.0,3.2);
      hMT_dPhimu[cj]->Sumw2();
      hMT_dPhitau[cj]= new TH2D ("MTtau_dPhi_"+nCut,"MTtau_dPhi "+cutName,40,0,200,64,0.0,3.2);
      hMT_dPhitau[cj]->Sumw2();
      hIso_sign[cj]= new TH2D ("Iso_sign_"+nCut,"Iso_sign "+cutName,100,0,5,2,0,2);
      hIso_sign[cj]->GetYaxis()->SetBinLabel(1,"SS");
      hIso_sign[cj]->GetYaxis()->SetBinLabel(2,"OS");
      hIso_sign[cj]->Sumw2();

    }
}

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



void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, double Chimass, double mintermediate,AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex){	

    AllJets_Lepton_noMet.clear();
    for (unsigned int i = 0; i <   JetsV.size(); i++) AllJets_Lepton_noMet.push_back (JetsV.at(i));

 if (mIndex<0 && tree_.muon_count > 0) mIndex=0;  
 if (eIndex<0 && tree_.electron_count > 0) eIndex=0;  
 if (tIndex<0 && tree_.tau_count > 0) tIndex=0;  

  TLorentzVector muV ; if (tree_.muon_count > 0) muV.SetPtEtaPhiM(tree_.muon_pt[mIndex], tree_.muon_eta[mIndex], tree_.muon_phi[mIndex], muonMass);
  TLorentzVector tauV; if (tree_.tau_count > 0) tauV.SetPtEtaPhiM(tree_.tau_pt[tIndex], tree_.tau_eta[tIndex], tree_.tau_phi[tIndex], tauMass);
  TLorentzVector elV;  if (tree_.electron_count > 0) elV.SetPtEtaPhiM(tree_.electron_pt[eIndex], tree_.electron_eta[eIndex], tree_.electron_phi[eIndex], tauMass);
		
if (Sel=="mutau" && mIndex !=-1 && MuV.size()>0)    AllJets_Lepton_noMet.push_back (muV);
if (Sel=="eltau" && eIndex !=-1 && ElV.size()>0)    AllJets_Lepton_noMet.push_back (elV);


	//cout<<" Init  "<<mIndex<<"  "<<eIndex<<"  "<<tIndex<<endl;

  double sumpT=0;

  double sumMuonpT=0;
  double sumElpT=0;
  double sumTaupT=0;  

  hnJet[CutIndex]->Fill(tree_.pfjet_count,EvWeight);
  hnMu[CutIndex]->Fill(tree_.muon_count,EvWeight);
  hnTau[CutIndex]->Fill(tree_.tau_count,EvWeight);
  hnEl[CutIndex]->Fill(tree_.electron_count,EvWeight);
  hnLep[CutIndex]->Fill(tree_.electron_count+tree_.muon_count+tree_.tau_count,EvWeight);
  hnpu[CutIndex]->Fill(tree_.numpileupinteractions,EvWeight);
  hnpv[CutIndex]->Fill(tree_.primvertex_count,EvWeight);
  hnrho[CutIndex]->Fill(tree_.rho,EvWeight);
  hWeights[CutIndex]->Fill(EvWeight);


   

		  
  //Float_T mT  = double Lester::mTBound(...)

  //     if (JetsV.size() > 0) h0JetpT[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);
  //  if(FillBJets){
  //      hnBJet[CutIndex]->Fill(Obj.nBJetGood,EvWeight);
  //  }
  //
  //cout << " Passing arguments befor" << JetsV.size()<<" muV "<<MuV.size()<<" tauV "<<TauV.size()<<" elV "<<ElV.size()<<" elI "<<eIndex<<" mI "<<mIndex<<"  tI "<<tIndex<<endl;
  //if (Sel=="mutau" && MuV.size()>0 && TauV.size()>0 && mIndex >-1 && tIndex >-1 ){

  if (Sel=="mutau" && tree_.muon_count > 0 && tree_.tau_count>0 ){

      	  double Mt2 = 0;
    Mt2 =asymm_mt2_lester_bisect::get_mT2(muonMass, tree_.muon_px[mIndex], tree_.muon_py[mIndex],tauMass,tree_.tau_px[tIndex],tree_.muon_py[tIndex],MetV.Px(),MetV.Py(),Chimass,Chimass,0);
    //   cout<<" Will call with Mt2 "<<Mt2<<"   "<<Chimass<<endl;
    hMt2mutau[CutIndex]->Fill(Mt2,EvWeight);
 
    double Centr = Centrality(AllJets_Lepton_noMet);
    hCentrality[CutIndex]->Fill(Centr,EvWeight);
    //cout<<"Centr "<<Centr<<"  "<< EvWeight <<endl;
    //double mTB = Lester::mTBound(muV.E(), muV.Px(), muV.Py(), muV.Pz(), tauV.E(),tauV.Px(),tauV.Py(),tauV.Pz(),MetV.Px(),MetV.Py(), mintermediate);
    //hTBoundmutau[CutIndex]->Fill(mTB);
    //cout<<"  "<<mTB<<endl;
    //	     	double Dr= deltaR(MuV.at(mIndex).Eta(), MuV.at(mIndex).Phi(),TauV.at(mIndex).Eta(),TauV.at(mIndex).Phi());
    //double Dr= deltaR(MuV.at(0).Eta(), MuV.at(0).Phi(),TauV.at(0).Eta(),TauV.at(0).Phi());
    //double Dr = MuV.at(mIndex).DeltaR(TauV.at(tIndex));
    double Dr = muV.DeltaR(tauV);
    //
    TLorentzVector DiL = muV  + tauV;
    double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
    double MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
    hMTmutau[CutIndex]->Fill(MT,EvWeight);
    hInvMassMuTau[CutIndex]->Fill(DiL.M(),EvWeight);
    hdR_mutau[CutIndex]->Fill(Dr,EvWeight);

    hdRmt_dPhi[CutIndex]->Fill(Dr,dPhi,EvWeight);



      TLorentzVector WBos = MetV + muV;
  

      dPhi=dPhiFrom2P( tree_.muon_px[mIndex], tree_.muon_py[mIndex], MetV.Px(),  MetV.Py() );
      hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      hdPhiMETMu[CutIndex]->Fill(dPhi,EvWeight);
      
      MT = TMath::Sqrt(2*tree_.muon_pt[mIndex]*MetV.Pt()*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hMTmu[CutIndex]->Fill(MT,EvWeight);

      hmet_dPhi[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);
      hmet_dPhimu[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);

      hmet_MT[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);
      hmet_MTmu[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);

      hMT_dPhi[CutIndex]->Fill(MT,dPhi,EvWeight);
      hMT_dPhimu[CutIndex]->Fill(MT,dPhi,EvWeight);
 
      

//	cout<<"  MUV SIZE ==================== "<<MuV.size()<<"  iCut "<<CutIndex<<"  mIndex  "<<mIndex<<" tIndex  "<<tIndex<<endl;

  }




  if (Sel=="eltau" && ElV.size()>0 && TauV.size()>0  && mIndex>-1 && eIndex > -1){
    //double Dr= deltaR(ElV.at(eIndex).Eta(), ElV.at(eIndex).Phi(),TauV.at(tIndex).Eta(),TauV.at(tIndex).Phi());
    double Dr = elV.DeltaR(tauV);
    hdR_eltau[CutIndex]->Fill(Dr,EvWeight);
     
    TLorentzVector DiL = elV + tauV;
    double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
    double MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
    hMTeltau[CutIndex]->Fill(MT,EvWeight);
    hInvMassElTau[CutIndex]->Fill(DiL.M(),EvWeight);
	
    hdRet_dPhi[CutIndex]->Fill(Dr,dPhi,EvWeight);
  }

  if (Sel=="muel" && MuV.size()>0 && ElV.size()>0  && mIndex>-1 && eIndex > -1 ){
    //Flaat_t Dr= deltaR(MuV.at(mIndex).Eta(), MuV.at(mIndex).Phi(),ElV.at(eIndex).Eta(),ElV.at(eIndex).Phi());
    double Dr = ElV.at(0).DeltaR(MuV.at(0));
    hdR_muel[CutIndex]->Fill(Dr,EvWeight);
    TLorentzVector DiL = ElV.at(0) + MuV.at(0);
    double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
    double MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
    hMTmuel[CutIndex]->Fill(MT,EvWeight);
    hInvMassMuEl[CutIndex]->Fill(DiL.M(),EvWeight);
    hdRme_dPhi[CutIndex]->Fill(Dr,dPhi,EvWeight);
  }
  //////////////// FIX THE INDEXES FOR TWO TAU CANDIDATES!!!!
  if (Sel=="tautau" &&  TauV.size()>1 && tIndex > -1){
    double Dr = TauV.at(0).DeltaR(TauV.at(0));
    //double Dr= deltaR(TauV.at(mIndex).Eta(), TauV.at(mIndex).Phi(),TauV.at(tIndex).Eta(),TauV.at(tIndex).Phi());
    hdR_tautau[CutIndex]->Fill(Dr,EvWeight);
    TLorentzVector DiL = TauV.at(0) + TauV.at(0);
    double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
    double MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
    hMTtautau[CutIndex]->Fill(MT,EvWeight);
    hdRtt_dPhi[CutIndex]->Fill(Dr,dPhi,EvWeight);
    if (mIndex >-1 && tIndex >-1){
      TLorentzVector dimass = TauV.at(mIndex) + TauV.at(tIndex);
      hInvMassTauTau[CutIndex]->Fill(dimass.M(),EvWeight);
    }
  }



  if (ElV.size() > 0 && eIndex>-1)
    {
	
        hElpt[CutIndex]->Fill(tree_.electron_pt[eIndex],EvWeight);
        hEleta[CutIndex]->Fill(tree_.electron_eta[eIndex],EvWeight);
        hel_dxy[CutIndex]->Fill(tree_.electron_dxy[eIndex],EvWeight);
        hel_dz[CutIndex]->Fill(tree_.electron_dz[eIndex],EvWeight);

      double dPhiJ=dPhiFrom2P( tree_.electron_px[eIndex], tree_.electron_py[eIndex], MetV.Px(),  MetV.Py() );
      hdPhiElMET[CutIndex]->Fill(dPhiJ,EvWeight);

      for (unsigned int ie=0;ie<ElV.size();ie++){

        sumElpT +=ElV.at(ie).Pt();
	hLeppt[CutIndex]->Fill(ElV.at(ie).Pt(),EvWeight);
	hLepeta[CutIndex]->Fill(ElV.at(ie).Eta(),EvWeight);

	hel_miniISO[CutIndex]->Fill(tree_.electron_miniISO[ie],EvWeight);
	double absIso= tree_.electron_r03_sumChargedHadronPt[ie]
	                    + max(tree_.electron_r03_sumNeutralHadronEt[ie] + tree_.electron_r03_sumPhotonEt[ie]
			    - 0.5 * tree_.electron_r03_sumPUPt[ie],0.0);

	double relIso = absIso/tree_.electron_pt[ie];

	hel_relISO[CutIndex]->Fill(relIso,EvWeight);
      }
	double absIso= tree_.electron_r03_sumChargedHadronPt[eIndex]
	                    + max(tree_.electron_r03_sumNeutralHadronEt[eIndex] + tree_.electron_r03_sumPhotonEt[eIndex]
			    - 0.5 * tree_.electron_r03_sumPUPt[eIndex],0.0);

	double relIsoL = absIso/tree_.electron_pt[eIndex];

      hel_relISOL[CutIndex]->Fill(relIsoL,EvWeight);
      hel_miniISOL[CutIndex]->Fill(tree_.electron_miniISO[mIndex],EvWeight);

	if (tIndex>-1){
      double q = tree_.electron_charge[eIndex] * tree_.muon_charge[mIndex];
      if (q>0) hIso_sign[CutIndex]->Fill(relIsoL,0.5,EvWeight);
      if (q<0) hIso_sign[CutIndex]->Fill(relIsoL,1.5,EvWeight);
	}

    }


  if (MuV.size() > 0 && tree_.muon_count > 0)
    {
	    if (mIndex<0) mIndex = 0;

        hMupt[CutIndex]->Fill(tree_.muon_pt[mIndex],EvWeight);
        hMueta[CutIndex]->Fill(tree_.muon_eta[mIndex],EvWeight);
        hmu_dxy[CutIndex]->Fill(tree_.muon_dxy[mIndex],EvWeight);
        hmu_dz[CutIndex]->Fill(tree_.muon_dz[mIndex],EvWeight);

      double dPhiJ=dPhiFrom2P( tree_.muon_px[mIndex], tree_.muon_py[mIndex], MetV.Px(),  MetV.Py() );
      hdPhiMuMET[CutIndex]->Fill(dPhiJ,EvWeight);

      for (unsigned int im=0;im<MuV.size();im++){

        sumMuonpT +=MuV.at(im).Pt();
	hLeppt[CutIndex]->Fill(MuV.at(im).Pt(),EvWeight);
	hLepeta[CutIndex]->Fill(MuV.at(im).Eta(),EvWeight);

	hmu_miniISO[CutIndex]->Fill(tree_.muon_miniISO[im],EvWeight);
	double absIso= tree_.muon_r03_sumChargedHadronPt[im]
	                    + max(tree_.muon_r03_sumNeutralHadronEt[im] + tree_.muon_r03_sumPhotonEt[im]
			    - 0.5 * tree_.muon_r03_sumPUPt[im],0.0);

	double relIso = absIso/tree_.muon_pt[im];

	hmu_relISO[CutIndex]->Fill(relIso,EvWeight);
      }

     double absIso= tree_.muon_r03_sumChargedHadronPt[mIndex]
	            + max(tree_.muon_r03_sumNeutralHadronEt[mIndex] + tree_.muon_r03_sumPhotonEt[mIndex]
	            - 0.5 * tree_.muon_r03_sumPUPt[mIndex],0.0);

      double relIsoL = absIso/tree_.muon_pt[mIndex];
      hmu_relISOL[CutIndex]->Fill(relIsoL,EvWeight);
      hmu_miniISOL[CutIndex]->Fill(tree_.muon_miniISO[mIndex],EvWeight);

      if (tIndex>-1 && tree_.tau_count > 0){
      double q = tree_.tau_charge[tIndex] * tree_.muon_charge[mIndex];
      if (q>0) hIso_sign[CutIndex]->Fill(relIsoL,0.5,EvWeight);
      if (q<0) hIso_sign[CutIndex]->Fill(relIsoL,1.5,EvWeight);
	}

    }


  if (TauV.size() > 0 && tree_.tau_count > 0)
    {
	    if (tIndex<0) tIndex = 0;
      hTaupt[CutIndex]->Fill(tree_.tau_pt[tIndex],EvWeight);
      hTaueta[CutIndex]->Fill(tree_.tau_eta[tIndex],EvWeight);
      htau_dxy[CutIndex]->Fill(tree_.tau_dxy[tIndex],EvWeight);
      htau_dz[CutIndex]->Fill(tree_.tau_dz[tIndex],EvWeight);

      htau_ISOL[CutIndex]->Fill(tree_.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex],EvWeight);
      double dPhiJ=dPhiFrom2P( tree_.tau_px[mIndex], tree_.tau_py[mIndex], MetV.Px(),  MetV.Py() );
      hdPhiTauMET[CutIndex]->Fill(dPhiJ,EvWeight);

      for (unsigned int it=0;it<TauV.size();it++){
      htau_ISO[CutIndex]->Fill(tree_.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex],EvWeight);
        sumTaupT +=TauV.at(it).Pt();
	hLeppt[CutIndex]->Fill(TauV.at(it).Pt(),EvWeight);
	hLepeta[CutIndex]->Fill(TauV.at(it).Eta(),EvWeight);

      }
    }


  int bjets=0;
  
  if (JetsV.size()>0){
  

    double dPhiJ=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), MetV.Px(),  MetV.Py() );
    hdPhiJ0MET[CutIndex]->Fill(dPhiJ,EvWeight);

    hPtJ0[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);

    if (JetsV.size()>1){
      hPtJ1[CutIndex]->Fill(JetsV.at(1).Pt(),EvWeight);
      hHT2[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(1).Px(), JetsV.at(1).Py(), MetV.Px(),  MetV.Py() );
      hdPhiJ1MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }
    if (JetsV.size()>2){
      hPtJ2[CutIndex]->Fill(JetsV.at(2).Pt(),EvWeight);
      hHT3[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(2).Px(), JetsV.at(2).Py(), MetV.Px(),  MetV.Py() );
      hdPhiJ2MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }
    if (JetsV.size()>3){
      hPtJ3[CutIndex]->Fill(JetsV.at(3).Pt(),EvWeight);
      hHT4[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt()+JetsV.at(3).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(3).Px(), JetsV.at(3).Py(), MetV.Px(),  MetV.Py() );
      hdPhiJ3MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }

    for (unsigned int ij=0;ij<JetsV.size();ij++){
      sumpT+=JetsV.at(ij).Pt();
      double dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), MetV.Px(),  MetV.Py() );
      hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);


      for (unsigned int ib = 0; ib <tree_.pfjet_count;ib++){
        //if (double(tree_.pfjet_pt[ib]) == double(JetsV.at(ij).Pt())  &&  tree_.pfjet_btag[ib][0]  > 0.89) {
         float comp1_ = tree_.pfjet_pt[ib];
	  float comp2_ =JetsV.at(ij).Pt();
	      if (    comp1_ == comp2_ && tree_.pfjet_btag[ib][0]  > 0.89)   {
	  bjets++;
	  //for (unsigned bb=0;bb<9;bb++)
	//   cout<<tree_.pfjet_pt[ib] <<"  "<<JetsV.at(ij).Pt()<<"  "<< tree_.pfjet_btag[ib][0]  <<"  "<<bjets<<"  "<<tree_.pfjet_pt[ib]/JetsV.at(ij).Pt()<<endl;
     	}
      }
    hnBJet[CutIndex]->Fill(bjets,EvWeight);
    }
      
    hHT[CutIndex]->Fill(sumpT,EvWeight);
    double HText =sumpT;
    if (Sel=="mutau" && mIndex !=-1)    HText += tree_.muon_pt[mIndex];
    if (Sel=="eltau" && eIndex !=-1)    HText += tree_.electron_pt[eIndex];
    hHText[CutIndex]->Fill(HText,EvWeight);
  if (Sel=="mutau" && mIndex !=-1)   hRht[CutIndex]->Fill(tree_.muon_pt[mIndex]/HText,EvWeight);
  if (Sel=="eltau" && eIndex !=-1)   hRht[CutIndex]->Fill(tree_.electron_pt[eIndex]/HText,EvWeight);

  }
  hMET[CutIndex]->Fill(MetV.Pt(),EvWeight);

  hMeffMuon[CutIndex]->Fill(sumMuonpT+sumpT,EvWeight);
  hMeffEl[CutIndex]->Fill(sumElpT+sumpT,EvWeight);
  hMeffTau[CutIndex]->Fill(sumTaupT+sumpT,EvWeight);
  double mett = MetV.Pt();
  if (MetV.Pt()>0 ) {
    hHTOsqrMET[CutIndex]->Fill(  sumpT/sqrt(MetV.Pt()),EvWeight);
    hMeffMuonOsqrMET[CutIndex]->Fill( (sumMuonpT+sumpT)/sqrt(MetV.Pt()),EvWeight);
    hMeffElOsqrMET[CutIndex]->Fill( (sumElpT+sumpT)/sqrt(MetV.Pt()),EvWeight);
    hMeffTauOsqrMET[CutIndex]->Fill( (sumTaupT+sumpT)/sqrt(MetV.Pt()),EvWeight);
  }

  else {
    hHTOsqrMET[CutIndex]->Fill(  -1,EvWeight);
    hMeffMuonOsqrMET[CutIndex]->Fill(-1,EvWeight);
    hMeffElOsqrMET[CutIndex]->Fill( -1,EvWeight);
    hMeffTauOsqrMET[CutIndex]->Fill( -1,EvWeight);
  }

  // hMT[CutIndex]->Fill(MT,EvWeight);
  //  hDZeta->Fill(DZeta,EvWeight);
}








