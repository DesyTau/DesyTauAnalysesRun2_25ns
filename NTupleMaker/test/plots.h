
#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/functionsSUSY.h"
#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/lester_mt2_bisect.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/Basic_Mt2_332_Calculator.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/Basic_Nt2_332_Calculator.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/mTBound.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/mTTrue.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/TMctLib.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/mctlib.h"
//#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/Basic_MPairProd_Calculator.h"
#include "TTree.h"
#include <iostream>

//#include <ScaleFactor.h>
#include <TLorentzVector.h>
#include <TVector2.h>
using namespace std;


const  int CutN=25;



unsigned int RunMin = 9999999;
unsigned int RunMax = 0;
     

//bool isData = false;


double ChiMass=0;
double mIntermediate = tauMass;
double sumpT = 0;
double XSec=-1;
double xs,fact,fact2;
 
vector<TLorentzVector> AllJets_Lepton_noMet;
vector<TLorentzVector> JetsMV;
vector<TLorentzVector>  ElMV;
vector<TLorentzVector>  MuMV;
vector<TLorentzVector>  TauMV;
vector<TLorentzVector>  LeptMV;




double HIP_SF(double pt,double eta)

{

	TString inputRootFile = "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/data/HIP_mediumID.root";
TFile * fileIn = new TFile(inputRootFile, "read");
	if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };


TIter nextkey (fileIn->GetListOfKeys ());
TKey *key = 0;
key = (TKey *) nextkey ();
TObject *obj = key->ReadObj ();
//cout << "2D histos name for SF = " << obj->GetName() << endl;
TH2D *histo = (TH2D*) obj;
//TH2D *histo = (TH2D*)fileIn->Get("histo2D");
double ptN = histo->GetXaxis()->FindBin(pt);
double etaN = histo->GetYaxis()->FindBin(eta);
double result = histo->GetBinContent(ptN,etaN);
delete histo;
delete fileIn;
return result;
}





double TauFakeRate(float pt,float eta, string sel,string working_point){
float SF = 1.;


//80x MVAid

if ( working_point == "MVA"){

/*
if ( sel =="mutau"){

if (  fabs(eta) < 0.9 )
        {
                if (pt>20 && pt<30) SF = 1.07968;
                if (pt>30 && pt<40) SF = 0.847355;
                if (pt>40 ) SF = 0.833833;
        }
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 )
        {

                if (pt>20 && pt<30) SF = 1.03273;
                if (pt>30 && pt<40) SF = 0.943387;
                if (pt>40 ) SF = 1.07113;
        }

if (  fabs(eta) > 1.2 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.08635;
                if (pt>30 && pt<40) SF = 1.13754;
                if (pt>40) SF = 1.0387;
        }
if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
        {

                if (pt>20 && pt<30) SF = 0.977868;
                if (pt>30 && pt<40) SF = 0.974665;
                if (pt>40) SF = 0.902185;
        }

}//mutau MVA
*/

if ( sel =="mutau"){

if (  fabs(eta) < 0.8 )
        {
                if (pt>20 && pt<30) SF = 1.06985;
                if (pt>30 && pt<40) SF = 0.916968;
                if (pt>40 ) SF = 0.811926;
        }
if (  fabs(eta) > 0.9 && fabs(eta) < 1.44 )
        {

                if (pt>20 && pt<30) SF = 1.10921;
                if (pt>30 && pt<40) SF = 0.88544;
                if (pt>40 ) SF = 1.08596;
        }

if (  fabs(eta) > 1.44 && fabs(eta) < 1.566 )
        {

                if (pt>20 && pt<30) SF = 1.28118;
                if (pt>30 && pt<40) SF = 1.1486;
                if (pt>40) SF = 0.91935;
        }
if (  fabs(eta) > 1.566 && fabs(eta) < 2.3 )
        {

                if (pt>20 && pt<30) SF = 0.991456;
                if (pt>30 && pt<40) SF = 1.06501;
                if (pt>40) SF = 0.976065;
        }

}//mutau MVA





if ( sel =="eltau"){

if (  fabs(eta) < 0.8 )
        {
                if (pt>20 && pt<30) SF = 1.06051;
                if (pt>30 && pt<40) SF = 0.920655;
                if (pt>40 ) SF = 0.804852;
        }
if (  fabs(eta) > 0.8 && fabs(eta) < 1.44 )
        {

                if (pt>20 && pt<30) SF = 1.03031;
                if (pt>30 && pt<40) SF = 1.00293;
                if (pt>40 ) SF = 0.87865;
        }

if (  fabs(eta) > 1.44 && fabs(eta) < 1.566 )
	{


		if (pt>20 && pt<30) SF = 1.2041;
		if (pt>30 && pt<40) SF = 0.516369;
		if (pt>40) SF = 0.695692;
	}
if (  fabs(eta) > 1.566 && fabs(eta) < 2.3 )
	{


		if (pt>20 && pt<30) SF = 0.932817;
		if (pt>30 && pt<40) SF = 1.0093;
		if (pt>40) SF = 0.871699;
	}
}//eltau MVA


/*
///mutau for 15invfb
if ( sel =="mutau"){

if (  fabs(eta) < 0.9 )
        {
                if (pt>20 && pt<30) SF = 1.17285;
                if (pt>30 && pt<40) SF = 1.12899;
                if (pt>40 ) SF = 1.19492;
        }
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 )
        {

                if (pt>20 && pt<30) SF = 1.11631;
                if (pt>30 && pt<40) SF = 1.04866;
                if (pt>40 ) SF = 1.36044;
        }

if (  fabs(eta) > 1.2 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.08863;
                if (pt>30 && pt<40) SF = 1.13749;
                if (pt>40) SF = 1.21236;
        }
if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
        {

                if (pt>20 && pt<30) SF = 0.986986;
                if (pt>30 && pt<40) SF = 0.982437;
                if (pt>40) SF = 0.995135;
        }

}//mutau MVA
//////////////
//


if ( sel =="eltau"){

if (  fabs(eta) < 0.8 )
        {
                if (pt>20 && pt<30) SF = 1.04111;
                if (pt>30 && pt<50) SF = 1.34817;
                if (pt>50 ) SF = 1.60973;
        }
if (  fabs(eta) > 0.8 && fabs(eta) < 1.44 )
        {

                if (pt>20 && pt<30) SF = 1.04158;
                if (pt>30 && pt<50) SF = 1.56632;
                if (pt>50 ) SF = 1.08038;
        }

if (  fabs(eta) > 1.44 && fabs(eta) < 1.566 )
	{


		if (pt>20 && pt<30) SF = 1.15617;
		if (pt>30 && pt<50) SF = 1.30554;
		if (pt>50) SF = 1.;
	}
if (  fabs(eta) > 1.566 && fabs(eta) < 2.1 )
	{


		if (pt>20 && pt<30) SF = 1.01365;
		if (pt>30 && pt<50) SF = 1.00284;
		if (pt>50) SF = 1.20248;
	}
}//eltau MVA

*/

/*
//eltau from 15invfb
if ( sel =="eltau"){

        {
                if (pt>20 && pt<30) SF = 1.08763;
                if (pt>30 && pt<50) SF = 1.52724;
                if (pt>50 ) SF = 1.70454;
        }
if (  fabs(eta) > 1.48 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.05136;
                if (pt>30 && pt<50) SF = 1.09525;
                if (pt>50 ) SF = 1.29595;
        }

if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
	{


		if (pt>20 && pt<30) SF = 0.752314;
		if (pt>30 && pt<50) SF = 1.41862;
		if (pt>50) SF = 1.26267;
	}
}//eltau MVA
if ( sel =="eltau"){

if (  fabs(eta) < 1.48 )
        {
                if (pt>20 && pt<30) SF = 1.08763;
                if (pt>30 && pt<50) SF = 1.52724;
                if (pt>50 ) SF = 1.70454;
        }
if (  fabs(eta) > 1.48 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.05136;
                if (pt>30 && pt<50) SF = 1.09525;
                if (pt>50 ) SF = 1.29595;
        }

if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
	{


		if (pt>20 && pt<30) SF = 0.752314;
		if (pt>30 && pt<50) SF = 1.41862;
		if (pt>50) SF = 1.26267;
	}
}//eltau MVA
if ( sel =="eltau"){

if (  fabs(eta) < 1.48 )
        {
                if (pt>20 && pt<30) SF = 1.08763;
                if (pt>30 && pt<50) SF = 1.52724;
                if (pt>50 ) SF = 1.70454;
        }
if (  fabs(eta) > 1.48 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.05136;
                if (pt>30 && pt<50) SF = 1.09525;
                if (pt>50 ) SF = 1.29595;
        }

if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
	{


		if (pt>20 && pt<30) SF = 0.752314;
		if (pt>30 && pt<50) SF = 1.41862;
		if (pt>50) SF = 1.26267;
	}
}//eltau MVA
if ( sel =="eltau"){

if (  fabs(eta) < 1.48 )
        {
                if (pt>20 && pt<30) SF = 1.08763;
                if (pt>30 && pt<50) SF = 1.52724;
                if (pt>50 ) SF = 1.70454;
        }
if (  fabs(eta) > 1.48 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.05136;
                if (pt>30 && pt<50) SF = 1.09525;
                if (pt>50 ) SF = 1.29595;
        }

if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
	{


		if (pt>20 && pt<30) SF = 0.752314;
		if (pt>30 && pt<50) SF = 1.41862;
		if (pt>50) SF = 1.26267;
	}
}//eltau MVA
*/
}




//CutBased
if (working_point == "CutBased" ){

	if (sel =="eltau"){

if (  fabs(eta) < 1.48 ) 
	{
		if (pt>20 && pt<30) SF = 1.125;
		if (pt>30 && pt<50) SF = 1.75138;
		if (pt>50 ) SF = 2.02406;
	}
if (  fabs(eta) > 1.48 && fabs(eta) < 2.1 ) 
	{

		if (pt>20 && pt<30) SF = 1.06778;
		if (pt>30 && pt<50) SF = 1.36354;
		if (pt>50 ) SF = 1.91359;
	}

if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 ) 
	{

		if (pt>20 && pt<30) SF = 0.856469;
		if (pt>30 && pt<50) SF = 1.28117;
		if (pt>60) SF = 1.80822;
	}
	
}//eltau CutBased


if (sel == "mutau"){

	if (  fabs(eta) < 0.9 )
        {
                if (pt>20 && pt<30) SF = 0.921993;
                if (pt>30 && pt<50) SF = 1.06514;
                if (pt>50 ) SF = 1.44677;
        }
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 )
        {

                if (pt>20 && pt<30) SF = 1.16935;
                if (pt>30 && pt<50) SF = 1.38829;
                if (pt>50 ) SF = 1.016;
        }

if (  fabs(eta) > 1.2 && fabs(eta) < 2.1 )
	{


		if (pt>20 && pt<30) SF = 1.0466;
		if (pt>30 && pt<50) SF = 1.10949;
		if (pt>50) SF = 1.265;
	}
if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
	{

		if (pt>20 && pt<30) SF = 1.1374;
		if (pt>30 && pt<50) SF = 0.859957;
		if (pt>50) SF = 0.584229;
	}

}//mutau CutBased


}//CutBased

//cout<<" working point "<<working_point<<"  "<<sel<<"  "<<pt<<"  "<<eta<<"  "<<SF<<endl;

return SF;

}

double TauiD(string sel,string working_point){

float SF = 1;

if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Loose" )
	SF = 0.81;

if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Medim" )
	SF = 0.79;
if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Tight" )
	SF = 0.79;


if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vloose" )
	SF = 0.83;
if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vloose" )
	SF= 0.84;
if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Medium" )
	SF = 0.84;

if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Tight" )
	SF = 0.83;

if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vtight" )
	SF = 0.80;



return SF;

}

float Dzeta(TLorentzVector LV, TLorentzVector muV, TLorentzVector MetV)
	{
				float LUnitX = LV.Px()/LV.Pt();
				float LUnitY = LV.Py()/LV.Pt();

				//	cout<<" CHECK =========== "<<tauV.Pt()<<"  "<<ta_pt[tIndex]<<endl;	
				float muonUnitX = muV.Px()/muV.Pt();
				float muonUnitY = muV.Py()/muV.Pt();

				float zetaX = LUnitX + muonUnitX;
				float zetaY = LUnitY + muonUnitY;

				float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

				zetaX = zetaX/normZeta;
				zetaY = zetaY/normZeta;

				float vectorX = MetV.Px() + muV.Px() + LV.Px();
				float vectorY = MetV.Py() + muV.Py() + LV.Py();

				float vectorVisX = muV.Px() + LV.Px();
				float vectorVisY = muV.Py() + LV.Py();

				// computation of DZeta variable
				// pfmet
				float pzetamiss = MetV.Px()*zetaX + MetV.Py()*zetaY;
				float pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
				float dzeta = pzetamiss - 0.85*pzetavis;


				return dzeta;
}




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
vector<string> var;
vector < string > vec;
double var_[1000];

vector<string> CutList;

TH1D * histRuns = new TH1D("histRuns","",6000,24000,30000);

TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);

TH1D * histTopPt = new TH1D("histTopPt","",1,-0.5,0.5);

TH1D * hWeights [CutN];
TH1D * hEventSign [CutN];

TH1D * hMudxy [CutN];
TH1D * hMudz [CutN];

TH1D * hMudxyerr [CutN];
TH1D * hMudzerr [CutN];
TH1D * hMuIPsigxy [CutN];
TH1D * hMuIPsigz [CutN];

TH1D * hLept1dxy [CutN];
TH1D * hLept1dz [CutN];

TH1D * hLept1dxyerr [CutN];
TH1D * hLept1dzerr [CutN];
TH1D * hLept1IPsigxy [CutN];
TH1D * hLept1IPsigz [CutN];

TH1D * hLept2dxy [CutN];
TH1D * hLept2dz [CutN];

TH1D * hLept2dxyerr [CutN];
TH1D * hLept2dzerr [CutN];
TH1D * hLept2IPsigxy [CutN];
TH1D * hLept2IPsigz [CutN];


TH1D * hEldxyerr [CutN];
TH1D * hEldzerr [CutN];

TH1D * hElIPsigxy [CutN];
TH1D * hElIPsigz [CutN];

TH1D *hElneutralHadIso[CutN]; 
TH1D *hElphotonIso[CutN];
TH1D *hElchargedHadIso[CutN];
TH1D *hElpuIso[CutN];
TH1D *hElneutralIso[CutN];
TH1D *hElabsIsoEl[CutN];
TH1D *hElrelIsoEl[CutN];


TH1D *hMuneutralHadIso[CutN]; 
TH1D *hMuphotonIso[CutN];
TH1D *hMuchargedHadIso[CutN];
TH1D *hMupuIso[CutN];
TH1D *hMuneutralIso[CutN];
TH1D *hMuabsIsoMu[CutN];
TH1D *hMurelIsoMu[CutN];

TH1D *hLept1neutralHadIso[CutN]; 
TH1D *hLept1photonIso[CutN];
TH1D *hLept1chargedHadIso[CutN];
TH1D *hLept1puIso[CutN];
TH1D *hLept1neutralIso[CutN];
TH1D *hLept1absIsoMu[CutN];
TH1D *hLept1relIsoMu[CutN];


TH1D *hLept2neutralHadIso[CutN]; 
TH1D *hLept2photonIso[CutN];
TH1D *hLept2chargedHadIso[CutN];
TH1D *hLept2puIso[CutN];
TH1D *hLept2neutralIso[CutN];
TH1D *hLept2absIsoMu[CutN];
TH1D *hLept2relIsoMu[CutN];



TH1D * hEldxy [CutN];
TH1D * hEldz [CutN];

TH1D * htau_dxy [CutN];
TH1D * htau_dz [CutN];

TH1D *hMeffMuon[CutN];
TH1D *hMeffEl[CutN];
TH1D *hMeff[CutN];
TH1D *hHTOsqrMET[CutN];
TH1D *hPtOHT[CutN];
TH1D *hMeffMuonOsqrMET[CutN];
TH1D *hMeffElOsqrMET[CutN];
TH1D *hMeffOsqrMET[CutN];

TH1D *hHT[CutN];
TH1D *hHText[CutN];
TH1D *hRht[CutN];
TH1D *hHT2[CutN];
TH1D *hHT3[CutN];
TH1D *hHT4[CutN];
//TH1D *hST[CutN];
//TH1D *h0JetpT[CutN];
TH1D *hnJet[CutN];
TH1D *hnpartons[CutN];
TH1D *hnBJet[CutN];

TH1D *hCentrality[CutN];

TH1D *hPtJ0[CutN];
TH1D *hPtJ1[CutN];
TH1D *hPtJ2[CutN];
TH1D *hPtJ3[CutN];

TH1D *hIsoMu[CutN];
TH1D *hIsoLept1[CutN];
TH1D *hIsoLept2[CutN];
TH1D *hIsoEl[CutN];
TH1D *hIsoTau[CutN];

TH1D *hInvMassMuTau[CutN];
TH1D *hInvMassMuEl[CutN];
TH1D *hInvMassTauTau[CutN];
TH1D *hInvMassElEl[CutN];
TH1D *hInvMassElTau[CutN];
TH1D *hInvMassMuMu[CutN];


TH1D *hnEl[CutN];
TH1D *hElpt[CutN];
TH1D *hEleta[CutN];


TH1D *hnMu[CutN];
TH1D *hMupt[CutN];
TH1D *hMueta[CutN];
TH1D *hLept1pt[CutN];
TH1D *hLept1eta[CutN];
TH1D *hLept2pt[CutN];
TH1D *hLept2eta[CutN];



TH1D *hnTau[CutN];
TH1D *hTaupt[CutN];
TH1D *hTaueta[CutN];


TH1D *hMET[CutN];
TH1D *hMETphi[CutN];
//TH1D *hnOver[CutN];
TH1D *hdPhiMETMu[CutN];
TH1D *hdPhiMETLept1[CutN];
TH1D *hdPhiMETLept2[CutN];
TH1D *hdPhiDil[CutN];
TH1D *hdEtaDil[CutN];
TH1D *hdEtaJ0J1[CutN];
TH1D *hdPhiMETEl[CutN];
TH1D *hdPhiMETTau[CutN];
TH1D *hdPhiJMET[CutN];
TH1D *hdPhiJ0J1[CutN];
TH1D *hdPhiJ0Tau[CutN];
TH1D *hdPhiJ0El[CutN];
TH1D *hdPhiJ0Mu[CutN];
TH1D *hdPhiJ0MET[CutN];
TH1D *hdPhiJ1MET[CutN];
TH1D *hdPhiJ2MET[CutN];
TH1D *hdPhiJ3MET[CutN];
TH1D *hdPhiMuMET[CutN];
TH1D *hdPhiLept1MET[CutN];
TH1D *hdPhiLept2MET[CutN];
TH1D *hdPhiElMET[CutN];
TH1D *hdPhiTauMET[CutN];
TH1D *hdPhiJ0Lept1[CutN];
TH1D *hdPhiJ0Lept2[CutN];

TH1D *hMT[CutN];
TH1D *hMTsum[CutN];
TH1D *hMTtot[CutN];
TH1D *hMTel[CutN];
TH1D *hMTmu[CutN];
TH1D *hMTtau[CutN];
TH1D *hDZeta[CutN];

TH1D *hMTmutau[CutN];
TH1D *hMTmuel[CutN];
TH1D *hMTeltau[CutN];
TH1D *hMTtautau[CutN];


TH1D *hMt2lesterDil[CutN];
TH1D *hMt2Dil[CutN];
TH1D *hMCTDil[CutN];
TH1D *hMCTbDil[CutN];
TH1D *hMTDil[CutN];
TH1D *hMTlept1[CutN];
TH1D *hMTlept2[CutN];
TH1D *hdR_Dil[CutN];
TH1D *hInvMassDil[CutN];
TH1D *hPtDil[CutN];
TH1D *hDiJetMass_J0J1[CutN];

TH1D *hMt2lestermutau[CutN];
TH1D *hMt2mutau[CutN];
TH1D *hMCTmutau[CutN];
TH1D *hMCTxmutau[CutN];
TH1D *hMCTymutau[CutN];
TH1D *hMCTbmutau[CutN];
TH1D *hMCTcor[CutN];


TH1D *hMt2lestermuel[CutN];
TH1D *hMt2muel[CutN];
TH1D *hMCTmuel[CutN];
TH1D *hMCTxmuel[CutN];
TH1D *hMCTymuel[CutN];
TH1D *hMCTbmuel[CutN];

TH1D *hMt2lestereltau[CutN];
TH1D *hMt2eltau[CutN];
TH1D *hMCTeltau[CutN];
TH1D *hMCTxeltau[CutN];
TH1D *hMCTyeltau[CutN];
TH1D *hMCTbeltau[CutN];

TH1D *hTBoundmutau[CutN];

TH1D *hdR_mutau[CutN];
TH1D *hdR_eltau[CutN];
TH1D *hdR_tautau[CutN];
TH1D *hdR_muel[CutN];
TH1D *hdR_taujet[CutN];


TH1D *hnpv[CutN];
TH1D *hnpu[CutN];
TH1D *hnrho[CutN];

TH2D *hmet_MT[CutN];
TH2D *hmet_MTsum[CutN];
TH2D *hmet_MTtot[CutN];
TH2D *hmet_MCTb[CutN];
TH2D *hmet_MT2lester[CutN];
TH2D *hmet_TauPt[CutN];
TH2D *hmet_DZeta[CutN];

TH2D *hMT_MTsum[CutN];
TH2D *hMT_MTtot[CutN];
TH2D *hMT_MCTb[CutN];
TH2D *hMT_MT2lester[CutN];
TH2D *hMT_TauPt[CutN];
TH2D *hMT_DZeta[CutN];

TH2D *hMTsum_MTtot[CutN];
TH2D *hMTsum_MCTb[CutN];
TH2D *hMTsum_MT2lester[CutN];
TH2D *hMTsum_TauPt[CutN];
TH2D *hMTsum_DZeta[CutN];

TH2D *hMTtot_MCTb[CutN];
TH2D *hMTtot_MT2lester[CutN];
TH2D *hMTtot_TauPt[CutN];
TH2D *hMTtot_DZeta[CutN];

TH2D *hMCTb_MT2lester[CutN];
TH2D *hMCTb_TauPt[CutN];
TH2D *hMCTb_DZeta[CutN];

TH2D *hMT2lester_TauPt[CutN];
TH2D *hMT2lester_DZeta[CutN];

TH2D *hTauPt_DZeta[CutN];


 
     



TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutN,1,CutN+1);
TH1D *CutFlowUnW= new TH1D("CutFlowUnW","Cut Flow",CutN,1,CutN+1);
TH1D *CutFlowUnWNorm= new TH1D("CutFlowUnWNorm","Cut Flow",CutN,1,CutN+1);

TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);


TH1D * muonPtAllH = new TH1D("muonPtAllH","",10,0,200);
TH1D * electronPtAllH = new TH1D("electronPtAllH","",10,0,200);
TH1D * TauPtAllH = new TH1D("TauPtAllH","",10,0,200);

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
/*
TLorentzVector ElV, MuV, TauV, JetsV, METV;

vector<TLorentzVector> AllJets_Lepton_noMet;
vector<TLorentzVector> JetsMV;
vector<TLorentzVector>  ElMV;
vector<TLorentzVector>  MuMV;
vector<TLorentzVector>  TauMV;
vector<TLorentzVector>  LeptMV;
*/

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



//string CutList[CutN];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNer){

for(int cj = 0; cj < CutNer; cj++)
    {
      CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
      CutFlowUnW->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
  //  }

//for(int cj = 0; cj < CutNer; cj++)
  //  {
      TString cutName=CutList[cj];
      TString nCut;
      nCut.Form("%i",cj);
      ///generic variables
      
      
	
      hMuneutralHadIso[cj] = new TH1D("mu_neutralHadIso_"+nCut,"mu_neutralHadIso "+cutName,150,0,300);
      hMuneutralHadIso[cj]->Sumw2();

      hMuphotonIso[cj] = new TH1D("mu_photonIso_"+nCut,"mu_photonIso "+cutName,150,0,300);
      hMuphotonIso[cj]->Sumw2();

      hMuchargedHadIso[cj] = new TH1D("mu_chargedHadIso_"+nCut,"mu_chargedHadIso "+cutName,150,0,300);
      hMuchargedHadIso[cj]->Sumw2();

      hMupuIso[cj] = new TH1D("mu_puIso_"+nCut,"mu_puIso "+cutName,200,0,200);
      hMupuIso[cj]->Sumw2();

      hMuneutralIso[cj] = new TH1D("mu_neutralIso_"+nCut,"mu_neutralIso "+cutName,150,0,300);
      hMuneutralIso[cj]->Sumw2();

      hMuabsIsoMu[cj] = new TH1D("mu_absIsoMu_"+nCut,"mu_absIsoMu "+cutName,150,0,300);
      hMuabsIsoMu[cj]->Sumw2();

      hMurelIsoMu[cj] = new TH1D("mu_relIsoMu_"+nCut,"mu_relIsoMu "+cutName,200,0,10);
      hMurelIsoMu[cj]->Sumw2();


      hLept1neutralHadIso[cj] = new TH1D("lept1_neutralHadIso_"+nCut,"lept1_neutralHadIso "+cutName,150,0,300);
      hLept1neutralHadIso[cj]->Sumw2();

      hLept1photonIso[cj] = new TH1D("lept1_photonIso_"+nCut,"lept1_photonIso "+cutName,150,0,300);
      hLept1photonIso[cj]->Sumw2();

      hLept1chargedHadIso[cj] = new TH1D("lept1_chargedHadIso_"+nCut,"lept1_chargedHadIso "+cutName,150,0,300);
      hLept1chargedHadIso[cj]->Sumw2();

      hLept1puIso[cj] = new TH1D("lept1_puIso_"+nCut,"lept1_puIso "+cutName,200,0,200);
      hLept1puIso[cj]->Sumw2();

      hLept1neutralIso[cj] = new TH1D("lept1_neutralIso_"+nCut,"lept1_neutralIso "+cutName,150,0,300);
      hLept1neutralIso[cj]->Sumw2();

      hLept1absIsoMu[cj] = new TH1D("lept1_absIsoMu_"+nCut,"lept1_absIsoMu "+cutName,150,0,300);
      hLept1absIsoMu[cj]->Sumw2();

      hLept1relIsoMu[cj] = new TH1D("lept1_relIsoMu_"+nCut,"lept1_relIsoMu "+cutName,200,0,10);
      hLept1relIsoMu[cj]->Sumw2();
/////////

      hLept2neutralHadIso[cj] = new TH1D("lept2_neutralHadIso_"+nCut,"lept2_neutralHadIso "+cutName,150,0,300);
      hLept2neutralHadIso[cj]->Sumw2();

      hLept2photonIso[cj] = new TH1D("lept2_photonIso_"+nCut,"lept2_photonIso "+cutName,150,0,300);
      hLept2photonIso[cj]->Sumw2();

      hLept2chargedHadIso[cj] = new TH1D("lept2_chargedHadIso_"+nCut,"lept2_chargedHadIso "+cutName,150,0,300);
      hLept2chargedHadIso[cj]->Sumw2();

      hLept2puIso[cj] = new TH1D("lept2_puIso_"+nCut,"lept2_puIso "+cutName,200,0,200);
      hLept2puIso[cj]->Sumw2();

      hLept2neutralIso[cj] = new TH1D("lept2_neutralIso_"+nCut,"lept2_neutralIso "+cutName,150,0,300);
      hLept2neutralIso[cj]->Sumw2();

      hLept2absIsoMu[cj] = new TH1D("lept2_absIsoMu_"+nCut,"lept2_absIsoMu "+cutName,150,0,300);
      hLept2absIsoMu[cj]->Sumw2();

      hLept2relIsoMu[cj] = new TH1D("lept2_relIsoMu_"+nCut,"lept2_relIsoMu "+cutName,200,0,10);
      hLept2relIsoMu[cj]->Sumw2();

///////////





      hElneutralHadIso[cj] = new TH1D("el_neutralHadIso_"+nCut,"el_neutralHadIso "+cutName,150,0,300);
      hElneutralHadIso[cj]->Sumw2();

      hElphotonIso[cj] = new TH1D("el_photonIso_"+nCut,"el_photonIso "+cutName,150,0,300);
      hElphotonIso[cj]->Sumw2();

      hElchargedHadIso[cj] = new TH1D("el_chargedHadIso_"+nCut,"el_chargedHadIso "+cutName,150,0,300);
      hElchargedHadIso[cj]->Sumw2();

      hElpuIso[cj] = new TH1D("el_puIso_"+nCut,"el_puIso "+cutName,200,0,200);
      hElpuIso[cj]->Sumw2();

      hElneutralIso[cj] = new TH1D("el_neutralIso_"+nCut,"el_neutralIso "+cutName,150,0,300);
      hElneutralIso[cj]->Sumw2();

      hElabsIsoEl[cj] = new TH1D("el_absIsoEl_"+nCut,"el_absIsoEl "+cutName,150,0,300);
      hElabsIsoEl[cj]->Sumw2();

      hElrelIsoEl[cj] = new TH1D("el_relIsoEl_"+nCut,"el_relIsoEl "+cutName,200,0,10);
      hElrelIsoEl[cj]->Sumw2();

      hMudxy[cj] = new TH1D ("mu_dxy_"+nCut,"mu_dxy "+cutName,20,-.2,.2);
      hMudxy[cj]->Sumw2();
      hMudxyerr[cj] = new TH1D ("mu_dxyerr_"+nCut,"mu_dxyerr "+cutName,20,-.2,.2);
      hMudxyerr[cj]->Sumw2();
      hLept1dxy[cj] = new TH1D ("lept1_dxy_"+nCut,"lept1_dxy "+cutName,20,-.2,.2);
      hLept1dxy[cj]->Sumw2();
      hLept1dxyerr[cj] = new TH1D ("lept1_dxyerr_"+nCut,"lept1_dxyerr "+cutName,20,-.2,.2);
      hLept1dxyerr[cj]->Sumw2();
      hLept2dxy[cj] = new TH1D ("lept2_dxy_"+nCut,"lept2_dxy "+cutName,20,-.2,.2);
      hLept2dxy[cj]->Sumw2();
      hLept2dxyerr[cj] = new TH1D ("lept2_dxyerr_"+nCut,"lept2_dxyerr "+cutName,20,-.2,.2);
      hLept2dxyerr[cj]->Sumw2();
 
      hEldxy[cj] = new TH1D ("el_dxy_"+nCut,"el_dxy "+cutName,20,-.2,.2);
      hEldxy[cj]->Sumw2();
      hEldxyerr[cj] = new TH1D ("el_dxyerr_"+nCut,"el_dxyerr "+cutName,20,-.2,.2);
      hEldxyerr[cj]->Sumw2();

      htau_dxy[cj] = new TH1D ("tau_dxy_"+nCut,"tau_dxy "+cutName,20,-.2,.2);
      htau_dxy[cj]->Sumw2();
      
      htau_dz[cj] = new TH1D ("tau_dz_"+nCut,"tau_dz "+cutName,20,-.2,.2);
      htau_dz[cj]->Sumw2();
      hMudz[cj] = new TH1D ("mu_dz_"+nCut,"mu_dz "+cutName,20,-.2,.2);
      hMudz[cj]->Sumw2();
      hLept1dz[cj] = new TH1D ("lept1_dz_"+nCut,"lept1_dz "+cutName,20,-.2,.2);
      hLept1dz[cj]->Sumw2();
      hLept2dz[cj] = new TH1D ("lept2_dz_"+nCut,"lept2_dz "+cutName,20,-.2,.2);
      hLept2dz[cj]->Sumw2();
      hEldz[cj] = new TH1D ("el_dz_"+nCut,"el_dz "+cutName,20,-.2,.2);
      hEldz[cj]->Sumw2();

      hMudzerr[cj] = new TH1D ("mu_dzerr_"+nCut,"mu_dzerr "+cutName,20,-.2,.2);
      hMudzerr[cj]->Sumw2();
      hLept1dzerr[cj] = new TH1D ("lept1_dzerr_"+nCut,"lept1_dzerr "+cutName,20,-.2,.2);
      hLept1dzerr[cj]->Sumw2();
      hLept2dzerr[cj] = new TH1D ("lept2_dzerr_"+nCut,"lept2_dzerr "+cutName,20,-.2,.2);
      hLept2dzerr[cj]->Sumw2();
      hEldzerr[cj] = new TH1D ("el_dzerr_"+nCut,"el_dzerr "+cutName,20,-.2,.2);
      hEldzerr[cj]->Sumw2();
	
      hMuIPsigxy[cj] = new TH1D ("muIPsigxy_"+nCut,"muIPsigxy "+cutName,50,0,10);
      hMuIPsigxy[cj]->Sumw2();
      hMuIPsigz[cj] = new TH1D ("muIPsigz_"+nCut,"muIPsigz "+cutName,50,0,10);
      hMuIPsigz[cj]->Sumw2();
      hLept1IPsigxy[cj] = new TH1D ("lept1IPsigxy_"+nCut,"lept1IPsigxy "+cutName,50,0,10);
      hLept1IPsigxy[cj]->Sumw2();
      hLept1IPsigz[cj] = new TH1D ("lept1IPsigz_"+nCut,"lept1IPsigz "+cutName,50,0,10);
      hLept1IPsigz[cj]->Sumw2();
      hLept2IPsigxy[cj] = new TH1D ("lept2IPsigxy_"+nCut,"lept2IPsigxy "+cutName,50,0,10);
      hLept2IPsigxy[cj]->Sumw2();
      hLept2IPsigz[cj] = new TH1D ("lept2IPsigz_"+nCut,"lept2IPsigz "+cutName,50,0,10);
      hLept2IPsigz[cj]->Sumw2();
 
      hElIPsigxy[cj] = new TH1D ("elIPsigxy_"+nCut,"elIPsigxy "+cutName,50,0,10);
      hElIPsigz[cj] = new TH1D ("elIPsigz_"+nCut,"elIPsigz "+cutName,50,0,10);

      hHTOsqrMET[cj] = new TH1D ("hHTOsqrMET_"+nCut,"hHTOsqrMET "+cutName,10,0.0,100.0);
      hHTOsqrMET[cj]->Sumw2();
      hPtOHT[cj] = new TH1D ("hPtOHT_"+nCut,"hPtOHT "+cutName,10,0.0,10.0);
      hPtOHT[cj]->Sumw2();
      
      hMeffMuonOsqrMET[cj] = new TH1D ("hMeffMuonOsqrMET_"+nCut,"hMeffMuonOsqrMET "+cutName,10,0.0,100.0);
      hMeffMuonOsqrMET[cj]->Sumw2();
      hMeffElOsqrMET[cj] = new TH1D ("hMeffElOsqrMET_"+nCut,"hMeffElOsqrMET "+cutName,10,0.0,100.0);
      hMeffElOsqrMET[cj]->Sumw2();
      hMeffOsqrMET[cj] = new TH1D ("hMeffOsqrMET_"+nCut,"hMeffOsqrMET "+cutName,10,0.0,100.0);
      hMeffOsqrMET[cj]->Sumw2();

      hMeffMuon[cj] = new TH1D ("hMeffMuon_"+nCut,"hMeffMuon "+cutName,50,0.0,500.0);
      hMeffMuon[cj]->Sumw2();
      hMeffEl[cj] = new TH1D ("hMeffEl_"+nCut,"hMeffEl "+cutName,50,0.0,500.0); 
      hMeffEl[cj]->Sumw2();
      hMeff[cj] = new TH1D ("hMeff_"+nCut,"hMeff "+cutName,50,0.0,500.0); 
      hMeff[cj]->Sumw2();

      hCentrality[cj]  = new TH1D ("hCentrality_"+nCut,"hCentrality "+cutName,20,0.,1.);
      hCentrality[cj]->Sumw2();

      hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,10,0.0,1000.0);
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
      hnpartons[cj] = new TH1D ("npartons_"+nCut,"npartons "+cutName,6,-0.5,5.5);
      hnpartons[cj]->Sumw2();
      hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,10,-0.5,9.5);
      hnBJet[cj]->Sumw2();

      hWeights[cj] = new TH1D ("hWeights_"+nCut,"hWeights "+cutName,10,-1,9);
      hWeights[cj]->Sumw2();
      hEventSign[cj] = new TH1D ("hEventSign_"+nCut,"hEventSign "+cutName,20,-10,10);
      hEventSign[cj]->Sumw2();
	
      hDZeta[cj] = new TH1D("hDZeta_"+nCut,"hDZeta"+cutName,60,-400,200);
      hDZeta[cj]->Sumw2();
      hInvMassMuTau[cj] = new TH1D ("hInvMassMuTau_"+nCut,"hInvMassMuTau "+cutName,50,0,500);
      hInvMassMuTau[cj]->Sumw2();
      hInvMassMuEl[cj] = new TH1D ("hInvMassMuEl_"+nCut,"hInvMassMuel "+cutName,50,0,500);
      hInvMassMuEl[cj]->Sumw2();
      hInvMassMuMu[cj] = new TH1D ("hInvMassMuMu_"+nCut,"hInvMassMuMu "+cutName,50,0,500);
      hInvMassMuMu[cj]->Sumw2();
      hInvMassElTau[cj] = new TH1D ("hInvMassElTau_"+nCut,"hInvMassElTau "+cutName,50,0,500);
      hInvMassElTau[cj]->Sumw2();
      hInvMassElEl[cj] = new TH1D ("hInvMassElEl_"+nCut,"hInvMassElEl "+cutName,50,0,500);
      hInvMassElEl[cj]->Sumw2();
      hInvMassTauTau[cj] = new TH1D ("hInvMassTauTau_"+nCut,"hInvMassTauTau "+cutName,50,0,500);
      hInvMassTauTau[cj]->Sumw2();
      hInvMassDil[cj] = new TH1D ("hInvMassDil_"+nCut,"hInvMassDil "+cutName,50,0,500);
      hInvMassDil[cj]->Sumw2();
      hDiJetMass_J0J1[cj] = new TH1D ("hDiJetMass_J0J1_"+nCut,"hDiJetMass_J0J1 "+cutName,50,0,500);
      hDiJetMass_J0J1[cj]->Sumw2();
        
      hPtDil[cj] = new TH1D ("PtDil_"+nCut,"PtDil "+cutName,80,0,800);
      hPtDil[cj]->Sumw2();
      //Leptons
      //
      //
        
      //Muons
      //
      //
      hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,-0.5,9.5);
      hnMu[cj]->Sumw2();
      hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,50,0,500);
      hMupt[cj]->Sumw2();
      hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,40,-4,4);
      hMueta[cj]->Sumw2();
      hLept1pt[cj] = new TH1D ("Lept1pT_"+nCut,"Lept1 pT "+cutName,50,0,500);
      hLept1pt[cj]->Sumw2();
      hLept1eta[cj] = new TH1D ("Lept1eta_"+nCut,"Lept1 eta "+cutName,40,-4,4);
      hLept1eta[cj]->Sumw2();
      hLept2pt[cj] = new TH1D ("Lept2pT_"+nCut,"Lept2 pT "+cutName,50,0,500);
      hLept2pt[cj]->Sumw2();
      hLept2eta[cj] = new TH1D ("Lept2eta_"+nCut,"Lept2 eta "+cutName,40,-4,4);
      hLept2eta[cj]->Sumw2();
        
      //Taus
      //
      //
      hnTau[cj] = new TH1D ("nTau_"+nCut,"nTau "+cutName,10,-0.5,9.5);
      hnTau[cj]->Sumw2();
      hTaupt[cj] = new TH1D ("TaupT_"+nCut,"Tau pT "+cutName,50,0,500);
      hTaupt[cj]->Sumw2();
      hTaueta[cj] = new TH1D ("Taueta_"+nCut,"Tau eta "+cutName,40,-4,4);
      hTaueta[cj]->Sumw2();
	
      //hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
      //Electrons
      //
      //
      hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,-0.5,9.5);
      hnEl[cj]->Sumw2();
      hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,50,0,500);
      hElpt[cj]->Sumw2();
      hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,40,-4,4);
      hEleta[cj]->Sumw2();
       
      hIsoMu[cj] = new TH1D ("IsoMu_"+nCut,"Mu Iso "+cutName,100,0,5);
      hIsoMu[cj]->Sumw2();
      hIsoLept1[cj] = new TH1D ("IsoLept1_"+nCut,"Lept1 Iso "+cutName,100,0,5);
      hIsoLept1[cj]->Sumw2();
      hIsoLept2[cj] = new TH1D ("IsoLept2_"+nCut,"Lept2 Iso "+cutName,100,0,5);
      hIsoLept2[cj]->Sumw2();
      hIsoEl[cj] = new TH1D ("IsoEl_"+nCut,"El Iso "+cutName,100,0,5);
      hIsoEl[cj]->Sumw2();
      hIsoTau[cj] = new TH1D ("IsoTau_"+nCut,"Tau Iso "+cutName,100,0,5);
      hIsoTau[cj]->Sumw2();
       
      hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,40.0,0.0,800.0);
      hMET[cj]->Sumw2();
	

      //dPhi
      //
      //
      hMETphi[cj] = new TH1D("METphi_"+nCut,"METphi "+cutName,64,0.0,3.2);
      hMETphi[cj]->Sumw2();

      hdEtaDil[cj] = new TH1D("dEtaDil_"+nCut,"dEtaDil "+cutName,100,-10.,10);
      hdEtaDil[cj]->Sumw2();
      hdEtaJ0J1[cj] = new TH1D("dEtaJ0J1_"+nCut,"dEtaJ0J1 "+cutName,100,-10.,10);
      hdEtaJ0J1[cj]->Sumw2();
    
      hdPhiDil[cj] = new TH1D("dPhiDil_"+nCut,"dPhiDil "+cutName,64,0.0,3.2);
      hdPhiDil[cj]->Sumw2();
      hdPhiJ0Lept1[cj] = new TH1D("dPhiJ0Lept1_"+nCut,"dPhiJ0Lept1 "+cutName,64,0.0,3.2);
      hdPhiJ0Lept1[cj]->Sumw2();
      hdPhiJ0Lept2[cj] = new TH1D("dPhiJ0Lept2_"+nCut,"dPhiJ0Lept2 "+cutName,64,0.0,3.2);
      hdPhiJ0Lept2[cj]->Sumw2();
      hdPhiMETMu[cj] = new TH1D("dPhiMETMu_"+nCut,"dPhiMETMu "+cutName,64,0.0,3.2);
      hdPhiMETMu[cj]->Sumw2();
      hdPhiMETEl[cj] = new TH1D("dPhiMETEl_"+nCut,"dPhiMETEl "+cutName,64,0.0,3.2);
      hdPhiMETEl[cj]->Sumw2();
      hdPhiMETTau[cj] = new TH1D("dPhiMETTau_"+nCut,"dPhiMETTau "+cutName,64,0.0,3.2);
      hdPhiMETTau[cj]->Sumw2();
      hdPhiMETLept1[cj] = new TH1D("dPhiMETLept1_"+nCut,"dPhiMETLept1 "+cutName,64,0.0,3.2);
      hdPhiMETLept1[cj]->Sumw2();
      hdPhiMETLept2[cj] = new TH1D("dPhiMETLept2_"+nCut,"dPhiMETLept2 "+cutName,64,0.0,3.2);
      hdPhiMETLept2[cj]->Sumw2();
    
      hdPhiJMET[cj] = new TH1D("dPhiJMET_"+nCut,"dPhiJMET "+cutName,64,0.0,3.2);
      hdPhiJMET[cj]->Sumw2();
      hdPhiJ0MET[cj] = new TH1D("dPhiJ0MET_"+nCut,"dPhiJ0MET "+cutName,64,0.0,3.2);
      hdPhiJ0MET[cj]->Sumw2();
      hdPhiJ0J1[cj] = new TH1D("dPhiJ0J1_"+nCut,"dPhiJ0J1 "+cutName,64,0.0,3.2);
      hdPhiJ0J1[cj]->Sumw2();
      hdPhiJ1MET[cj] = new TH1D("dPhiJ1MET_"+nCut,"dPhiJ1MET "+cutName,64,0.0,3.2);
      hdPhiJ1MET[cj]->Sumw2();
      hdPhiJ2MET[cj] = new TH1D("dPhiJ2MET_"+nCut,"dPhiJ2MET "+cutName,64,0.0,3.2);
      hdPhiJ2MET[cj]->Sumw2();
      hdPhiJ3MET[cj] = new TH1D("dPhiJ3MET_"+nCut,"dPhiJ3MET "+cutName,64,0.0,3.2);
      hdPhiJ3MET[cj]->Sumw2();

      hdPhiJ0Tau[cj] = new TH1D("dPhiJ0Tau_"+nCut,"dPhiJ0Tau "+cutName,64,0.0,3.2);
      hdPhiJ0Tau[cj]->Sumw2();
      hdPhiJ0Mu[cj] = new TH1D("dPhiJ0Mu_"+nCut,"dPhiJ0Mu "+cutName,64,0.0,3.2);
      hdPhiJ0Mu[cj]->Sumw2();
      hdPhiJ0El[cj] = new TH1D("dPhiJ0El_"+nCut,"dPhiJ0El "+cutName,64,0.0,3.2);
      hdPhiJ0El[cj]->Sumw2();

      hdPhiMuMET[cj] = new TH1D("dPhiMuMET_"+nCut,"dPhiMuMET "+cutName,64,0.0,3.2);
      hdPhiMuMET[cj]->Sumw2();
      hdPhiLept1MET[cj] = new TH1D("dPhiLept1MET_"+nCut,"dPhiLept1MET "+cutName,64,0.0,3.2);
      hdPhiLept1MET[cj]->Sumw2();
      hdPhiLept2MET[cj] = new TH1D("dPhiLept2MET_"+nCut,"dPhiLept2MET "+cutName,64,0.0,3.2);
      hdPhiLept2MET[cj]->Sumw2();

      hdPhiElMET[cj] = new TH1D("dPhiElMET_"+nCut,"dPhiElMET "+cutName,64,0.0,3.2);
      hdPhiElMET[cj]->Sumw2();
   
      hdPhiTauMET[cj] = new TH1D("dPhiTauMET_"+nCut,"dPhiTauMET "+cutName,64,0.0,3.2);
      hdPhiTauMET[cj]->Sumw2();
   
      //MT
      //
      //
      hMTsum[cj] = new TH1D ("MTsum_"+nCut,"MTsum "+cutName,30,0,300);
      hMTsum[cj]->Sumw2();
      hMTtot[cj] = new TH1D ("MTtot_"+nCut,"MTtot "+cutName,30,0,300);
      hMTtot[cj]->Sumw2();
      hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,30,0,300);
      hMT[cj]->Sumw2();
      hMTel[cj] = new TH1D ("MTel_"+nCut,"MTel "+cutName,30,0,300);
      hMTel[cj]->Sumw2();
      hMTmu[cj] = new TH1D ("MTmu_"+nCut,"MTmu "+cutName,30,0,300);
      hMTmu[cj]->Sumw2();

      hMTtau[cj] = new TH1D ("MTtau_"+nCut,"MTtau "+cutName,30,0,300);
      hMTtau[cj]->Sumw2();

      hMTmutau[cj] = new TH1D ("MTmutau_"+nCut,"MTmutau "+cutName,30,0,300);
      hMTmutau[cj]->Sumw2();
     
      hMTmuel[cj] = new TH1D ("MTmuel_"+nCut,"MTmuel "+cutName,30,0,300);
      hMTmuel[cj]->Sumw2();
       
      hMTeltau[cj] = new TH1D ("MTeltau_"+nCut,"MTeltau "+cutName,30,0,300);
      hMTeltau[cj]->Sumw2();
      
      hMTtautau[cj] = new TH1D ("MTtautau_"+nCut,"MTtautau "+cutName,50,0,1000);
      hMTtautau[cj]->Sumw2();


      hMTDil[cj] = new TH1D ("MTDil_"+nCut,"MTDil "+cutName,30,0,300);
      hMTDil[cj]->Sumw2();
      hMTlept1[cj] = new TH1D ("MTlept1_"+nCut,"MTlept1 "+cutName,30,0,300);
      hMTlept1[cj]->Sumw2();
      hMTlept2[cj] = new TH1D ("MTlept2_"+nCut,"MTlept2 "+cutName,30,0,300);
      hMTlept2[cj]->Sumw2();






      hMt2mutau[cj] = new TH1D ("Mt2mutau_"+nCut,"Mt2mutau "+cutName,50,0,1000);
      hMt2mutau[cj]->Sumw2();
      
      hMt2lestermutau[cj] = new TH1D ("Mt2lestermutau_"+nCut,"Mt2lestermutau "+cutName,30,0,600);
      hMt2lestermutau[cj]->Sumw2();

      hMCTcor[cj] = new TH1D ("MCTcor_"+nCut,"MCTcor "+cutName,50,0,500);
      hMCTcor[cj]->Sumw2();
      hMCTmutau[cj] = new TH1D ("MCTmutau_"+nCut,"MCTmutau "+cutName,50,0,500);
      hMCTmutau[cj]->Sumw2();
      hMCTxmutau[cj] = new TH1D ("MCTxmutau_"+nCut,"MCTxmutau "+cutName,50,0,500);
      hMCTxmutau[cj]->Sumw2();
      hMCTymutau[cj] = new TH1D ("MCTymutau_"+nCut,"MCTymutau "+cutName,50,0,500);
      hMCTymutau[cj]->Sumw2();
      
      hMCTbmutau[cj] = new TH1D ("MCTbmutau_"+nCut,"MCTbmutau "+cutName,50,0,500);
      hMCTbmutau[cj]->Sumw2();

      hMt2muel[cj] = new TH1D ("Mt2muel_"+nCut,"Mt2muel "+cutName,50,0,1000);
      hMt2muel[cj]->Sumw2();
      

      hMt2lesterDil[cj] = new TH1D ("Mt2lesterDil_"+nCut,"Mt2lesterDil "+cutName,30,0,600);
      hMt2lesterDil[cj]->Sumw2();

      hMCTDil[cj] = new TH1D ("MCTDil_"+nCut,"MCTDil "+cutName,50,0,500);
      hMCTDil[cj]->Sumw2();
      
      hMCTbDil[cj] = new TH1D ("MCTbDil_"+nCut,"MCTbDil "+cutName,50,0,500);
      hMCTbDil[cj]->Sumw2();

      hMt2Dil[cj] = new TH1D ("Mt2Dil_"+nCut,"Mt2Dil "+cutName,50,0,1000);
      hMt2Dil[cj]->Sumw2();




      hMt2lestermuel[cj] = new TH1D ("Mt2lestermuel_"+nCut,"Mt2lestermuel "+cutName,30,0,600);
      hMt2lestermuel[cj]->Sumw2();

      hMCTmuel[cj] = new TH1D ("MCTmuel_"+nCut,"MCTmuel "+cutName,50,0,500);
      hMCTmuel[cj]->Sumw2();
      hMCTxmuel[cj] = new TH1D ("MCTxmuel_"+nCut,"MCTxmuel "+cutName,50,0,500);
      hMCTxmuel[cj]->Sumw2();
      hMCTymuel[cj] = new TH1D ("MCTymuel_"+nCut,"MCTymuel "+cutName,50,0,500);
      hMCTymuel[cj]->Sumw2();
      
      hMCTbmuel[cj] = new TH1D ("MCTbmuel_"+nCut,"MCTbmuel "+cutName,50,0,500);
      hMCTbmuel[cj]->Sumw2();



      hTBoundmutau[cj] = new TH1D ("TBoundmutau_"+nCut,"TBoundmutau "+cutName,100,0,2000);
      hTBoundmutau[cj] ->Sumw2();
      hMt2eltau[cj] = new TH1D ("Mt2eltau_"+nCut,"Mt2eltau "+cutName,50,0,500);
      hMt2eltau[cj]->Sumw2();
      hMt2lestereltau[cj] = new TH1D ("Mt2lestereltau_"+nCut,"Mt2lestereltau "+cutName,30,0,600);
      hMt2lestereltau[cj]->Sumw2();
 
      hMCTeltau[cj] = new TH1D ("MCTeltau_"+nCut,"MCTeltau "+cutName,50,0,500);
      hMCTeltau[cj]->Sumw2();
      hMCTxeltau[cj] = new TH1D ("MCTxeltau_"+nCut,"MCTxeltau "+cutName,50,0,500);
      hMCTxeltau[cj]->Sumw2();
      hMCTyeltau[cj] = new TH1D ("MCTyeltau_"+nCut,"MCTyeltau "+cutName,50,0,500);
      hMCTyeltau[cj]->Sumw2();
      
      hMCTbeltau[cj] = new TH1D ("MCTbeltau_"+nCut,"MCTbeltau "+cutName,50,0,500);
      hMCTbeltau[cj]->Sumw2();

 
      hdR_eltau[cj]= new TH1D ("dR_eltau_"+nCut,"dR_eltau "+cutName,60,0,6);;
      hdR_eltau[cj]->Sumw2();
        
      hdR_mutau[cj]= new TH1D ("dR_mutau_"+nCut,"dR_mutau "+cutName,60,0,6);;
      hdR_mutau[cj]->Sumw2();

      hdR_taujet[cj]= new TH1D ("dR_taujet_"+nCut,"dR_taujet "+cutName,60,0,6);;
      hdR_taujet[cj]->Sumw2();

      hdR_tautau[cj]= new TH1D ("dR_tautau_"+nCut,"dR_tautau "+cutName,60,0,6);;
      hdR_tautau[cj]->Sumw2();
	
      hdR_muel[cj]= new TH1D ("dR_muel_"+nCut,"dR_muel "+cutName,60,0,6);;
      hdR_muel[cj]->Sumw2();
      hdR_Dil[cj]= new TH1D ("dR_Dil_"+nCut,"dR_Dil "+cutName,60,0,6);;
      hdR_Dil[cj]->Sumw2();

      hnpv[cj]= new TH1D ("npv_"+nCut,"npv "+cutName,100,-0.5,99.5);;
      hnpv[cj]->Sumw2();
      hnpu[cj]= new TH1D ("npu_"+nCut,"npu "+cutName,100,-0.5,99.5);;
      hnpu[cj]->Sumw2();
      hnrho[cj]= new TH1D ("nrho_"+nCut,"nrho "+cutName,100,-0.5,99.5);;
      hnrho[cj]->Sumw2();
 

   int nBinsTauPt = 4;
   double binsTauPt[5] = {0, 40, 80,120,1000};

   int nBinsmet = 4;
   double binsmet[5] = {0, 40, 80,120,1000};

   int nBinsMTsum = 4;
   double binsMTsum[5] = {0, 40,120,260,1000};

   int nBinsMTtot = 4;
   double binsMTtot[5] = {0, 40,120,200,1000};

   int nBinsMCTb = 3;
   double binsMCTb[4] = {0, 50,100,1000};

   int nBinsMT = 4;
   double binsMT[5] = {0, 40,120,160,1000};

   int nBinsDZeta = 4;
   double binsDZeta[5] = {-500, -150,-100,50,1000};

   int nBinsMT2lester = 3;
   double binsMT2lester[4] = {0, 40,80,1000};


   //met 


      hmet_MT[cj] = new TH2D ("met_MT_"+nCut,"met_MT "+cutName, nBinsmet, binsmet ,nBinsMT,binsMT);
      hmet_MT[cj]->Sumw2();

      hmet_MTsum[cj] = new TH2D ("met_MTsum_"+nCut,"met_MTsum "+cutName, nBinsmet, binsmet , nBinsMTsum,binsMTsum);
      hmet_MTsum[cj]->Sumw2();

      hmet_MTtot[cj] = new TH2D ("met_MTtot_"+nCut,"met_MTtot "+cutName, nBinsmet, binsmet , nBinsMTtot,binsMTtot);
      hmet_MTtot[cj]->Sumw2();

      hmet_MCTb[cj] = new TH2D ("met_MCTb_"+nCut,"met_MCTb "+cutName, nBinsmet, binsmet , nBinsMCTb,binsMCTb);
      hmet_MCTb[cj] ->Sumw2();

      hmet_MT2lester[cj] = new TH2D ("met_MT2lester_"+nCut,"met_MT2lester "+cutName, nBinsmet, binsmet , nBinsMT2lester,binsMT2lester);
      hmet_MT2lester[cj] ->Sumw2();

      hmet_TauPt[cj] = new TH2D ("met_TauPt_"+nCut,"met_TauPt "+cutName, nBinsmet, binsmet ,nBinsTauPt,binsTauPt);
      hmet_TauPt[cj]->Sumw2();
      
      hmet_DZeta[cj] = new TH2D ("met_DZeta_"+nCut,"met_DZeta "+cutName, nBinsmet, binsmet , nBinsDZeta,binsDZeta);
      hmet_DZeta[cj] -> Sumw2();
     

///////////MT variables
      hMT_MTsum[cj] = new TH2D ("MT_MTsum_"+nCut,"MT_MTsum "+cutName, nBinsMT,binsMT,nBinsMTsum,binsMTsum);
      hMT_MTsum[cj]->Sumw2();

      hMT_MTtot[cj] = new TH2D ("MT_MTtot_"+nCut,"MT_MTtot "+cutName, nBinsMT,binsMT,nBinsMTtot,binsMTtot);
      hMT_MTtot[cj]->Sumw2();

      hMT_MCTb[cj] = new TH2D ("MT_MCTb_"+nCut,"MT_MCTb "+cutName, nBinsMT,binsMT,nBinsMCTb,binsMCTb);
      hMT_MCTb[cj]->Sumw2();

      hMT_MT2lester[cj] = new TH2D ("MT_MT2lester_"+nCut,"MT_MT2lester "+cutName, nBinsMT,binsMT,nBinsMT2lester,binsMT2lester);
      hMT_MT2lester[cj]->Sumw2();

      hMT_TauPt[cj] = new TH2D ("MT_TauPt_"+nCut,"MT_TauPt "+cutName, nBinsMT,binsMT,nBinsTauPt,binsTauPt);
      hMT_TauPt[cj]->Sumw2();
      
      hMT_DZeta[cj] = new TH2D ("MT_DZeta_"+nCut,"MT_DZeta "+cutName, nBinsMT,binsMT,nBinsDZeta,binsDZeta);
      hMT_DZeta[cj]->Sumw2();

/////////////////MTsum
      hMTsum_MTtot[cj] = new TH2D ("MTsum_MTtot_"+nCut,"MTsum_MTtot "+cutName, nBinsMTsum, binsMTsum , nBinsMTtot,binsMTtot);
      hMTsum_MTtot[cj] ->Sumw2();

      hMTsum_MCTb[cj] = new TH2D ("MTsum_MCTb_"+nCut,"MTsum_MCTb "+cutName, nBinsMTsum, binsMTsum , nBinsMCTb,binsMCTb);
      hMTsum_MCTb[cj] ->Sumw2();

      hMTsum_MT2lester[cj] = new TH2D ("MTsum_MT2lester_"+nCut,"MTsum_MT2lester "+cutName, nBinsMTsum, binsMTsum , nBinsMT2lester,binsMT2lester);
      hMTsum_MT2lester[cj] ->Sumw2();


      hMTsum_TauPt[cj] = new TH2D ("MTsum_TauPt_"+nCut,"MTsum_TauPt "+cutName, nBinsMTsum, binsMTsum , nBinsTauPt,binsTauPt);
      hMTsum_TauPt[cj] ->Sumw2();


      hMTsum_DZeta[cj] = new TH2D ("MTsum_DZeta_"+nCut,"MTsum_DZeta "+cutName, nBinsMTsum, binsMTsum , nBinsDZeta,binsDZeta);
      hMTsum_DZeta[cj] ->Sumw2();

//////////MTtot
      hMTtot_MCTb[cj] = new TH2D ("MTtot_MCTb_"+nCut,"MTtot_MCTb "+cutName, nBinsMTtot, binsMTtot , nBinsMCTb,binsMCTb);
      hMTtot_MCTb[cj] ->Sumw2();

      hMTtot_MT2lester[cj] = new TH2D ("MTtot_MT2lester_"+nCut,"MTtot_MT2lester "+cutName, nBinsMTtot, binsMTtot , nBinsMT2lester,binsMT2lester);
      hMTtot_MT2lester[cj] ->Sumw2();

      hMTtot_TauPt[cj] = new TH2D ("MTtot_TauPt_"+nCut,"MTtot_TauPt "+cutName, nBinsMTtot, binsMTtot , nBinsTauPt,binsTauPt);
      hMTtot_TauPt[cj] ->Sumw2();

      hMTtot_DZeta[cj] = new TH2D ("MTtot_DZeta_"+nCut,"MTtot_DZeta "+cutName, nBinsMTtot, binsMTtot , nBinsDZeta,binsDZeta);
      hMTtot_DZeta[cj] ->Sumw2();



/////////MCTb

      hMCTb_MT2lester[cj] = new TH2D ("MCTb_MT2lester_"+nCut,"MCTb_MT2lester "+cutName, nBinsMCTb,binsMCTb ,nBinsMT2lester,binsMT2lester);
      hMCTb_MT2lester[cj]->Sumw2();

      hMCTb_TauPt[cj] = new TH2D ("MCTb_TauPt_"+nCut,"MCTb_TauPt "+cutName, nBinsMCTb,binsMCTb ,nBinsTauPt,binsTauPt);
      hMCTb_TauPt[cj]->Sumw2();

      hMCTb_DZeta[cj] = new TH2D ("MCTb_DZeta_"+nCut,"MCTb_DZeta "+cutName, nBinsMCTb,binsDZeta ,nBinsDZeta,binsDZeta);
      hMCTb_DZeta[cj]->Sumw2();

////////////// MT2lester

       hMT2lester_DZeta[cj] = new TH2D ("MT2lester_DZeta_"+nCut,"MT2lester_DZeta "+cutName, nBinsMT2lester, binsMT2lester , nBinsDZeta,binsDZeta);
	hMT2lester_DZeta[cj]->Sumw2();

       hMT2lester_TauPt[cj] = new TH2D ("MT2lester_TauPt_"+nCut,"MT2lester_TauPt "+cutName, nBinsMT2lester, binsMT2lester , nBinsTauPt,binsTauPt);
	hMT2lester_TauPt[cj]->Sumw2();



     
//////DZeta
      hTauPt_DZeta[cj] = new TH2D ("TauPt_DZeta_"+nCut,"TauPt_DZeta "+cutName, nBinsTauPt,binsTauPt,nBinsDZeta, binsDZeta);
      hTauPt_DZeta[cj] ->Sumw2();


////////////////////////////////////////////////////////////


 

    }

  /*
      char arg[100];
     for (unsigned int i = 0; i < vec.size (); i++)
    {
      var_[i] = -8888.;
      //string name = vec[i].c_str()+"_"+ds;
      sprintf (arg, "%s/F", vec[i].c_str ());
      T->Branch (vec[i].c_str (), &var_[i], arg);

      //cout << " creating the TTree. " << vec[i].c_str()<< "  "<<i<<endl;
    }
*/
}


//void FillHists(int CutIndex, Double_t EvWeight, TLorentzVector  &ElV, TLorentzVector  &MuV, TLorentzVector  &TauV, vector<TLorentzVector>  &JetsV, TLorentzVector  &MetV, double Chimass, double mintermediate, string & Sel, int  mIndex, int eIndex, int  tIndex){}

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

