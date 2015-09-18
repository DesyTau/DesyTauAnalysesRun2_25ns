#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
const float MuMass = 0.105658367;

int binNumber(float x, int nbins, float * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double PtoPhi(double Px, double Py) {
  return TMath::ATan2(Py,Px);
}

double PtoPt(double Px, double Py) {
  return TMath::Sqrt(Px*Px+Py*Py);
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

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

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;

}
bool electronMvaIdTight(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.73) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.57) passed = true;
  }
  else {
    if (mva>0.05) passed = true;
  }

  return passed;

}

struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject;

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // Run range
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // **** end of configuration

  // file name and tree name
	char ff[100];

	sprintf(ff,"%s/%s",argv[3],argv[2]);

  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  string SelectionSign="mumu";

  TString era=argv[3];
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  //TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  TFile * file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());
  TH1D * hxsec = new TH1D("xsec","",1,0,0);
  TH1D * histWeights = new TH1D("histWeights","",1,0,0);
  TH1D * histWeights2 = new TH1D("histWeights2","",2,0,0);
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);

  TH1D * JPsiMassAllMuonsH = new TH1D("JPsiMassAllMuonsH","",200,2,4);
  TH1D * JPsiMassAllMuonsDRCutH = new TH1D("JPsiMassAllMuonsDRCutH","",200,2,4);
  TH1D * JPsiMassIdMuonsH = new TH1D("JPsiMassIdMuonsH","",200,2,4);
  TH1D * JPsiMassIdMuonsDRCutH = new TH1D("JPsiMassIdMuonsDRCutH","",200,2,4);

  TH1D * YpsilonMassAllMuonsH = new TH1D("YpsilonMassAllMuonsH","",400,8,12);
  TH1D * YpsilonMassAllMuonsDRCutH = new TH1D("YpsilonMassAllMuonsDRCutH","",400,8,12);
  TH1D * YpsilonMassIdMuonsH = new TH1D("YpsilonMassIdMuonsH","",400,8,12);
  TH1D * YpsilonMassIdMuonsDRCutH = new TH1D("YpsilonMassIdMuonsDRCutH","",400,8,12);

  TH1D * ZMassAllMuonsH = new TH1D("ZMassAllMuonsH","",60,60,120);
  TH1D * ZMassAllMuonsDRCutH = new TH1D("ZMassAllMuonsDRCutH","",400,60,120);
  TH1D * ZMassIdMuonsH = new TH1D("ZMassIdMuonsH","",60,60,120);
  TH1D * ZMassIdMuonsDRCutH = new TH1D("ZMassIdMuonsDRCutH","",60,60,120);
  TH1D * ZMassIsoMuonsH = new TH1D("ZMassIsoMuonsH","",60,60,120);
  TH1D * ZMassIsoMuonsDRCutH = new TH1D("ZMassIsoMuonsDRCutH","",60,60,120);

  Float_t XSec=-1;
  Float_t xs,fact,fact2;
  int nFiles = 0;
  int nEvents = 0;
  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;
 
  ifstream ifs("xsecs");
  string line;

  while(std::getline(ifs, line)) // read one line from ifs
    {
		
      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt;
      iss >> dt >> xs >> fact >> fact2;
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      //if (dt==argv[2]) {
      //if (std::string::npos != dt.find(argv[2])) {
      if (  dt == argv[2]) {
	XSec= xs*fact*fact2;
	cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
      }
        
    }

  if (XSec<0) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}


  std::vector<unsigned int> allRuns; allRuns.clear();
	//nTotalFiles=50;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;



    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      Float_t weight = 1;
      
	bool isData= false;
    	bool lumi=false;

      if (XSec == 1)  isData = true;
      if (!isData && XSec !=1 )  { weight *=analysisTree.genweight;   lumi=true;} 
   
      if( !isData && analysisTree.genweight){
     histWeights->Fill(1,weight); 
     histWeights2->Fill(weight); 
      }  
      else histWeights->Fill(1); 

      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

	
      //if (analysisTree.event_run != 251251) continue;

//if (nEvents%1000==0) cout<<" run "<<analysisTree.event_run<<"  "<<analysisTree.event_luminosityblock<<endl;

    
      if (isData){
      if (analysisTree.event_run<RunRangeMin) continue;
      if (analysisTree.event_run>RunRangeMax) continue;
      

      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;

            //std::cout << " Run : " << analysisTree.event_run << std::endl;

      bool isNewRun = true;
      if (allRuns.size()>0) {
	for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	  if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }

      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);

   
   std::vector<Period> periods;
    
    std::fstream inputFileStream("temp", std::ios::in);
    for(std::string s; std::getline(inputFileStream, s); )
    {
        periods.push_back(Period());
        std::stringstream ss(s);
        ss >> periods.back();
    }
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;

    std::string num = std::to_string(n);
    std::string lnum = std::to_string(lum);
    for(const auto& a : periods)
    {
        
         if ( num.c_str() ==  a.name ) {
        //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
       //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

         for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

	//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	}
        auto last = std::prev(a.ranges.end());
        //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
   	if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	 }
	
}
    
}
//	if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
    //std::remove("myinputfile");
	if (!lumi) continue;


      // vertex cuts

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;

      

      // muon selection

      vector<unsigned int> allMuons; allMuons.clear();
      vector<unsigned int> idMuons; idMuons.clear();
      vector<unsigned int> isoMuons; isoMuons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	allMuons.push_back(im);
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (!analysisTree.muon_isMedium[im]) continue;
	idMuons.push_back(im);
	float absIso = analysisTree.muon_chargedHadIso[im];
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut) continue;
	isoMuons.push_back(im);
      }
/*
       std::cout << "allMuons : " << allMuons.size() << std::endl;
       std::cout << "idMuons  : " << idMuons.size() << std::endl;
       std::cout << "isoMuons : " << isoMuons.size() << std::endl;
*/
      //      continue;

      bool isAllMuonsPair = false;
      if (allMuons.size()>1) {
	//	std::cout << "allMuons : " << allMuons.size() << std::endl;
	for (unsigned int im1=0; im1<allMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<allMuons.size(); ++im2) {
	    unsigned int index1 = allMuons[im1];
	    unsigned int index2 = allMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {
	      isAllMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      JPsiMassAllMuonsH->Fill(mass,weight);
	      YpsilonMassAllMuonsH->Fill(mass,weight);
	      ZMassAllMuonsH->Fill(mass,weight);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassAllMuonsDRCutH->Fill(mass,weight);
		YpsilonMassAllMuonsDRCutH->Fill(mass,weight);
		ZMassAllMuonsDRCutH->Fill(mass,weight);
	      }

	    }
	  }
	}
	//	std::cout << "allMuons : " << allMuons.size() << "  :  OK" << std::endl;
      }

      //      continue;

      bool isIdMuonsPair = false;
      if (idMuons.size()>1) {
	for (unsigned int im1=0; im1<idMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<idMuons.size(); ++im2) {
	    unsigned int index1 = idMuons[im1];
	    unsigned int index2 = idMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {
	      isIdMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      	  //    std::cout << "Mass = " << mass << std::endl;
	      JPsiMassIdMuonsH->Fill(mass,weight);
	      YpsilonMassIdMuonsH->Fill(mass,weight);
	      ZMassIdMuonsH->Fill(mass,weight);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassIdMuonsDRCutH->Fill(mass,weight);
		YpsilonMassIdMuonsDRCutH->Fill(mass,weight);
		ZMassIdMuonsDRCutH->Fill(mass,weight);
	      }

	    }
	  }
	}
      }

      bool isIsoMuonsPair = false;
      if (isoMuons.size()>1) {
	for (unsigned int im1=0; im1<isoMuons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
	    unsigned int index1 = isoMuons[im1];
	    unsigned int index2 = isoMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    if (q1*q2<0) {
	      isIsoMuonsPair = true;
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[index2],
						  analysisTree.muon_py[index2],
						  analysisTree.muon_pz[index2],
						  muonMass);
	      TLorentzVector dimuon = muon1 + muon2;
	      float mass = dimuon.M();
	      ZMassIsoMuonsH->Fill(mass,weight);
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	      if (dR>dRleptonsCut) {
		ZMassIsoMuonsDRCutH->Fill(mass,weight);
	      }

	    }
	  }
	}
      }

      if (isAllMuonsPair) selEventsAllMuons++;
      if (isIdMuonsPair)  selEventsIdMuons++;
      if (isIsoMuonsPair) selEventsIsoMuons++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                     = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events (muon pairs)     = " << selEventsAllMuons << std::endl;
  std::cout << "Total number of selected events (id muon pairs)  = " << selEventsIdMuons << std::endl;
  std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;
  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;

  file->cd(SelectionSign.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  histWeights->Write();
  histWeights2->Write();
  file->Write();
  file->Close();
  delete file;
  
}



