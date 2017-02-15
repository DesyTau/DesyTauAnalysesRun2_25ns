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
#include "TChain.h"
#include "TMath.h"
#include "Riostream.h"

#include "TRandom.h"
#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "AnalysisMacro.h"

using namespace std;





int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);
  string SelectionSign=argv[3];

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  
  // kinematic cuts on Jets
  const float etaJetCut   = cfg.get<float>("etaJetCut");
  const float ptJetCut   = cfg.get<float>("ptJetCut");
  
  
  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZetaCut       = cfg.get<float>("dZetaCut");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  
  const float Lumi   = cfg.get<float>("Lumi");

  const float bTag 	     = cfg.get<float>("bTag");
  const Float_t metcut         = cfg.get<Float_t>("metcut");
  // **** end of configuration
 
   //TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   //dir.ReplaceAll("basic.C","");
   //dir.ReplaceAll("/./","/");
   //ifstream in;
 CutList.clear();
 CutList.push_back("No cut");
 CutList.push_back("Trigger");
 CutList.push_back("2l dR > "+to_string(dRleptonsCut));
 CutList.push_back("b-Veto ");
 CutList.push_back("lep SumpT> 0");
 CutList.push_back("MET $>$ 50");
 CutList.push_back("MET $>$ 100");
 CutList.push_back("dPhi > 1");


int CutNumb = int(CutList.size());

         xs=1;fact=1;fact2=1;
        
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

	
if (SelectionSign !="OS" && SelectionSign !="SS") {
       cout <<" Wrong selection...you should use OS or SS  as input "<<endl;
   //    SelectionSign="2l";
       return 1;
}       



//CutList[CutNumb]=CutListt[CutNumb];

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());




  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
  int iCut=0;
    double CFCounter[CutNumb];
    double statUnc[CutNumb];
    int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; i++){
          CFCounter[i] = 0;
         iCFCounter[i] = 0;
         statUnc[i] =0;
        }

  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
 
  SetupHists(CutNumb); 
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
    //NE=1000;

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
     //numberOfEntries = 1000;
    
     std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;

     
     numberOfEntries=100;

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
     
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
       
      JetsMV.clear();
      ElMV.clear();
      MuMV.clear();
      LeptMV.clear();



      Float_t MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      
      METV.SetPx(analysisTree.pfmet_ex);	      
      METV.SetPy(analysisTree.pfmet_ey);
 
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
 

      for (unsigned int ij = 0; ij<analysisTree.pfjet_count; ++ij) {
	     JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      JetsMV.push_back(JetsV);
      } 




      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	      MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
       MuMV.push_back(MuV);

      }

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	      ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
       ElMV.push_back(ElV);
      }


      float weight = 1;
       iCut = 0;
      
      Double_t EvWeight = 1.0;
      EvWeight *= weight ;
      
      FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
     // hnJets[->Fill(pfjet_count);

      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      // int indexW  = -1;
      // int indexNu = -1; 
      // int indexMu = -1;
      // int indexE  = -1;
      // int nGenMuons = 0;
      // int nGenElectrons = 0;
      // for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      // 	float pxGen = analysisTree.genparticles_px[igen];
      // 	float pyGen = analysisTree.genparticles_py[igen];
      // 	float pzGen = analysisTree.genparticles_pz[igen];
      // 	float etaGen = PtoEta(pxGen,pyGen,pzGen);
      // 	float ptGen  = PtoPt(pxGen,pyGen);

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
      // 	  indexW = igen;
      // 	if ((fabs(analysisTree.genparticles_pdgid[igen])==12 
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==14
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==16) 
      // 	    && analysisTree.genparticles_info[igen]== (1<<1) )
      // 	  indexNu = igen;

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==13) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexMu = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenMuons++;
      // 	  }
      // 	}
      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==11) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexE = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenElectrons++;
      // 	  }
      // 	}
      // }

      // trigger selection
 
      //selecTable.Fill(1,0, weight );      
      bool trigAccept = false;

      for (int i=0; i<kMaxhltriggerresults; ++i) {
	if ((i==5||i==6)&&analysisTree.hltriggerresults_second[i]==1) {
	  //	  std::cout << analysisTree.run_hltnames->at(i) << " : " << analysisTree.hltriggerresults_second[i] << std::endl;
	  trigAccept = true;
	}
      }

      if (!trigAccept) continue;
      //Trigger
      FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      /////now clear the Mu.El.Jets again to fill them again after cleaning
      MuMV.clear();
      ElMV.clear();
      // electron selection

      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	electronPtAllH->Fill(analysisTree.electron_pt[ie],weight);
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ie];
	if (relIso>isoElectronHighCut) continue;
	if (relIso<isoElectronLowCut) continue;
	bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[ie],
						analysisTree.electron_mva_id_nontrigPhys14[ie]);
	if (!electronMvaId&&applyElectronId) continue;
	electrons.push_back(ie);
        ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
        ElMV.push_back(ElV);
	LeptMV.push_back(ElV);
        hel_miniISO[1]->Fill(analysisTree.electron_miniISO[ie],weight);
      }

      // muon selection


      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	muonPtAllH->Fill(analysisTree.muon_pt[im],weight);
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	float neutralIso = 
	  analysisTree.muon_neutralHadIso[im] + 
	  analysisTree.muon_photonIso[im] - 
	  0.5*analysisTree.muon_puIso[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.muon_chargedHadIso[im] + neutralIso;
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonHighCut) continue;
	if (relIso<isoMuonLowCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	muons.push_back(im);
	MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
        MuMV.push_back(MuV);
	LeptMV.push_back(MuV);
        hmu_miniISO[1]->Fill(analysisTree.muon_miniISO[im],weight);
      }


      ///sort leptons vector
      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      
	//for (unsigned int j=0;j<LeptMV.size();++j) cout<<" j "<<j<<"  "<<LeptMV.at(j).Pt()<<endl;
	//cout<<""<<endl;
      ////////jets cleaning 
      Float_t DRmax=0.4;
      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;
      
      
      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      for (unsigned int il = 0; il<LeptMV.size(); ++il) {
      
	 for (unsigned int ij = 0; ij<JetsMV.size(); ++ij) {
        
		 if(fabs(JetsMV.at(ij).Eta())>etaJetCut) continue;
                 if(fabs(JetsMV.at(ij).Pt())<ptJetCut) continue;
      
       Float_t Dr= deltaR(LeptMV.at(il).Eta(), LeptMV.at(il).Phi(),JetsMV.at(ij).Eta(),JetsMV.at(ij).Phi());

     if (  Dr  < DRmax) {
	     
	     JetsMV.erase (JetsMV.begin()+ij);
    		 }	
		       
	 }
      }
      

      if (electrons.size()==0) continue;
      if (muons.size()==0) continue;
      //if (jets.size()<4) continue;

      // selecting muon and electron pair (OS or SS);
      float ptScalarSum = -1;
      float dRleptons = -1;
      int electronIndex = -1;
      int muonIndex = -1;
      



      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	int eIndex = electrons[ie];
	    
	if ( analysisTree.electron_pt[eIndex]< ptElectronLowCut ) ElMV.erase(ElMV.begin()+ie);
	
			for (unsigned int im=0; im<muons.size(); ++im) {
	  int mIndex = muons[im];
	  float qProd = analysisTree.electron_charge[eIndex]*analysisTree.muon_charge[mIndex];
	
	  if ( analysisTree.muon_pt[mIndex]<ptMuonLowCut ) MuMV.erase(MuMV.begin()+im);

	  //// This is for the SS to be true
	  if (SelectionSign == "SS") {
	  if (oppositeSign || qProd<0) continue;
	  //if (!oppositeSign && qProd<0) continue;
	  
	  }
	  if (SelectionSign == "OS") {
	  if (!oppositeSign || qProd>0) continue;
	  //if (oppositeSign &&  qProd<0) continue;
	  }
           
	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;
          //CFCounter[iCut]+= weight;
          //iCFCounter[iCut]++;
          //iCut++;

		
           Float_t leadLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] > analysisTree.muon_pt[mIndex]	: analysisTree.muon_pt[mIndex];
           Float_t trailLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] < analysisTree.muon_pt[mIndex] : analysisTree.muon_pt[mIndex];

 
	   //cout<<"  Leading "<<leadLep<<" trailing "<<trailLep<<endl;
	   //cout<<"analysisTree.electron_pt[eIndex] "<<analysisTree.electron_pt[eIndex]<<"  nalysisTree.muon_pt[mIndex] "<<analysisTree.muon_pt[mIndex]<<endl;
           //cout<<""<<endl;

	  bool kinematics = 
            (analysisTree.electron_pt[eIndex]> ptElectronLowCut && analysisTree.muon_pt[mIndex]>ptMuonHighCut) || 
	    (analysisTree.electron_pt[eIndex]>ptElectronHighCut && analysisTree.muon_pt[mIndex]>ptMuonLowCut);

	  if (!kinematics) continue;

	  float sumPt = analysisTree.electron_pt[eIndex] + analysisTree.muon_pt[mIndex];
	  if (sumPt>ptScalarSum) {
	    ptScalarSum = sumPt;
	    dRleptons = dR;
	    electronIndex = ie;
	    muonIndex = im;
	  	}
	  
		}
      	}/// dR cut

          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	
	  //if (JetsMV.size()<3) continue;


 	  bool btagged= false;
	  for (unsigned int ib = 0; ib <analysisTree.pfjet_count;ib++){
            if (analysisTree.pfjet_btag[ib][6]  > bTag) btagged = true;
  		  //cout<<" pfjet_b "<<ib<<"  "<<analysisTree.pfjet_btag[ib][6]<<endl;
	  }
	  //if (btagged) continue;

          // Jets
	  FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
          // pt Scalar
	  if (SelectionSign == "OS" || SelectionSign == "SS") {
    	  if (ptScalarSum<0  ) continue;
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	  }
      // computations of kinematic variables

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
						    analysisTree.electron_py[electronIndex],
						    analysisTree.electron_pz[electronIndex],
						    electronMass);
      

      TLorentzVector dileptonLV = muonLV + electronLV;
      float dileptonMass = dileptonLV.M();
      float dileptonPt = dileptonLV.Pt();
      float dileptonEta = dileptonLV.Eta();

      float ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // bisector of electron and muon transverse momenta
      float electronUnitX = electronLV.Px()/electronLV.Pt();
      float electronUnitY = electronLV.Py()/electronLV.Pt();
	
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();

      float zetaX = electronUnitX + muonUnitX;
      float zetaY = electronUnitY + muonUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + electronLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + electronLV.Py();
      
      float vectorVisX = muonLV.Px() + electronLV.Px();
      float vectorVisY = muonLV.Py() + electronLV.Py();

      // computation of DZeta variable
      float PZeta = vectorX*zetaX + vectorY*zetaY;
      float PVisZeta = vectorVisX*zetaX + vectorVisY*zetaY;
      float DZeta = PZeta - 1.85*PVisZeta;

      // computation of MT variable
      float dPhi=-999; 
      float MT = -999;

	   	dPhi=dPhiFrom2P( dileptonLV.Px(), dileptonLV.Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      		//MT=TMath::Sqrt(2*dileptonPt*ETmiss*(1-TMath::Cos(dPhi)));
	




      // filling histograms after dilepton selection

      electronPtH->Fill(electronLV.Pt(),weight);
      electronEtaH->Fill(electronLV.Eta(),weight);
      
      muonPtH->Fill(muonLV.Pt(),weight);
      muonEtaH->Fill(muonLV.Eta(),weight);
      
      dileptonMassH->Fill(dileptonMass,weight);
      dileptonPtH->Fill(dileptonPt,weight);
      dileptonEtaH->Fill(dileptonEta,weight);
      dileptondRH->Fill(dRleptons,weight);
      
      ETmissH->Fill(ETmiss,weight);
      MtH->Fill(MT,weight);
      DZetaH->Fill(DZeta,weight);

      if (ETmiss < metcut) continue;
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
       
       if (ETmiss < 2*metcut) continue;
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
      // topological cut
      //if (DZeta<dZetaCut) continue;
       if (dPhi<1) continue; 
       
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV, analysisTree);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	  

      electronPtSelH->Fill(electronLV.Pt(),weight);
      electronEtaSelH->Fill(electronLV.Eta(),weight);
      
      muonPtSelH->Fill(muonLV.Pt(),weight);
      muonEtaSelH->Fill(muonLV.Eta(),weight);
      
      dileptonMassSelH->Fill(dileptonMass,weight);
      dileptonPtSelH->Fill(dileptonPt,weight);
      dileptonEtaSelH->Fill(dileptonEta,weight);
      dileptondRSelH->Fill(dRleptons,weight);
      
      ETmissSelH->Fill(ETmiss,weight);
      MtSelH->Fill(MT,weight);
      DZetaSelH->Fill(DZeta,weight);


      //      std::cout << std::endl;
      
      selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    
    file_->Close();
    delete file_;
}
 
    cout << endl << "Finished event loop" << endl;
    for (int i=0;i<CutNumb;++i){
	         CFCounter[i] *= float(XSec*Lumi/inputEventsH->GetSum());
                 if (iCFCounter[i] <0.2) statUnc[i] =0;
                else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
        }

    //write out cutflow
    ofstream tfile;
   // TString outname = argv[argc-1];
    TString outname=argv[2];
    TString textfilename = "cutflow_"+outname+"_"+SelectionSign+".txt";
    tfile.open(textfilename);
    tfile << "########################################" << endl;
 //   tfile << "Cut efficiency numbers:" << endl;

 //       tfile << " Cut "<<"\t & \t"<<"#Evnts for "<<Lumi/1000<<" fb-1 & \t"<<" Uncertainty \t"<<" cnt\t"<<endl;
       for(int ci = 0; ci < CutNumb; ci++)
        {
                tfile << CutList[ci]<<"\t & \t"
                      << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
                CutFlow->SetBinContent(1+ci,CFCounter[ci]);
        }

    tfile.close();
    //ofstream tfile1;
    //TString textfile_Con = "CMG_cutflow_Con_Mu_"+outname+".txt";
    //tfile1.open(textfile_Con);
    //tfile1 << "########################################" << endl;
    //tfile1 << "RCS:" << endl;






  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd(SelectionSign.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  CutFlow->Write();

  muonPtAllH->Write();
  electronPtAllH->Write();

  // histograms (dilepton selection)
  electronPtH->Write();  
  electronEtaH ->Write();
  muonPtH ->Write();
  muonEtaH ->Write();

  dileptonMassH ->Write();
  dileptonPtH ->Write();
  dileptonEtaH ->Write();
  dileptondRH ->Write();
  ETmissH ->Write();
  MtH ->Write();
  DZetaH ->Write();

  // histograms (dilepton selection + DZeta cut DZeta)
  electronPtSelH ->Write();
  electronEtaSelH ->Write();
  muonPtSelH  ->Write();
  muonEtaSelH ->Write();

  dileptonMassSelH ->Write();
  dileptonPtSelH ->Write();
  dileptonEtaSelH ->Write();
  dileptondRSelH ->Write();
  ETmissSelH ->Write();
  MtSelH ->Write();
  DZetaSelH ->Write();


  file->Write();
  file->Close();
  
  delete file;
  
}



