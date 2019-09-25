#include "HttStylesNew.cc"
#include "HtoH.h"
#include "FitRecoil.C"

void FitRecoilFromTree( TString met_type = "pfmet" 
                   // can be changed e.g. to puppimet 
                 ) 
{

  SetStyle();
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);

  // location of root files 
  TString dir("./../final/");

  TString fileName_Data(dir+"SingleMuon_Run2016.root");
  TString fileName_DrellYan(dir+"DYJetsToLL_M-50_13TeV_madgraphMLM_pythia8.root");

  TString samples[9] = {
    "WW_13TeV-pythia8",                
    "WZ_13TeV-pythia8",                
    "ZZ_13TeV-pythia8",                
    "WJetsToLNu_13TeV-madgraphMLM-pythia8", 
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8", 
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8",
    "ST_t-channel_top_4f_InclusiveDecays_13TeV-powheg-pythia8", 
    "ST_t-channel_antitop_4f_InclusiveDecays_13TeV-powheg-pythia8",
    "TT_13TeV-powheg-pythia8"
  };

  float xsec[9] = {118.7,
                    27.68,
                    12.19,
                    52760*1.166,
                    35.85,
                    35.85,
                    44.33,
                    26.38,
                    831.8
                 };    

  float xsecDY = 5345*1.079; // Drell Yan cross section
  float lumi = 35900;  // integrated luminosity

  TString baseString("recoilZ");

  TString Proj[2] = {baseString+"Perp",
		             baseString+"Paral"};
  
  TString xtit[2] = {"U_{2} [GeV]","U_{1} [GeV]"};
  TString ytit[2] = {"Events","Events"};

  TString weights = "(mc_weight*pu_weight*zpt_weight*toppt_weight*idiso_weight_1*idiso_weight_2*trig_weight*track_weight_1*track_weight_2)";
  TString selection = "(iso_2<0.15 && m_ll > 70 && m_ll < 110)"; 
  TString os_cut = "(os>0.5)";
  TString ss_cut = "(os<0.5)";

  // extrapolation factor for QCD estimation 
  float ss_os_ratio = 2.;   

  // set the name of the variables for TTree->Draw
  TString variable[2];
  if (met_type == "pfmet" ) {
    variable[0] = "met_recoil_perp";
    variable[1] = "met_recoil_paral"; 
  }
  else if (met_type == "puppimet") {
    variable[0] = "puppimet_recoil_perp";
    variable[1] = "puppimet_recoil_paral"; 
  }
    
  // histograms binning
  // final range overwritten later by xrange
  int nbins = 200;
  float xmin = -400;
  float xmax = 400;

  // njets binning 
  int nJetBins = 3;
  float jetBinsX[4] = {-0.5,0.5,1.5,2.5};
  TString jetBins[3] = {"NJet0",
			            "NJet1",
			            "NJetGe2"};

  TString jetBinCut[3] = {"(njets30==0)", "(njets30==1)", "(njets30>1)"}; 


  // Z Pt binning 
  int nZPtBins = 5;
  float zPtBinsX[6] = {0,10,20,30,50,1000};
  TString ptBins[5] = {"Pt0to10",
		               "Pt10to20",
		               "Pt20to30",
		               "Pt30to50",
		               "PtGt50"};

  TString ZPtCut[5] = {"(pt_ll>=0 && pt_ll<10)",
                       "(pt_ll>=10 && pt_ll<20)",
                       "(pt_ll>20 && pt_ll<30)",
                       "(pt_ll>=30 && pt_ll<50)",
                       "(pt_ll>=50)"
                      };

  TFile * fileOutput = new TFile(baseString+".root","recreate");
  // write the used binning in the output file 
  fileOutput->cd("");
  TH1D * projH = new  TH1D("projH","",2,-1,1);
  TH1D * zPtBinsH = new TH1D("ZPtBinsH","",nZPtBins,zPtBinsX);
  TH1D * nJetBinsH = new TH1D("nJetBinsH","",nJetBins,jetBinsX);
  for (int i=0; i<2; ++i)
    projH->GetXaxis()->SetBinLabel(i+1,Proj[i]);
  for (int i=0; i<nJetBins; ++i)
    nJetBinsH->GetXaxis()->SetBinLabel(i+1,jetBins[i]);
  for (int i=0; i<nZPtBins; ++i)
    zPtBinsH->GetXaxis()->SetBinLabel(i+1,ptBins[i]);
  
  projH->Write("projH");
  zPtBinsH->Write("ZPtBinsH");
  nJetBinsH->Write("nJetBinsH");

  TFile * fileData   = new TFile(fileName_Data);
  TTree * treeData = (TTree*)(fileData->Get("ZMuMu"));

  TFile * fileDrellYan     = new TFile(fileName_DrellYan);
  TTree * treeDrellYan = (TTree*)(fileDrellYan->Get("ZMuMu"));

  TH1D * nEventsDY = (TH1D*)fileDrellYan->Get("inputEventsH");
  float normDY = xsecDY*lumi/nEventsDY->GetSumOfWeights();

  // create histograms for data and Drell Yan
  // then, subtract all bkgs. from data to obtain Z->mumu events  
  // for each recoil projection, jet bin and Z Pt bin

  for (int iProj=0; iProj<2; ++iProj) {
    for (int iJet=0; iJet<nJetBins; ++iJet) {
      for (int iPt=0; iPt<nZPtBins; ++iPt)    {

        TString histName = Proj[iProj] + "_" + jetBins[iJet] + ptBins[iPt];
        TString full_selection = selection + "&&" + jetBinCut[iJet] + "&&" + ZPtCut[iPt] ; 
        cout << histName << endl;

        TH1D * histData = new TH1D(histName+"_data", "", nbins, xmin, xmax);  
        treeData->Draw(variable[iProj]+">>"+histName+"_data", full_selection + "&&" + os_cut);   

        TH1D * histMC = new TH1D(histName+"_mc", "", nbins, xmin, xmax); 
        treeDrellYan->Draw(variable[iProj]+">>"+histName+"_mc", weights + "* (" + full_selection + "&&" + os_cut + ")");   
        histMC->Scale(normDY);

        for (int iS=0; iS<9; ++iS) {

          TFile * fileSample = new TFile(dir+samples[iS]+".root");
          TTree * treeSample = (TTree*)(fileSample->Get("ZMuMu"));
          TH1D * inputEventsH = (TH1D*)fileSample->Get("inputEventsH");
          double normSample;

          if (inputEventsH->GetSumOfWeights() > 0)
             normSample = lumi*xsec[iS]/inputEventsH->GetSumOfWeights();
          else
            std::cout << "Sample " << dir << samples[iS] <<".root : inputEventsH-> GetSumOfWeights() = 0 !" << std::endl;   
      
          TH1D * histSample =  new TH1D(histName+"_"+samples[iS], "", nbins, -xmin, xmax);
          treeSample->Draw(variable[iProj]+">>"+histName+"_"+samples[iS], weights + "* (" + full_selection + "&&" + os_cut + ")");
          histSample->Scale(normSample); 
          // subtract this bkg from data
          histData->Add(histData,histSample,1,-1);
          delete fileSample;
        }
        // subtract the QCD background. Evaluated from ss data, scaled by SS/OS ratio. 
        TH1D * histQCD =  new TH1D(histName+"_qcd", "", nbins, xmin, xmax); 
        treeData->Draw(variable[iProj]+">>"+histName+"_qcd", full_selection + "&&" + ss_cut ); 
        histData->Add(histData, histQCD, 1,-1*ss_os_ratio);  

        bool rebin = false;
	    float xrange = 120; 
	    bool logY = false;
	    bool asymGauss = true;
	    //	if (iJet==2) rebin = true;
	    if (iProj==0) asymGauss = false;
	    if (iJet==0) xrange = 120;
	    FitRecoil(histData,histMC,fileOutput,histName,-xrange,xrange,xtit[iProj],ytit[iProj],-180,180,asymGauss,rebin,logY);
      }
    }
  }

  fileOutput->Close();

}
