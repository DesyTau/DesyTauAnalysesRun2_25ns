//----------------------Version 2.0-------------------------//
//Plotting Macro for ele -> Tau Fake Rates study
//Author: Yiwen Wen & Andrea Cardini
//DESY
//-----------------------------------------------------------//
#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <math.h>
void ProduceDatacardInputs_ETauFR(
				 TString Variable = "m_vis",
                                 TString etaSuffix = "Lt1p460",
                                 TString wp = "Medium",
				 int nBins = 12,
                                 float xmin = 60,
                                 float xmax = 120,
				 TString Year="2018",
                                 bool passProbe = true,
				 TString wpIso = "Tight",
                                 bool DeepTau = true
                                 ){
	
        bool applyPU = true;
	TString tauIso = "tauby"+wpIso+"IsolationMVArun2v1DBoldDMwLT";  //Medium wpIso not in the Ntuple for MVA discriminator
 	TString againstMu = "tauagainstMuonLoose3";
	TString DeepTauString ="MVA";

	if(DeepTau>0.5){
	   DeepTauString ="DeepTau";
	   tauIso = "tauby"+wpIso+"DeepTau2017v2VSjet";
	   againstMu = "tauagainstMuLooseDeepTau";
	}
  
	
using namespace std;
	SetStyle();
	TString directory="./";

	TString MCweight="effweight*mcweight*";
	const int nSamples = 20;
				 
	TString Cut="(iso_1<0.1&&mt_1<30)";
	
	TString samples[nSamples] =
	  {
	    "EGamma_Run2018",//(0)data
	    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",// (1)Drell-Yan Z->EE
	    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (2)Drell-Yan ZJ
	    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (3)Drell-Yan ZTT(tau -> lepton)
	    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (4)Drell-Yan ZTT(hadronic tau)
	    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",//(5)TTbar leptonic
	    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",//(6) semileptonic
	    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",//(7) hadronic
	    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
	    "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (9)WJets
	    "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (10)WJets
	    "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (11)WJets
	    "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (12)WJets
	    "WW_TuneCP5_13TeV-pythia8",// (13)WW
	    "WZ_TuneCP5_13TeV-pythia8",// (14)WZ
	    "ZZ_TuneCP5_13TeV-pythia8",// (15)ZZ
	    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (16) SingleTop tW tbar
	    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (17) SingleTop tW t
	    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",// (18) SingleTop t antitop
	    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"// (19) SingleTop t top
	  };
	float WJetsweight=0.8369;
	float DYweight=1.;
	
	double xsec[nSamples] = {
	  1, // (0)data
	  DYweight,// (1)Drell-Yan Z->EE
	  DYweight,// (2)Drell-Yan ZJ
	  DYweight,// (3)Drell-Yan ZTT(tau->lepton)
	  DYweight,// (4)Drell-Yan ZTT(hadronic tau)
	  88.29, // (5)TTPowHeg leptonic
	  377.96, // (6)TTBar semileptonic
	  365.35, // (7)TTBar hadronic
	  WJetsweight, // (8)WJetsToLNu_MG
	  WJetsweight, // (9)WJetsToLNu_MG
	  WJetsweight, // (10)WJetsToLNu_MG
	  WJetsweight, // (11)WJetsToLNu_MG
	  WJetsweight, // (12)WJetsToLNu_MG
	  63.21,  // WW (13)
	  22.82,   // WZ (14)
	  10.32,  // ZZ (15)
	  38.09,   // ST_tW_antitop_5f_inclusiveDecays (16)
	  38.09,   // ST_tW_top_5f_inclusiveDecays (17)
	  80.95, // (18) SingleTop t antitop
	  136.02// (19) SingleTop t top
	};
	float lumi;
	if(Year=="2018") lumi = 59970;
	else if(Year=="2017") lumi = 41860;
	else if(Year=="2016") lumi = 36773;
	else exit(EXIT_FAILURE);
	//float lumi = 41860;
	
	TString suffix("");
	if(passProbe)
	  suffix="_pass";
	else
	  suffix="_fail";
	
	
	TH1D * hist[nSamples];
	TH1D * histSS[nSamples];
	
	//inistiating cuts
	TString cuts[nSamples];
	TString cutsSS[nSamples];
   	
 	TString qcdweight ="1.06*";
    
    //0.892518 +/- 0.00446623
	for (int i=0; i<nSamples; ++i)
	{
	  cuts[i] = MCweight+Cut+"*(os>0.5) ";
        cutsSS[i] = qcdweight+MCweight+Cut+"*(os<0.5)";
  	}
	TString ZEEcut    = MCweight+Cut+"* (gen_match_1 == 1 && gen_match_2 == 1) ";
	TString ZJcut     = MCweight+Cut+"* (gen_match_2 == 6) ";
	TString ZTT_elcut = MCweight+Cut+"* (gen_match_1 == 3 && (gen_match_2 == 3 || gen_match_2 == 4)) ";
	TString ZTT_etcut = MCweight+Cut+"* (gen_match_1 == 3 && gen_match_2 == 5)";
  	cuts[0] = Cut       +"*(os>0.5)";
	cuts[1] = ZEEcut    +"*(os>0.5)";
	cuts[2] = ZJcut     +"*(os>0.5)";
	cuts[3] = ZTT_elcut +"*(os>0.5)";
	cuts[4] = ZTT_etcut +"*(os>0.5)";

	if (DeepTau){
		ZTT_etcut += "*0.83";
	}else{
		ZTT_etcut += "*0.9";
	}
   
	cutsSS[0] = qcdweight+Cut       +"*(os<0.5)";
	cutsSS[1] = qcdweight+ZEEcut    +"*(os<0.5)";
	cutsSS[2] = qcdweight+ZJcut     +"*(os<0.5)";
	cutsSS[3] = qcdweight+ZTT_elcut +"*(os<0.5)";
	cutsSS[4] = qcdweight+ZTT_etcut +"*(os<0.5)";
  
	Cut+="*("+againstMu+">0.5&&"+tauIso+">0.5)";

	if (DeepTau){
		Cut="(DecayModeProbe!=5&&DecayModeProbe!=6)*"+Cut;
			
	}else{
		Cut="(tau_decayModeFinding>0.5)*"+Cut;
		
	}	  

	if(passProbe)	{
        
	for (int i=0; i<nSamples; ++i)
	  {
            cuts[i] ="(tauagainstEle"+wp+DeepTauString+">0.5&&"+tauIso+">0.5&&"+againstMu+">0.5)*"+cuts[i];
            cutsSS[i] ="(tauagainstEle"+wp+DeepTauString+">0.5&&"+tauIso+">0.5&&"+againstMu+">0.5)*"+cutsSS[i];
        }
        
	}
	else{
        
	for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(tauagainstEle"+wp+DeepTauString+"<0.5&&"+tauIso+">0.5&&tauagainstEleVVVLooseDeepTau>0.5&&"+againstMu+">0.5)*"+cuts[i];
            cutsSS[i] ="(tauagainstEle"+wp+DeepTauString+"<0.5&&"+tauIso+">0.5&&tauagainstEleVVVLooseDeepTau>0.5&&"+againstMu+">0.5)*"+cutsSS[i];
        }
    }


    if(applyPU)
	{
		for (int i=1; i<nSamples; ++i)
		{
			cuts[i] ="puweight*"+cuts[i];
			cutsSS[i] ="puweight*"+cutsSS[i];
		}
	}
    
    if(etaSuffix=="Lt1p460")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] = "(fabs(EtaProbe)<1.460)*"+cuts[i];
            cutsSS[i] = "(fabs(EtaProbe)<1.460)*"+cutsSS[i];
        }
    }
    
    if(etaSuffix=="Gt1p558")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(fabs(EtaProbe)>1.558)*"+cuts[i];
            cutsSS[i] ="(fabs(EtaProbe)>1.558)*"+cutsSS[i];
        }
    }
    
	for (int i=0; i<nSamples; ++i)
	{
		std::cout <<i<<":"+samples[i]<<std::endl;
		TFile * file = new TFile(directory+samples[i]+".root");
		TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
		TTree * tree = (TTree*)file->Get("ETauFR");
		double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
		if(samples[i].Contains("Jets")){
		  norm=xsec[i]*lumi;
		  cuts[i]+="*stitchweight";
		  cutsSS[i]+="*stitchweight";
		}
		TString histName = samples[i] + "_"+Variable;
		TString histNameSS = samples[i] + "_"+Variable+"_ss";
		hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
		hist[i]->Sumw2();
		histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
		histSS[i]->Sumw2();  
		if( (Variable.Contains("scale")) && (i<1 || i>4)) continue;

		tree->Draw(Variable+">>"+histName,cuts[i]);
		tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
		if(i > 0)
		  {		
		    hist[i]->Scale(norm);
		    histSS[i]->Scale(norm);
		    
		  }
		cout << samples[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;

	}


  // ********************************************
  // ***** Stitching of Drell-Yan samples *******
  // ********************************************
    
	TH1D * histZEE[4];
	TH1D * histZJ[4];
	TH1D * histZTT_EL[4];
	TH1D * histZTT_ET[4];
	TH1D * histSSZEE[4];
	TH1D * histSSZJ[4];
	TH1D * histSSZTT_EL[4];
	TH1D * histSSZTT_ET[4];
	TString refSamplesDY=samples[1];

	
	
	// filling histograms for DY samples
	for (int i=0; i<4; i++) { // run over samples
	  TString name=refSamplesDY;
	  TString njets=TString::Itoa(i+1,10);
	  njets+="Jets";
	  TString filename=name.ReplaceAll("Jets",njets);
	  cout<<filename<<endl;
	
	  TFile * file = new TFile(directory+filename+".root");
	  TTree * tree = (TTree*)file->Get("ETauFR");
	  double norm = lumi;
	  
	  TString histNameZEE   = filename + Variable + "_zee_os";
	  TString histNameSSZEE = filename + Variable + "_zee_ss";
	  TString histNameZJ   = filename + Variable + "_zj_os";
	  TString histNameSSZJ = filename + Variable + "_zj_ss";
	  TString histNameZTT_EL   = filename + Variable + "_ztt_el_os";
	  TString histNameSSZTT_EL = filename + Variable + "_ztt_el_ss";
	  TString histNameZTT_ET   = filename + Variable + "_ztt_et_os";
	  TString histNameSSZTT_ET = filename + Variable + "_ztt_et_ss";
	  histZEE[i]   = new TH1D(histNameZEE,"",nBins,xmin,xmax);
	  histSSZEE[i] = new TH1D(histNameSSZEE,"",nBins,xmin,xmax);
	  histZJ[i]   = new TH1D(histNameZJ,"",nBins,xmin,xmax);
	  histSSZJ[i] = new TH1D(histNameSSZJ,"",nBins,xmin,xmax);
	  histZTT_EL[i]   = new TH1D(histNameZTT_EL,"",nBins,xmin,xmax);
	  histSSZTT_EL[i] = new TH1D(histNameSSZTT_EL,"",nBins,xmin,xmax);
	  histZTT_ET[i]   = new TH1D(histNameZTT_ET,"",nBins,xmin,xmax);
	  histSSZTT_ET[i] = new TH1D(histNameSSZTT_ET,"",nBins,xmin,xmax);
	  

	  tree -> Draw(Variable+">>"+histNameZEE,     cuts[1]);
	  tree -> Draw(Variable+">>"+histNameSSZEE,   cutsSS[1]);  
	  tree -> Draw(Variable+">>"+histNameZJ,      cuts[2]);
	  tree -> Draw(Variable+">>"+histNameSSZJ,    cutsSS[2]);
	  tree -> Draw(Variable+">>"+histNameZTT_EL,  cuts[3]);
	  tree -> Draw(Variable+">>"+histNameSSZTT_EL,cutsSS[3]);
	  tree -> Draw(Variable+">>"+histNameZTT_ET,  cuts[4]);
	  tree -> Draw(Variable+">>"+histNameSSZTT_ET,cutsSS[4]);
	
	  histZEE[i]   ->Scale(norm); 
	  histSSZEE[i] ->Scale(norm); 
	  histZJ[i]   ->Scale(norm); 
	  histSSZJ[i] ->Scale(norm); 
	  histZTT_EL[i]   ->Scale(norm);
	  histSSZTT_EL[i] ->Scale(norm);
	  histZTT_ET[i]   ->Scale(norm);
	  histSSZTT_ET[i] ->Scale(norm);
	  
  }


	for (int iDY=0; iDY<4; ++iDY) {
	  hist[1]  -> Add(hist[1],histZEE[iDY]);
	  hist[2]  -> Add(hist[2],histZJ[iDY]);
	  hist[3]  -> Add(hist[3],histZTT_EL[iDY]);
	  hist[4]  -> Add(hist[4],histZTT_ET[iDY]);
	  histSS[1]-> Add(histSS[1],histSSZEE[iDY]);
	  histSS[2]-> Add(histSS[2],histSSZJ[iDY]);
	  histSS[3]  -> Add(histSS[3],histSSZTT_EL[iDY]);
	  histSS[4]  -> Add(histSS[4],histSSZTT_ET[iDY]);
	}
	
    //adding up ttbar backgrounds
    for (int iH=6; iH<8; ++iH)
    {
        hist[5]->Add(hist[5],hist[iH]);
    }

    //adding up WJets backgrounds
                                 
    for (int iH=9; iH<13; ++iH)
    {
        hist[8]->Add(hist[8],hist[iH]);
    }
     

    //  adding up VV backgrounds
    for (int iH=14; iH<nSamples; ++iH)
    {
        hist[13]->Add(hist[9],hist[iH]);
    }

    for (int iH=1; iH<nSamples; ++iH)
    {
        histSS[0]->Add(histSS[0],histSS[iH],1,-1);

    }
    
    //Getting rid of negative bin in QCD
    for (int iB=1; iB<=nBins; ++iB)
    {
        float ySS = histSS[0]->GetBinContent(iB);
        if (ySS<0)
        {
            histSS[0]->SetBinContent(iB,0.);
            histSS[0]->SetBinError(iB,0.);
        }
        
    }
    TFile * file = new TFile("./shapes/ETauFRAgainstEle"+wp+etaSuffix+suffix+"_"+Variable+".root","recreate");
    file->mkdir("ETauFR"+suffix);
    file->cd("ETauFR"+suffix);

    TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
    TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
    //if(!passProbe)
        //QCD = (TH1D*)hist[16]->Clone("QCD");
    TH1D * ZEE = (TH1D*)hist[1]->Clone("ZEE");

    TH1D * ZJ = (TH1D*)hist[2]->Clone("ZJ");
    TH1D * ZTT_el = (TH1D*)hist[3]->Clone("ZTT_el");
    TH1D * ZTT_et = (TH1D*)hist[4]->Clone("ZTT_et");
      
    TH1D * W = (TH1D*)hist[8]->Clone("W");
    TH1D * TT = (TH1D*)hist[5]->Clone("TT");
    TH1D * VV = (TH1D*)hist[13]->Clone("VV");
    
    float totData = data_obs->GetSumOfWeights();
    float totZEE = ZEE->GetSumOfWeights();
    float totZJ = ZJ->GetSumOfWeights();
    float totZTT_el = ZTT_el->GetSumOfWeights();
    float totZTT_et = ZTT_et->GetSumOfWeights();
    float totQCD = QCD->GetSumOfWeights();
    float totTT = TT->GetSumOfWeights();
    float totVV = VV->GetSumOfWeights();
    float totW = W->GetSumOfWeights();
    float totMC = totZEE+totZJ+totZTT_el+totZTT_et+totTT+totW+totVV+totQCD;
    
    std::cout << "data : " << totData << std::endl;
    std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
    std::cout << "VV  : " << VV->GetSumOfWeights() << std::endl;
    std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
    std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
    std::cout << "ZEE : " << ZEE->GetSumOfWeights() << std::endl;
    std::cout << "ZJ : " << totZJ << std::endl;
    std::cout << "ZTT_el : " << totZTT_el << std::endl;
    std::cout << "ZTT_et : " << totZTT_et << std::endl;
    std::cout << "MC : " << totMC << std::endl;
    if(!(Variable.Contains("scale"))){
      data_obs->Write("data_obs");
      QCD->Write("QCD");
      ZJ ->Write("ZJ");
      ZTT_el->Write("ZTT_el");
      W->Write("W");
      TT->Write("TT");
      VV->Write("VV");
      ZEE->Write("ZEE");
      ZTT_et->Write("ZTT_et");
    }else{
      TString suffix="_";
      if(Variable.Contains("probe")){
	suffix+="probe";
	if(Variable.Contains("tau"))suffix+="tau";
	else suffix+="ele";
      }else if(Variable.Contains("tagele"))suffix+="tagele";
      else suffix+="reso";

      if(Variable.Contains("scale_up"))suffix+="_Up";
      else if(Variable.Contains("scale_down"))suffix+="_Down";
      else{
	cout << "ERROR: scale variation not accounted for, please check!"<<endl;
	return;
      }
      ZEE->Write("ZEE"+suffix);
      ZTT_et->Write("ZTT_et"+suffix);
    }

    //file->Write();
    file->Close();


}

