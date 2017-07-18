
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//
//    This code uses the root files produced by the TagAndProbe code and computes the efficiencies fitting the Z peak 
//		with the fitting tool. In this code the eta and pt bins are defined. Root files with will be produced. 
//		The plots with the various fits will be stored in dedicated directories.
//		This code is for electrons.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "DesyTauAnalyses/NTupleMaker/test/TP_2016/FitPassAndFail.C"

void TP_eff_e(TString fileName = "SingleElectron_Run2016B_TP", // RooT file with tag-&-probe histograms 

 	   TString what = "IdIso", //what do you want to evaluated
 	   float iso_max = 0.1, //isolation cut to be used
	   float iso_min = 0.0,
	   float norm = 1 // luminosity normalization factor (1 for data) 
	   ) 
{

	gErrorIgnoreLevel = kFatal;
	TH1::SetDefaultSumw2(kTRUE);

  // output inizialization 
    TString lepton = "Electron";
	TString OutFileName = fileName + "_" + lepton + "_" + what + "_IsoGt" + Form("%.2f", iso_min)  + "_IsoLt" + Form("%.2f", iso_max) + "_eff_full2016";
	TFile * outputFile = new TFile(OutFileName+".root","recreate");

	// Title of axis in plots  
	TString yTitle = "Efficiency";
	TString xTitle = lepton+"  p_{T}[GeV]";
	TString xtit; 
	xtit = "m_{ee}[GeV]"; 

	//names of final graphs - suffix
	bool isData=false;
	if (fileName.Contains("SingleElectron")) isData = true;
  TString SampleName("_MC");
  if (isData) SampleName = "_Data";

  //open input file
	TFile * file = new TFile(fileName+".root");
  file->cd();
  TTree *t = (TTree*)(file->Get("TagProbe"));


  //binning inizialization

//for triggers up to 2.5

  int nEtaBins = 3;
  float etaBins[4] = {0,1.48,2.1,2.5};

  TString EtaBins[3] = {"EtaLt1p48","Eta1p48to2p1",
				"EtaGt2p1"};
  float ptBins_def[11] = {10,15,20,25,30,40,50,60,70,100,1000};
  TString PtBins_def[10] = {"Pt10to15",
       "Pt15to20",
       "Pt20to25",
       "Pt25to30",
       "Pt30to40",
       "Pt40to50",
       "Pt50to60",
       "Pt60to70",
       "Pt70to100",
       "PtGt100"};
/*
  float ptBinsTrig_def[26] = {10,
			  13,
			  16,
              17, 
              18,
			  19,
              20,
              21,
			  22, 
			  23,
			  24,
			  25,
			  26,
			  27,
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
			  200,
			  1000};

	TString PtBinsTrig_def[25] = {"Pt10to13",
		    "Pt13to16","Pt16to17","Pt17to18",
		    "Pt18to19","Pt19to20","Pt20to21",
		    "Pt21to22","Pt22to23","Pt23to24",
		    "Pt24to25","Pt25to26","Pt26to27",
		    "Pt27to28",
		    "Pt28to31",
		    "Pt31to34",
		    "Pt34to37",
		    "Pt37to40",
		    "Pt40to45",
		    "Pt45to50",
		    "Pt50to60",
		    "Pt60to70",
		    "Pt70to100",
			"Pt100to200",
			"PtGt200"};

	int nPtBins = 25; if(what == "IdIso") nPtBins = 10;
*/

// start from pt = 24
  float ptBinsTrig_def[16] = {
			  24,
			  25,
			  26,
			  27,
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
			  200,
			  1000};

	TString PtBinsTrig_def[15] = {
		    "Pt24to25","Pt25to26","Pt26to27",
		    "Pt27to28",
		    "Pt28to31",
		    "Pt31to34",
		    "Pt34to37",
		    "Pt37to40",
		    "Pt40to45",
		    "Pt45to50",
		    "Pt50to60",
		    "Pt60to70",
		    "Pt70to100",
			"Pt100to200",
			"PtGt200"};

	int nPtBins = 15; if(what == "IdIso") nPtBins = 10;



/*
/// for anti-isolated regions
// up to Eta2p1

  int nEtaBins = 2;
  float etaBins[3] = {0, 1.48, 2.1};

  TString EtaBins[2] = {"EtaLt1p48","Eta1p48to2p1"};

  float ptBins_def[8] = {10,20,25,30,40,50,100,1000};
//pt binning for the triggers: 10., 20., 22., 24., 26., 28., 30., 40., 50., 100., 200., 1000.


  TString PtBins_def[7] = {"Pt10to20",
       "Pt20to25",
       "Pt25to30",
       "Pt30to40",
       "Pt40to50",
       "Pt50to100",
       "PtGt100"};

  float ptBinsTrig_def[11] = {10,
              20,
			  22, 
			  24,
			  26,
			  28,
			  30,
			  40,
			  50,
			  100,
			  1000};

	TString PtBinsTrig_def[10] = {"Pt10to20",
			"Pt20to22","Pt22to24",
		    "Pt24to26","Pt26to28",
		    "Pt28to30",
		    "Pt30to40",
		    "Pt40to50",
		    "Pt50to100",
			"PtGt100"};

	int nPtBins = 10; if(what == "IdIso") nPtBins = 7;*/



/*  int nEtaBins = 5;
  float etaBins[6] = {0, 1.0, 1.4442, 1.56, 2.1, 2.5};

  TString EtaBins[5] = {"EtaLt1p0","Eta1p0to1p44", "Eta1p44to1p56", "Eta1p56to2p1","EtaGt2p1"};

  float ptBins_def[9] = {10,20,25,30,40,50,100,200,1000};
//pt binning for the triggers: 10., 20., 22., 24., 26., 28., 30., 40., 50., 100., 200., 1000.


  TString PtBins_def[8] = {"Pt10to20",
       "Pt20to25",
       "Pt25to30",
       "Pt30to40",
       "Pt40to50",
       "Pt50to100",
       "Pt100to200",
       "PtGt200"};

  float ptBinsTrig_def[13] = {10,
              20,
			  22, 
			  24,
			  26,
			  27,
			  28,
			  30,
			  40,
			  50,
			  100,
			  200,
			  1000};

	TString PtBinsTrig_def[12] = {"Pt10to20",
			"Pt20to22","Pt22to24",
		    "Pt24to26","Pt26to27",
		    "Pt27to28",
		    "Pt28to30",
		    "Pt30to40",
		    "Pt40to50",
		    "Pt50to100",
			"Pt100to200",
			"PtGt200"};

	int nPtBins = 12; if(what == "IdIso") nPtBins = 8;*/


	float * ptBins = new float[nPtBins+1];
	TString * PtBins = new TString[nPtBins];

	if(what == "IdIso"){
		for(int i=0; i<nPtBins; ++i){
			ptBins[i] = ptBins_def[i];
			PtBins[i] = PtBins_def [i];
		}
		ptBins[nPtBins] = ptBins_def[nPtBins];
	} else {
		for(int i=0; i<nPtBins; ++i){
			ptBins[i] = ptBinsTrig_def[i];
			PtBins[i] = PtBinsTrig_def[i];
		}
		ptBins[nPtBins] = ptBinsTrig_def[nPtBins];
	}


  //	create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //

	TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  etaBinsH->Draw();


	//	create pt histogram_s with pt ranges associated to their names (eg. Pt10to13, ..)   ***** //

	TH1D * ptBinsH =  new TH1D("ptBinsH", "ptBinsH", nPtBins, ptBins);
  ptBinsH->Draw();
  ptBinsH->GetXaxis()->Set(nPtBins, ptBins);
  for (int i=0; i<nPtBins; i++){ ptBinsH->GetXaxis()->SetBinLabel(i+1, PtBins[i]);}
  ptBinsH->Draw();

  float ptBins_edges[nPtBins+1];

	for (int i=0; i<nPtBins; i++) { ptBins_edges[i]=ptBinsH->GetBinLowEdge(i+1); }
	ptBins_edges[nPtBins]= ptBinsH->GetBinLowEdge(nPtBins+1); 

  
	// define if in the fit of failing probes
  // the FSR component will be used in the 
  // signal function 
  bool fitWithFSR[nPtBins];

  for (int i=0; i<nPtBins; i++)  fitWithFSR[i] = false;

  if(what == "IdIso"){
   	fitWithFSR[2]=true;
    fitWithFSR[3]=true;
  } 

	// building the histogram base name
  TString prefix = "ZMass";
  TString which = what; if (what == "IdIso") which = "";
  which = "";
  TString histBaseName; 


  TCut cut_flag_idiso_pass = Form("id_probe == 1 && iso_probe < %f && iso_probe >=%f", iso_max, iso_min );
  
  TCut cut_flag_hlt_pass, cut_flag_hlt_fail, cut_pt, cut_eta;


  if(what == "hlt_1") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_1_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_1_probe == 0"; }
  if(what == "hlt_2") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_2_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_2_probe == 0 "; }
  //if(what == "hlt_3") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_3_probe == 1  && run<280385"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_3_probe == 0  && run<280385"; }
  // Ele25Tigth OR Ele27Loose
  if(what == "hlt_3") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && (hlt_3_probe == 1 || hlt_1_probe ==1)"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_3_probe == 0  && hlt_1_probe==0"; }
  if(what == "hlt_4") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_4_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_4_probe == 0"; }
  if(what == "hlt_5") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_5_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_5_probe == 0"; }
  if(what == "hlt_6") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_6_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_6_probe == 0"; }
  if(what == "hlt_7") {cut_flag_hlt_pass = "fabs(eta_tag)<2 && hlt_7_probe == 1"; cut_flag_hlt_fail = "fabs(eta_tag)<2 && hlt_7_probe == 0"; }
// for Ele24 leg of the e+tau trigger
  if(what == "hlt_8") {cut_flag_hlt_pass = "(fabs(eta_tag)<2 && hlt_tag_match ==1 && hlt_8_probe == 1 && levelone_match_probe ==1 && fabs(eta_probe)<2.1 && trigobjpt_probe>24)"; cut_flag_hlt_fail = "(fabs(eta_tag)<2 && (hlt_8_probe == 0 || levelone_match_probe == 0 || trigobjpt_probe<=24 || fabs(eta_probe)>=2.1) )"; }
  if(what == "hlt_9") {cut_flag_hlt_pass = "hlt_9_probe == 1"; cut_flag_hlt_fail = "hlt_9_probe == 0"; }
  if(what == "hlt_10") {cut_flag_hlt_pass = "hlt_10_probe == 1"; cut_flag_hlt_fail = "hlt_10_probe == 0"; }
  if(what == "hlt_11") {cut_flag_hlt_pass = "hlt_11_probe == 1"; cut_flag_hlt_fail = "hlt_11_probe == 0"; }
  if(what == "hlt_12") {cut_flag_hlt_pass = "hlt_12_probe == 1"; cut_flag_hlt_fail = "hlt_12_probe == 0"; }
  if(what == "hlt_13") {cut_flag_hlt_pass = "hlt_13_probe == 1"; cut_flag_hlt_fail = "hlt_13_probe == 0"; }
  if(what == "hlt_14") {cut_flag_hlt_pass = "hlt_14_probe == 1"; cut_flag_hlt_fail = "hlt_14_probe == 0"; }
  if(what == "hlt_15") {cut_flag_hlt_pass = "hlt_15_probe == 1"; cut_flag_hlt_fail = "hlt_15_probe == 0"; }
  if(what == "hlt_16") {cut_flag_hlt_pass = "hlt_16_probe == 1"; cut_flag_hlt_fail = "hlt_16_probe == 0"; }
  if(what == "hlt_17") {cut_flag_hlt_pass = "hlt_17_probe == 1"; cut_flag_hlt_fail = "hlt_17_probe == 0"; }
  if(what == "hlt_18") {cut_flag_hlt_pass = "hlt_18_probe == 1"; cut_flag_hlt_fail = "hlt_18_probe == 0"; }
  //if(what == "hlt_19") {cut_flag_hlt_pass = "hlt_19_probe == 1"; cut_flag_hlt_fail = "hlt_19_probe == 0"; }

  if(what == "hlt_19") {cut_flag_hlt_pass = "(fabs(eta_tag)<2 && hlt_19_probe == 1)"; cut_flag_hlt_fail = "(fabs(eta_tag)<2 && hlt_19_probe == 0 )"; }
  if(what == "hlt_20") {cut_flag_hlt_pass = "(fabs(eta_tag)<2 && hlt_20_probe == 1 )"; cut_flag_hlt_fail = "(fabs(eta_tag)<2 && hlt_20_probe == 0 )"; }

	TString dir_name = "Electron_";
	dir_name += what;
	//dir_name += Form("%.2f", iso); // "_IsoGt" + Form("%.2f", iso_min)  + "_IsoLt" + Form("%.2f", iso_max) 
    dir_name += Form("_Iso%.2fto%.2f", iso_min, iso_max);
	if (!isData) dir_name += "_MC";
	dir_name += "_eff";
	gSystem->mkdir(dir_name, kTRUE);

	for (int iEta = 0; iEta < nEtaBins; iEta++) {

		histBaseName = prefix+which+EtaBins[iEta];

		cut_eta = Form("abs(eta_probe)>= %f && abs(eta_probe)< %f", etaBins[iEta], etaBins[iEta+1]);

		TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins_edges);
		TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins_edges);
		TH1F * ratioH   = new TH1F("ratioH","",nPtBins,ptBins_edges);
	
	  for (int iPt=0; iPt<nPtBins; ++iPt) {

	  	cut_pt = Form("pt_probe >= %f && pt_probe < %f", ptBins[iPt], ptBins[iPt+1]);

	  	TH1F * histPassOld = new TH1F("histPassOld","",250,50,300);
	  	TH1F * histFailOld = new TH1F("histFailOld","",250,50,300);
	  	
	  	if (what == "IdIso") {
		  	t->Draw("m_vis>>histPassOld",  "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_idiso_pass));
		  	t->Draw("m_vis>>histFailOld",  "pu_weight*mcweight" + (cut_eta && cut_pt && !cut_flag_idiso_pass));
		  }else{
		  	t->Draw("m_vis>>histPassOld",  "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_hlt_pass && cut_flag_idiso_pass));
		  	t->Draw("m_vis>>histFailOld",  "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_hlt_fail && cut_flag_idiso_pass));
		  }

	  	int nBinsX = histPassOld->GetNbinsX();

	    for (int iB=1;iB<=nBinsX;++iB) {
	      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
	      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
	      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
	      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
	    }

	    float output[2];
	    TCanvas * c1 = new TCanvas("c1","",700,600);
	    TCanvas * c2 = new TCanvas("c2","",700,600);
	    bool fitPass = true; 
	    bool fitFail = true;
	    bool rebinPass = false;
	    bool rebinFail = false;
	    if(what != "IdIso") {fitPass= false; fitFail = false;}

	    FitPassAndFail(fileName,
	    	histBaseName+PtBins[iPt],
	    	xtit,
	    	histPassOld,
	    	histFailOld,
	    	fitPass,
	    	fitFail,
	    	fitWithFSR[iPt],
	    	rebinPass,
	    	rebinFail,
	    	c1,
	    	c2,
	    	output,
	    	dir_name);

	    c1->cd();
	    c1->Update();
	    c2->cd();
	    c2->Update();
	    numeratorH->SetBinContent(iPt+1,output[0]);
	    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);

		if (output[0] == 0. and output[1] ==0.) ratioH->SetBinContent(iPt+1,1.);
		else ratioH->SetBinContent(iPt+1,output[0]/(output[0]+output[1]));

		std::cout << "iPt+1 : " << iPt+1 << " num: " << output[0] << " bin content : " << numeratorH->GetBinContent(iPt+1) << std::endl;
		std::cout << "iPt+1 : " << iPt+1 << " den: " << output[0]+output[1] << " bin content : " << denominatorH->GetBinContent(iPt+1) << std::endl;
		std::cout << "ratio : " <<   output[0]/(output[0]+output[1]) << std::endl;
	  }

	  outputFile->cd();


	  //produce efficiencies plot
	  TGraphAsymmErrors * eff = new TGraphAsymmErrors(ratioH);
	  
	  //eff->Divide(numeratorH, denominatorH,"v");
	  std:cout << "eff->GetN() : " << eff->GetN() << std::endl;

	  //TGraphAsymmErrors * eff = new TGraphAsymmErrors();
	  //eff->Divide(numeratorH,denominatorH);

	  eff->GetXaxis()->SetTitle(xTitle);
	  //  eff->GetXaxis()->SetRangeUser(10.01,59.99);
	  eff->GetYaxis()->SetRangeUser(0,1.0);
	  eff->GetXaxis()->SetRangeUser(0,99.99);
	  eff->GetYaxis()->SetTitle(yTitle);
	  eff->GetXaxis()->SetTitleOffset(1.1);
	  eff->GetXaxis()->SetNdivisions(510);
	  eff->GetYaxis()->SetTitleOffset(1.1);
	  eff->SetMarkerStyle(21);
	  eff->SetMarkerSize(1);
	  eff->SetMarkerColor(kBlue);
	  eff->SetLineWidth(2);
	  eff->SetLineColor(kBlue);
	  if(!isData){
	  	eff->SetMarkerColor(kRed);
	  	eff->SetLineColor(kRed);
	  }


	  TCanvas * canv = new TCanvas("canv","",700,600);
	  eff->Draw("APE");
	  canv->SetGridx();
	  canv->SetGridy();
	  canv->Update();

	  canv->SaveAs(dir_name + "/" + fileName+"_" + histBaseName + ".png");
	  eff->Write(histBaseName+SampleName);

 	  std::cout << histBaseName+SampleName<< " : eff->GetN() " << eff->GetN() << std::endl;
//	  for(int ip=0; ip<nPtBins; ++ip){
//		cout<<"PtBins "<<ip<<" content: "<<numeratorH->GetBinContent(ip)/ denominatorH->GetBinContent(ip)<<endl;
//		}

	}



	outputFile->cd(); 
  etaBinsH->Write();
  outputFile->Close();

}
