//////////////////////////////////////////////////INSTRUCTIONS///////////////////////////////////////////////////////
//                                                                                                                    
//This macro loop over all the (nested)(sub)TDirectories inside a root file which name is passed as argument,       
//when a directory with a TH1 inside is found the function plot_hist() plots the histograms found in the directory, 
//then the macro will continue to loop over the other directories     
//
// Arguments: filename, (mandatory, name of the root file with the histogram)
//            channel name, (optional, passed as a Latex string e.g. "#mu#tau#", default: "#mu#tau#"")                                              
//            output directory name (optional, if passed all the output directories will be put in a new directory with whis name)
//                                                                                                                   
//WARNING: the directories with the TH1s have to be the most nested ones, i.e. with no other subdirectories.        
//                                                                                                                    
//Auxiliary functions definitions are at the bottom
//                                                           
//OUTPUT: the macro will create directories named "output_"+ /name of the directories with inside the variables 
//subdirectories/ and then put the plots and the latex tables inside these directories. 
//There is an overwrite control with the following options:
//1) Overwrite file
//2) Overwrite ALL
//3) Rename old file
//4) Rename ALL
//5) Abort
//This question is asked to the user only when the macro already find an output file.
//
//Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//                                                                                                                    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include <iostream>
#include <iomanip>

int g_over_flag;
bool g_over_all_flag = kFALSE;
TString g_ch = "#mu#mu";
TString g_dir_name = "";
TTimeStamp *time_stamp = new TTimeStamp(); 

void  SetAxis (TString var, TString *xtitle, TString *UdM);
void  SetBinErrorZero(TH1 *h);
float CountError (TH1 *h, float errLumi, float errSist);
void  plot_hist();
void  directory_loop(TDirectory *target);



void Plot_zmm(
	TString DataFile = "", //file name
  TString ch = "#mu#mu",
  TString dir_name_input = ""
	) {

	TFile *file = new TFile(DataFile);  //opening the TFile
	if (file->IsZombie()) {	
		std::cout << "File "<<DataFile<<" is Zombie" << endl ;
		exit(1);
	}
  g_ch = ch;
  g_dir_name = dir_name_input;
  if(g_dir_name != "") gSystem->mkdir(g_dir_name, kTRUE);
	directory_loop(file);

}


void directory_loop(TDirectory *target) {

	target->cd();

	TDirectory *current_sourcedir = gDirectory;		//inizializing the key and the iter
	TIter nextkey( current_sourcedir->GetListOfKeys() );
	TKey *key, *oldkey=0;

	while ( (key = (TKey*)nextkey())) {

		if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

		TObject *obj = key->ReadObj();

		if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      //if a TH1 is found the actual plotting function is called
			plot_hist();
			break;

		} else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
      //reiteration over the subdirectory
      directory_loop( (TDirectory*)obj );	
		}

	}
}

            ////////////MAIN PLOTTING FUNCTION/////////////

void plot_hist() {

  std::cout << " enter plottting function" << std::endl;
  TString Variable = gDirectory->GetName();
  TString path( (char*)strstr( gDirectory->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  TString dir_name = path( 0 , path.Last('/'));


 	SetStyle(); gStyle->SetOptTitle(0);
 	TH1::SetDefaultSumw2(kTRUE);

//	TString Variable = DataFile(DataFile.Last('-')+1, DataFile.Last('.')-DataFile.Last('-')-1);

  	TH1 * histData = (TH1*)gDirectory->Get("data_obs");
  	histData->SetStats(0); 	


    //binning and axis titles
    int nBins = histData ->GetNbinsX(); //-underflowbin -overflowbin
  	int xmin = histData->GetXaxis()->GetXmin();
    int xmax = histData->GetXaxis()->GetXmax();
  	float binWidth = (xmax-xmin)/nBins;
  	TString xtitle = "x title";  TString ytitle = "y title"; TString UdM = ""; 

//    SetAxis(Variable, &xtitle, &UdM);

    ytitle = Form("Events/ %.2f ", binWidth); ytitle += UdM; 

    TH1 * Z = (TH1*)gDirectory->Get("Z");
    TH1 * TT  = (TH1*)gDirectory->Get("TT");
    TH1 * W   = (TH1*)gDirectory->Get("W");
    TH1 * VV  = (TH1*)gDirectory->Get("VV");
    TH1 * QCD = (TH1*)gDirectory->Get("QCD");


    // merging all electroweak process to the W histogram
    W->Add(VV);

    //number of events for samples
    double nQCD	= QCD->Integral();
    std::cout << "QCD: " << nQCD;
    double nVV 	= VV->Integral();
    std::cout << "VV: " << nVV;
    double nW 	= W->Integral();
    std::cout << "W: " << nW;
    double nTT 	= TT->Integral();
    std::cout << "TT: " << nTT;
    double nZ	= Z->Integral();
    std::cout << "Z: " << nZ;

    //errors values

    float errLumi = 0.0;
    float errQCD = 0.0;
    float errVV = 0.0;
    float errW = 0.0;
    float errTT = 0.0;
    float errZ = 0.0;


    TH1 * dummy = (TH1*)Z->Clone("dummy");

    for (int iB=1; iB<=nBins; ++iB) {
      float eQCD = errQCD*QCD->GetBinContent(iB);
      float eW = errW*W->GetBinContent(iB);
      float eTT = errTT*TT->GetBinContent(iB);
      float eZ = errZ*Z->GetBinContent(iB);
      float err2 = eQCD*eQCD + eW*eW + eTT*eTT + eZ*eZ;
      float errTot = TMath::Sqrt(err2);
      dummy->SetBinError(iB,errTot);
    }

    float etotQCD = CountError(QCD, errLumi, errQCD);
    float etotW = CountError(W, errLumi, errW);
    float etotTT = CountError(TT, errLumi, errTT);
    float etotZ = CountError(Z, errLumi, errZ);

    //hist sum
    W->Add(W,QCD);
    TT->Add(TT,W);
    Z->Add(Z,TT);

    TH1 * bkgdErr = (TH1*)Z->Clone("bkgdErr");
    bkgdErr->SetFillStyle(3013);
    bkgdErr->SetFillColor(1);
    bkgdErr->SetMarkerStyle(21);
    bkgdErr->SetMarkerSize(0); 

    for (int iB=1; iB<=nBins; ++iB) {
      float eStat =  bkgdErr->GetBinError(iB);
      float X = bkgdErr->GetBinContent(iB);
      float eLumi = errLumi * X;    
      float eBkg = dummy->GetBinError(iB);
      float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg);
      bkgdErr->SetBinError(iB, Err);
    }

    SetBinErrorZero(QCD);
    SetBinErrorZero(W);
    SetBinErrorZero(TT);
    SetBinErrorZero(Z);
    //graphics
    InitData(histData);
    InitHist(QCD,"","",TColor::GetColor("#FFCCFF"),1001);
    InitHist(W,"","",TColor::GetColor("#DE5A6A"),1001);
    InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
    InitHist(Z,"","",TColor::GetColor("#4496C8"),1001);

    histData->GetYaxis()->SetTitleOffset(1.3);
    histData->GetYaxis()->SetTitleSize(0.06);
    float yUpper = histData->GetMaximum();
    if(histData->GetMaximum() < Z->GetMaximum()) yUpper = Z->GetMaximum();
    histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);

   	histData->SetMarkerSize(1.2);
    histData->GetXaxis()->SetLabelSize(0);
    histData->GetYaxis()->SetLabelSize(0.06);
    histData->GetYaxis()->SetTitleOffset(1.6);

    TCanvas * c1 = new TCanvas("c1", "", 700, 800);
    	
   	TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  	upper->Draw();
    upper->cd();
    upper->SetFillColor(0);
    upper->SetBorderMode(0);
    upper->SetBorderSize(10);
    upper->SetTickx(1);
    upper->SetTicky(1); upper->SetLeftMargin(0.20);
    upper->SetRightMargin(0.05);
    upper->SetBottomMargin(0.02);
    upper->SetFrameFillStyle(0);
    upper->SetFrameLineStyle(0);
    upper->SetFrameLineWidth(2);
    upper->SetFrameBorderMode(0);
    upper->SetFrameBorderSize(10);
    upper->SetFrameFillStyle(0);
    upper->SetFrameLineStyle(0);
    upper->SetFrameLineWidth(2);
    upper->SetFrameBorderMode(0);
    upper->SetFrameBorderSize(10);

    //DRAWING
  	histData->Draw("e1");
  	Z->Draw("sameh");
   	TT->Draw("sameh");
    W->Draw("sameh");
  	QCD->Draw("sameh");
  	histData->Draw("e1same");  
    bkgdErr->Draw("e2same");

  	//CHI2
  	int dof=0;
    float chi2 = 0;
    for (int iB=1; iB<=nBins; ++iB) {
      	float xData = histData->GetBinContent(iB);
      	float xMC = Z->GetBinContent(iB);
      	float xErr = bkgdErr->GetBinError(iB);
      	if (xMC>1e-1) {
          dof++;
      		float diff2 = pow(xData-xMC,2);
      		chi2 += diff2/(xErr*xErr);
      	}
    }


  	//LEGEND
  	TLegend * leg = new TLegend(0.6,0.4,0.82,0.78);
    SetLegendStyle(leg);
    leg->SetTextSize(0.05);
    leg->AddEntry(histData,"Observed","lp");
    leg->AddEntry(Z,"Z#rightarrow #mu#mu","f");
    leg->AddEntry(TT,"t#bar{t}","f");
    leg->AddEntry(W,"Electroweak","f");
    leg->AddEntry(QCD,"QCD","f");

    writeExtraText = false;
    //extraText = "Work in Progress      ";
    CMS_lumi(upper,4,33); 

    plotchannel(g_ch);


    TLatex t; 
    t.SetNDC();

    // Drawing legend and show chi2
    bool plotLegend = false;
	bool plotChi2 = false;
	bool LogYscale = true;

	/*if ((Variable != "eta1")&&(Variable != "eta2")) {
		leg->Draw();
		t.DrawLatex(.6,.3, Form("#bf{#chi^{2}/d.o.f. = %.2f}", chi2/dof) );
	}*/

	if (plotLegend) leg->Draw();
	if (plotChi2) t.DrawLatex(.6,.3, Form("#bf{#chi^{2}/d.o.f. = %.2f}", chi2/dof) );

    upper->Draw("SAME");

	if(LogYscale) {histData->SetMinimum(10e-3); upper->SetLogy();}

    upper->RedrawAxis();
    upper->Modified();
    upper->Update();
    c1->cd();

    //RATIO

    TH1 * ratioH = (TH1*)histData->Clone("ratioH");
    TH1 * ratioErrH = (TH1*)bkgdErr->Clone("ratioErrH");
    ratioH->SetMarkerColor(1);
    ratioH->SetMarkerStyle(20);
    ratioH->SetMarkerSize(1.2);
    ratioH->SetLineColor(1);
    ratioH->GetYaxis()->SetRangeUser(0.4,1.6);
    ratioH->GetYaxis()->SetNdivisions(505);
    ratioH->GetXaxis()->SetLabelFont(42);
    ratioH->GetXaxis()->SetLabelOffset(0.04);
    ratioH->GetXaxis()->SetLabelSize(0.14);
    ratioH->GetXaxis()->SetTitleSize(0.13);
    ratioH->GetXaxis()->SetTitleOffset(1.2);
    ratioH->GetYaxis()->SetTitle("obs/exp");
    ratioH->GetYaxis()->SetLabelFont(42);
    ratioH->GetYaxis()->SetLabelOffset(0.015);
    ratioH->GetYaxis()->SetLabelSize(0.13);
    ratioH->GetYaxis()->SetTitleSize(0.14);
    ratioH->GetYaxis()->SetTitleOffset(0.5);
    ratioH->GetXaxis()->SetTickLength(0.07);
    ratioH->GetYaxis()->SetTickLength(0.04);
    ratioH->GetYaxis()->SetLabelOffset(0.01);

    for (int iB=1; iB<=nBins; ++iB) {
      float x1 = histData->GetBinContent(iB);   float x2 = Z->GetBinContent(iB);
      ratioErrH->SetBinContent(iB,1.0);
      ratioErrH->SetBinError(iB,0.0);
      float xBkg = bkgdErr->GetBinContent(iB);
      float errBkg = bkgdErr->GetBinError(iB);
      if (xBkg>0) {
        float relErr = errBkg/xBkg;
        ratioErrH->SetBinError(iB,relErr);
      }
      if (x1>0&&x2>0) {
        float e1 = histData->GetBinError(iB);
        float ratio = x1/x2;
        float eratio = e1/x2;
        ratioH->SetBinContent(iB,ratio);
        ratioH->SetBinError(iB,eratio);
      } else {ratioH->SetBinContent(iB,1000);}
  	}

  	TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
    	lower->Draw();
    	lower->cd();
    	lower->SetFillColor(0);
    	lower->SetBorderMode(0);
    	lower->SetBorderSize(10);
    	lower->SetGridy();
    	lower->SetTickx(1);
    	lower->SetTicky(1);
    	lower->SetLeftMargin(0.20);
    	lower->SetRightMargin(0.05);
    	lower->SetTopMargin(0.026);
    	lower->SetBottomMargin(0.35);
    	lower->SetFrameFillStyle(0);
    	lower->SetFrameLineStyle(0);
    	lower->SetFrameLineWidth(2);
    	lower->SetFrameBorderMode(0);
    	lower->SetFrameBorderSize(10);
    	lower->SetFrameFillStyle(0);
    	lower->SetFrameLineStyle(0);
    	lower->SetFrameLineWidth(2);
    	lower->SetFrameBorderMode(0);
    	lower->SetFrameBorderSize(10);

    	ratioH->Draw("e1");
    	ratioErrH->Draw("e2same");

    	lower->Modified();
    	lower->RedrawAxis();

    	c1->cd();
    	c1->Modified();
    	c1->cd();
    	c1->SetSelected(c1);
      
      TString plot_name = "plot_" + Variable + ".png";
      TString plot_name_copy = plot_name;

      TString output_path;

      if(g_dir_name != "") {
        gSystem->mkdir(g_dir_name + "/" + "output_" + dir_name, kTRUE);
        output_path = g_dir_name + "/" + "output_" + dir_name + "/";
      } else {
        gSystem->mkdir("output_" + dir_name, kTRUE);
        output_path = "output_" + dir_name + "/";
      } 
      
      //controls if outputfile already exists
      const char* fileExists = gSystem->FindFile((const char*)output_path, plot_name_copy);
      if ( !(fileExists == 0)){ 
        if (g_over_all_flag) { 

          if(g_over_flag == 2) c1->SaveAs(output_path + plot_name); //if ALL flag =  true then always overwrite
          if(g_over_flag == 4) {
              gSystem->Rename(output_path + plot_name, output_path + "plot_" + Variable + "_"+ Form("%i", time_stamp->GetTime()) +".png" );
              c1->SaveAs(output_path + plot_name);
          } 
        
        } else {
          cout<<endl<<"Output file "<<plot_name<<" in path"<<endl<<output_path<<endl<<"already exists."<<endl<<" What do you want to do?"<<endl;
          cout<<"1 Overwrite"<<endl<<"2 Overwrite ALL"<<endl<<"3 Rename old file (a timestamp will be added) (latex table will be renamed)"<<endl<<"4 Rename old files ALL"<<endl<<"5 Abort process"<<endl<<endl;
          cout<<"Enter your choice: ";
          cin>>g_over_flag;
          cout<<endl;
          
          switch(g_over_flag){

            case(1):
              cout<<"You choose to overwrite "<<plot_name<<" in "<<output_path<<endl;
              c1->SaveAs(output_path + plot_name);
              break;

            case(2):
              cout<<"You chose to overwrite all the existent plot."<<endl;
              c1->SaveAs(output_path + plot_name);
              g_over_all_flag = kTRUE;
              break;

            case(3):
              cout<<"You chose to rename the existent plot and save the new one."<<endl;
              gSystem->Rename(output_path + plot_name, output_path + "plot_" + Variable + "_"+ Form("%i", time_stamp->GetTime(kFALSE)) +".png" );
              c1->SaveAs(output_path + plot_name);
              break;

            case(4):
              cout<<"You chose to rename all the existent plots."<<endl;
              gSystem->Rename(output_path + plot_name, output_path + "plot_" + Variable + "_"+ Form("%i", time_stamp->GetTime(kFALSE)) +".png" );
              c1->SaveAs(output_path + plot_name);
              ::g_over_all_flag = kTRUE;
              break;

            case(5):
              cout<<"Process aborted."<<endl;
              exit(0);
              break;

            default: 
              cout<<"Choice not valid. Process aborted."<<endl;
              exit(1);

          }
        }


      } else {c1->SaveAs(output_path + plot_name); }



    // LEGEND ON SEPARATE CANVAS 
    TString legend_name = "legend";
    if( (g_over_flag == 3) || (g_over_flag == 4) ){ 
      gSystem->Rename(output_path + legend_name, output_path + "legend" + Form("%i", time_stamp->GetTime(kFALSE)) +".png" );
    }
    TCanvas * legendCanvas = new TCanvas("legendCanvas", "", 700, 800);
	legendCanvas->cd();
	leg->Draw();
    legendCanvas->SaveAs(output_path+legend_name+".png");
}



												///////////FUNCTIONS//////////

void SetAxis (TString var, TString *xtitle, TString *UdM) {

  	if (var == "mvis") {*xtitle= "m_{vis}"; *UdM = "GeV";}
  	if (var == "pt_1") {*xtitle = "muon p_{T}"; *UdM = "GeV";}
  	if (var == "pt_2")  {*xtitle = "tau p_{T}"; *UdM = "GeV";}
  	if (var == "eta_1") {*xtitle = "muon #eta"; *UdM = "";} 
  	if (var == "eta_2") {*xtitle = "tau #eta"; *UdM = "";} 
  	if (var == "met") {*xtitle = "MET"; *UdM = "GeV";}
  	if (var == "puppimet") {*xtitle = "Puppi MET"; *UdM = "GeV";}
  	if (var == "mvamet") {*xtitle = "MVA MET"; *UdM = "GeV";}
  	if (var == "mt_1"){*xtitle= "m_{T,1}"; *UdM = "GeV";}
  	if (var == "pfmt_1"){*xtitle= "PF m_{T,1}"; *UdM = "GeV";}
  	if (var == "puppimt_1"){*xtitle= "Puppi m_{T,1}"; *UdM = "GeV";}
  	if (var == "npv"){*xtitle= "Number of primary vertices"; *UdM = "";}

} 


void SetBinErrorZero(TH1 *h) { //accepts a TH1 and set all the Bin Errors to 0

	int nBins = h->GetNbinsX();
	for(int iB=1; iB<=nBins; ++iB){
		h->SetBinError(iB, 0);
	}

}


float CountError (TH1 *h, float errLumi, float errSist) { //calculate the integral error for a histogram
	
	float errTot = 0;
	float errSistTot=0;
	float errLumiTot=0;
	float errStatTot=0;
	int nBins = h->GetNbinsX();
	for(int iB=1; iB<=nBins; ++iB){
		float eStat =  h->GetBinError(iB);
		float X = h->GetBinContent(iB);
	    float eLumi = errLumi * X;
	    float eh = errSist * X;

	    errSistTot += eh;
	    errLumiTot += eLumi;
	    errStatTot += eStat*eStat;
	}
	errStatTot=TMath::Sqrt(errStatTot);
	errTot=TMath::Sqrt(errSistTot*errSistTot+errLumiTot*errLumiTot+errStatTot*errStatTot);
	return errTot;
}
