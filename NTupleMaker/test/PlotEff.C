#include "tdrstyle.C" 
void PlotEff(TString fileName = "DYJetsEff",
	     bool posTag = true,
	     bool isMuonId = true) {

  setTDRStyle();

  TString mass("massEEId_");
  if (isMuonId) 
    mass = "massMuMuId_";

  TString TagQ("Neg");
  if (posTag)
    TagQ = "Pos";

  int nEtaBins = 3;
  TString EtaBinName[3] = {"0To0p8",
			   "0p8To1p5",
			   "1p5To2p3"};
  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};
  TString PtBinName[7] = {"pt10to15",
			  "pt15to20",
			  "pt20to30",
			  "pt30to40",
			  "pt40to50",
			  "pt50to70",
			  "pt70to100"};
  
  TString passName[2] = {"Pass","Fail"};

  TFile * file = new TFile(fileName+".root");

  TH1F * histos[3][2];
  for (int iEta=0; iEta<3; ++iEta) {
    for (int iPass=0; iPass<2; ++iPass) {
      histos[iEta][iPass] = new TH1F(passName[iPass]+EtaBinName[iEta],"",nPtBins,ptBins);
      for (int iPt=0; iPt<7; ++iPt) {
	TH1F * hist = (TH1F*)file->Get(mass+TagQ+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass]);
	float yield = hist->GetSum();
	histos[iEta][iPass]->SetBinContent(iPt+1,yield);
	// std::cout << mass << " " << TagQ << " " 
	// 	  << " " << EtaBinName[iEta] << " " << PtBinName[iPt] 
	// 	  << " " << passName[iPass] << " : " << yield << std::endl;
      }
    }
  }

  int color[3] = {1,2,4};
  int symbol[3] = {20,21,22};

  TGraphAsymmErrors * eff[3];
  for (int iEta=0; iEta<3; ++iEta) {
    eff[iEta] = new TGraphAsymmErrors();
    eff[iEta]->SetLineColor(color[iEta]);
    eff[iEta]->SetMarkerColor(color[iEta]);
    eff[iEta]->SetLineWidth(2);
    eff[iEta]->SetMarkerSize(1.6);
    eff[iEta]->SetMarkerStyle(symbol[iEta]);
    histos[iEta][1]->Add(histos[iEta][0],histos[iEta][1]);
    eff[iEta]->Divide(histos[iEta][0],histos[iEta][1],"Pois");
  }  

 
  TH2F * frame = new TH2F("frame","",2,0,100,2,0,1);
  TCanvas * canv = new TCanvas("canv","",700,700);
  frame->Draw();
  eff[0]->Draw("PS");
  eff[1]->Draw("PS");
  eff[2]->Draw("PS");
  canv->Update();
  


}
