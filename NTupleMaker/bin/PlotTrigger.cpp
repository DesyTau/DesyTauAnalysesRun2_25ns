#include "HttStylesNew.cc"
#include <TString.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include "Functions.h"

void PlotCrossTrigger(TString WP,
		      TTree * tree,
		      int iEta,
		      int iWP,
		      bool embedded,
		      int nBins,
		      float * bins,
		      int DM) {

  SetStyle();

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TString decayMode = "inclusive";
  TString dmCut = "";
  TString DMLabel = "_incl";
  TString Sample = "Embedded";
  TString SampleLabel = "_emb";
  if (DM==0) {
    decayMode = "1-prong";
    dmCut = "&&tauDecay==0";
    DMLabel = "_1pr";
  }
  if (DM==1) {
    decayMode = "1-prong+#pi^{0}s";
    dmCut = "&&(tauDecay==1||tauDecay==2)";
    DMLabel = "_1pr1pi0";
  }
  if (DM==2) {
    decayMode = "3-prong";
    dmCut = "&&tauDecay>=10";
    DMLabel = "_3pr";
  }
  if (!embedded) {
    Sample = "MC SUSY";
    SampleLabel = "_mc";
  }

  TString cutsFilterSingle("tauSinglePFTau180Trk50");
  TString cutsFilterDouble("(tauMuTauTrigger1||tauMuTauTrigger2||tauMuTauTrigger3)");
  TString cutsFilter = cutsFilterSingle + "&&" + cutsFilterDouble;

  TString cutsBase("ngentaus==2&&deltaR>0.5&&genTauDecay>=0&&genTauDecay<8&&genTauFoundReco&&tauNewDM>0.5&&tauby"+WP+"DeepTau2017v2p1VSjet&&tauPt>80&&TMath::Abs(tauEta)<2.1");
 
  cutsBase = cutsBase + "&&" + etaCut[iEta] + "&&" + WPCut[iWP];

  TString cutsMuTau = ("&&genTau1Decay==8&&genTau1FoundReco&&tau1Pt>28&&tau1Iso<0.15&&tau1SingleMuon&&TMath::Abs(tau1Eta)<2.1");

  cutsBase += dmCut;
  cutsBase += cutsMuTau;
  
  TString cutsPass  = cutsBase + "&&("+cutsFilter+")";
  TString cutsFail  = cutsBase + "&&!("+cutsFilter+")";

  TString cutsPassSingle  = cutsBase + "&&("+cutsFilterSingle+")";
  TString cutsFailSingle  = cutsBase + "&&!("+cutsFilterSingle+")";

  TString cutsPassDouble  = cutsBase + "&&("+cutsFilterDouble+")";
  TString cutsFailDouble  = cutsBase + "&&!("+cutsFilterDouble+")";

  TString weight("puWeight*");
  if (embedded)
    weight = "genWeight*embWeight*"; 

  TString var("tauPt");
  TCanvas * dummy = new TCanvas("dummy","",600,600);

  std::cout << "OK both" << std::endl;
  TH1D * passH = Make1DHistoFromTree(tree,var,"passH",nBins,bins,weight+"("+cutsPass+")");
  std::cout << "OK1 both" << std::endl;
  TH1D * failH = Make1DHistoFromTree(tree,var,"failH",nBins,bins,weight+"("+cutsFail+")");
  std::cout << "OK2 both" << std::endl;
  TH1D * passSingleH = Make1DHistoFromTree(tree,var,"passSingleH",nBins,bins,weight+"("+cutsPassSingle+")");
  std::cout << "OK3 both" << std::endl;
  TH1D * failSingleH = Make1DHistoFromTree(tree,var,"failSingleH",nBins,bins,weight+"("+cutsFailSingle+")");
  std::cout << "OK4 both" << std::endl;
  TH1D * passDoubleH = Make1DHistoFromTree(tree,var,"passDoubleH",nBins,bins,weight+"("+cutsPassDouble+")");
  std::cout << "OK5 both" << std::endl;
  TH1D * failDoubleH = Make1DHistoFromTree(tree,var,"failDoubleH",nBins,bins,weight+"("+cutsFailDouble+")");
  std::cout << "OK6 both" << std::endl;


  delete dummy;

  TH1D * effH = make1DEfficiency(passH,failH,"effH");
  TH1D * effSingleH = make1DEfficiency(passSingleH,failSingleH,"effSingleH");
  TH1D * effDoubleH = make1DEfficiency(passDoubleH,failDoubleH,"effDoubleH");

  effH->SetLineColor(1);
  effH->SetMarkerColor(1);
  effH->SetMarkerSize(0);

  effSingleH->SetLineColor(2);
  effSingleH->SetMarkerColor(2);
  effSingleH->SetMarkerSize(0);

  effDoubleH->SetLineColor(4);
  effDoubleH->SetMarkerColor(4);
  effDoubleH->SetMarkerSize(0);

  TH2D * frame = new TH2D("frame","",2,30,2000,2,0,1.1);
  frame->GetXaxis()->SetTitle("#tau p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("#epsilon(trigger)");
  
  TCanvas * canv = new TCanvas("canv","",700,700);  
  frame->Draw();
  effH->Draw("e1same");
  effDoubleH->Draw("e1same");
  effSingleH->Draw("e1same");

  TLegend * leg = new TLegend(0.4,0.2,0.9,0.4);
  leg->SetHeader(Sample + "  " + decayMode);
  leg->AddEntry(effH,"single-#tau AND di-#tau ","elp");
  leg->AddEntry(effSingleH,"only single-#tau ","elp");
  leg->AddEntry(effDoubleH,"only di-#tau ","elp");
  //  writeExtraText = true;
  //  extraText = "Simulation";
  //  CMS_lumi(canv,4,33); 
  canv->SetLogx(true);
  canv->SetGridx(true);
  canv->SetGridy(true);
  leg->Draw();
  canv->Update();
  canv->Print("figures/tauCrossTrigger_"+SampleLabel+DMLabel+etaLabel[iEta]+WPLabel[iWP]+".png");

  TFile * fileOutput = new TFile("eff/crossTrg_"+SampleLabel+DMLabel+etaLabel[iEta]+WPLabel[iWP]+".root","recreate");
  fileOutput->cd("");
  effH->Write("eff_sd"+SampleLabel+DMLabel+etaLabel[iEta]+WPLabel[iWP]);
  effDoubleH->Write("eff_d"+SampleLabel+DMLabel+etaLabel[iEta]+WPLabel[iWP]);
  fileOutput->Close();
  delete fileOutput;
  delete canv;

}

void PlotTriggerWP(TString WP,
		   bool fullReco,
		   bool embedded,
		   TTree * treeMuTau,
		   TTree * treeElTau,
		   int nBins,
		   float * bins,
		   int DM) {

  SetStyle();

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TString decayMode = "inclusive";
  TString dmCut = "";
  TString DMLabel = "_incl";
  TString Sample = "Embedded";
  TString SampleLabel = "_emb";
  if (DM==0) {
    decayMode = "1-prong";
    dmCut = "&&tauDecay==0";
    DMLabel = "_1pr";
  }
  if (DM==1) {
    decayMode = "1-prong+#pi^{0}s";
    dmCut = "&&(tauDecay==1||tauDecay==2)";
    DMLabel = "_1pr1pi0";
  }
  if (DM==2) {
    decayMode = "3-prong";
    dmCut = "&&tauDecay>=10";
    DMLabel = "_3pr";
  }
  if (!embedded) {
    Sample = "MC SUSY";
    SampleLabel = "_mc";
  }

  TString cutsFilter("tauSinglePFTau180Trk50");

  TString cutsMuTau("&&genTau1Decay==8");
  TString cutsElTau("&&genTau1Decay==9");

  if (fullReco) {
    cutsMuTau = "&&genTau1Decay==8&&genTau1FoundReco&&tau1Pt>25&&tau1Iso<0.15&&tau1SingleMuon";
    cutsElTau = "&&genTau1Decay==9&&genTau1FoundReco&&tau1Pt>30&&tau1Iso<0.10&&tau1SingleElectron";
  }

  TString cutsBase("ngentaus==2&&deltaR>0.5&&genTauDecay>=0&&genTauDecay<8&&genTauFoundReco&&tauNewDM>0.5&&tauby"+WP+"DeepTau2017v2p1VSjet&&tauPt>100&&TMath::Abs(tauEta)<2.1");

  cutsBase = cutsBase + dmCut;
  TString cutsBaseMuTau = cutsBase + cutsMuTau;
  TString cutsBaseElTau = cutsBase + cutsElTau;

  int color[3] = {1,2,4};
  TString cutsMuTauPass[3][3];
  TString cutsMuTauFail[3][3];
  TString cutsElTauPass[3][3];  
  TString cutsElTauFail[3][3];
  for (unsigned int i = 0; i<3; ++i) {
    for (unsigned int j = 0; j<3; ++j) {
      cutsMuTauPass[i][j] = cutsBaseMuTau+"&&("+cutsFilter+")&&"+WPCut[j]+"&&"+etaCut[i];
      cutsMuTauFail[i][j] = cutsBaseMuTau+"&&!("+cutsFilter+")&&"+WPCut[j]+"&&"+etaCut[i];
      cutsElTauPass[i][j] = cutsBaseElTau+"&&("+cutsFilter+")&&"+WPCut[j]+"&&"+etaCut[i];
      cutsElTauFail[i][j] = cutsBaseElTau+"&&!("+cutsFilter+")&&"+WPCut[j]+"&&"+etaCut[i];
    }
  }

  TString var("tauPt");

  TString weight("puWeight*");
  if (embedded)
    weight = "genWeight*embWeight*"; 

  TCanvas * dummy = new TCanvas("dummy","",600,600);

  TH1D * passMuTauH[3][3];
  TH1D * failMuTauH[3][3];
  TH1D * passElTauH[3][3];
  TH1D * failElTauH[3][3];
  TH1D * effH[3][3];

  for (unsigned int i = 0; i<1; ++i) {
    for (unsigned int j = 0; j<1; ++j) {
      passMuTauH[i][j] = Make1DHistoFromTree(treeMuTau,var,"pass"+etaLabel[i]+WPLabel[j]+"_MuTauH",nBins,bins,weight+"("+cutsMuTauPass[i][j]+")");
      failMuTauH[i][j] = Make1DHistoFromTree(treeMuTau,var,"fail"+etaLabel[i]+WPLabel[j]+"_MuTauH",nBins,bins,weight+"("+cutsMuTauFail[i][j]+")");
      passElTauH[i][j] = Make1DHistoFromTree(treeElTau,var,"pass"+etaLabel[i]+WPLabel[j]+"_ElTauH",nBins,bins,weight+"("+cutsElTauPass[i][j]+")");
      failElTauH[i][j] = Make1DHistoFromTree(treeElTau,var,"fail"+etaLabel[i]+WPLabel[j]+"_ElTauH",nBins,bins,weight+"("+cutsElTauFail[i][j]+")");

      passMuTauH[i][j]->Add(passMuTauH[i][j],passElTauH[i][j]);
      failMuTauH[i][j]->Add(failMuTauH[i][j],failElTauH[i][j]);
      effH[i][j] = make1DEfficiency(passMuTauH[i][j],failMuTauH[i][j],"eff"+etaLabel[i]+WPLabel[i]);
    //    if (embedded) {
    //      for (int iB=8; iB<=9; ++iB) {
    //    	double x = effH[i]->GetBinContent(iB);
    //    	x = 0.5*(1.0+x);
    //    	effH[i]->SetBinContent(iB,x);
    //      }
    }
  }
  TH2D * frame = new TH2D("frame","",2,bins[0],bins[nBins],2,0,1.25);
  frame->GetXaxis()->SetTitle("#tau p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("single-#tau trigger efficiency");

  for (int i=0; i<3; ++i) {
    effH[0][i]->SetLineColor(color[i]);
    effH[0][i]->SetMarkerColor(color[i]);
  
  }
  effH[0][0]->SetMarkerSize(1.6);

  TCanvas * canv = new TCanvas("canv","",700,700);  
  frame->Draw();
  effH[0][0]->Draw("e1same");
  effH[0][1]->Draw("e1same");
  effH[0][2]->Draw("e1same");

  TLegend * leg = new TLegend(0.4,0.2,0.9,0.4);
  leg->SetHeader(Sample + " : " + decayMode);
  leg->AddEntry(effH[0][0],"#tau+#tau WP","elp");
  leg->AddEntry(effH[0][1],"#mu+#tau WP","elp");
  leg->AddEntry(effH[0][2],"e+#tau WP","elp");
  //  writeExtraText = true;
  //  extraText = "Simulation";
  //  CMS_lumi(canv,4,33); 
  canv->SetLogx(true);
  canv->SetGridx(true);
  canv->SetGridy(true);
  leg->Draw();
  canv->Update();
  canv->Print("figures/SingleTauTrigger_"+Sample+DMLabel+"_wp.png");

  // Second plot

  TCanvas * canv1 = new TCanvas("canv1","",700,700);  

  for (int i=1; i<3; ++i) {
    effH[i][0]->SetLineColor(color[i]);
    effH[i][0]->SetMarkerColor(color[i]);
    effH[i][0]->SetMarkerSize(0.);
    
  }
  effH[1][0]->SetMarkerSize(1.6);
  
  frame->Draw();
  effH[1][0]->Draw("e1same");
  effH[2][0]->Draw("e1same");
  TLegend * leg1 = new TLegend(0.4,0.2,0.9,0.4);
  leg1->SetHeader(Sample + " : " + decayMode);
  leg1->AddEntry(effH[1][0],"Barrel","elp");
  leg1->AddEntry(effH[2][0],"Endcap","elp");
  canv1->SetLogx(true);
  canv1->SetGridx(true);
  canv1->SetGridy(true);
  leg1->Draw();
  canv1->Update();
  canv1->Print("figures/SingleTauTrigger_"+SampleLabel+DMLabel+"_eta.png");

  delete canv;
  delete canv1;
  
  TFile * fileOutput = new TFile("eff/eff_singleTau"+SampleLabel+DMLabel+".root","recreate");
  fileOutput->cd();
  for (unsigned int i=0; i<3; ++i) { 
    for (unsigned int j=0; j<3; ++j) {
      effH[i][j]->Write("eff_s"+SampleLabel+DMLabel+etaLabel[i]+WPLabel[j]);
    }
  }
  fileOutput->Close();
  delete fileOutput;
}  

int main(int argc, char * argv[]) {

  TString era(argv[1]);
  
  TString fileNameMuTau = "EmbeddingRun"+era+"_MuTau.root";
  TString fileNameElTau = "EmbeddingRun"+era+"_ElTau.root";
  TString fileNameMC("SUSYGluGluToHToTauTau_M_200-3200.root");

  TFile * fileMuTau = new TFile(fileNameMuTau);
  TFile * fileElTau = new TFile(fileNameElTau);
  TFile * fileMC = new TFile(fileNameMC);

  TTree * treeMuTau = (TTree*)fileMuTau->Get("NTuple");
  TTree * treeElTau = (TTree*)fileElTau->Get("NTuple");
  TTree * treeMC = (TTree*)fileMC->Get("NTuple");
  
  int nBinsCross = 11;
  float binsCross[12] = {80,90,100,120,140,160,180,200,225,250,300,2000};

  int nBins = 12;
  float bins[13] = {80,100,120,140,160,180,200,225,250,300,400,700,2000};

  for (int iDM=-1; iDM<3; ++iDM) {
    for (int iEta=0; iEta<1; ++iEta) {
      for (int iWP=0; iWP<1; ++iWP) {
	PlotCrossTrigger("Medium",treeMC,iEta,iWP,false,nBinsCross,binsCross,iDM);
	std::cout << "Simulation : iDM = " << iDM << " iEta = " << iEta << " iWP = " << iWP << std::endl;
	PlotCrossTrigger("Medium",treeMuTau,iEta,iWP,true,nBinsCross,binsCross,iDM);
	std::cout << "Embedded   : iDM = " << iDM << " iEta = " << iEta << " iWP = " << iWP << std::endl;
      }
    }
  }


}
