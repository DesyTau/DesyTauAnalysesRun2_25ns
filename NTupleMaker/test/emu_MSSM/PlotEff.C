#include "HttStylesNew.cc"

double getEff(TFile * file, TString Cuts, TString Weight) {
  TH1D * hist = new TH1D("hist","",100,0,4000);
  TTree * tree = (TTree*)file->Get("TauCheck");
  TH1D * histWeightsH = (TH1D*)file->Get("nWeightedEvents");
  double nevents = histWeightsH->GetSumOfWeights(); 
  tree->Draw("mt_tot_puppi>>hist",Weight+"("+Cuts+")");
  double eff = 100.0*hist->GetSumOfWeights()/nevents;
  delete hist;
  return eff;
}

void PlotEff(TString era ="2018",
	     TString sig="bbH", 
	     bool useSF = false) {

  TString Sample = "SUSYGluGluToHToTauTau_M-";
  if (sig=="bbH")
    Sample = "SUSYGluGluToBBHToTauTau_M-";

  SetStyle();
  unsigned int nPoints = masses.size();
  double x[50];
  double sig_single[50];
  double sig_emu[50];
  double sig_tot[50];
  double sig_singleand[50];
  double sig_and[50];
  double rat[50];
  double ratioTot[50];

  TString dir = "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan1/"+era;

  TString Weight("puweight*mcweight*prefiringweight*");
  TString WeightSingle = Weight + "effweight*";
  TString WeightEMu = Weight + "effweightExcl*";
  TString WeightTot = Weight + "effweightExcl*";
  if (!useSF) {
    WeightSingle = "(" +Weight + "effweight/trigweight)*";
    WeightEMu = "(" +Weight + "effweightExcl/trigweightExcl)*";
    WeightTot = "(" +Weight + "effweightExcl/trigweightExcl)*";
  }
  
  TString CutsEMu("((trg_muhigh_elow>0.5&&pt_2>24.0&&pt_1>13.0)||(trg_ehigh_mulow>0.5&&pt_1>24.0&&pt_2>10.0))");

  TString CutsSingle = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>33.))";
  TString CutsSingleAnd = "((trg_singlemuon>0.5&&pt_2>25.)&&(trg_singleelectron>0.5&&pt_1>33.))";
  if (era=="2017") {
    CutsSingle = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>28.))";
    CutsSingleAnd = "((trg_singlemuon>0.5&&pt_2>25.)&&(trg_singleelectron>0.5&&pt_1>28.))";
  }
  if (era=="2016") {
    CutsSingle = "((trg_singlemuon>0.5&&pt_2>23.)||(trg_singleelectron>0.5&&pt_1>26.))";
    CutsSingleAnd = "((trg_singlemuon>0.5&&pt_2>23.)&&(trg_singleelectron>0.5&&pt_1>26.))";
  }

  TString CutsTot = "(" + CutsEMu + "||" + CutsSingle + ")";
  TString CutsAnd = "((" + CutsEMu + ")&&(" + CutsSingleAnd + "))";

  TString Cuts = TString("&&iso_1<0.15&&iso_2<0.20&&extraelec_veto<0.5&&extramuon_veto<0.5&&dr_tt>0.3&&pt_1>15.0&&pt_2>15.&&puppipzeta>-40.0");

  CutsEMu += Cuts;
  CutsSingle += Cuts;
  CutsTot += Cuts;
  CutsAnd += Cuts;
  CutsSingleAnd += Cuts;
  std::cout << " Mass  single   e+mu    total   and_single  and" << std::endl;
             //"  140   0.796   0.735   1.147   0.731   0.678   1.042
  for (unsigned int i=0; i<nPoints; ++i) {

    TFile * file = new TFile(dir+"/"+Sample+masses.at(i)+".root");

    x[i] = Masses[masses.at(i)];

    sig_single[i] = getEff(file,CutsSingle,WeightSingle);
    sig_emu[i] = getEff(file,CutsEMu,WeightEMu);
    sig_tot[i] = getEff(file,CutsTot,WeightTot);
    sig_singleand[i] = getEff(file,CutsSingleAnd,WeightTot);
    sig_and[i] = getEff(file,CutsAnd,WeightTot);

    rat[i] = sig_single[i]/sig_emu[i];
    ratioTot[i] = sig_tot[i]/sig_emu[i];

    //    rat[i] = sig_singleand[i]/sig_and[i];
    //    ratioTot[i] = sig_[i]/sig_emu[i];

    printf("%4.0f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",
	   x[i],
	   sig_single[i],
	   sig_emu[i],
	   sig_tot[i],
	   sig_singleand[i],
	   sig_and[i]);
  }

  
  TGraph * gr1 = new TGraph(nPoints,x,sig_single);
  TGraph * gr2 = new TGraph(nPoints,x,sig_emu);
  TGraph * gr3 = new TGraph(nPoints,x,sig_tot);

  gr1->SetLineColor(1);
  gr1->SetLineWidth(3);

  gr2->SetLineColor(2);
  gr2->SetLineWidth(3);

  gr3->SetLineColor(3);
  gr3->SetLineWidth(3);

  TGraph * gr_ratio = new TGraph(nPoints,x,rat);
  gr_ratio->SetLineColor(1);
  gr_ratio->SetLineWidth(3);

  TGraph * gr_ratio_tot = new TGraph(nPoints,x,ratioTot);
  gr_ratio_tot->SetLineColor(3);
  gr_ratio_tot->SetLineWidth(3);

  TH2D * frame = new TH2D("frame","",2,140,3200,2,0,6);
  frame->GetYaxis()->SetTitle("#epsilon(bb#phi)");
  if (sig=="ggH")
    frame->GetYaxis()->SetTitle("#epsilon(gg#rightarrow #phi)");
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleSize(0.07);

  TH2D * frameRatio = new TH2D("frameRatio","",2,140,3200,2,0.5,1.5);
  frameRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  frameRatio->GetYaxis()->SetNdivisions(505);
  frameRatio->GetXaxis()->SetLabelFont(42);
  frameRatio->GetXaxis()->SetLabelOffset(0.04);
  frameRatio->GetXaxis()->SetLabelSize(0.14);
  frameRatio->GetXaxis()->SetTitleSize(0.13);
  frameRatio->GetXaxis()->SetTitleOffset(1.2);
  frameRatio->GetYaxis()->SetTitle("ratio");
  frameRatio->GetYaxis()->SetLabelFont(42);
  frameRatio->GetYaxis()->SetLabelOffset(0.015);
  frameRatio->GetYaxis()->SetLabelSize(0.13);
  frameRatio->GetYaxis()->SetTitleSize(0.14);
  frameRatio->GetYaxis()->SetTitleOffset(0.5);
  frameRatio->GetXaxis()->SetTickLength(0.07);
  frameRatio->GetYaxis()->SetTickLength(0.04);
  frameRatio->GetYaxis()->SetLabelOffset(0.01);
  frameRatio->GetXaxis()->SetTitle("m_{#phi} [GeV]");

  TCanvas * canv = new TCanvas("canv","",700,600);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
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

  frame->Draw();
  gr1->Draw("ls");
  gr2->Draw("ls");
  gr3->Draw("ls");
  //  gr4->Draw("ls");
  //  gr5->Draw("ls");

  //  upper->SetLogy(true);
  upper->SetLogx(true);
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  TLegend * leg = new TLegend(0.2,0.7,0.6,0.9);
  SetLegendStyle(leg);
  leg->SetTextSize(0.05);
  //  leg->AddEntry(gr1,"single-e && single-mu","l");
  //  leg->AddEntry(gr2,"e-mu && single-e && single-mu ","l");
  leg->AddEntry(gr1,"single-lep","l");
  leg->AddEntry(gr2,"e+mu","l");
  leg->AddEntry(gr3,"combination","l");
  leg->Draw();
  canv->cd();
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
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

  frameRatio->Draw();
  gr_ratio->Draw("ls");
  gr_ratio_tot->Draw("ls");
  //  gr_r4->Draw("ls");
  //  gr_r5->Draw("ls");

  lower->Modified();
  lower->RedrawAxis();
  lower->SetLogx(true);
  canv->cd();
  canv->Modified();
  canv->cd();
  canv->SetSelected(canv);
  canv->Update();
  canv->Print("eff_tot_"+sig+".png");
  



}

