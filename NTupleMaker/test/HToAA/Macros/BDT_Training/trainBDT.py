#!/usr/bin/env python
 
import ROOT
 
# in order to start TMVA
ROOT.TMVA.Tools.Instance()

############################################################################################
###################************* Switch of input variables ************#####################
############################################################################################
MuLMuT_Mass_isIN            = False
MuLMuT_DR_isIN              = True
MuLMuT_DPhi_isIN            = False
MuLTrk_Mass_isIN            = True
MuLTrk_Pt_isIN              = True
MuLTrk_DR_isIN              = True
MuLTrkMET_DPhi_isIN         = False
MuTTrk_Mass_isIN            = True
MuTTrk_Pt_isIN              = True
MuTTrk_DR_isIN              = True
MuTTrkMET_DPhi_isIN         = True
MuLTrkMuTTrk_Mass_isIN      = True
MuLTrkMuTTrk_Pt_isIN        = False
MuLTrkMuTTrkMET_Mass_isIN   = False
MET_Pt_isIN                 = True


############################################################################################
#######################*************** Cross Sections ************##########################
############################################################################################
mass_ma_str='herereplacemasspointstring'
mass_ma_float=float(mass_ma_str)

massTau = 1.777
massMu  = 0.106
massRatio = (massMu*massMu)/(massTau*massTau)
aF = 2*massTau/mass_ma_float
SF = 2*massRatio/ROOT.TMath.Sqrt(1-aF*aF)


xsecGGH = 48.52
xsecVBF = 3.779
xsecVH  = (1.369 + 0.8824)
xsecTTH = 0.5065
xsecMMTT = (xsecGGH+xsecVBF+xsecVH+xsecTTH) * SF

############################################################################################
#######################*************** Files & Trees *************##########################
############################################################################################
fileS_ggh = ROOT.TFile("SUSYGluGluToHToAA_AToTauTau_SCDoMuFilter_M-"+ mass_ma_str +".root")
fileS_vbf = ROOT.TFile("SUSYVBFToHToAA_AToTauTau_M-"+ mass_ma_str +".root")
fileS_vh = ROOT.TFile("SUSYVH_HToAA_AToTauTau_M-"+ mass_ma_str +".root")
fileS_tth = ROOT.TFile("SUSYttH_HToAA_AToTauTau_M-"+ mass_ma_str +".root")
fileS_mmtt = ROOT.TFile("SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-"+ mass_ma_str +".root")

fileB = ROOT.TFile("DoubleMuon_Run2016.root")

tree_s_ggh = fileS_ggh.Get("tree_Sel")
tree_s_vbf = fileS_vbf.Get("tree_Sel")
tree_s_vh = fileS_vh.Get("tree_Sel")
tree_s_tth = fileS_tth.Get("tree_Sel")
tree_s_mmtt = fileS_mmtt.Get("tree_Sel")

tree_b = fileB.Get("tree_SemiIso")

############################################################################################
###################*************** Weights for Training ************########################
############################################################################################
histWeightsGGH = fileS_ggh.Get('histWeightsH')
nGenGGH = histWeightsGGH.GetSumOfWeights()

histWeightsVBF = fileS_vbf.Get('histWeightsH')
nGenVBF = histWeightsVBF.GetSumOfWeights()

histWeightsVH = fileS_vh.Get('histWeightsH')
nGenVH = histWeightsVH.GetSumOfWeights()

histWeightsttH = fileS_tth.Get('histWeightsH')
nGenttH = histWeightsttH.GetSumOfWeights()

histWeightsMMTT = fileS_mmtt.Get('histWeightsH')
nGenMMTT = histWeightsMMTT.GetSumOfWeights()

SCDoMuFilter_Eff = 0.0034

Weight_GGH = 1. * SCDoMuFilter_Eff
Weight_VBF = (xsecVBF*tree_s_vbf.GetEntries()/nGenVBF)/(xsecGGH*tree_s_ggh.GetEntries()/nGenGGH)
Weight_VH = (xsecVH*tree_s_vh.GetEntries()/nGenVH)/(xsecGGH*tree_s_ggh.GetEntries()/nGenGGH)
Weight_TTH = (xsecTTH*tree_s_tth.GetEntries()/nGenttH)/(xsecGGH*tree_s_ggh.GetEntries()/nGenGGH)
Weight_MMTT = (xsecMMTT*tree_s_mmtt.GetEntries()/nGenMMTT)/(xsecGGH*tree_s_ggh.GetEntries()/nGenGGH)

Weight_b = 1.


############################################################################################
###################**************** BDT Configuration **************########################
############################################################################################

###################**************** Output File **************########################
fout = ROOT.TFile("test.root","RECREATE")

###################**************** General Configuration **************########################
factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([    "!V",
                                          "!Silent",
                                          "Color",
                                          "DrawProgressBar",
                                          "Transformations=I;D;P;G,D",
                                          "AnalysisType=Classification"]
                                     ))

###################**************** Input Variables **************########################
if MuLMuT_Mass_isIN:
                             factory.AddVariable("MuLMuT_Mass","F")
if MuLMuT_DR_isIN:                    
                             factory.AddVariable("MuLMuT_DR","F")
if MuLMuT_DPhi_isIN:
                             factory.AddVariable("MuLMuT_DPhi","F")
if MuLTrk_Mass_isIN:
                             factory.AddVariable("MuLTrk_Mass","F")
if MuLTrk_Pt_isIN:
                             factory.AddVariable("MuLTrk_Pt","F")
if MuLTrk_DR_isIN:
                             factory.AddVariable("MuLTrk_DR","F")
if MuLTrkMET_DPhi_isIN:
                             factory.AddVariable("MuLTrkMET_DPhi","F")
if MuTTrk_Mass_isIN:
                             factory.AddVariable("MuTTrk_Mass","F")
if MuTTrk_Pt_isIN:
                             factory.AddVariable("MuTTrk_Pt","F")
if MuTTrk_DR_isIN:
                             factory.AddVariable("MuTTrk_DR","F")
if MuTTrkMET_DPhi_isIN:
                             factory.AddVariable("MuTTrkMET_DPhi","F")
if MuLTrkMuTTrk_Mass_isIN:
                             factory.AddVariable("MuLTrkMuTTrk_Mass","F")
if MuLTrkMuTTrk_Pt_isIN:
                             factory.AddVariable("MuLTrkMuTTrk_Pt","F")
if MuLTrkMuTTrkMET_Mass_isIN:
                             factory.AddVariable("MuLTrkMuTTrkMET_Mass","F")
if MET_Pt_isIN:
                             factory.AddVariable("MET_Pt","F")

###################**************** Addding Trees with weights **************########################
factory.AddSignalTree(tree_s_ggh,Weight_GGH)
factory.AddSignalTree(tree_s_vbf,Weight_VBF)
factory.AddSignalTree(tree_s_vh,Weight_VH)
factory.AddSignalTree(tree_s_tth,Weight_TTH)
factory.AddSignalTree(tree_s_mmtt,Weight_MMTT)

factory.AddBackgroundTree(tree_b,Weight_b)


###################**************** Variables Cuts **************########################
sigCut = ROOT.TCut("MuLTrk_Mass<22. && MuTTrk_Mass<22. && MuLTrkMuTTrk_Mass<125.")
bgCut = ROOT.TCut("MuLTrk_Mass<22. && MuTTrk_Mass<22. && MuLTrkMuTTrk_Mass<125.")
 
###################**************** General Options for Training **************########################
factory.PrepareTrainingAndTestTree(sigCut, 
                                   bgCut, 
                                   ":".join(["nTrain_Signal=0",
                                             "nTrain_Background=0",
                                             "nTest_Signal=0",
                                             "nTest_Background=0",
                                             "SplitMode=Random",
                                             "NormMode=NumEvents",
                                             "!V"
                                             ]))
 
###################**************** BDT Options **************######################## 
method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
                            ":".join([ "!H",
                                       "!V",
                                       "NTrees=300",
                                       "MaxDepth=3",
                                       "BoostType=AdaBoost",
                                       "AdaBoostBeta=0.5",
                                       "SeparationType=GiniIndex",
                                       "nCuts=20",
                                       "MinNodeSize=5",
                                       ]))
 
###################**************** Self-Explanatory **************########################
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
