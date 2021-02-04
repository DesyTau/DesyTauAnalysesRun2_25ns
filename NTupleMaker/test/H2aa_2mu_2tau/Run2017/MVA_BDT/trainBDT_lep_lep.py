#!/usr/bin/env python
 
import ROOT
 
# in order to start TMVA
ROOT.TMVA.Tools.Instance()

############################################################################################
###################************* Switch of input variables ************#####################
############################################################################################
MuMu_Mass_isIN            = True
MuMu_Pt_isIN              = False
MuMu_DR_isIN              = True
MuMuMET_DPhi_isIN         = False
TrkTrk_Mass_isIN          = False
TrkTrk_Pt_isIN            = False
TrkTrk_DR_isIN            = True
TrkTrkMET_DPhi_isIN       = False
MuMuTrkTrk_Mass_isIN      = False
MuMuTrkTrk_Pt_isIN        = False
MuMuTauTau_Mass_isIN      = False
MuMuTauTau_Pt_isIN        = False
MuMuTrkTrkMET_Mass_isIN   = True
MET_Pt_isIN               = False
 
############################################################################################
#######################*************** Mass Points ***************##########################
############################################################################################
mass_str='herereplacemasspointstring'


############################################################################################
#######################*************** Files & Trees *************##########################
############################################################################################
fileS = ROOT.TFile('../../Workspace_Interpolation.root')
fileB = ROOT.TFile('../MergeTrees/tree_b_file.root')

ws = fileS.Get('ws_lep_lep_'+mass_str)
SignalModel = ws.pdf('model') 

Dataset_Old = SignalModel.generate(ws.set('observables'),300000)
ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)
Dataset = ROOT.RooDataSet('Dataset', 'Dataset', Dataset_Old, Dataset_Old.get())


tree_s = Dataset.tree().Clone()

tree_b_lep_lep = fileB.Get('tree_b_lep_lep')
tree_b_lep_had = fileB.Get('tree_b_lep_had')
tree_b_had_had = fileB.Get('tree_b_had_had')

############################################################################################
###################*************** Weights for Training ************########################
############################################################################################

Weight_s = 1.

Weight_b_lep_lep = 5.
Weight_b_lep_had = 2.
Weight_b_had_had = 1.


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
if MuMu_Mass_isIN:
                             factory.AddVariable("MuMu_Mass","F")
if MuMu_Pt_isIN:                    
                             factory.AddVariable("MuMu_Pt","F")
if MuMu_DR_isIN:
                             factory.AddVariable("MuMu_DR","F")
if MuMuMET_DPhi_isIN:
                             factory.AddVariable("MuMuMET_DPhi","F")
if TrkTrk_Mass_isIN:
                             factory.AddVariable("TrkTrk_Mass","F")
if TrkTrk_Pt_isIN:
                             factory.AddVariable("TrkTrk_Pt","F")
if TrkTrk_DR_isIN:
                             factory.AddVariable("TrkTrk_DR","F")
if TrkTrkMET_DPhi_isIN:
                             factory.AddVariable("TrkTrkMET_DPhi","F")
if MuMuTrkTrk_Mass_isIN:
                             factory.AddVariable("MuMuTrkTrk_Mass","F")
if MuMuTrkTrk_Pt_isIN:
                             factory.AddVariable("MuMuTrkTrk_Pt","F")
if MuMuTauTau_Mass_isIN:
                             factory.AddVariable("MuMuTauTau_Mass","F")
if MuMuTauTau_Pt_isIN:
                             factory.AddVariable("MuMuTauTau_Pt","F")
if MuMuTrkTrkMET_Mass_isIN:
                             factory.AddVariable("MuMuTrkTrkMET_Mass","F")
if MET_Pt_isIN:
                             factory.AddVariable("MET_Pt","F")


###################**************** Addding Trees with weights **************########################
factory.AddSignalTree(tree_s,Weight_s)

factory.AddBackgroundTree(tree_b_lep_lep,Weight_b_lep_lep)
factory.AddBackgroundTree(tree_b_lep_had,Weight_b_lep_had)
factory.AddBackgroundTree(tree_b_had_had,Weight_b_had_had)

###################**************** Variables Cuts **************########################
sigCut = ROOT.TCut("1")
bgCut = ROOT.TCut("1")
 
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
method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT_lep_lep",
                            ":".join([ "!H",
                                       "!V",
                                       "NTrees=500",
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
