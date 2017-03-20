######################################################################################
# Script to plot the central, up and down distributions								 #
# for each MC process and syst uncertainty. 										 #
# The lower pad plots the ratios. 													 #
# It is meant to run on ROOT files produced by the datacardProducer.py script		 #
# 																					 #
# Please adjust																		 #
# - Input directory (indir)															 #	
# - Output directory where to save plots (plotOutDir)								 #
# - Dict with input ROOT files and variable name (varDict)							 #
# - The name of the categories and systematics (cat, systName)						 #
#																					 #	
# Input argument: channel (et or mt)												 #
# To run, e.g.:																		 #
# python ratioPlotSyst.py mt														 #		
# 																					 #
######################################################################################



import os, sys, ROOT
from ROOT import *

def setOutputName(dirName, catName, mcName, systName, varName):
	return dirName+"/"+catName+"_"+mcName+"_"+systName+"_"+varName

def makeRatioPlot(hcentr, hup, hdown, dirName, catName, mcName, systName, varName):

	c = TCanvas("c"+varName+"_"+catName+"_"+mcName+"_"+systName, "c"+varName+"_"+catName+"_"+mcName+"_"+systName, 600,700)
	topPad = TPad("top","top",0, 0.3, 1,1)
	topPad.SetBottomMargin(0)
	topPad.Draw()
	topPad.cd()

	hcentr.SetTitle(varName + " in " +catName+": "+mcName+" "+systName)
	hcentr.Draw()
	hup.Draw("same")
	hdown.Draw("same")
	hcentr.SetLineColor(kBlack)
	hcentr.SetMarkerStyle(20)
	hcentr.SetMarkerSize(0.7)
	hcentr.SetMarkerColor(kBlack)

	hup.SetLineColor(kBlue)
	hup.SetMarkerStyle(22)
	hup.SetMarkerSize(0.7)
	hup.SetMarkerColor(kBlue)

	hdown.SetLineColor(kRed)
	hdown.SetMarkerStyle(23)
	hdown.SetMarkerSize(0.7)
	hdown.SetMarkerColor(kRed)

	hup.SetStats(0)
	hdown.SetStats(0)
	hcentr.SetStats(0)

	leg = TLegend(0.6, 0.7, 0.9, 0.9)
	leg.AddEntry(hcentr, "central", "lp")
	leg.AddEntry(hup, "syst Up","lp")
	leg.AddEntry(hdown, "syst Down","lp") 	
	leg.Draw()

	c.cd()
	downPad = TPad("down","down",0,0,1,0.3)
	downPad.SetTopMargin(0)
	downPad.Draw()
	downPad.cd()

	#making clones before dividing
	hdown_c = hdown.Clone("")
	hup_c  = hup.Clone("")
	hcentr_c = hcentr.Clone("")
	
	hup_c.SetStats(0)
	hdown_c.SetStats(0)
	hcentr_c.SetStats(0)

	hup_c.Divide(hcentr_c)
	hup_c.SetAxisRange(0.,2.,"Y")
	hup_c.Draw("e")
	hdown_c.Divide(hcentr)
	hdown_c.Draw("eSAME")
	hcentr_c.Divide(hcentr)
	hcentr_c.Draw("eSAME")
	hup_c.SetTitle("")

	c.Update()

	oname = setOutputName(dirName, catName, mcName, systName, varName)
	c.SaveAs(oname+".png")

	#delete objects (maybe not needed)
	hdown_c.Delete()
	hup_c.Delete()
	hcentr_c.Delete()



#######################
# beginning of "main" #
#######################
ROOT.TH1.SetDefaultSumw2(kTRUE)

channel = sys.argv[1]

# directory with input ROOT file 
indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/DatacardProducer/"

#specify here (varDict) the name of the input ROOT file (datacard file with systematic variations) and name of the variable (eg: "mvis", "msv", ..)
varDict = {"bin": "cuts_mt_nominal_mt_nominal_htt_mt.inputs-sm-13TeV-mvis.root"}

#mcProc = ["ZTT","ZL","ZJ","TTJ","TTT","W","VV","QCD"]
mcProc = ["ZL"]

#directory for plots
plotOutDir = "."

#systematics to include and cathegories for each channel
#Must match the names in the input ROOT file
if channel == "et":
	cat = ["et_inclusive", "et_inclusivemt40"]
	systName = ["_CMS_scale_t_13TeV","_topPtWeight","_CMS_scale_eEB_13TeV", "_CMS_scale_eEE_13TeV"]
if channel == "mt":
	#cat = ["mt_boosted", "mt_0jet", "mt_vbf"]
	cat = ["mt_0jet"]
	#systName = ["_CMS_scale_t_1prong_13TeV","_CMS_scale_t_3prong_13TeV", "_CMS_scale_t_1prong1pizero_13TeV"]
	systName = ["_CMS_mFakeTau_1prong_13TeV", "_CMS_mFakeTau_1prong1pizero_13TeV", "_CMS_ZLShape_mt_0jet_1prong_13TeV", "_CMS_ZLShape_mt_0jet_1prong1pizero_13TeV"]
else:
	print "invalid channel choice"


for varName, varFile in varDict.iteritems(): 	
	inFile = TFile(indir+varFile,"read")
	
	#loop on file folders
	for aCat in cat: 

		print "Category ", aCat
		#define totMC histos
		h_totMC_centr = TH1F()
		h_totMC_up = TH1F()
		h_totMC_down = TH1F()
		
		#loop on MC processes
		for MCname in mcProc: 
			#take histos
			print "Getting ", aCat+"/"+MCname
			hCentral = inFile.Get(aCat+"/"+MCname)
			h_totMC_centr.Add(hCentral)
			for aSyst in systName:
				print "Getting ", aCat+"/"+MCname+aSyst+"Down"
				hDown = inFile.Get(aCat+"/"+MCname+aSyst+"Down")
				h_totMC_down.Add(hDown)
				print "Getting ", aCat+"/"+MCname+aSyst+"Up"
				hUp = inFile.Get(aCat+"/"+MCname+aSyst+"Up")
				h_totMC_up.Add(hUp)
				#make ratios for each MC process, each systematic
				makeRatioPlot(hCentral, hUp, hDown, plotOutDir, aCat, MCname, aSyst, varName)

		#make ratio plot of total MC, all syst up and all down
		makeRatioPlot(h_totMC_centr, h_totMC_up, h_totMC_down, plotOutDir, aCat, "allMC", "allSyst",varName)

	inFile.Close()


	

