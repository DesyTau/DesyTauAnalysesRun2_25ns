###################################################################################################################
#
# This is a datacard producer 
# Author: Valeria Botta
#
#
# Arguments ( MUST be in THIS order )
# 1. channel (zmm)
# 2. json with selection cuts
# 3. json with weights
# 4. json with hitogram settings
# 5. mode ("histos" or "datacard")
# 6. name of the systematic (optional, only in datacard mode)
#
# ---- IMPORTANT -----
# - Please adjust the total lumi (1/pb) and the input directory
# - Stick to conventional names for root files, as in "define file lists" section
#
# ---- HOW TO RUN -----
#
# -- histos mode (for control plots)
# > python datacard_producer_zmumu.py zmm cuts_zmm.json weights.json histos.json histos
#
# -- datacard mode, without syst variations
# > python datacard_producer_zmumu.py zmm cuts_zmm.json weights.json histos.json datacard 
#
# -- datacard mode, with syst variation (additional argument is the suffix of the syst variation)
# > python datacard_producer_ztt2015.py mt cuts_mt.json weights.json histos.json datacard _CMS_scale_t_13TeVUp
# (run many times with all syst variations you want to include, then you can hadd the root files)
#
#################################################################################################################

import os, sys, json, ROOT
from ROOT import TFile, TH1F, TH1D, TCanvas, TPad, gStyle, kRed, kBlack, TLine, TLegend, TROOT, TChain

def makeDatacard_ssos(treeName, files, cut, weight, histo, histo_ss):
	print "opposite sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 < 0.", weight, "_os", histo)
	print "same sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 > 0.", weight, "_ss", histo_ss)

def makeDatacard_sample_inc(treeName, sampleKey,sampleDict, catName, weight, cut, additionalCut):
	hsample = TH1F(sampleKey+catName, sampleKey+catName, nbins, xmin, xmax)
	for fileName in sampleDict["files"]:
		print "file name", fileName
		f=TFile(indir+"/"+fileName+".root", "read")
		t = f.Get(treeName)
		h = TH1F(fileName+catName, fileName+catName, nbins, xmin, xmax)
		t.Draw(var+" >> "+fileName+catName, weight + " * ("+ cut + " && "+additionalCut+")")
		xs = sampleDict["xsec"]
		hevents=f.Get("inputEventsH")
		nevents=hevents.GetBinContent(1)
		h.Scale(xs*lumi/nevents)
		hsample.Add(h)
	return hsample

def makeDatacard_sample_bin(treeName, sampleKey, sampleDict, catName, weight, cut, additionalCut):
	t = TChain(treeName)
	for fileName in sampleDict["files"]:
		t.Add(indir+"/"+fileName+".root")
	h = TH1F(sampleKey+catName, sampleKey+catName, nbins, xmin, xmax)
	binString = makeBinString(sampleDict["weights"],sampleDict["binVar"])
	
	print "Drawing Entries : ", t.Draw(var+" >> "+sampleKey+catName, binString + " * "+ weight + " * ("+ cut + " && "+additionalCut+")")
	h.Scale(lumi)
	return h

def makeDatacard_category(treeName, catDict, cut, additionalCut, weight, catName, ohisto): #ohisto is output histogram for that category
	for sampleKey,sampleDict in catDict.iteritems():
		if (isBinned(sampleDict)):
			h = makeDatacard_sample_bin(treeName, sampleKey, sampleDict, catName, weight, cut, additionalCut)
		else: 
			h = makeDatacard_sample_inc(treeName, sampleKey, sampleDict, catName, weight, cut, additionalCut)
		ohisto.Add(h)


def isBinned(sampleDict):
		if 'weights' in sampleDict:
			return True 
		else:
			return False

def makeBinString(weightsDict, binVar):
	htwstr = "("
	for htwkey, htwvalue in weightsDict.iteritems():
		htwstr += "( "+ binVar +">="+ str(htwkey[0])
		htwstr += " && "+ binVar+ " < " + str(htwkey[1]) + ")"
		htwstr += " * " + str(htwvalue) + " + "
	htwstr = htwstr[0:-2]
	htwstr += ")"
	return htwstr


##########################
# read arguments 
##########################

channel = sys.argv[1] 
cuts = json.loads(open(sys.argv[2]).read())
weighting = json.loads(open(sys.argv[3]).read())
histos = json.loads(open(sys.argv[4]).read())
mode = sys.argv[5]

systName = ""
if mode == "datacard":
	if len(sys.argv)>6:
		systName = sys.argv[6] #name of the tree with the systematic

#eg:
#"_CMS_scale_t_13TeVUp"
#"_topPtWeightUp"


##########################################
# adjust lumi and input directory  here  #
##########################################

lumi = 35900 #/pb
indir = "./"

#directory with input root files (synch ntuples) with conventional names!

#if channel == "mt":
#	indir = "/nfs/dust/cms/user/fcost/HIGGS/newMVAMET/met0p6/CMSSW_7_6_3_patch2/src/DesyTauAnalyses/NTupleMaker/mutau/final"
#if channel == "et":
#	indir = "/nfs/dust/cms/user/fcost/HIGGS/newMVAMET/met0p6/CMSSW_7_6_3_patch2/src/DesyTauAnalyses/NTupleMaker/etau/final"	

ROOT.TH1.SetDefaultSumw2(True)

##########################
# get data tree 
##########################

treeName = "ZMuMu"
fdata =  TFile(indir+"/"+"SingleMuon_Run2016.root","read")
tdata = fdata.Get(treeName)

#redefine tree name for MC
treeName = "ZMuMu"+systName

##########################
# define file lists
##########################


lvv={"ST_t-channel_top_4f_leptonDecays" : {"xsec": 44.33, "isBinned": False, "binVar": "", "files": ["ST_t-channel_top_4f_InclusiveDecays_13TeV-powheg-pythia8"]},
     "ST_t-channel_antitop_4f_leptonDecays" : {"xsec": 26.38, "isBinned": False, "binVar": "", "files": ["ST_t-channel_antitop_4f_InclusiveDecays_13TeV-powheg-pythia8"]},
     "ST_tW_antitop_5f_inclusiveDecays" : {"xsec": 35.85, "isBinned": False, "binVar": "", "files": ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8"]},
     "ST_tW_top_5f_inclusiveDecays" : {"xsec": 35.85, "isBinned": False, "binVar": "", "files": ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8"]},
     "ZZ" : {"xsec": 12.19, "isBinned": False, "binVar": "", "files": ["ZZ_13TeV-pythia8"]},
     "WW" : {"xsec": 118.7, "isBinned": False, "binVar": "", "files": ["WW_13TeV-pythia8"]},
     "WZ" : {"xsec": 27.68, "isBinned": False, "binVar": "", "files": ["WZ_13TeV-pythia8"]}
}


ltt = {"TTPowHeg": {"xsec": 831.8, "files": ["TT_13TeV-powheg-pythia8"]}}

ldy = {"DYJetsToLL_M-50": {"xsec": 5345*1.079, "files": ["DYJetsToLL_M-50_13TeV_madgraphMLM_pythia8_RC"]} }

lw  = {"WJetsToLNu": {"xsec":52760*1.166, "files":["WJetsToLNu_13TeV-madgraphMLM-pythia8"]}}



# Example for a binned sample:
#
# ldy = {"DYJetsToLL_M-50": { "isBinned": True, "binVar": "gen_noutgoing", 
#                            "files": ["DYJetsToLL_M-50_MG_ext1", "DYJetsToLL_M-50_MG_ext2", 
#                                      "DY1JetsToLL_M-50_MG",
#                                      "DY2JetsToLL_M-50_MG",
#                                      "DY3JetsToLL_M-50_MG",
#                                      "DY4JetsToLL_M-50_MG"], 
#                            "weights":{ (-0.5, 0.5) : 0.0000395423,
#                                        (0.5, 1.5) :  0.0000128925,
#                                        (1.5, 2.5) :  0.0000130801,
#                                        (2.5, 3.5) :  0.0000130062,
#                                        (3.5, 100.5) : 0.0000105677	
#                              }
#                    },
#       "DYJetsToLL_M-10to50" : {"xsec": 18610, "isBinned": False, "binVar": "", "files": ["DYJetsToLL_M-10to50_MG"]}
#}
#

if mode == "datacard": 
	for hkey, hvalue in histos.iteritems():
		for wkey, wvalue in weighting.iteritems():
			filename = wkey+"_htt_"+channel+".inputs-sm-13TeV-"+hkey+systName+".root"
			#remove old files with same name as the one that will be created, if any
			os.system("if [ -f "+filename+" ]; then rm "+filename+";fi")

for wkey, wvalue in weighting.iteritems():
	print "weight : ", wkey
	if mode == "histos":
		suffix = os.path.splitext(os.path.basename(sys.argv[4]))[0] #add name of the (histos)json file as suffix
		outfile=TFile( wkey+"_htt_"+channel+".inputs-sm-13TeV-histos-"+suffix+".root", "recreate")	

	for ckey, cvalue in cuts.iteritems():
		print "cut : ", ckey
		
		for hkey, hvalue in histos.iteritems():

			var = hvalue["var"]	
			nbins = hvalue["nbins"]
			xmin = hvalue["xmin"]
			xmax = hvalue["xmax"]
			xtit = hvalue["xtit"]
			ytit = hvalue["ytit"]

			print "plotting ", hkey, ": var ", var, " | bins ", nbins, " | xmin ", xmin, " | xmax ", xmax
			#build string for axis titles
			axisTit = ";"+xtit+";"+ytit 
			# data
			hdata = TH1F("hdata", "hdata"+axisTit, nbins, xmin, xmax)
			hdata_ss = TH1F("hdata_ss", "hdata_ss"+axisTit,  nbins, xmin, xmax) 
			tdata.Draw(var+" >> hdata", cvalue + " && q_1 * q_2 < 0.")
			tdata.Draw(var+" >> hdata_ss", cvalue + " && q_1 * q_2 > 0.")
			print "Data, done"

			#VV section
			hvv = TH1F("VV", "VV"+axisTit, nbins, xmin, xmax)
			hvv_ss = TH1F("VV_ss", "VV_ss"+axisTit, nbins, xmin, xmax)
			hvv.SetDirectory(0)
			hvv_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lvv, cvalue, wvalue, hvv, hvv_ss)
			print "VV, done"
			#TT section
			htt = TH1F("TT", "TT"+axisTit, nbins, xmin, xmax)
			htt_ss = TH1F("TT_ss", "TT_ss"+axisTit, nbins, xmin, xmax)
			htt.SetDirectory(0)
			htt_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ltt, cvalue, wvalue, htt, htt_ss)
			print "TT, done"
			#Wjets
			hw =  TH1F("W", "W"+axisTit,nbins, xmin, xmax)
			hw_ss =  TH1F("W_ss","W_ss"+axisTit,nbins, xmin, xmax) 
			hw.SetDirectory(0)
			hw_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lw, cvalue, wvalue, hw, hw_ss)
			print "W, done"
			#DY
			hdy = TH1F("DY", "DY"+axisTit, nbins, xmin, xmax)
			hdy_ss = TH1F("DY_ss", "DY_ss"+axisTit, nbins, xmin, xmax) 
			hdy.SetDirectory(0)
			hdy_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ldy, cvalue, wvalue, hdy, hdy_ss)
			print "DY, done"
			##########################
			# QCD estimate
			##########################

			print "qcd calculation..."
			
			SSOSratio = 1.06
		
			#histograms to put the *total* MC yield 
			hMC = TH1F("hMC", "hMC"+axisTit, nbins, xmin, xmax)
			hMC_ss = TH1F("hMC_ss", "hMC_ss"+axisTit, nbins, xmin, xmax)
		
			hMC.Add(hvv)
			hMC.Add(htt)
			hMC.Add(hdy)
			hMC.Add(hw)

			hMC_ss.Add(hvv_ss)
			hMC_ss.Add(htt_ss)		
			hMC_ss.Add(hdy_ss)
			hMC_ss.Add(hw_ss)

			hqcd = TH1F("QCD", "QCD"+axisTit, nbins, xmin, xmax)
			hqcd_ss = TH1F("QCD_ss", "QCD_ss"+axisTit, nbins, xmin, xmax)

			hqcd_ss.Add(hdata_ss)
			hqcd_ss.Add(hMC_ss,-1)
			hqcd.Add(hqcd_ss, SSOSratio)

			print " ... done"

			##########################
			#writing to output file
			##########################

			
			if mode == "datacard":
				outfile=TFile( wkey+"_htt_"+channel+".inputs-sm-13TeV-"+hkey+systName+".root", "update")

				outfile.cd()
				outfile.mkdir(channel+"_"+ckey)
				outfile.cd(channel+"_"+ckey)

			elif mode == "histos":
				outfile.cd()
				outfile.mkdir(channel+"_"+ckey+"/"+hkey)
				outfile.cd(channel+"_"+ckey+"/"+hkey)

			print "writing to out file..."


			if (mode == "histos"):
				hdata.Write("data_obs")

			#write data histogram only once, in the file w/o syst variations
			if (mode == "datacards" and systName == ""): 
				hdata.Write("data_obs")

			hdy.Write("Z"+systName)
			htt.Write("TT"+systName)
			hw.Write("W"+systName)
			hvv.Write("VV"+systName)
			hqcd.Write("QCD"+systName)
			print "... done"

			if mode == "datacard":
				outfile.Close()

	if mode == "histos":
		outfile.Close()


