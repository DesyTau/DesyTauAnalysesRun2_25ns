###################################################################################################################
#
# This is a datacard producer 
# Author: Valeria Botta
#
# Please read the instructions in instructions.txt before running it
#
# Arguments ( MUST be in THIS order )
# 1. channel (et, mt)
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
# -- histos mode
# > python datacard_producer_htt2016.py mt cuts_mt.json weights.json histos.json histos
#
# -- datacard mode, without syst variations
# > python datacard_producer_htt2016.py mt cuts_mt.json weights.json histos.json datacard 
#
# -- datacard mode, with syst variation (additional argument is the suffix of the syst variation)
# > python datacard_producer_htt2016.py mt cuts_mt.json weights.json histos.json datacard _CMS_scale_t_13TeVUp
# (run many times with all syst variations you want to include, then you can hadd the root files)
#
#################################################################################################################

import os, sys, json, ROOT
from ROOT import TFile, TH1F, TH1D, TH2D, TH2F, TCanvas, TPad, gStyle, kRed, kBlack, TLine, TLegend, TROOT, TChain
from array import array

def applyRelaxedSelection(categoryName):
	if categoryName in ["vbf", "boosted"]:
		return True
	else: 
		return False

def defineBlinding(blind,var,category):
	blinding = {
	"m_vis":"(m_vis<60 || m_vis>120)", 
	"m_sv":"(m_sv<100 || m_sv>150)",
	"mt_tot":"(mt_tot<80 || mt_tot>130)"
	}
	blinding_2Dfit = {
	"0jet":"(m_vis<65 || m_vis>110)", 
	"boosted":"(m_sv<100 || m_sv>150)",
	"vbf":"(m_sv<95 || m_sv > 155)"
	}	
	if (blind and (var in blinding) and ("inclusive" not in category)):
		if (datacardFor2Dfit):
			blindingCut = blinding_2Dfit[ckey]
		else:
			blindingCut = blinding[var]
		print "applying blinding for " , var , " : " , blindingCut
	else:
		blindingCut = " 1 "
	return blindingCut
	

def makeDatacard_ssos(treeName, files, cut, weight, histo, histo_ss):
	#print "opposite sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 < 0.", weight, "_os", histo)
	#print "same sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 > 0.", weight, "_ss", histo_ss)

def makeDatacard_sample_inc(treeName, sampleKey,sampleDict, catName, weight, cut, additionalCut):
	hsample = hdummy.Clone(sampleKey+catName)
	for fileName in sampleDict["files"]:
		#print "file name", fileName
		f=TFile(indir+"/"+fileName+".root", "read")
		t = f.Get(treeName)
		h = hdummy.Clone(fileName+catName)
		t.Draw(var+" >> "+fileName+catName, weight + " * ("+ cut + " && "+additionalCut+")")
		xs = sampleDict["xsec"]
		hevents=f.Get("nWeightedEvents")
		nevents=hevents.GetBinContent(1)
		try:
			h.Scale(xs*lumi/nevents)
		except ZeroDivisionError:
			pass
		hsample.Add(h)
	return hsample

def makeDatacard_sample_bin(treeName, sampleKey, sampleDict, catName, weight, cut, additionalCut):
	t = TChain(treeName)
	for fileName in sampleDict["files"]:
		t.Add(indir+"/"+fileName+".root")
	h = hdummy.Clone(sampleKey+catName)
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

def unroll(histo2D):			
	if histo2D is not None:
		htitle = histo2D.GetTitle()
		hname = histo2D.GetName()
		n_bins_x = (histo2D.GetNbinsX())*(histo2D.GetNbinsY())
		histo1D = TH1D(hname, htitle, n_bins_x, 0, n_bins_x)
		k =1 
		for i in xrange(1, histo2D.GetNbinsX()+1):
			for j in xrange(1, histo2D.GetNbinsY()+1):
				histo1D.SetBinContent(k, histo2D.GetBinContent(histo2D.GetBin(i,j)))
				histo1D.SetBinError(k, histo2D.GetBinError(histo2D.GetBin(i,j)))
				k+=1
		histo1D.GetXaxis().SetTitle("2D bin number")
		histo1D.GetYaxis().SetTitle("Events / bin width")
		return histo1D
	else:
		pass

#histo_list should contain tuples of (histogram, factor)
def add_all(result, histo_list):
	for h in histo_list:
		result.Add(h[0],h[1])


def Wjets_dataMC_sf(SSOSratio, hw): #hw contains the result of the calculation, should be the w histogram as it comes from MC, nominal selection

	print "Wjets data-driven normalisation from high mt region"
	print "using SS/OS QCD  = ", SSOSratio

	hdata_os_highmt = hdummy.Clone("hdata_os_highmt")
	hdata_ss_highmt = hdummy.Clone("hdata_ss_highmt")
	tdata.Draw(var+" >> hdata_os_highmt", cvalue+"&&"+highMTcut+" && q_1 * q_2 < 0. ")
	tdata.Draw(var + " >> hdata_ss_highmt", cvalue+"&&"+highMTcut+" && q_1 * q_2 > 0. " )

	hvv_os_highmt = hdummy.Clone("hvv_os_highmt")
	hvv_ss_highmt = hdummy.Clone("hvv_ss_highmt")
	makeDatacard_ssos(treeName, lvv, cvalue+"&&"+highMTcut, wvalue, hvv_os_highmt, hvv_ss_highmt)
	hewkz_os_highmt = hdummy.Clone("hewkz_os_highmt")
	hewkz_ss_highmt = hdummy.Clone("hewkz_ss_highmt")
	makeDatacard_ssos(treeName, lewkz, cvalue+"&&"+highMTcut, wvalue, hewkz_os_highmt, hewkz_ss_highmt)
	htt_os_highmt = hdummy.Clone("htt_os_highmt")
	htt_ss_highmt = hdummy.Clone("htt_ss_highmt")
	makeDatacard_ssos(treeName, ltt, cvalue+"&&"+highMTcut, wvalue, htt_os_highmt, htt_ss_highmt)
	hdy_os_highmt = hdummy.Clone("hdy_os_highmt")
	hdy_ss_highmt = hdummy.Clone("hdy_ss_highmt")
	makeDatacard_ssos(treeName, ldy, cvalue+"&&"+highMTcut, wvalue, hdy_os_highmt, hdy_ss_highmt)
	hw_os_highmt = hdummy.Clone("hw_os_highmt")
	hw_ss_highmt = hdummy.Clone("hw_ss_highmt")		
	makeDatacard_ssos(treeName, lw, cvalue+"&&"+highMTcut, wvalue, hw_os_highmt, hw_ss_highmt)

	hqcd_ss_highmt = hdummy.Clone("hqcd_ss_highmt")
	add_all(hqcd_ss_highmt, [(hdata_ss_highmt,1),(hvv_ss_highmt,-1),(htt_ss_highmt,-1),(hdy_ss_highmt,-1),(hw_ss_highmt,-1),(hewkz_ss_highmt,-1)])	
			
	#if relaxed selection, 
	# 1. define and fill histograms with ss + highmt + relaxed selection
	# 2. calculate qcd ss relaxed by subtracting data ss - all other MC
	# 3. normalise qcd ss relaxed to the yield of qcd ss w/o relaxed selection.
	
	#1

	if (relaxedSel is not None):
		print "relaxed sel in QCD in Wjets:" , relaxedSel+"&&"+highMTcut+"&& q_1*q_2>0"

		hdata_ss_highmt_rel = hdummy.Clone("hdata_ss_highmt_rel")
		tdata.Draw(var+">> hdata_ss_highmt_rel", relaxedSel+"&&"+highMTcut+"&& q_1*q_2>0")
		hvv_ss_highmt_rel = hdummy.Clone("hvv_ss_highmt_rel")
		makeDatacard_category(treeName, lvv, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", hvv_ss_highmt_rel )
		hewkz_ss_highmt_rel = hdummy.Clone("hewkz_ss_highmt_rel")
		makeDatacard_category(treeName, lewkz, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", hewkz_ss_highmt_rel )
		htt_ss_highmt_rel = hdummy.Clone("htt_ss_highmt_rel")
		makeDatacard_category(treeName, ltt, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", htt_ss_highmt_rel )
		hdy_ss_highmt_rel = hdummy.Clone("hdy_ss_highmt_rel")
		makeDatacard_category(treeName, ldy, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", hdy_ss_highmt_rel )
		hw_ss_highmt_rel = hdummy.Clone("hw_ss_highmt_rel")
		makeDatacard_category(treeName, lw, relaxedSel+"&&"+highMTcut, "q_1*q_2>0",  wvalue, "_highmt_ss_rel", hw_ss_highmt_rel )
	#2
		hqcd_ss_highmt_rel = hdummy.Clone("hqcd_ss_highmt_rel")
		add_all(hqcd_ss_highmt_rel, [(hdata_ss_highmt_rel,1),(hvv_ss_highmt_rel,-1),(htt_ss_highmt_rel,-1),(hdy_ss_highmt_rel,-1),(hw_ss_highmt_rel,-1),(hewkz_ss_highmt_rel,-1)])	
	#3
		hqcd_ss_highmt_rel.Scale(hqcd_ss_highmt.Integral()/(hqcd_ss_highmt_rel.Integral()))
		hqcd_ss_highmt = hqcd_ss_highmt_rel

	hw_os_highmt_datadriven = hdummy.Clone("hw_os_highmt_datadriven") 
	add_all(hw_os_highmt_datadriven, [(hdata_os_highmt,1),(hvv_os_highmt,-1),(htt_os_highmt,-1),(hdy_os_highmt,-1),(hewkz_os_highmt,-1),(hqcd_ss_highmt, -1*SSOSratio)] )

	# apply some relaxed selection. check which. For now taking the usual one. 
	relaxed_cuts_w_MC = {
	"mt_boosted":"pt_2>30 &&(njets==1 || (njets==2 && mjj<300) || njets>2) && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
	"mt_vbf":"pt_2>30 && njets==2 && mjj>300  && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
	"mt_0jet": "pt_2 > 30 && njets==0 && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
	"et_boosted":"pt_2>30 && (njets==1 || (njets==2 && mjj<300) || njets>2) && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
	"et_vbf":"pt_2>30 && njets==2 && mjj>300 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0" , 
	"et_0jet":"pt_2>30 && njets==0 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0"		
	}

	hw_os_MC = hdummy.Clone("hw_os_MC")
	hw_os_highmt_MC = hdummy.Clone("hw_os_highmt_MC")
	makeDatacard_category(treeName, lw, relaxed_cuts_w_MC[channel+"_"+ckey], "q_1*q_2<0"+lowMTcut,       wvalue, "hw_os_MC", hw_os_MC )
	makeDatacard_category(treeName, lw, relaxed_cuts_w_MC[channel+"_"+ckey], "q_1*q_2<0"+"&&"+highMTcut, wvalue, "hw_os_highmt_MC", hw_os_highmt_MC)

	# W scale factor 
	W_sf =  hw_os_highmt_datadriven.Integral()/hw_os_highmt_MC.Integral()
	hw.Scale(hw_os_MC.Integral()/hw.Integral()) 
	#adjust yield to norm_w = W_sf*hw_os_MC.Integral()
	print " --- W SF = ", W_sf
	return W_sf



def QCD_datadriven(SSOSratio, W_sf, hqcd): #hqcd contains the result of the QCD calculation
	hqcd_ss = hdummy.Clone("QCD_ss")
	add_all(hqcd_ss, [(hdata_ss,1),(hvv_ss,-1),(htt_ss,-1),(hdy_ss,-1),(hw_ss,-1*W_sf),(hewkz_ss,-1)])
	hqcd.Add(hqcd_ss, SSOSratio)

	if (relaxedSel is not None): 
		print "QCD relaxed selection"
		hqcd_rel_ss = hdummy.Clone("QCD_rel_ss")
		add_all(hqcd_rel_ss, [(hdata_rel_ss,1),(hvv_rel_ss,-1),(htt_rel_ss,-1),(hdy_rel_ss,-1),(hw_rel_ss,-1*W_sf),(hewkz_rel_ss,-1)])
		norm_nom = hqcd_ss.Integral()
		norm_rel = hqcd_rel_ss.Integral()
		hqcd_rel_ss.Scale(SSOSratio*norm_nom/norm_rel)
		hqcd = hqcd_rel_ss


##########################
# read arguments 
##########################

channel = sys.argv[1] #et, mt
cuts = json.loads(open(sys.argv[2]).read())
weighting = json.loads(open(sys.argv[3]).read())
histos = json.loads(open(sys.argv[4]).read())
mode = sys.argv[5]

systName = ""
if mode == "datacard":
	if len(sys.argv)>6:
		systName = sys.argv[6] #name of the tree with the systematic

#eg:
#"_CMS_shape_t_13TeVUp" Down --> _CMS_shape_t_mt_13TeVUp
#"_topPtWeightUp" Down --> CMS_shape_dyShape_13TeVUP (Down)
#"_CMS_shape_dyShape_13TeVUp" Down


##########################################
# adjust lumi and input directory  here  #
##########################################
lumi = 12900

#directory with input root files (synch ntuples) with conventional names!

if channel == "mt":
	#indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Nov/mutau/final"
	indir = "/afs/cern.ch/work/v/vbotta/NTuple_Nov/mutau/final"
	indirData = indir

if channel == "et":
	indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Nov/etau/final"
	indirData = indir

isBlinded = False 
datacardFor2Dfit = True

ROOT.TH1.SetDefaultSumw2(True)

##########################
# get data tree 
##########################

treeName = "TauCheck"

if channel == "mt":
	fdata =  TFile(indirData+"/"+"SingleMuon_Run2016BCD.root","read")
if channel == "et":
	fdata = TFile(indirData+"/"+"SingleElectron_Run2016BCD.root","read")

tdata = fdata.Get(treeName)


#redefine tree name for MC
treeName = "TauCheck"+systName


####################################################
# scale factors and up/down variations
# other selection cuts
# NB: more selection cuts in relaxed selections
####################################################

ssosvalues = {"et_0jet": 0.74, "et_boosted": 1.00, "et_vbf": 1.15, 
			  "mt_0jet": 1.02, "mt_boosted": 1.22, "mt_vbf": 1.13}

SSOS_variations = {"mt_0jet":{"up":1.10, "down":0.90}, "mt_boosted":{"up":1.20, "down":0.80}, "mt_vbf":{"up":1.20, "down":0.80}, 
				   "et_0jet":{"up":1.15, "down":0.85}, "et_boosted":{"up":1.35, "down":0.65}, "et_vbf":{"up":1.40, "down":0.60}}

Wsf_variations = {"0jet":{"up":1.05, "down":0.95}, "boosted":{"up":1.05, "down":0.95}, "vbf":{"up":1.25, "down":0.75}}

renorm_weights = {"et_0jet":"(0.9730+0.00034*pt_2)", "et_boosted":"(0.9857-0.000027787*pt_sv)", "et_vbf":"(0.9708+0.0003266*mjj)", 
				  "mt_0jet":"(0.9289+0.00017022*pt_2)", "mt_boosted":"(0.9195+0.0010055*pt_sv)", "mt_vbf":"(1.0258+0.00006596*mjj)"}

lowMTcut = "&& mt_1<50"
highMTcut = "mt_1>80"

##########################
# define file lists
##########################

lvv={"ST_t-channel_top_4f_leptonDecays" : {"xsec": 136.95*3*0.108, "isBinned": False, "binVar": "", "files": ["ST_t-channel_top_4f_leptonDecays"]},
     "ST_t-channel_antitop_4f_leptonDecays" : {"xsec": 80.95*3*0.108, "isBinned": False, "binVar": "", "files": ["ST_t-channel_antitop_4f_leptonDecays"]},
     "ST_tW_antitop_5f_inclusiveDecays" : {"xsec": 35.6, "isBinned": False, "binVar": "", "files": ["ST_tW_antitop_5f_inclusiveDecays"]},
     "ST_tW_top_5f_inclusiveDecays" : {"xsec": 35.6, "isBinned": False, "binVar": "", "files": ["ST_tW_top_5f_inclusiveDecays"]},
     "ZZTo2L2Q" : {"xsec": 3.22, "isBinned": False, "binVar": "", "files": ["ZZTo2L2Q"]},
     "ZZTo4L" : {"xsec": 1.212, "isBinned": False, "binVar": "", "files": ["ZZTo4L"]},
     "WZTo1L3Nu" : {"xsec": 3.05, "isBinned": False, "binVar": "", "files": ["WZTo1L3Nu"]},
     #"WZJTo3L1Nu" : {"xsec": 4.708, "isBinned": False, "binVar": "", "files": ["WZJTo3L1Nu"]},
     "WWTo1L1Nu2Q" : {"xsec": 49.997, "isBinned": False, "binVar": "", "files": ["WWTo1L1Nu2Q"]},
     "WZTo1L1Nu2Q" : {"xsec": 10.71, "isBinned": False, "binVar": "", "files": ["WZTo1L1Nu2Q"]},
     "VVTo2L2Nu" : {"xsec": 11.95, "isBinned": False, "binVar": "", "files": ["VVTo2L2Nu"]},
     "WZTo2L2Q" : {"xsec": 5.595, "isBinned": False, "binVar": "", "files": ["WZTo2L2Q"]}
}

ltt = {"TTPowHeg": {"xsec": 831.76, "files": ["TTpowheg"]}}

ldy = {"DYJetsToLL_M-50": { "isBinned": True, "binVar": "gen_noutgoing", 
                            "files": ["DYJetsToLL_M-50_MG",
                                      "DY1JetsToLL_M-50_MG",
                                      "DY2JetsToLL_M-50_MG",
                                      "DY3JetsToLL_M-50_MG",
                                      "DY4JetsToLL_M-50_MG"], 
                            "weights":{ (-0.5, 0.5) : 0.000118442,
                                        (0.5, 1.5) :  0.000015831,
                                        (1.5, 2.5) :  0.0000169775,
                                        (2.5, 3.5) :  0.0000169091,
                                        (3.5, 100.5) : 0.0000131252	
                              }
                    },
       "DYJetsToLL_M-10to50" : {"xsec": 18610, "isBinned": False, "binVar": "", "files": ["DYJetsToLL_M-10to50_MG"]}
}

lw = {"WJetsToLNu_MG": {"isBinned": True, "binVar": "gen_noutgoing", 
                              "files": ["WJetsToLNu_MG",
                                        "W1JetsToLNu_MG",
                                        "W2JetsToLNu_MG",
                                        "W3JetsToLNu_MG",
                                        "W4JetsToLNu_MG"], 
                              "weights":{ (-0.5, 0.5) : 0.00241933,
                                          (0.5, 1.5) : 0.000258246, 
                                          (1.5, 2.5) : 0.000120034,
                                          (2.5, 3.5) : 0.000056226,
                                          (3.5, 100.5) : 0.000067417
                                  }
                   },
	  "EWKWMinus2Jets_WToLNu_M-50":{"xsec": 20.25, "isBinned":False, "files":["EWKWMinus2Jets_WToLNu_M-50"]},
	  "EWKWPlus2Jets_WToLNu_M-50":{"xsec": 25.62, "isBinned":False, "files":["EWKWPlus2Jets_WToLNu_M-50"]}
}

lewkz = {"EWKZ2Jets_ZToLL_M-50": {"xsec": 3.987, "isBinned":False, "files":["EWKZ2Jets_ZToLL_M-50"]},
		 "EWKZ2Jets_ZToNuNu":    {"xsec": 10.01, "isBinned":False, "files":["EWKZ2Jets_ZToNuNu"]}
}

#Higgs to tau tau samples
lvbfhtt = {"VBFHTT_M125" : {"xsec": 0.22692, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M125"]}}

lwhtt = {"WmHTT_M125" : {"xsec": 0.0504, "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M125"]},
         "WpHTT_M125" : {"xsec": 0.031968, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M125"]} 
}

lzhtt = {"ZHTT_M125" : {"xsec": 0.053034, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M125"]}}

lgghtt = {"GluGluHTT_M125" : {"xsec": 2.9148, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M125"]}}


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

		#########################################
		# relaxed selection for some categories #
		#########################################

		relaxed_cuts = {
			"mt_boosted":"pt_2>30 &&(njets==1 || (njets==2 && mjj<300) || njets>2) && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
			"mt_vbf":"pt_2>30 && njets==2 && mjj>300  && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
			"et_boosted":"pt_2>30 && (njets==1 || (njets==2 && mjj<300) || njets>2) && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
			"et_vbf":"pt_2>30 && njets==2 && mjj>300 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0"
			}
		if (applyRelaxedSelection(ckey)):
			relaxedSel =  relaxed_cuts[channel+"_"+ckey]
		else: 
			relaxedSel = None
		print " applyRelaxedSel : " , applyRelaxedSelection(ckey)

		for hkey, hvalue in histos.iteritems():

			var = hvalue["var"]	
			nbins = hvalue["nbins"]
			xmin = hvalue["xmin"]
			xmax = hvalue["xmax"]
			xtit = hvalue["xtit"]
			ytit = hvalue["ytit"]

			#insert blinding to update cvalue in OS region
			blindingCut = defineBlinding(isBlinded,var,ckey)

			#here : update var to var1::var2 for 2D fit 
			# define binning for var2 according to category 

			#build string for axis titles
			axisTit = ";"+xtit+";"+ytit
			hdummy = 0
			if (datacardFor2Dfit==False):	
				print "Entered in datacardFor2Dfit==False "		
				hdummy = TH1F("hdummy","hdummy"+axisTit, nbins, xmin, xmax)
				print "plotting ", hkey, ": var ", var, " | bins ", nbins, " | xmin ", xmin, " | xmax ", xmax
			if (datacardFor2Dfit==True):
				print "Entered in datacardFor2Dfit==True "		
				bins_2Dfit = {
				"0jet": {"varX":"pt_2","binX":[30,35,40,45,50,55,300], "axistitX":"tau p_{T}", 
						 "varY":"m_vis","binY":[0,60,65,70,75,80,85,90,95,100,105,110,400], "axistitY":"m_{vis}"},
				"boosted": {"varX":"pt_sv","binX":[0,100,150,200,250,300,500], "axistitX":"Higgs p_{T}", 
							"varY":"m_sv","binY":[0,80,90,100,110,120,130,140,150,160,300], "axistitY":"m_{#tau#tau}"},
				"vbf": {"varX":"mjj","binX":[300,700,1100,1500,3000], "axistitX":"m_{jj}",
						"varY":"m_sv","binY":[0,95,115,135,155,400], "axistitY":"m_{#tau#tau}"}
				}

				varX = bins_2Dfit[ckey]["varX"]
				binsX = array('d',bins_2Dfit[ckey]["binX"])
				nbinsX = len(binsX)-1
				titX = bins_2Dfit[ckey]["axistitX"]

				varY = bins_2Dfit[ckey]["varY"]
				binsY = array('d',bins_2Dfit[ckey]["binY"])
				nbinsY = len(binsY)-1
				titY =  bins_2Dfit[ckey]["axistitY"]
				titZ = "Events / bin width"
				axisTit = ";"+titX+";"+titY+";"+titZ
				print "used for hdummy: nbinsX = ", nbinsX, " |  binsX = " , binsX, " | nbinsY = ", nbinsY, " | binsY ", binsY 
				hdummy = TH2F("hdummy","hdummy"+axisTit,nbinsX,binsX,nbinsY,binsY)
				var = varY+":"+varX

				blindingCut = defineBlinding(isBlinded,varY,ckey)

				print "plotting ", var	
			# data
			hdata = hdummy.Clone("hdata")
			hdata_ss = hdummy.Clone("hdata_ss")
			tdata.Draw(var+" >> hdata_ss", cvalue+ lowMTcut + " && q_1 * q_2 > 0.")
			if (relaxedSel is not None):
				hdata_rel = hdummy.Clone("hdata_rel")
				hdata_rel_ss = hdummy.Clone("hdata_rel_ss")
				tdata.Draw(var+" >> hdata_rel",  relaxedSel + lowMTcut +" && q_1 * q_2 < 0. ")
				tdata.Draw(var+" >> hdata_rel_ss",  relaxedSel + lowMTcut +" && q_1 * q_2 > 0. ")

			print "cut data: ", cvalue +lowMTcut + "&&" , blindingCut , "&& q_1 * q_2 < 0."
			print "Data OS: ", tdata.Draw(var+" >> hdata", cvalue +lowMTcut + "&&" + blindingCut + "&& q_1 * q_2 < 0.")
			#print "Data, done"

			### SIGNAL ###

			#ggH
			hggHtt = hdummy.Clone("ggH125")
			hggHtt.SetDirectory(0)
			makeDatacard_category(treeName, lgghtt, cvalue + lowMTcut , "q_1 * q_2 < 0.", wvalue, "_os", hggHtt)
			#VBF 
			hvbfHtt = hdummy.Clone("qqH125")
			hvbfHtt.SetDirectory(0)
			makeDatacard_category(treeName, lvbfhtt, cvalue + lowMTcut , "q_1 * q_2 < 0.", wvalue, "_os", hvbfHtt)
			#WH
			hWHtt = hdummy.Clone("hWHtt")
			hWHtt.SetDirectory(0)
			makeDatacard_category(treeName, lwhtt, cvalue + lowMTcut , "q_1 * q_2 < 0.", wvalue, "_os", hWHtt)
			#ZH
			hZHtt = hdummy.Clone("hZHtt")
			hZHtt.SetDirectory(0)
			makeDatacard_category(treeName, lzhtt, cvalue + lowMTcut , "q_1 * q_2 < 0.", wvalue, "_os", hZHtt)
			#no need for SS for signal

			#estimating renormalization uncertainty for ggH process. done when running the nominal datacard (systName = "")
			hggHtt_Up = None
			hggHtt_Down = None
			if (mode == "datacard" and systName == ""):
				print "calculating renormalisation uncertainty"
				renorm_weight_up = renorm_weights[channel+"_"+ckey]
				renorm_weight_down = "(1./ "+renorm_weight_up+")"
				hggHtt_Up = hdummy.Clone("hggHtt_Up")
				hggHtt_Down = hdummy.Clone("hggHtt_Down")
				hggHtt_Up.SetDirectory(0)
				hggHtt_Down.SetDirectory(0)
				makeDatacard_category(treeName, lgghtt, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_up, "_os", hggHtt_Up )
				makeDatacard_category(treeName, lgghtt, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_down, "_os", hggHtt_Down )




			### BACKGROUNDS ###

			#VV section
			hvv = hdummy.Clone("VV")
			hvv_ss = hdummy.Clone("VV_ss")
			hvv.SetDirectory(0)
			hvv_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lvv, cvalue+lowMTcut, wvalue, hvv, hvv_ss)
			if (relaxedSel is not None):
				hvv_rel = hdummy.Clone("VV_rel")
				hvv_rel_ss = hdummy.Clone("VV_rel_ss")
				makeDatacard_ssos(treeName, lvv, relaxedSel+lowMTcut, wvalue, hvv_rel, hvv_rel_ss)
			print "VV, done"

			#EWKZ
			hewkz = hdummy.Clone("EWKZ")
			hewkz_ss = hdummy.Clone("EWKZ_ss")
			hewkz.SetDirectory(0)
			hewkz_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lewkz, cvalue+lowMTcut, wvalue, hewkz, hewkz_ss)
			if (relaxedSel is not None):
				hewkz_rel = hdummy.Clone("EWKZ_rel")
				hewkz_rel_ss = hdummy.Clone("EWKZ_rel_ss")
				makeDatacard_ssos(treeName, lewkz, relaxedSel+lowMTcut, wvalue, hewkz_rel, hewkz_rel_ss)
			print "EWKZ, done"

			#TT section. Split in TTT and TTJ
			httt = hdummy.Clone("TTT")
			httj = hdummy.Clone("TTJ")
			htt = hdummy.Clone("TT")
			htt_ss = hdummy.Clone("TT_ss")
			httt.SetDirectory(0)
			httj.SetDirectory(0)
			htt_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ltt, cvalue+lowMTcut, wvalue, htt, htt_ss)
			makeDatacard_category(treeName, ltt, cvalue+lowMTcut, "gen_match_2 == 5 && q_1*q_2 < 0", wvalue, "TTT", httt)
			makeDatacard_category(treeName, ltt, cvalue+lowMTcut, "gen_match_2 != 5 && q_1*q_2 < 0", wvalue, "TTJ", httj)
			if (relaxedSel is not None):
				htt_rel = hdummy.Clone("TT_rel")
				htt_rel_ss = hdummy.Clone("TT_rel_ss")
				makeDatacard_ssos(treeName, ltt, relaxedSel+lowMTcut, wvalue, htt_rel, htt_rel_ss)
			print "TT, done"

			#DY
			hdy = hdummy.Clone("DY")
			hdy_ss = hdummy.Clone("DY_ss")
			hdy.SetDirectory(0)
			hdy_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ldy, cvalue+lowMTcut, wvalue, hdy, hdy_ss)
			if (relaxedSel is not None):
				hdy_rel = hdummy.Clone("DY_rel")
				hdy_rel_ss = hdummy.Clone("DY_rel_ss")				
				makeDatacard_ssos(treeName, ldy, relaxedSel+lowMTcut, wvalue, hdy_rel, hdy_rel_ss)

			#ZTT , ZL, ZJ
			hztt = hdummy.Clone("ztt")
			hzl = hdummy.Clone("zl")
			hzj = hdummy.Clone("zj")

			hztt.SetDirectory(0)
			hzl.SetDirectory(0)		
			hzj.SetDirectory(0)		

			ZTTcut = "gen_match_2 == 5 && q_1 * q_2 < 0." 
			ZLcut = "gen_match_2 < 5 && q_1 * q_2 < 0." 
			ZJcut = "gen_match_2 == 6 && q_1 * q_2 < 0." 

			makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut, wvalue, "ZTT", hztt)
			makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut, wvalue, "ZL", hzl)
			makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut, wvalue, "ZJ", hzj)

			print "DY all done"

			#Wjets
			hw = hdummy.Clone("W")
			hw_ss = hdummy.Clone("WW_ss")
			hw.SetDirectory(0)
			hw_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lw, cvalue+lowMTcut, wvalue, hw, hw_ss)
			if (relaxedSel is not None):
				hw_rel = hdummy.Clone("W_rel")
				hw_rel_ss = hdummy.Clone("W_rel_ss")
				makeDatacard_ssos(treeName, lw, relaxedSel+lowMTcut, wvalue, hw_rel, hw_rel_ss)

			#Wjets shape from relaxed selection
			if (relaxedSel is not None):
				print "Wjets relaxed selection "
				norm_nom = hw.Integral()
				norm_rel = hw_rel.Integral()
				hw_rel.Scale(norm_nom/norm_rel)
				hw = hw_rel 
			

			####### Wjets data-driven normalisation from high mt region

			SSOSratio = ssosvalues.get(channel+"_"+ckey, 1.06)
			hw_corrected = hw.Clone("hw_corrected")
			W_sf = Wjets_dataMC_sf(SSOSratio, hw_corrected)
			hw_corrected.Scale(W_sf)
			
			##########################
			# QCD estimate
			##########################

			print "qcd calculation..."

			hqcd = hdummy.Clone("QCD")
			QCD_datadriven(SSOSratio, W_sf, hqcd)

			print " ... done"

			####################################
			# shape uncertainties on QCD and W
			####################################

			if (mode=="datacard" and systName ==""):
				W_sf_up = Wsf_variations.get(ckey, 1).get("up",1)		
				W_sf_down = Wsf_variations.get(ckey, 1).get("down",1)			
				hw_wsf_up = hw.Clone("hw_wsf_up")
				hw_wsf_up.Scale(W_sf_up)
				hw_wsf_down = hw.Clone("hw_wsf_down")
				hw_wsf_down.Scale(W_sf_down)	
				hqcd_wsf_up = hdummy.Clone("hqcd_wsf_up")
				hqcd_wsf_down = hdummy.Clone("hqcd_wsf_down")
				QCD_datadriven(SSOSratio, W_sf*W_sf_up, hqcd_wsf_up)
				QCD_datadriven(SSOSratio, W_sf*W_sf_down, hqcd_wsf_down)

			if (mode=="datacard" and systName ==""):				
				SSOS_up = SSOS_variations.get(channel+"_"+ckey, 1).get("up",1)
				SSOS_down = SSOS_variations.get(channel+"_"+ckey, 1).get("down",1)
				hqcd_ssos_up = hdummy.Clone("hqcd_ssos_up")
				hqcd_ssos_down = hdummy.Clone("hqcd_ssos_down")
				QCD_datadriven(SSOSratio*SSOS_up, W_sf, hqcd_ssos_up)
				QCD_datadriven(SSOSratio*SSOS_down, W_sf, hqcd_ssos_down)
				hw_ssos_up = hw.Clone("hw_ssos_up")
				hw_ssos_down = hw.Clone("hw_ssos_down")
				W_sf = Wjets_dataMC_sf(SSOSratio*SSOS_up, hw_ssos_up)
				hw_ssos_up.Scale(W_sf)
				W_sf = Wjets_dataMC_sf(SSOSratio*SSOS_down, hw_ssos_down)
				hw_ssos_down.Scale(W_sf)

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

			#implement unrolling here! if 2d fit ...
			
			if (datacardFor2Dfit):
				hdata = unroll(hdata)	
				hggHtt = unroll(hggHtt)
				hvbfHtt = unroll(hvbfHtt)
				hWHtt = unroll(hWHtt)
				hZHtt = unroll(hZHtt)
				hztt = unroll(hztt)
				hzj = unroll(hzj)
				hzl = unroll(hzl)
				httt = unroll(httt)
				httj = unroll(httj)
				hw = unroll(hw)
				hw_corrected = unroll(hw_corrected)
				hewkz = unroll(hewkz)
				hvv = unroll(hvv)
				hqcd = unroll(hqcd)
				if (mode=="datacard" and systName == ""):
					hggHtt_Up = unroll(hggHtt_Up)
					hggHtt_Down = unroll(hggHtt_Down)
					hw_wsf_up = unroll(hw_wsf_up)
					hw_wsf_down = unroll(hw_wsf_down)
					hqcd_wsf_up = unroll(hqcd_wsf_up)
					hqcd_wsf_down = unroll(hqcd_wsf_down)
					hqcd_ssos_up = unroll(hqcd_ssos_up)
					hqcd_ssos_down = unroll(hqcd_ssos_down)
					hw_ssos_up = unroll(hw_ssos_up)
					hw_ssos_down = unroll(hw_ssos_down)


			
			if (mode == "histos"):
				hdata.Write("data_obs")

			#write data histogram only once, in the file w/o syst variations
			if (mode == "datacard" and systName == ""): 
				hdata.Write("data_obs")
				hggHtt_Up.Write("ggH125_CMS_scale_gg_13TeVUp")
				hggHtt_Down.Write("ggH125_CMS_scale_gg_13TeVDown")
				hw_wsf_up.Write("W_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hw_wsf_down.Write("W_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hqcd_wsf_up.Write("QCD_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hqcd_wsf_down.Write("QCD_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hqcd_ssos_up.Write("QCD_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hqcd_ssos_down.Write("QCD_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hw_ssos_up.Write("W_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hw_ssos_down.Write("W_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				
			print "writing to out file..."

			hggHtt.Write("ggH125"+systName)
			hvbfHtt.Write("qqH125"+systName)
			hWHtt.Write("WH125"+systName)
			hZHtt.Write("ZH125"+systName)

			hztt.Write("ZTT"+systName)
			hzl.Write("ZL"+systName)
			hzj.Write("ZJ"+systName)
			httt.Write("TTT"+systName)
			httj.Write("TTJ"+systName)
			hw.Write("W_noWSF"+systName)
			hw_corrected.Write("W"+systName)
			hewkz.Write("EWKZ"+systName)
			hvv.Write("VV"+systName)
			hqcd.Write("QCD"+systName)
			print "... done"

			if mode == "datacard":
				outfile.Close()

	if mode == "histos":
		outfile.Close()


