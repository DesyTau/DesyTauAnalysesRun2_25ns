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
	if (blind and (var in blinding) and ("inclusive" not in category)):
		blindingCut = blinding[var]
		print "applying blinding for " , var , " : " , blindingCut
	else:
		blindingCut = " 1 "
	return blindingCut
	

def makeDatacard_ssos(treeName, files, cut, weight, histo, histo_ss):
	print "opposite sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 < 0.", weight, "_os", histo)
	print "same sign"
	makeDatacard_category(treeName, files, cut, "q_1 * q_2 > 0.", weight, "_ss", histo_ss)

def makeDatacard_sample_inc(treeName, sampleKey,sampleDict, catName, weight, cut, additionalCut):
	hsample = hdummy.Clone(sampleKey+catName)
	for fileName in sampleDict["files"]:
		print "file name", fileName
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
'''
def unroll_all(histo_list, out_list):
	#out_list = []
	for h in histo_list:
		out_list.append(unroll(h))
'''
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
#"_CMS_scale_t_13TeVUp"
#"_topPtWeightUp"


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

isBlinded = True 
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
		#applyRelaxedSel = False
		#print "ckey : ", ckey
		#if ckey in ["1jet_high","vbf_low","vbf_high","vbf","1jet"]:
		#	applyRelaxedSel = True
		#	relaxedSel = "njets==2 && mjj>500 && (mjj<800 || pt_sv<100) && mt_1<50 && pt_1>23 && pt_2>20 && iso_1 < 0.15 && byTightIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 "
		#print " applyRelaxedSel : " , applyRelaxedSel

		if (applyRelaxedSelection(ckey)):
			relaxed_cuts = {
			"boosted":"pt_2>30 && mt_1<50 &&(njets==1 || (njets==2 && mjj<300) || njets>2) && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
			"vbf":"pt_2>30 && mt_1<50 && njets==2 && mjj>300  && pt_1>23 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 "
			}
			relaxedSel =  relaxed_cuts[ckey]
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

			print "plotting ", hkey, ": var ", var, " | bins ", nbins, " | xmin ", xmin, " | xmax ", xmax

			# data
			hdata = hdummy.Clone("hdata")
			hdata_ss = hdummy.Clone("hdata_ss")
			tdata.Draw(var+" >> hdata_ss", cvalue + " && q_1 * q_2 > 0.")
			if (relaxedSel is not None):
				hdata_rel = hdummy.Clone("hdata_rel")
				hdata_rel_ss = hdummy.Clone("hdata_rel_ss")
				tdata.Draw(var+" >> hdata_rel",  relaxedSel + " && q_1 * q_2 < 0. ")
				tdata.Draw(var+" >> hdata_rel_ss",  relaxedSel + " && q_1 * q_2 > 0. ")

			print "cut data: ", cvalue, "&&" , blindingCut , "&& q_1 * q_2 < 0."
			print "Data OS: ", tdata.Draw(var+" >> hdata", cvalue + "&&" + blindingCut + "&& q_1 * q_2 < 0.")
			print "Data, done"

			### SIGNAL ###

			#ggH
			hggHtt = hdummy.Clone("ggH125")
			hggHtt.SetDirectory(0)
			makeDatacard_category(treeName, lgghtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hggHtt)
			#VBF 
			hvbfHtt = hdummy.Clone("qqH125")
			hvbfHtt.SetDirectory(0)
			makeDatacard_category(treeName, lvbfhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hvbfHtt)
			#WH
			hWHtt = hdummy.Clone("hWHtt")
			hWHtt.SetDirectory(0)
			makeDatacard_category(treeName, lwhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hWHtt)
			#ZH
			hZHtt = hdummy.Clone("hZHtt")
			hZHtt.SetDirectory(0)
			makeDatacard_category(treeName, lzhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hZHtt)
			#no need for SS for signal

			#estimating renormalization uncertainty for ggH process. done when running the nominal datacard (systName = "")
			hggHtt_Up = None
			hggHtt_Down = None
			if (mode == "datacard" and systName == ""):
				print "calculating renormalisation uncertainty"
				renorm_weights = {"et_0jet":"(0.9730+0.00034*pt_2)", "et_boosted":"(0.9857-0.000027787*pt_sv)", "et_vbf":"(0.9708+0.0003266*mjj)", 
								  "mt_0jet":"(0.9289+0.00017022*pt_2)", "mt_boosted":"(0.9195+0.0010055*pt_sv)", "mt_vbf":"(1.0258+0.00006596*mjj)"}
				renorm_weight_up = renorm_weights[channel+"_"+ckey]
				renorm_weight_down = "(1./ "+renorm_weight_up+")"
				hggHtt_Up = hdummy.Clone("hggHtt_Up")
				hggHtt_Down = hdummy.Clone("hggHtt_Down")
				hggHtt_Up.SetDirectory(0)
				hggHtt_Down.SetDirectory(0)
				makeDatacard_category(treeName, lgghtt, cvalue,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_up, "_os", hggHtt_Up )
				makeDatacard_category(treeName, lgghtt, cvalue,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_down, "_os", hggHtt_Down )




			### BACKGROUNDS ###

			#VV section
			hvv = hdummy.Clone("VV")
			hvv_ss = hdummy.Clone("VV_ss")
			hvv.SetDirectory(0)
			hvv_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lvv, cvalue, wvalue, hvv, hvv_ss)
			if (relaxedSel is not None):
				hvv_rel = hdummy.Clone("VV_rel")
				hvv_rel_ss = hdummy.Clone("VV_rel_ss")
				makeDatacard_ssos(treeName, lvv, relaxedSel, wvalue, hvv_rel, hvv_rel_ss)
			print "VV, done"

			#EWKZ
			hewkz = hdummy.Clone("EWKZ")
			hewkz_ss = hdummy.Clone("EWKZ_ss")
			hewkz.SetDirectory(0)
			hewkz_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lewkz, cvalue, wvalue, hewkz, hewkz_ss)
			if (relaxedSel is not None):
				hewkz_rel = hdummy.Clone("EWKZ_rel")
				hewkz_rel_ss = hdummy.Clone("EWKZ_rel_ss")
				makeDatacard_ssos(treeName, lewkz, relaxedSel, wvalue, hewkz_rel, hewkz_rel_ss)
			print "EWKZ, done"

			#TT section. Split in TTT and TTJ
			httt = hdummy.Clone("TTT")
			httj = hdummy.Clone("TTJ")
			htt = hdummy.Clone("TT")
			htt_ss = hdummy.Clone("TT_ss")
			httt.SetDirectory(0)
			httj.SetDirectory(0)
			htt_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ltt, cvalue, wvalue, htt, htt_ss)
			makeDatacard_category(treeName, ltt, cvalue, "gen_match_2 == 5 && q_1*q_2 < 0", wvalue, "TTT", httt)
			makeDatacard_category(treeName, ltt, cvalue, "gen_match_2 == 5 && q_1*q_2 < 0", wvalue, "TTJ", httj)
			if (relaxedSel is not None):
				htt_rel = hdummy.Clone("TT_rel")
				htt_rel_ss = hdummy.Clone("TT_rel_ss")
				makeDatacard_ssos(treeName, ltt, relaxedSel, wvalue, htt_rel, htt_rel_ss)
			print "TT, done"

			#Wjets
			hw = hdummy.Clone("W")
			hw_ss = hdummy.Clone("WW_ss")
			hw.SetDirectory(0)
			hw_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lw, cvalue, wvalue, hw, hw_ss)
			if (relaxedSel is not None):
				hw_rel = hdummy.Clone("W_rel")
				hw_rel_ss = hdummy.Clone("W_rel_ss")
				makeDatacard_ssos(treeName, lw, relaxedSel, wvalue, hw_rel, hw_rel_ss)

			#DY
			hdy = hdummy.Clone("DY")
			hdy_ss = hdummy.Clone("DY_ss")
			hdy.SetDirectory(0)
			hdy_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, ldy, cvalue, wvalue, hdy, hdy_ss)
			if (relaxedSel is not None):
				hdy_rel = hdummy.Clone("DY_rel")
				hdy_rel_ss = hdummy.Clone("DY_rel_ss")				
				makeDatacard_ssos(treeName, ldy, relaxedSel, wvalue, hdy_rel, hdy_rel_ss)

			#ZTT , ZL, ZJ, ZLL 
			hztt = hdummy.Clone("ztt")
			hzl = hdummy.Clone("zl")
			hzj = hdummy.Clone("zj")
			#hzll = hdummy.Clone("zll")

			hztt.SetDirectory(0)
			hzl.SetDirectory(0)		
			hzj.SetDirectory(0)		
			#hzll.SetDirectory(0)

			#looping on different cuts that define categories 
			#for now, done by hand here with additional cut string
			ZTTcut = "gen_match_2 == 5 && q_1 * q_2 < 0." 
			ZLcut = "gen_match_2 < 5 && q_1 * q_2 < 0." 
			ZJcut = "gen_match_2 == 6 && q_1 * q_2 < 0." 
			#ZLLcut = "gen_match_2 != 5 && q_1 * q_2 < 0."

			makeDatacard_category(treeName, ldy, cvalue, ZTTcut, wvalue, "ZTT", hztt)
			makeDatacard_category(treeName, ldy, cvalue, ZLcut, wvalue, "ZL", hzl)
			makeDatacard_category(treeName, ldy, cvalue, ZJcut, wvalue, "ZJ", hzj)
			#makeDatacard_category(treeName, ldy, cvalue, ZLLcut, wvalue, "ZLL", hzll)

			print "DY all done"
			
			##########################
			# QCD estimate
			##########################

			print "qcd calculation..."

			ssosvalues = {"et_0jet": 0.74, "et_boosted": 1.00 , "et_vbf": 1.15, 
						  "mt_0jet": 1.02, "mt_boosted": 1.22, "mt_vbf": 1.13}
			SSOSratio = ssosvalues.get(channel+"_"+ckey, 1.06)
			print "using SS/OS QCD  = ", SSOSratio
		
			#histograms to put the *total* MC yield 
			#hMC = TH1F("hMC", "hMC"+axisTit, nbins, xmin, xmax)
			#hMC_ss = TH1F("hMC_ss", "hMC_ss"+axisTit, nbins, xmin, xmax)
			hMC = hdummy.Clone("hMC")
			hMC_ss = hdummy.Clone("hMC_ss")	
	
			hMC.Add(hvv)
			hMC.Add(htt)
			hMC.Add(hdy)
			hMC.Add(hw)
			hMC.Add(hewkz)

			hMC_ss.Add(hvv_ss)
			hMC_ss.Add(htt_ss)		
			hMC_ss.Add(hdy_ss)
			hMC_ss.Add(hw_ss)
			hMC_ss.Add(hewkz_ss)

			#hqcd = TH1F("QCD", "QCD"+axisTit, nbins, xmin, xmax)
			#hqcd_ss = TH1F("QCD_ss", "QCD_ss"+axisTit, nbins, xmin, xmax)
			hqcd = hdummy.Clone("QCD")
			hqcd_ss = hdummy.Clone("QCD_ss")
		
			hqcd_ss.Add(hdata_ss)
			hqcd_ss.Add(hMC_ss,-1)
			hqcd.Add(hqcd_ss, SSOSratio)

			print " ... done"

			if (relaxedSel is not None): 
				print "QCD relaxed selection"
				hqcd_rel_ss = hdummy.Clone("QCD_rel_ss")
				hqcd_rel_ss.Add(hdata_rel_ss)
				hqcd_rel_ss.Add(hvv_rel_ss,-1)
				hqcd_rel_ss.Add(htt_rel_ss,-1)
				hqcd_rel_ss.Add(hdy_rel_ss,-1)
				hqcd_rel_ss.Add(hw_rel_ss,-1)
				hqcd_rel_ss.Add(hewkz_rel_ss,-1)
				norm_nom = hqcd_ss.Integral()
				norm_rel = hqcd_rel_ss.Integral()
				hqcd_rel_ss.Scale(SSOSratio*norm_nom/norm_rel)
				hqcd = hqcd_rel_ss

			#Wjets from control region 
			#if (relaxedSel!=None):
			if (relaxedSel is not None):
				print "Wjets relaxed selection "
				norm_nom = hw.Integral()
				norm_rel = hw_rel.Integral()
				hw_rel.Scale(norm_nom/norm_rel)
				hw = hw_rel 




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
				hggHtt_Up = unroll(hggHtt_Up)
				hggHtt_Down = unroll(hggHtt_Down)
				hvbfHtt = unroll(hvbfHtt)
				hWHtt = unroll(hWHtt)
				hZHtt = unroll(hZHtt)
				hztt = unroll(hztt)
				hzj = unroll(hzj)
				hzl = unroll(hzl)
				httt = unroll(httt)
				httj = unroll(httj)
				hw = unroll(hw)
				hewkz = unroll(hewkz)
				hvv = unroll(hvv)
				hqcd = unroll(hqcd)
			
			if (mode == "histos"):
				hdata.Write("data_obs")

			#write data histogram only once, in the file w/o syst variations
			if (mode == "datacard" and systName == ""): 
				hdata.Write("data_obs")
				hggHtt_Up.Write("ggH125_CMS_scale_gg_13TeVUp")
				hggHtt_Down.Write("ggH125_CMS_scale_gg_13TeVDown")
				
			print "writing to out file..."

			hggHtt.Write("ggH125"+systName)
			hvbfHtt.Write("qqH125"+systName)
			hWHtt.Write("WH125"+systName)
			hZHtt.Write("ZH125"+systName)

			hztt.Write("ZTT"+systName)
			hzl.Write("ZL"+systName)
			hzj.Write("ZJ"+systName)
			#hzll.Write("ZLL"+systName)
			httt.Write("TTT"+systName)
			httj.Write("TTJ"+systName)
			hw.Write("W"+systName)
			hewkz.Write("EWKZ"+systName)
			hvv.Write("VV"+systName)
			hqcd.Write("QCD"+systName)
			print "... done"

			if mode == "datacard":
				outfile.Close()

	if mode == "histos":
		outfile.Close()


