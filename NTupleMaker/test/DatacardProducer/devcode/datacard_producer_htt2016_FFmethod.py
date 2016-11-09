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
from ROOT import TFile, TH1F, TH1D, TCanvas, TPad, gStyle, kRed, kBlack, TLine, TLegend, TROOT, TChain

#def makeDatacard_fakerest(treeName, files, cut, weight, histo_fake, histo_rest):
#	print "fake"
#	makeDatacard_category(treeName, files, cut, "q_1 * q_2 < 0. && gen_match_2==6", weight, "_fake", histo_fake)
#	print "rest"
#	makeDatacard_category(treeName, files, cut, "q_1 * q_2 < 0. && gen_match_2!=6", weight, "_rest", histo_rest)

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
		#t = f.Get("TauCheck")
		t = f.Get(treeName)
		h = TH1F(fileName+catName, fileName+catName, nbins, xmin, xmax)
		t.Draw(var+" >> "+fileName+catName, weight + " * ("+ cut + " && "+additionalCut+")")
		xs = sampleDict["xsec"]
		hevents=f.Get("nWeightedEvents")
		nevents=hevents.GetBinContent(1)
		if (nevents == 0):
			print "nevents : ", nevents
		h.Scale(xs*lumi/nevents)
		hsample.Add(h)
	return hsample

def makeDatacard_sample_bin(treeName, sampleKey, sampleDict, catName, weight, cut, additionalCut):
	#t = TChain("TauCheck")
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
#lumi = 6260
#lumi = 3990
#directory with input root files (synch ntuples) with conventional names!

if channel == "mt":
	indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/NTuple_10Oct/mutau/final"
	indirData = indir

if channel == "et":
	indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/NTuple_10Oct/etau/final"
	indirData = indir

isBlinded = True 

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
                   }
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
			tdata.Draw(var+" >> hdata_ss", cvalue + " && q_1 * q_2 > 0.")

			#insert blinding to update cvalue in OS region
			blinding = {
			"m_vis":"(m_vis<60 || m_vis>120)", 
			"m_sv":"(m_sv<100 || m_sv>150)",
			"mt_tot":"(mt_tot<80 || mt_tot>130)"
			}
			if (isBlinded and (var in blinding) and ("inclusive" not in ckey)):
				blindingCut = blinding[var]
				print "applying blinding for " , var , " : " , blindingCut
			else:
				blindingCut = " 1 "

			print "cut data: ", cvalue, "&&" , blindingCut , "&& q_1 * q_2 < 0."
			print "Data OS: ", tdata.Draw(var+" >> hdata", cvalue + "&&" + blindingCut + "&& q_1 * q_2 < 0.")
			print "Data, done"

			#HTT section (signal)
			#hHtt = TH1F("HTT", "HTT"+axisTit, nbins, xmin, xmax)
			#hHtt.SetDirectory(0)
			#makeDatacard_category(treeName, lhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hHtt)
			#no need for SS for signal

			### SIGNAL ###

			#ggH
			hggHtt = TH1F("ggH125", "ggH125"+axisTit, nbins, xmin, xmax)
			hggHtt.SetDirectory(0)
			makeDatacard_category(treeName, lgghtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hggHtt)
			#VBF 
			hvbfHtt = TH1F("qqH125", "qqH125"+axisTit, nbins, xmin, xmax)
			hvbfHtt.SetDirectory(0)
			makeDatacard_category(treeName, lvbfhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hvbfHtt)
			#WH
			hWHtt = TH1F("WH125", "WH125"+axisTit, nbins, xmin, xmax)
			hWHtt.SetDirectory(0)
			makeDatacard_category(treeName, lwhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hWHtt)
			#ZH
			hZHtt = TH1F("ZH125", "ZH125"+axisTit, nbins, xmin, xmax)
			hZHtt.SetDirectory(0)
			makeDatacard_category(treeName, lzhtt, cvalue, "q_1 * q_2 < 0.", wvalue, "_os", hZHtt)
			#no need for SS for signal

			### BACKGROUNDS ###
			notFakeTauCut = "gen_match_2 != 6 "
			#VV section
			hvv = TH1F("VV", "VV"+axisTit, nbins, xmin, xmax)
			hvv.SetDirectory(0)
			makeDatacard_category(treeName, lvv, cvalue, "q_1 * q_2 < 0. &&"+notFakeTauIsoCut, wvalue, "_rest",hvv)
			print "VV, done"
			#TT section
			htt = TH1F("TT", "TT"+axisTit, nbins, xmin, xmax)
			htt.SetDirectory(0)
			makeDatacard_category(treeName, ltt, cvalue, "q_1 * q_2 < 0. &&"+notFakeTauIsoCut, wvalue, "_rest", htt)
			print "TT, done"
			#Wjets
			hw =  TH1F("W", "W"+axisTit,nbins, xmin, xmax)
			hw.SetDirectory(0)
			makeDatacard_category(treeName, lw, cvalue, "q_1 * q_2 < 0. &&"+notFakeTauIsoCut, wvalue, "_rest", hw)
			print "W, done"
			#ZTT, ZL, ZLL 
			hztt = TH1F("ztt", "ztt"+axisTit, nbins, xmin, xmax)
			hzl = TH1F("zl", "zl"+axisTit, nbins, xmin, xmax)
			hzll = TH1F("zll", "zll"+axisTit, nbins, xmin, xmax)
			
			hztt.SetDirectory(0)
			hzl.SetDirectory(0)		
			hzll.SetDirectory(0)

			#looping on different cuts that define categories 
			#for now, done by hand here with additional cut string
			ZTTcut = "gen_match_2 == 5" 
			ZLcut = "gen_match_2 < 5" 
			ZLLcut = "gen_match_2 != 5"

			makeDatacard_category(treeName, ldy, cvalue, ZTTcut+"&& q_1 * q_2 < 0. &&"+notFakeTauIsoCut , wvalue, "ZTT", hztt)
			makeDatacard_category(treeName, ldy, cvalue, ZLcut+"&& q_1 * q_2 < 0. &&"+notFakeTauIsoCut, wvalue, "ZL", hzl)
			makeDatacard_category(treeName, ldy, cvalue, ZLLcut+"&& q_1 * q_2 < 0. &&"+notFakeTauIsoCut, wvalue, "ZLL", hzll)

			print "DY all done"

			##########################
			# FF estimate
			##########################

			print "FF estimate ... "
			#define fakeweight
			fakeweight_dict = {"inclusive":"fakeweight_incl", 
							   "inclusivemt50":"fakeweight_incl", 
							   "0jet_low":"fakeweight_0jetlow",
							   "0jet_high":"fakeweight_0jethigh",
							   "1jet_low":"fakeweight_1jetlow",
							   "1jet_high":"fakeweight_1jethigh",
							   "vbf_low":"fakeweight_vbflow",
							   "vbf_high":"fakeweight_vbfhigh"} 
			try:
				fakeweight = fakeweight_dict[ckey]
			except KeyError:
				print "the category ", ckey, "does not have a corresponding fake factor weight in the dictionary "

			#data in anti-isolated region estimated with FF weights
			antiIsoCut = "byVLooseIsolationMVArun2v1DBoldDMwLT_2 ==1 && byTightIsolationMVArun2v1DBoldDMwLT_2 == 0"
			hdata_aSR = TH1F("hdata_aSR", "hdata_aSR"+axisTit, nbins, xmin, xmax)
			tdata.Draw(var+" >> hdata_aSR", fakeweight+"*("+cvalue + "  && q_1 * q_2 < 0. &&"+antiIsoCut+")") 

			aSRrestCut = "gen_match_2 != 6 && q_1*q_2<0 &&"+antiIsoCut
			#VV
			hvv_FFrestSR = TH1F("VV_FFrestSR", "VV_FFrestSR"+axisTit, nbins, xmin, xmax)
			hvv_FFrestSR.SetDirectory(0)
			makeDatacard_category(treeName,lvv,cvalue,aSRrestCut,wvalue*fakeweight,"VV_FFrestSR",hvv_FFrestSR)

			#TT
			htt_FFrestSR = TH1F("TT_FFrestSR", "TT_FFrestSR"+axisTit, nbins, xmin, xmax)
			htt_FFrestSR.SetDirectory(0)
			makeDatacard_category(treeName,ltt,cvalue,aSRrestCut,wvalue*fakeweight,"TT_FFrestSR",htt_FFrestSR)

			#W
			hw_FFrestSR = TH1F("W_FFrestSR", "W_FFrestSR"+axisTit, nbins, xmin, xmax)
			hw_FFrestSR.SetDirectory(0)
			makeDatacard_category(treeName,lw,cvalue,aSRrestCut,wvalue*fakeweight,"W_FFrestSR", hw_FFrestSR)

			#DY
			hdy_FFrestSR = TH1F("DY_FFrestSR", "DY_FFrestSR"+axisTit, nbins, xmin, xmax)
			hdy_FFrestSR.SetDirectory(0)
			makeDatacard_category(treeName,ldy,cvalue,aSRrestCut,wvalue*fakeweight,"DY_FFrestSR",hdy_FFrestSR)

			hff = TH1F("FF", "FF"+axisTit, nbins, xmin, xmax)
			hff.Add(hdata_aSR)

			hff.Add(hvv_FFrestSR,-1)
			hff.Add(htt_FFrestSR,-1)
			hff.Add(hw_FFrestSR,-1)
			hff.Add(hdy_FFrestSR,-1)

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

			hggHtt.Write("ggH125"+systName)
			hvbfHtt.Write("qqH125"+systName)
			hWHtt.Write("WH125"+systName)
			hZHtt.Write("ZH125"+systName)

			hztt.Write("ZTT"+systName)
			hzl.Write("ZL"+systName)
			hzll.Write("ZLL"+systName)
			htt.Write("TT"+systName)
			hw.Write("W"+systName)
			hvv.Write("VV"+systName)
			hff.Write("FF"+systName)
			print "... done"

			if mode == "datacard":
				outfile.Close()

	if mode == "histos":
		outfile.Close()


