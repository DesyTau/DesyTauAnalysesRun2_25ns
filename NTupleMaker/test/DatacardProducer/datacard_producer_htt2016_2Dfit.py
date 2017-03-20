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
	if categoryName in ["vbf", "boosted", "antiiso_vbf_cr", "antiiso_boosted_cr", "wjets_vbf_cr", "wjets_boosted_cr"]:
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
			#print " #######################"
			#print "Sample " , fileName
			#print "xs : " , xs, "lumi: ", lumi, "nevents: ", nevents , "total : ", xs*lumi/nevents
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
	global var
	#print "var= " , var
	tdata.Draw(var+" >> hdata_os_highmt", cvalue+"&&"+highMTcut+" && q_1 * q_2 < 0. ")
	tdata.Draw(var + " >> hdata_ss_highmt", cvalue+"&&"+highMTcut+" && q_1 * q_2 > 0. " )

	hvv_os_highmt = hdummy.Clone("hvv_os_highmt")
	hvv_ss_highmt = hdummy.Clone("hvv_ss_highmt")
	makeDatacard_ssos(treeName, lvv, cvalue+"&&"+highMTcut, wvalue, hvv_os_highmt, hvv_ss_highmt)
	hewkz_os_highmt = hdummy.Clone("hewkz_os_highmt")
	hewkz_ss_highmt = hdummy.Clone("hewkz_ss_highmt")
	makeDatacard_ssos(treeName, lewkz, cvalue+"&&"+highMTcut, wvalue, hewkz_os_highmt, hewkz_ss_highmt)

	hvbfhww_os_highmt = hdummy.Clone("hvbfhww_os_highmt")
	hvbfhww_ss_highmt = hdummy.Clone("hvbfhww_ss_highmt")
	makeDatacard_ssos(treeName, lvbfhww, cvalue+"&&"+highMTcut, wvalue, hvbfhww_os_highmt, hvbfhww_ss_highmt)

	hgghww_os_highmt = hdummy.Clone("hgghww_os_highmt")
	hgghww_ss_highmt = hdummy.Clone("hgghww_ss_highmt")
	makeDatacard_ssos(treeName, lgghww, cvalue+"&&"+highMTcut, wvalue, hgghww_os_highmt, hgghww_ss_highmt)

	htt_os_highmt = hdummy.Clone("htt_os_highmt")
	htt_ss_highmt = hdummy.Clone("htt_ss_highmt")
	makeDatacard_ssos(treeName, ltt, cvalue+"&&"+highMTcut, wvalue, htt_os_highmt, htt_ss_highmt)
	
	#Drell Yan components
	hztt_os_highmt = hdummy.Clone("hztt_os_highmt")
	hztt_ss_highmt = hdummy.Clone("hztt_ss_highmt")
	makeDatacard_ssos(treeName, ldy, cvalue+"&&"+highMTcut+"&&"+ZTTcut, wvalue+"*"+dy_weight, hztt_os_highmt, hztt_ss_highmt)
	hzl_os_highmt = hdummy.Clone("hzl_os_highmt")
	hzl_ss_highmt = hdummy.Clone("hzl_ss_highmt")
	#if (("0jet" in ckey) and datacardFor2Dfit):
	#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
	makeDatacard_ssos(treeName, ldy, cvalue+"&&"+highMTcut+"&&"+ZLcut, wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), hzl_os_highmt, hzl_ss_highmt)
	#var = varY+":"+varX
	hzj_os_highmt = hdummy.Clone("hzj_os_highmt")
	hzj_ss_highmt = hdummy.Clone("hzj_ss_highmt")
	makeDatacard_ssos(treeName, ldy, cvalue+"&&"+highMTcut+"&&"+ZJcut, wvalue+"*"+dy_weight, hzj_os_highmt, hzj_ss_highmt)
	# add Drell Yan together
	hdy_os_highmt = hdummy.Clone("hdy_os_highmt")
	hdy_os_highmt.Add(hztt_os_highmt)
	hdy_os_highmt.Add(hzl_os_highmt)
	hdy_os_highmt.Add(hzj_os_highmt)
	hdy_ss_highmt = hdummy.Clone("hdy_ss_highmt")	
	hdy_ss_highmt.Add(hztt_ss_highmt)
	hdy_ss_highmt.Add(hzl_ss_highmt)
	hdy_ss_highmt.Add(hzj_ss_highmt)
	#makeDatacard_ssos(treeName, ldy, cvalue+"&&"+highMTcut, wvalue+"*"+dy_weight, hdy_os_highmt, hdy_ss_highmt)

	hw_os_highmt = hdummy.Clone("hw_os_highmt")
	hw_ss_highmt = hdummy.Clone("hw_ss_highmt")		
	makeDatacard_ssos(treeName, lw, cvalue+"&&"+highMTcut, wvalue+"*"+dy_weight, hw_os_highmt, hw_ss_highmt)
	
	hqcd_ss_highmt = hdummy.Clone("hqcd_ss_highmt")
	add_all(hqcd_ss_highmt, [(hdata_ss_highmt,1),(hvv_ss_highmt,-1),(htt_ss_highmt,-1),(hdy_ss_highmt,-1),(hw_ss_highmt,-1),(hewkz_ss_highmt,-1), (hvbfhww_ss_highmt,-1), (hgghww_ss_highmt,-1)])	
			
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
		makeDatacard_category(treeName, lewkz, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue+"*"+dy_weight, "_highmt_ss_rel", hewkz_ss_highmt_rel )
		htt_ss_highmt_rel = hdummy.Clone("htt_ss_highmt_rel")
		makeDatacard_category(treeName, ltt, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", htt_ss_highmt_rel )

		hztt_ss_highmt_rel = hdummy.Clone("hztt_ss_highmt_rel")
		hzl_ss_highmt_rel = hdummy.Clone("hzl_ss_highmt_rel")
		hzj_ss_highmt_rel = hdummy.Clone("hzj_ss_highmt_rel")
		makeDatacard_category(treeName, ldy, relaxedSel+"&&"+highMTcut+"&&"+ZTTcut, "q_1*q_2>0", wvalue+"*"+dy_weight, "_ztt_highmt_ss_rel", hztt_ss_highmt_rel )
		#if (("0jet" in ckey) and datacardFor2Dfit):
		#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
		makeDatacard_category(treeName, ldy, relaxedSel+"&&"+highMTcut+"&&"+ZLcut, "q_1*q_2>0", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "_zl_highmt_ss_rel", hzl_ss_highmt_rel )
		#var = varY+":"+varX
		makeDatacard_category(treeName, ldy, relaxedSel+"&&"+highMTcut+"&&"+ZJcut, "q_1*q_2>0", wvalue+"*"+dy_weight, "_zj_highmt_ss_rel", hzj_ss_highmt_rel )
		hdy_ss_highmt_rel = hdummy.Clone("hdy_ss_highmt_rel")
		hdy_ss_highmt_rel.Add(hztt_ss_highmt_rel)
		hdy_ss_highmt_rel.Add(hzl_ss_highmt_rel)
		hdy_ss_highmt_rel.Add(hzj_ss_highmt_rel)		
		#makeDatacard_category(treeName, ldy, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue+"*"+dy_weight, "_highmt_ss_rel", hdy_ss_highmt_rel )

		hw_ss_highmt_rel = hdummy.Clone("hw_ss_highmt_rel")
		makeDatacard_category(treeName, lw, relaxedSel+"&&"+highMTcut, "q_1*q_2>0",  wvalue, "_highmt_ss_rel", hw_ss_highmt_rel )
		hvbfhww_ss_highmt_rel = hdummy.Clone("hvbfhww_ss_highmt_rel")
		makeDatacard_category(treeName, lvbfhww, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", hvbfhww_ss_highmt_rel )
		hgghww_ss_highmt_rel = hdummy.Clone("hgghww_ss_highmt_rel")
		makeDatacard_category(treeName, lgghww, relaxedSel+"&&"+highMTcut, "q_1*q_2>0", wvalue, "_highmt_ss_rel", hgghww_ss_highmt_rel )
	#2
		hqcd_ss_highmt_rel = hdummy.Clone("hqcd_ss_highmt_rel")
		add_all(hqcd_ss_highmt_rel, [(hdata_ss_highmt_rel,1),(hvv_ss_highmt_rel,-1),(htt_ss_highmt_rel,-1),(hdy_ss_highmt_rel,-1),(hw_ss_highmt_rel,-1),(hewkz_ss_highmt_rel,-1), (hvbfhww_ss_highmt_rel,-1), (hgghww_ss_highmt_rel,-1)])	
	#3
		hqcd_ss_highmt_rel.Scale(hqcd_ss_highmt.Integral()/(hqcd_ss_highmt_rel.Integral()))
		hqcd_ss_highmt = hqcd_ss_highmt_rel

	hw_os_highmt_datadriven = hdummy.Clone("hw_os_highmt_datadriven") 
	add_all(hw_os_highmt_datadriven, [(hdata_os_highmt,1),(hvv_os_highmt,-1),(htt_os_highmt,-1),(hdy_os_highmt,-1),(hewkz_os_highmt,-1), (hvbfhww_os_highmt,-1), (hgghww_os_highmt,-1), (hqcd_ss_highmt, -1*SSOSratio)] )

	# apply some relaxed selection. check which. For now taking the usual one. 
	relaxed_cuts_w_MC = {
	"mt_inclusive": "pt_2 > 30 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
	"mt_0jet": "pt_2 > 30 && njets==0 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
	"mt_boosted":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 ) ) ) && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
	"mt_vbf":"pt_2>40 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
    "mt_wjets_boosted_cr":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 ) ) ) && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
	"mt_wjets_vbf_cr":"pt_2>40 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
	"mt_antiiso_boosted_cr":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 )) ) && pt_1>20 && iso_1> 0.15 && iso_1<0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
    "mt_antiiso_vbf_cr":"pt_2>40 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>20 && iso_1> 0.15 && iso_1<0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
	"et_inclusive":"pt_2>30 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
	"et_0jet":"pt_2>30 && njets==0 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
	"et_boosted":"pt_2>30 &&  (njets==1 || (njets>=2 &&(mjj<300 || pt_sv<=50) ) )  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
	"et_vbf":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0" , 
	"et_wjets_boosted_cr":"pt_2>30 &&  (njets==1 || (njets>=2 &&(mjj<300 || pt_sv<=50) ) )  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
	"et_wjets_vbf_cr":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0" , 
    "et_antiiso_boosted_cr":"pt_2>30 && (njets==1 || (njets>=2 &&(mjj<300 ||pt_sv<=50) ) ) && pt_1>26 && fabs(eta_1)<2.1 && iso_1 > 0.1 && iso_1<0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ",
   "et_antiiso_vbf_cr":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 > 0.1 && iso_1<0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 "		
	}



	hw_os_MC = hdummy.Clone("hw_os_MC")
	hw_os_highmt_MC = hdummy.Clone("hw_os_highmt_MC")
	makeDatacard_category(treeName, lw, relaxed_cuts_w_MC.get(channel+"_"+ckey, cvalue), "q_1*q_2<0"+lowMTcut,       wvalue, "hw_os_MC", hw_os_MC )
	makeDatacard_category(treeName, lw, relaxed_cuts_w_MC.get(channel+"_"+ckey, cvalue), "q_1*q_2<0"+"&&"+highMTcut, wvalue, "hw_os_highmt_MC", hw_os_highmt_MC)

	# W scale factor 
	W_sf =  hw_os_highmt_datadriven.Integral()/hw_os_highmt_MC.Integral()
	hw.Scale(hw_os_MC.Integral()/hw.Integral()) 
	#adjust yield to norm_w = W_sf*hw_os_MC.Integral()
	print " --- W SF = ", W_sf
	return W_sf



def QCD_datadriven(SSOSratio, W_sf, hqcd): #hqcd contains the result of the QCD calculation
	hqcd_ss = hdummy.Clone("QCD_ss")
	add_all(hqcd_ss, [(hdata_ss,1),(hvv_ss,-1),(htt_ss,-1),(hdy_ss,-1),(hw_ss,-1*W_sf),(hewkz_ss,-1), (hvbfhww_ss,-1), (hgghww_ss,-1)])
	hqcd.Add(hqcd_ss, SSOSratio)

	if (relaxedSel is not None): 
		print "QCD relaxed selection"
		hqcd_rel_ss = hdummy.Clone("QCD_rel_ss")
		add_all(hqcd_rel_ss, [(hdata_rel_ss,1),(hvv_rel_ss,-1),(htt_rel_ss,-1),(hdy_rel_ss,-1),(hw_rel_ss,-1*W_sf),(hewkz_rel_ss,-1),(hvbfhww_rel_ss,-1),(hgghww_rel_ss,-1) ])
		norm_nom = hqcd_ss.Integral()
		norm_rel = hqcd_rel_ss.Integral()
		hqcd_rel_ss.Scale(SSOSratio*norm_nom/norm_rel)
		hqcd = hqcd_rel_ss

def make_signal_histo(histo, dataset_dict):
	histo.SetDirectory(0)
	makeDatacard_category(treeName, dataset_dict, cvalue + lowMTcut , "q_1 * q_2 < 0.", wvalue, "_os", histo)

	
def make_renormscale_shape(histo_up, histo_down, dataset_dict): 
	'''
	For gg fusion signal, estimating renormalization uncertainty for ggH process. 
	'''
	print "calculating renormalisation uncertainty"
	renorm_weight_up = renorm_weights.get(channel+"_"+ckey,"1")
	renorm_weight_down = "(2-( "+renorm_weight_up+"))"
	histo_up.SetDirectory(0)
	histo_down.SetDirectory(0)
	makeDatacard_category(treeName, dataset_dict, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_up, "_gg_up", histo_up)
	makeDatacard_category(treeName, dataset_dict, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_down, "_gg_down", histo_down )

def MTcut(region, category):
	mtcut = {"low":"&& mt_1<50" , "high":"mt_1>80"}
	if "wjets_" in category:
		mtcut["low"]= "&& mt_1>80"
	return mtcut[region]

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
lumi = 35870

#directory with input root files (synch ntuples) with conventional names!

if channel == "mt":
	indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Mar/mutau/final"
	indirData = "/nfs/dust/cms/user/bottav/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Mar/mutau/final"

if channel == "et":
	indir = "/nfs/dust/cms/user/bottav/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Mar/etau/final"
	indirData = "/nfs/dust/cms/user/bottav/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Mar/etau/final"

isBlinded = False 
datacardFor2Dfit = True

ROOT.TH1.SetDefaultSumw2(True)

##########################
# get data tree 
##########################

treeName = "TauCheck"

if channel == "mt":
	fdata =  TFile(indirData+"/"+"SingleMuon_Run2016BtoH.root","read")
if channel == "et":
	fdata = TFile(indirData+"/"+"SingleElectron_Run2016BtoH.root","read")

tdata = fdata.Get(treeName)


#redefine tree name for MC
treeName = "TauCheck"+systName


####################################################
# scale factors and up/down variations
# other selection cuts
# NB: more selection cuts in relaxed selections
####################################################
'''
ssosvalues = {"et_0jet": 1.00, "et_boosted": 1.15, "et_vbf": 1.20, 
			  "mt_0jet": 1.00, "mt_boosted": 1.15, "mt_vbf": 1.20,
			  "et_wjets_0jet_cr": 1.00, "et_wjets_boosted_cr": 1.15, "et_wjets_vbf_cr": 1.20, 
			  "mt_wjets_0jet_cr": 1.00, "mt_wjets_boosted_cr": 1.15, "mt_wjets_vbf_cr": 1.20,
			  "et_antiiso_0jet_cr": 1.00, "et_antiiso_boosted_cr": 1.15, "et_antiiso_vbf_cr": 1.20, 
			  "mt_antiiso_0jet_cr": 1.00, "mt_antiiso_boosted_cr": 1.15, "mt_antiiso_vbf_cr": 1.20}


SSOS_variations = {"mt_0jet":{"up":1.15, "down":0.85}, "mt_boosted":{"up":1.15, "down":0.85}, "mt_vbf":{"up":1.30, "down":0.70}, 
				   "et_0jet":{"up":1.15, "down":0.85}, "et_boosted":{"up":1.15, "down":0.85}, "et_vbf":{"up":1.30, "down":0.70},
				   "et_wjets_0jet_cr":{"up":1.15, "down":0.85}, "et_wjets_boosted_cr":{"up":1.15, "down":0.85}, "et_wjets_vbf_cr":{"up":1.30, "down":0.70},
				   "mt_wjets_0jet_cr":{"up":1.15, "down":0.85}, "mt_wjets_boosted_cr":{"up":1.15, "down":0.85}, "mt_wjets_vbf_cr":{"up":1.30, "down":0.70},
				   "et_antiiso_0jet_cr":{"up":1.15, "down":0.85}, "et_antiiso_boosted_cr":{"up":1.15, "down":0.85}, "et_antiiso_vbf_cr":{"up":1.30, "down":0.70},
				   "mt_antiiso_0jet_cr":{"up":1.15, "down":0.85}, "mt_antiiso_boosted_cr":{"up":1.15, "down":0.85}, "mt_aniiso_vbf_cr":{"up":1.30, "down":0.70}}
'''
#new values for full 2016 dataset. Up and Down are NOT multiplied by the central value. 
ssosvalues = {"et_0jet": 1.00, "et_boosted": 1.28, "et_vbf": 1.00, 
			  "mt_0jet": 1.07, "mt_boosted": 1.06, "mt_vbf": 1.00,
			  "et_wjets_0jet_cr": 1.00, "et_wjets_boosted_cr": 1.28, "et_wjets_vbf_cr": 1.00, 
			  "mt_wjets_0jet_cr": 1.07, "mt_wjets_boosted_cr": 1.06, "mt_wjets_vbf_cr": 1.00,
			  "et_antiiso_0jet_cr": 1.00, "et_antiiso_boosted_cr": 1.28, "et_antiiso_vbf_cr": 1.00, 
			  "mt_antiiso_0jet_cr": 1.07, "mt_antiiso_boosted_cr": 1.06, "mt_antiiso_vbf_cr": 1.00}


SSOS_variations = {"mt_0jet":{"up":1.11, "down":1.03}, "mt_boosted":{"up":1.12, "down":1.00}, "mt_vbf":{"up":1.20, "down":0.80}, 
				   "et_0jet":{"up":1.07, "down":0.93}, "et_boosted":{"up":1.41, "down":1.15}, "et_vbf":{"up":1.20, "down":0.80},
				   "et_wjets_0jet_cr":{"up":1.07, "down":0.93}, "et_wjets_boosted_cr":{"up":1.41, "down":1.15}, "et_wjets_vbf_cr":{"up":1.20, "down":0.80},
				   "mt_wjets_0jet_cr":{"up":1.11, "down":1.03}, "mt_wjets_boosted_cr":{"up":1.12, "down":1.00}, "mt_wjets_vbf_cr":{"up":1.20, "down":0.80},
				   "et_antiiso_0jet_cr":{"up":1.07, "down":0.93}, "et_antiiso_boosted_cr":{"up":1.41, "down":1.15}, "et_antiiso_vbf_cr":{"up":1.20, "down":0.80},
				   "mt_antiiso_0jet_cr":{"up":1.11, "down":1.03}, "mt_antiiso_boosted_cr":{"up":1.12, "down":1.00}, "mt_aniiso_vbf_cr":{"up":1.20, "down":0.80}}


Wsf_variations = {"0jet":{"up":1.10, "down":0.90}, "boosted":{"up":1.10, "down":0.90}, "vbf":{"up":1.30, "down":0.70},
                  "wjets_0jet_cr":{"up":1.10, "down":0.90}, "wjets_boosted_cr":{"up":1.10, "down":0.90}, "wjets_vbf_cr":{"up":1.30, "down":0.70},
                  "antiiso_0jet_cr":{"up":1.10, "down":0.90}, "antiiso_boosted_cr":{"up":1.10, "down":0.90}, "antiiso_vbf_cr":{"up":1.30, "down":0.70}}

renorm_weights = {"et_0jet":"(0.929+0.0001702*pt_2)", "et_boosted":"(0.919+0.0010055*pt_sv)", "et_vbf":"(1.026+0.000066*mjj)", 
				  "mt_0jet":"(0.929+0.0001702*pt_2)", "mt_boosted":"(0.919+0.0010055*pt_sv)", "mt_vbf":"(1.026+0.000066*mjj)",
                  "et_wjets_0jet_cr":"(0.929+0.0001702*pt_2)", "et_wjets_boosted_cr":"(0.919+0.0010055*pt_sv)", "et_wjets_vbf_cr":"(1.026+0.000066*mjj)", 
				  "mt_wjets_0jet_cr":"(0.929+0.0001702*pt_2)", "mt_wjets_boosted_cr":"(0.919+0.0010055*pt_sv)", "mt_wjets_vbf_cr":"(1.026+0.000066*mjj)",
                  "et_antiiso_0jet_cr":"(0.929+0.0001702*pt_2)", "et_antiiso_boosted_cr":"(0.919+0.0010055*pt_sv)", "et_antiiso_vbf_cr":"(1.026+0.000066*mjj)", 
				  "mt_antiiso_0jet_cr":"(0.929+0.0001702*pt_2)", "mt_antiiso_boosted_cr":"(0.919+0.0010055*pt_sv)", "mt_antiiso_vbf_cr":"(1.026+0.000066*mjj)"}


zmumu_CR_weights = {"vbf":               {"binVar":"mjj", "weights":{(300.,700.):1.070, (700.,1000.):1.090, (1100.,1500.):1.055, (1500.,14000):1.015}}, 
					"wjets_vbf_cr":      {"binVar":"mjj", "weights":{(300.,700.):1.070, (700.,1000.):1.090, (1100.,1500.):1.055, (1500.,14000):1.015}}, 
					"antiiso_vbf_cr":    {"binVar":"mjj", "weights":{(300.,700.):1.070, (700.,1000.):1.090, (1100.,1500.):1.055, (1500.,14000):1.015}},
					"down_shape":        {"binVar":"mjj", "weights":{(300.,700.):1.1449, (700.,1000.):1.1881, (1100.,1500.):1.113, (1500.,14000):1.030}}}

# correct mvis shape for ZL in 0jet. Applied to nominal value. 
'''
mvis_zl_corr = {"mt": "(1+(0.01*(tau_decay_mode_2==0)))", 
				"mt_1prong_up": "(1+(0.013*(tau_decay_mode_2==0)))", #1.01*1.003 is the up shift 
				"mt_1prong_down": "(1+(0.007*(tau_decay_mode_2==0)))", #1.01*0.997 is the down shift
				"mt_1prong1pi_up": "(1+(0.01*(tau_decay_mode_2==0))+(0.004*(tau_decay_mode_2==1)))", # 1.004
				"mt_1prong1pi_down": "(1+(0.01*(tau_decay_mode_2==0))+(-0.004*(tau_decay_mode_2==1)))", #0.996
				"et": "(1+((0.017*(tau_decay_mode_2==0))+(0.03*(tau_decay_mode_2==1))))",
				"et_1prong_up": "(1+((0.022*(tau_decay_mode_2==0))+(0.03*(tau_decay_mode_2==1))))", #1.017*1.005
				"et_1prong_down": "(1+((0.0119*(tau_decay_mode_2==0))+(0.03*(tau_decay_mode_2==1))))", #1.017*0.995
				"et_1prong1pi_up": "(1+((0.017*(tau_decay_mode_2==0))+(0.035*(tau_decay_mode_2==1))))", #1.03*1.005
				"et_1prong1pi_down": "(1+((0.017*(tau_decay_mode_2==0))+(0.02485*(tau_decay_mode_2==1))))" #1.03*0.995
				}

'''
ZL_yield_SF = {"mt": "((0.74*(tau_decay_mode_2==0))+(1*(tau_decay_mode_2!=0)))",
				"mt_1prong_up":"((0.74*1.25*(tau_decay_mode_2==0))+(1*(tau_decay_mode_2!=0)))", 
				"mt_1prong_down":"((0.74*0.75*(tau_decay_mode_2==0))+(1*(tau_decay_mode_2!=0)))",
				"mt_1prong1pi_up":"((0.74*(tau_decay_mode_2==0))+(1.25*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))", 
				"mt_1prong1pi_down":"((0.74*(tau_decay_mode_2==0))+(0.75*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))",
			   "et": "((0.98*(tau_decay_mode_2==0))+(1.20*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))",
				"et_1prong_up":"((0.98*1.12*(tau_decay_mode_2==0))+(1.20*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))", 
				"et_1prong_down":"((0.98*0.88*(tau_decay_mode_2==0))+(1.20*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))",
				"et_1prong1pi_up":"((0.98*(tau_decay_mode_2==0))+(1.20*1.12*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))", 
				"et_1prong1pi_down":"((0.98*(tau_decay_mode_2==0))+(1.20*0.88*(tau_decay_mode_2==1))+(1*(tau_decay_mode_2==10)))"
			  }

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
     "WZJTo3L1Nu" : {"xsec": 4.708, "isBinned": False, "binVar": "", "files": ["WZJTo3L1Nu"]},
     "WWTo1L1Nu2Q" : {"xsec": 49.997, "isBinned": False, "binVar": "", "files": ["WWTo1L1Nu2Q"]},
     "WZTo1L1Nu2Q" : {"xsec": 10.71, "isBinned": False, "binVar": "", "files": ["WZTo1L1Nu2Q"]},
     "VVTo2L2Nu" : {"xsec": 11.95, "isBinned": False, "binVar": "", "files": ["VVTo2L2Nu"]},
     "WZTo2L2Q" : {"xsec": 5.595, "isBinned": False, "binVar": "", "files": ["WZTo2L2Q"]}
}

ltt = {"TTPowHeg": {"xsec": 831.76, "files": ["TTpowheg"]}}

ldy = {"DYJetsToLL_M-50": { "isBinned": True, "binVar": "gen_noutgoing", 
                            "files": ["DYJetsToLL_M-50_MG_ext1", "DYJetsToLL_M-50_MG_ext2", 
                                      "DY1JetsToLL_M-50_MG",
                                      "DY2JetsToLL_M-50_MG",
                                      "DY3JetsToLL_M-50_MG",
                                      "DY4JetsToLL_M-50_MG"], 
                            "weights":{ (-0.5, 0.5) : 0.0000395423,
                                        (0.5, 1.5) :  0.0000128925,
                                        (1.5, 2.5) :  0.0000130801,
                                        (2.5, 3.5) :  0.0000130062,
                                        (3.5, 100.5) : 0.0000105677	
                              }
                    },
       #"DYJetsToLL_M-10to50" : {"xsec": 18610, "isBinned": False, "binVar": "", "files": ["DYJetsToLL_M-10to50_MG"]}
}

lw = {"WJetsToLNu_MG": {"isBinned": True, "binVar": "gen_noutgoing", 
                              "files": ["WJetsToLNu_MG",
                                        "W1JetsToLNu_MG",
                                        "W2JetsToLNu_MG", "W2JetsToLNu_MG_ext1",
                                        "W3JetsToLNu_MG", "W3JetsToLNu_MG_ext1",
                                        "W4JetsToLNu_MG", "W4JetsToLNu_MG_ext1", "W4JetsToLNu_MG_ext2"], 
                              "weights":{ (-0.5, 0.5) : 0.0007093903,
                                          (0.5, 1.5) : 0.0001870249, 
                                          (1.5, 2.5) : 0.0000577854,
                                          (2.5, 3.5) : 0.0000190373,
                                          (3.5, 100.5) : 0.0000205964
                                  }
                   },
	  "EWKWMinus2Jets_WToLNu_M-50":{"xsec": 20.25, "isBinned":False, "files":["EWKWMinus2Jets_WToLNu_M-50_all"]},
	  "EWKWPlus2Jets_WToLNu_M-50":{"xsec": 25.62, "isBinned":False, "files":["EWKWPlus2Jets_WToLNu_M-50_all"]}
}

lewkz = {"EWKZ2Jets_ZToLL_M-50": {"xsec": 3.987, "isBinned":False, "files":["EWKZ2Jets_ZToLL_M-50_all"]},
		 "EWKZ2Jets_ZToNuNu":    {"xsec": 10.01, "isBinned":False, "files":["EWKZ2Jets_ZToNuNu_ext1"]}
}


lvbfhww = {"VBFHToWWTo2L2Nu_M125": {"xsec": 3.782*0.2137*0.324*0.324, "isBinned":False, "files":["VBFHToWWTo2L2Nu_M125"]}}

lgghww = {"GluGluHToWWTo2L2Nu_M125": {"xsec": 48.58*0.2137*0.324*0.324, "isBinned":False, "files":["GluGluHToWWTo2L2Nu_M125"]}}


#Higgs to tau tau samples
#MH=125
lvbfhtt = {"VBFHTT_M125" : {"xsec": 3.782*0.0627, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M125"]}}
lwhtt = {"WmHTT_M125" : {"xsec": 0.5328*0.0627, "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M125"]},
         "WpHTT_M125" : {"xsec": 0.8400*0.0627, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M125"]} }
lzhtt = {"ZHTT_M125" : {"xsec": 0.8839*0.0627, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M125"]}}
lgghtt = {"GluGluHTT_M125" : {"xsec": 48.58*0.0627, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M125"]}}

#MH= 110
lvbfhtt_m110 = {"VBFHTT_M110" : {"xsec":4.434*0.0791, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M110"]}}
lwhtt_m110 = {"WmHTT_M110" : {"xsec":0.8587*0.0791, "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M110"]},
              "WpHTT_M110" : {"xsec":1.335*0.0791, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M110"]}}
lzhtt_m110 = {"ZHTT_M110" : {"xsec":1.309*0.0791, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M110"]}}
lgghtt_m110 = {"GluGluHTT_M110" : {"xsec": 57.90*0.0791, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M110"]}}

#MH= 120
lvbfhtt_m120 = {"VBFHTT_M120" : {"xsec":3.935*0.0698, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M120"]}}
lwhtt_m120 = {"WmHTT_M120" : {"xsec":0.6092*0.0698, "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M120"]},
              "WpHTT_M120" : {"xsec":0.9558*0.0698, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M120"]}}
lzhtt_m120 = {"ZHTT_M120" : {"xsec":0.9939*0.0698, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M120"]}}
lgghtt_m120 = {"GluGluHTT_M120" : {"xsec": 52.22*0.0698, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M120"]}}

#MH= 130
lvbfhtt_m130 = {"VBFHTT_M130" : {"xsec": 3.637*0.0541, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M130"]}}
lwhtt_m130 = {"WmHTT_M130" : {"xsec": 0.4676*0.0541 , "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M130"]},
              "WpHTT_M130" : {"xsec": 0.7414*0.0541, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M130"]}}
lzhtt_m130 = {"ZHTT_M130" : {"xsec":0.7899*0.0541, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M130"]}}
lgghtt_m130 = {"GluGluHTT_M130" : {"xsec":45.31*0.0541, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M130"]}}

#MH= 140
lvbfhtt_m140 = {"VBFHTT_M140" : {"xsec":3.492*0.0360, "isBinned": False, "binVar": "", "files": ["VBFHToTauTau_M140"]}}
lwhtt_m140 = {"WmHTT_M140" : {"xsec":0.3940*0.0360, "isBinned": False, "binVar": "", "files": ["WminusHToTauTau_M140"]},
              "WpHTT_M140" : {"xsec":0.6308*0.0360, "isBinned": False, "binVar": "", "files": ["WplusHToTauTau_M140"]}
              }
lzhtt_m140 = {"ZHTT_M140" : {"xsec":0.6514*0.0360, "isBinned": False, "binVar": "", "files": ["ZHToTauTau_M140"]}}
lgghtt_m140 = {"GluGluHTT_M140" : {"xsec": 36.00*0.0360, "isBinned": False, "binVar": "", "files": ["GluGluHToTauTau_M140"]}}


if mode == "datacard": 
	for hkey, hvalue in histos.iteritems():
		for wkey, wvalue in weighting.iteritems():
			selection_name = sys.argv[2].replace(".json","")
			filename = selection_name+"_"+wkey+"_htt_"+channel+".inputs-sm-13TeV-"+hkey+systName+".root"
			#remove old files with same name as the one that will be created, if any
			os.system("if [ -f "+filename+" ]; then rm "+filename+";fi")

for wkey, wvalue in weighting.iteritems():
	print "weight : ", wkey
	if mode == "histos":
		suffix = os.path.splitext(os.path.basename(sys.argv[4]))[0] #add name of the (histos)json file as suffix
		outfile=TFile( wkey+"_htt_"+channel+".inputs-sm-13TeV-histos-"+suffix+".root", "recreate")	

	for ckey, cvalue in cuts.iteritems():
		print "cut : ", ckey

		lowMTcut = MTcut("low",ckey)
		highMTcut = MTcut("high",ckey)

		#########################################
		# relaxed selection for some categories #
		#########################################

		relaxed_cuts = {
			"mt_boosted":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 ) ) ) && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
			"mt_vbf":"pt_2>40 && njets>=2 && mjj>300 &&  pt_sv>50 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
            "mt_wjets_boosted_cr":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 ) ) ) && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0", 
			"mt_wjets_vbf_cr":"pt_2>40 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>20 && iso_1 < 0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
			"mt_antiiso_boosted_cr":"pt_2>30 &&(njets==1 || (njets>=2 && (mjj<300 || pt_2<=40 || pt_sv<=50 )) ) && pt_1>20 && iso_1> 0.15 && iso_1<0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
    		"mt_antiiso_vbf_cr":"pt_2>40 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>20 && iso_1> 0.15 && iso_1<0.3 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && againstElectronVLooseMVA6_2>0.5 && againstMuonTight3_2>0.5  && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
			"et_boosted":"pt_2>30 &&  (njets==1 || (njets>=2 &&(mjj<300|| pt_sv<=50) ) ) && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
			"et_vbf":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0",
			"et_wjets_boosted_cr":"pt_2>30 &&  (njets==1 || (njets>=2 &&(mjj<300 || pt_sv<=50) ) )  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ", 
			"et_wjets_vbf_cr":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50  && pt_1>26 && fabs(eta_1)<2.1 && iso_1 < 0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0" , 
    		"et_antiiso_boosted_cr":"pt_2>30 && (njets==1 || (njets>=2 &&(mjj<300 ||pt_sv<=50) ) ) && pt_1>26 && fabs(eta_1)<2.1 && iso_1 > 0.1 && iso_1<0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 ",
   			"et_antiiso_vbf_cr":"pt_2>30 && njets>=2 && mjj>300 && pt_sv>50 && pt_1>26 && fabs(eta_1)<2.1 && iso_1 > 0.1 && iso_1<0.3 && againstMuonLoose3_2>0.5 && againstElectronTightMVA6_2>0.5 && byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && dilepton_veto == 0 && extraelec_veto == 0 && extramuon_veto == 0 "		
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
				"inclusive": {"varX":"pt_2","binX":[30,35,40,45,50,55,300], "axistitX":"tau p_{T}", 
						      "varY":"m_vis","binY":[0,60,65,70,75,80,85,90,95,100,105,110,400], "axistitY":"m_{vis}"},
				"0jet": {#"varX":"pt_2","binX":[30,35,40,45,50,55,300], "axistitX":"tau p_{T}", 
						 "varX":"tau_decay_mode_2","binX":[0,0.1,1.1,10.1], "axistitX":"tau DM", 
						 "varY":"m_vis","binY":[0,60,65,70,75,80,85,90,95,100,105,110,400], "axistitY":"m_{vis}"},

				"boosted": {"varX":"pt_sv","binX":[0,100,150,200,250,300,5000], "axistitX":"Higgs p_{T}", 
							"varY":"m_sv","binY":[0,80,90,100,110,120,130,140,150,160,300], "axistitY":"m_{#tau#tau}"},
				"vbf": {"varX":"mjj","binX":[300,700,1100,1500,10000], "axistitX":"m_{jj}",
						"varY":"m_sv","binY":[0,95,115,135,155,400], "axistitY":"m_{#tau#tau}"},
				#"boosted": {"varX":"pt_tt","binX":[0,100,150,200,250,300,5000], "axistitX":"Higgs p_{T}", 
				#			"varY":"m_vis","binY":[0,80,90,100,110,120,130,140,150,160,300], "axistitY":"m_{#tau#tau}"},
				#"vbf": {"varX":"mjj","binX":[300,700,1100,1500,10000], "axistitX":"m_{jj}",
				#		"varY":"m_vis","binY":[0,95,115,135,155,400], "axistitY":"m_{#tau#tau}"},
				"antiiso_0jet_cr": {"varX":"pt_2","binX":[0,14000], "axistitX":"tau p_{T}", 
						            "varY":"m_vis","binY":[40,80,120,160,200], "axistitY":"m_{vis}"},
 									###"varY":"m_vis","binY":[0,40,80,120,160,200,240,280,320], "axistitY":"m_{vis}"},
				"antiiso_boosted_cr": {"varX":"pt_sv","binX":[0, 14000], "axistitX":"Higgs p_{T}", 
							           "varY":"m_sv","binY":[40,80,120,160,200], "axistitY":"m_{#tau#tau}"},
							           ###"varY":"m_sv","binY":[0,40,80,120,160,200,240,280,320], "axistitY":"m_{#tau#tau}"},
				"antiiso_vbf_cr": {"varX":"mjj","binX":[0,10000], "axistitX":"m_{jj}",
						           "varY":"m_sv","binY":[40,80,120,160,200], "axistitY":"m_{#tau#tau}"},
							      ###"varY":"m_sv","binY":[0,40,80,120,160,200,240,280,320], "axistitY":"m_{#tau#tau}"},
				"wjets_0jet_cr":{"varX":"pt_2","binX":[0,14000], "axistitX":"tau p_{T}", 
						         "varY":"mt_1","binY":[80,200], "axistitY":"m_{T,1}"},
				"wjets_boosted_cr":{"varX":"pt_2","binX":[0,14000], "axistitX":"tau p_{T}", 
						            "varY":"mt_1","binY":[80,200], "axistitY":"m_{T,1}"},
				"wjets_vbf_cr":{"varX":"pt_2","binX":[0,14000], "axistitX":"tau p_{T}", 
						        "varY":"mt_1","binY":[80,200], "axistitY":"m_{T,1}"},
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
			
			#MH=125
			hggHtt = hdummy.Clone("ggH125")
			hvbfHtt = hdummy.Clone("qqH125")
			hWHtt = hdummy.Clone("hWH125")
			hZHtt = hdummy.Clone("hZH125")
			make_signal_histo(hggHtt, lgghtt)
			make_signal_histo(hvbfHtt, lvbfhtt)
			make_signal_histo(hWHtt, lwhtt)
			make_signal_histo(hZHtt, lzhtt)
			#MH=130
			hggHtt_m130 = hdummy.Clone("ggH130")
			hvbfHtt_m130 = hdummy.Clone("qqH130")
			hWHtt_m130 = hdummy.Clone("hWH130")
			hZHtt_m130 = hdummy.Clone("hZH130")
			make_signal_histo(hggHtt_m130, lgghtt_m130)
			make_signal_histo(hvbfHtt_m130, lvbfhtt_m130)
			make_signal_histo(hWHtt_m130, lwhtt_m130)
			make_signal_histo(hZHtt_m130, lzhtt_m130)
			#MH=120
			hggHtt_m120 = hdummy.Clone("ggH120")
			hvbfHtt_m120 = hdummy.Clone("qqH120")
			hWHtt_m120 = hdummy.Clone("hWH120")
			hZHtt_m120 = hdummy.Clone("hZH120")
			make_signal_histo(hggHtt_m120, lgghtt_m120)
			make_signal_histo(hvbfHtt_m120,lvbfhtt_m120)
			make_signal_histo(hWHtt_m120, lwhtt_m120)
			make_signal_histo(hZHtt_m120, lzhtt_m120)
			#MH=110
			hggHtt_m110 = hdummy.Clone("ggH110")
			hvbfHtt_m110 = hdummy.Clone("qqH110")
			hWHtt_m110 = hdummy.Clone("hWH110")
			hZHtt_m110 = hdummy.Clone("hZH110")
			make_signal_histo(hggHtt_m110, lgghtt_m110)
			make_signal_histo(hvbfHtt_m110,lvbfhtt_m110)
			make_signal_histo(hWHtt_m110, lwhtt_m110)
			make_signal_histo(hZHtt_m110, lzhtt_m110)
			#MH=140
			hggHtt_m140 = hdummy.Clone("ggH140")
			hvbfHtt_m140 = hdummy.Clone("qqH140")
			hWHtt_m140 = hdummy.Clone("hWH140")
			hZHtt_m140 = hdummy.Clone("hZH140")
			make_signal_histo(hggHtt_m140, lgghtt_m140)
			make_signal_histo(hvbfHtt_m140,lvbfhtt_m140)
			make_signal_histo(hWHtt_m140, lwhtt_m140)
			make_signal_histo(hZHtt_m140, lzhtt_m140)
			
			#estimating renormalization uncertainty for ggH process. done when running the nominal datacard (systName = "")
			
			hggHtt_Up = None
			hggHtt_Down = None
			if (mode == "datacard" and systName == ""):
				
				print "calculating renormalisation uncertainty"
				renorm_weight_up = renorm_weights.get(channel+"_"+ckey,"1")
				renorm_weight_down = "(2-( "+renorm_weight_up+"))"
				hggHtt_Up = hdummy.Clone("hggHtt_Up")
				hggHtt_Down = hdummy.Clone("hggHtt_Down")
				hggHtt_Up.SetDirectory(0)
				hggHtt_Down.SetDirectory(0)
				makeDatacard_category(treeName, lgghtt, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_up, "_os", hggHtt_Up )
				makeDatacard_category(treeName, lgghtt, cvalue+lowMTcut,  "q_1 * q_2 < 0.", wvalue+"*"+renorm_weight_down, "_os", hggHtt_Down )
				
				hggHtt_Up = hdummy.Clone("hggHtt_Up")
				hggHtt_Down = hdummy.Clone("hggHtt_Down")
				make_renormscale_shape(hggHtt_Up, hggHtt_Down, lgghtt)

				hggHtt_Up_m120 = hdummy.Clone("hggHtt_Up_m120")
				hggHtt_Down_m120 = hdummy.Clone("hggHtt_Down_m120")
				make_renormscale_shape(hggHtt_Up_m120, hggHtt_Down_m120, lgghtt_m120)

				hggHtt_Up_m130 = hdummy.Clone("hggHtt_Up_m130")
				hggHtt_Down_m130 = hdummy.Clone("hggHtt_Down_m130")
				make_renormscale_shape(hggHtt_Up_m130, hggHtt_Down_m130, lgghtt_m130)

				hggHtt_Up_m110 = hdummy.Clone("hggHtt_Up_m110")
				hggHtt_Down_m110 = hdummy.Clone("hggHtt_Down_m110")
				make_renormscale_shape(hggHtt_Up_m110, hggHtt_Down_m110, lgghtt_m110)

				hggHtt_Up_m140 = hdummy.Clone("hggHtt_Up_m140")
				hggHtt_Down_m140 = hdummy.Clone("hggHtt_Down_m140")
				make_renormscale_shape(hggHtt_Up_m140, hggHtt_Down_m140, lgghtt_m140)

			

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


			#vbfHWW
			hvbfhww = hdummy.Clone("qqHWW")
			hvbfhww_ss = hdummy.Clone("qqHWW_ss")
			hvbfhww.SetDirectory(0)
			hvbfhww_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lvbfhww, cvalue+lowMTcut, wvalue, hvbfhww, hvbfhww_ss)
			if (relaxedSel is not None):
				hvbfhww_rel = hdummy.Clone("qqHWW_rel")
				hvbfhww_rel_ss = hdummy.Clone("qqHWW_rel_ss")
				makeDatacard_ssos(treeName, lvbfhww, relaxedSel+lowMTcut, wvalue, hvbfhww_rel, hvbfhww_rel_ss)
			print "qqHWW, done"


			#gghww
			hgghww = hdummy.Clone("ggHWW")
			hgghww_ss = hdummy.Clone("ggHWW_ss")
			hgghww.SetDirectory(0)
			hgghww_ss.SetDirectory(0)
			makeDatacard_ssos(treeName, lgghww, cvalue+lowMTcut, wvalue, hgghww, hgghww_ss)
			if (relaxedSel is not None):
				hgghww_rel = hdummy.Clone("ggHWW_rel")
				hgghww_rel_ss = hdummy.Clone("ggHWW_rel_ss")
				makeDatacard_ssos(treeName, lgghww, relaxedSel+lowMTcut, wvalue, hgghww_rel, hgghww_rel_ss)
			print "ggHWW, done"



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

			# weights for DY from Zmumu CR
			zmumu_weight = zmumu_CR_weights.get(ckey,None)
			dy_weight = "1.00"
			if zmumu_weight is not None:
				dy_weight += "*" + makeBinString(zmumu_weight["weights"],zmumu_weight["binVar"])

			#ZTT , ZL, ZJ
			hztt = hdummy.Clone("ztt")
			hzl = hdummy.Clone("zl")
			hzj = hdummy.Clone("zj")

			hztt.SetDirectory(0)
			hzl.SetDirectory(0)		
			hzj.SetDirectory(0)

			hztt_ss = hdummy.Clone("ztt_ss")
			hzl_ss = hdummy.Clone("zl_ss")
			hzj_ss = hdummy.Clone("zj_ss")

			hztt_ss.SetDirectory(0)
			hzl_ss.SetDirectory(0)		
			hzj_ss.SetDirectory(0)

			ZTTcut = "gen_match_2 == 5"# && q_1 * q_2 < 0." 
			ZLcut = "gen_match_2 < 5 "#&& q_1 * q_2 < 0." 
			ZJcut = "gen_match_2 == 6 "#&& q_1 * q_2 < 0." 

			makeDatacard_ssos(treeName, ldy, cvalue+lowMTcut+"&&"+ZTTcut, wvalue+"*"+dy_weight, hztt, hztt_ss)
			makeDatacard_ssos(treeName, ldy, cvalue+lowMTcut+"&&"+ZJcut, wvalue+"*"+dy_weight, hzj, hzj_ss)
			makeDatacard_ssos(treeName, ldy, cvalue+lowMTcut+"&&"+ZLcut, wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), hzl, hzl_ss)
			# correct m_vis for ZL in 0jet
			#if (("0jet" in ckey) and datacardFor2Dfit):
			#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
			#	print "#################"
			#	print " using var = ", var , "for ZL"
			# back to normal "var"
 			#var = varY+":"+varX

			hdy = hdummy.Clone("DY")
			hdy.SetDirectory(0)
			hdy.Add(hztt)
			hdy.Add(hzl)
			hdy.Add(hzj)

			hdy_ss = hdummy.Clone("DY_ss")
			hdy_ss.SetDirectory(0)
			hdy_ss.Add(hztt_ss)
			hdy_ss.Add(hzl_ss)
			hdy_ss.Add(hzj_ss)

			if (relaxedSel is not None):
				hztt_rel = hdummy.Clone("ZTT_rel")
				hzl_rel = hdummy.Clone("ZL_rel")
				hzj_rel = hdummy.Clone("ZJ_rel")	
				hztt_rel.SetDirectory(0)
				hzl_rel.SetDirectory(0)
				hzj_rel.SetDirectory(0)
				hztt_rel_ss = hdummy.Clone("ZTT_rel_ss")
				hzl_rel_ss = hdummy.Clone("ZL_rel_ss")
				hzj_rel_ss = hdummy.Clone("ZJ_rel_ss")			
				hztt_rel_ss.SetDirectory(0)
				hzl_rel_ss.SetDirectory(0)
				hzj_rel_ss.SetDirectory(0)

				makeDatacard_ssos(treeName, ldy, relaxedSel+lowMTcut+"&&"+ZTTcut, wvalue+"*"+dy_weight, hztt_rel, hztt_rel_ss)
				makeDatacard_ssos(treeName, ldy, relaxedSel+lowMTcut+"&&"+ZJcut, wvalue+"*"+dy_weight, hzj_rel, hzj_rel_ss)
				makeDatacard_ssos(treeName, ldy, relaxedSel+lowMTcut+"&&"+ZLcut, wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), hzl_rel, hzl_rel_ss)
				#if (("0jet" in ckey) and datacardFor2Dfit):
				#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
 				#var = varY+":"+varX

				hdy_rel = hdummy.Clone("DY_rel")
				hdy_rel.SetDirectory(0)
				hdy_rel.Add(hztt_rel)
				hdy_rel.Add(hzj_rel)
				hdy_rel.Add(hzl_rel)
				hdy_rel_ss = hdummy.Clone("DY_rel_ss")	
				hdy_rel_ss.SetDirectory(0)
				hdy_rel_ss.Add(hztt_rel_ss)
				hdy_rel_ss.Add(hzj_rel_ss)
				hdy_rel_ss.Add(hzl_rel_ss)

			print "DY all done"

			#Wjets
			hw = hdummy.Clone("W")
			hw_ss = hdummy.Clone("WW_ss")
			hw.SetDirectory(0)
			hw_ss.SetDirectory(0)
			print "cvalue+lowMTcut : ", cvalue+lowMTcut
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
				W_sf_up = Wsf_variations.get(ckey, {}).get("up",1)		
				W_sf_down = Wsf_variations.get(ckey, {}).get("down",1)			
				hw_wsf_up = hw.Clone("hw_wsf_up")
				hw_wsf_up.Scale(W_sf_up)
				hw_wsf_down = hw.Clone("hw_wsf_down")
				hw_wsf_down.Scale(W_sf_down)	
				hqcd_wsf_up = hdummy.Clone("hqcd_wsf_up")
				hqcd_wsf_down = hdummy.Clone("hqcd_wsf_down")
				QCD_datadriven(SSOSratio, W_sf*W_sf_up, hqcd_wsf_up)
				QCD_datadriven(SSOSratio, W_sf*W_sf_down, hqcd_wsf_down)
				#rescale to nominal yield
				hqcd_wsf_up.Scale(hqcd.Integral()/hqcd_wsf_up.Integral())
				hqcd_wsf_down.Scale(hqcd.Integral()/hqcd_wsf_down.Integral())

			if (mode=="datacard" and systName ==""):				
				SSOS_up = SSOS_variations.get(channel+"_"+ckey, {}).get("up",1)
				SSOS_down = SSOS_variations.get(channel+"_"+ckey, {}).get("down",1)
				hqcd_ssos_up = hdummy.Clone("hqcd_ssos_up")
				hqcd_ssos_down = hdummy.Clone("hqcd_ssos_down")
				#QCD_datadriven(SSOSratio*SSOS_up, W_sf, hqcd_ssos_up)
				#QCD_datadriven(SSOSratio*SSOS_down, W_sf, hqcd_ssos_down)
				QCD_datadriven(SSOS_up, W_sf, hqcd_ssos_up)
				QCD_datadriven(SSOS_down, W_sf, hqcd_ssos_down)
				hw_ssos_up = hw.Clone("hw_ssos_up")
				hw_ssos_down = hw.Clone("hw_ssos_down")
				#W_sf = Wjets_dataMC_sf(SSOSratio*SSOS_up, hw_ssos_up)
				W_sf = Wjets_dataMC_sf(SSOS_up, hw_ssos_up)
				hw_ssos_up.Scale(W_sf)
				#W_sf = Wjets_dataMC_sf(SSOSratio*SSOS_down, hw_ssos_down)
				W_sf = Wjets_dataMC_sf(SSOS_down, hw_ssos_down)
				hw_ssos_down.Scale(W_sf)


			##########################################
			# shape uncertainties for jet->tau fakes #
			##########################################
			#additional weight depending on tau pt. Propagated to W, ZJ, TTJ. Renormalise to nominal yield only for W. 
			jet_tau_fake_weights = {"up":"( ( (1-0.2*(pt_2/100))*(pt_2<=200 && gen_match_2==6))+(0.6*(pt_2>200&&gen_match_2==6))+(1*(gen_match_2 != 6)) )", 
									"down": "( ( ((1+0.2*(pt_2/100))*(pt_2<=200 && gen_match_2==6))+(1.4*(pt_2>200 && gen_match_2==6))+(1*(gen_match_2 != 6))) )"}

			if (mode=="datacard" and systName ==""):		
				# up variation		
				hw_jettaufake_up = hdummy.Clone("W_jettaufake_up")
				hzj_jettaufake_up =  hdummy.Clone("ZJ_jettaufake_up")
				httj_jettaufake_up =  hdummy.Clone("TTJ_jettaufake_up")
				weight_up = wvalue+"*"+jet_tau_fake_weights["up"]
				#calculate shift and rescale to nominal yield
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1 * q_2 < 0.", weight_up+"*"+dy_weight, "ZJ_jettaufake_up", hzj_jettaufake_up)
				#if (hzj_jettaufake_up.Integral()!=0):
				#	hzj_jettaufake_up.Scale(hzj.Integral()/hzj_jettaufake_up.Integral())

				makeDatacard_category(treeName, ltt, cvalue+lowMTcut, "gen_match_2 != 5 && q_1*q_2 < 0", weight_up, "TTJ_jettaufake_up", httj_jettaufake_up)
				#if (httj_jettaufake_up.Integral()!=0):
				#	httj_jettaufake_up.Scale(httj.Integral()/httj_jettaufake_up.Integral())
				if (relaxedSel is not None):
					makeDatacard_category(treeName, lw, relaxedSel+lowMTcut,"q_1*q_2 < 0", weight_up, "W_jettaufake_up", hw_jettaufake_up)
				else:
					makeDatacard_category(treeName, lw, cvalue+lowMTcut, "q_1*q_2 < 0", weight_up, "W_jettaufake_up", hw_jettaufake_up)
				if ( hw_jettaufake_up.Integral() != 0 ):
					hw_jettaufake_up.Scale(hw_corrected.Integral()/hw_jettaufake_up.Integral())				

				# down variation				
				hw_jettaufake_down =  hdummy.Clone("W_jettaufake_down")
				hzj_jettaufake_down =  hdummy.Clone("ZJ_jettaufake_down")
				httj_jettaufake_down =  hdummy.Clone("TTJ_jettaufake_down")
				weight_down = wvalue+"*"+jet_tau_fake_weights["down"]
				#calculate shift and rescale to nominal yield
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1 * q_2 < 0.", weight_down+"*"+dy_weight, "ZJ_jettaufake_down", hzj_jettaufake_down)
				#if (hzj_jettaufake_down.Integral()!=0):
				#	hzj_jettaufake_down.Scale(hzj.Integral()/hzj_jettaufake_down.Integral())

				makeDatacard_category(treeName, ltt, cvalue+lowMTcut, "gen_match_2 != 5 && q_1*q_2 < 0", weight_down, "TTJ_jettaufake_down", httj_jettaufake_down)
				#if (httj_jettaufake_down.Integral()!=0):
				#	httj_jettaufake_down.Scale(httj.Integral()/httj_jettaufake_down.Integral())

				if (relaxedSel is not None):
					makeDatacard_category(treeName, lw, relaxedSel+lowMTcut,"q_1*q_2 < 0", weight_down, "W_jettaufake_down", hw_jettaufake_down)
				else:
					makeDatacard_category(treeName, lw, cvalue+lowMTcut, "q_1*q_2 < 0", weight_down, "W_jettaufake_down", hw_jettaufake_down)

				if (hw_jettaufake_down.Integral()!=0):
					hw_jettaufake_down.Scale(hw_corrected.Integral()/hw_jettaufake_down.Integral())			
	

			##########################################
			# shape uncertainties for DY weights     #
			##########################################
			
			if (mode=="datacard" and systName =="" and ("vbf" in ckey) ):
				dy_weight_up = "1.0"
				hztt_dyweight_up = hdummy.Clone("hztt_dyweight_up")
				hzl_dyweight_up = hdummy.Clone("hzl_dyweight_up")
				hzj_dyweight_up = hdummy.Clone("hzj_dyweight_up")
				hewkz_dyweight_up = hdummy.Clone("hewkz_dyweight_up") 
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_up, "ZTT_dyweight_up", hztt_dyweight_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_up+"*"+ZL_yield_SF.get(channel,"1"), "ZL_dyweight_up", hzl_dyweight_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_up, "ZJ_dyweight_up", hzj_dyweight_up)
				makeDatacard_category(treeName, lewkz, cvalue+lowMTcut, "q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_up, "EWKZ_dyweight_up", hewkz_dyweight_up)


				dy_weight_down_dict = zmumu_CR_weights.get("down_shape",None)
				dy_weight_down = makeBinString(dy_weight_down_dict["weights"],dy_weight_down_dict["binVar"])
				hztt_dyweight_down = hdummy.Clone("hztt_dyweight_down")
				hzl_dyweight_down = hdummy.Clone("hzl_dyweight_down")
				hzj_dyweight_down = hdummy.Clone("hzj_dyweight_down")
				hewkz_dyweight_down = hdummy.Clone("hewkz_dyweight_down")
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_down, "ZTT_dyweight_down", hztt_dyweight_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_down+"*"+ZL_yield_SF.get(channel,"1"), "ZL_dyweight_down", hzl_dyweight_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_down, "ZJ_dyweight_down", hzj_dyweight_down)
				makeDatacard_category(treeName, lewkz, cvalue+lowMTcut,"q_1*q_2 < 0", wvalue+"*"+dy_weight+"*"+dy_weight_down, "EWKZ_dyweight_down", hewkz_dyweight_down)

			


			##########################################
			# shape uncertainties for tau DM         #
			##########################################
			
			if (mode=="datacard" and systName ==""):
				#apply 3% shift for ZTT, to allow for DM migration. Normalise to ZTT yield (per DM)
				hztt_taudm_1prong_up = hdummy.Clone("ZTT_tauDM_1prong_up")
				hztt_taudm_3prong_up = hdummy.Clone("ZTT_tauDM_3prong_up")
				hztt_taudm_1prong1pi_up = hdummy.Clone("ZTT_tauDM_1prong1pi_up")

				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1+(0.03*(tau_decay_mode_2==0)))", "ZTT_tauDM_1prong_up", hztt_taudm_1prong_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1+(0.03*(tau_decay_mode_2==10)))",  "ZTT_tauDM_3prong_up", hztt_taudm_3prong_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1+(0.03*(tau_decay_mode_2==1)))",  "ZTT_tauDM_1prong1pi_up", hztt_taudm_1prong1pi_up)
			
				hztt_taudm_1prong_down = hdummy.Clone("ZTT_tauDM_1prong_down")
				hztt_taudm_3prong_down = hdummy.Clone("ZTT_tauDM_3prong_down")
				hztt_taudm_1prong1pi_down = hdummy.Clone("ZTT_tauDM_1prong1pi_down")

				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1-(0.03*(tau_decay_mode_2==0)))", "ZTT_tauDM_1prong_down", hztt_taudm_1prong_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1-(0.03*(tau_decay_mode_2==10)))",  "ZTT_tauDM_3prong_down", hztt_taudm_3prong_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+"(1-(0.03*(tau_decay_mode_2==1)))",  "ZTT_tauDM_1prong1pi_down", hztt_taudm_1prong1pi_down)	

				#get the nominal ZTT yield per DM and renormalise. 
				hztt_taudm_1prong_up.Scale(hztt.Integral()/hztt_taudm_1prong_up.Integral())
				hztt_taudm_1prong_down.Scale(hztt.Integral()/hztt_taudm_1prong_down.Integral())
				hztt_taudm_3prong_up.Scale(hztt.Integral()/hztt_taudm_3prong_up.Integral())
				hztt_taudm_3prong_down.Scale(hztt.Integral()/hztt_taudm_3prong_down.Integral())
				hztt_taudm_1prong1pi_up.Scale(hztt.Integral()/hztt_taudm_1prong1pi_up.Integral())
				hztt_taudm_1prong1pi_down.Scale(hztt.Integral()/hztt_taudm_1prong1pi_down.Integral()) 


			########################################################
			# shape uncertainties for lepton-> tau fakes in 0 jet  #
			########################################################
			'''
			if (mode=="datacard" and systName =="" and ("0jet" in ckey)):
				#apply 80% shift for ZL, 1 prong only. Normalised to ZL yield. 
				hzl_leptaufake_up = hdummy.Clone("ZL_leptaufake_up")
				hzl_leptaufake_down = hdummy.Clone("ZL_leptaufake_down")
				var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1")+"*"+"(1-(0.8*(tau_decay_mode_2==0)))", "ZL_leptaufake_down", hzl_leptaufake_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1")+"*"+"(1+(0.8*(tau_decay_mode_2==0)))", "ZL_leptaufake_up", hzl_leptaufake_up)

				hzl_leptaufake_up.Scale(hzl.Integral()/hzl_leptaufake_up.Integral())
				hzl_leptaufake_down.Scale(hzl.Integral()/hzl_leptaufake_down.Integral())
				var = varY+":"+varX
			'''

			#uncertainty on the ZL scale factor per DM
			if (mode=="datacard" and systName ==""):
				#if (("0jet") in ckey):
				#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX

				hzl_sf_1prong_up = hdummy.Clone("ZL_SF_1prong_up")
				hzl_sf_1prong_down = hdummy.Clone("ZL_SF_1prong_down")

				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel+"_1prong_up","1"), "ZL_SF_1prong_up", hzl_sf_1prong_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel+"_1prong_down","1"), "ZL_SF_1prong_down", hzl_sf_1prong_down)								

				hzl_sf_1prong1pi_up = hdummy.Clone("ZL_SF_1prong1pi_up")
				hzl_sf_1prong1pi_down = hdummy.Clone("ZL_SF_1prong1pi_down")
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel+"_1prong1pi_up","1"), "ZL_SF_1prong1pi_up", hzl_sf_1prong1pi_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel+"_1prong1pi_down","1"), "ZL_SF_1prong1pi_down", hzl_sf_1prong1pi_down)								


			########################################################
			# shape uncertainties for ZL m_vis shape               #
			########################################################
			'''				
			if (mode=="datacard" and systName =="" and ("0jet" in ckey) and datacardFor2Dfit):
				# shift the visible mass for ZL, for 1 prong and 1prong 1pizero.
				hzl_mvis_1prong_up = hdummy.Clone("ZL_mvis_1prong_up")
				hzl_mvis_1prong_down = hdummy.Clone("ZL_mvis_1prong_down")
				# 1 prong up
				#mvis_shift_1prong_up = "(1+(0.03*(tau_decay_mode_2==0)))"
				#var = varY+"*"+mvis_zl_corr.get(channel,"1")+"*"+mvis_shift_1prong_up+":"+varX

				var = varY+"*"+mvis_zl_corr.get(channel+"_1prong_up","1")+":"+varX
				print " 1 prong up --- var: ", var
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mvis_1prong_up", hzl_mvis_1prong_up)
				# 1 prong down
				#mvis_shift_1prong_down = "(1-(0.03*(tau_decay_mode_2==0)))"
				#var = varY+"*"+mvis_zl_corr.get(channel,"1")+"*"+mvis_shift_1prong_down+":"+varX

				var = varY+"*"+mvis_zl_corr.get(channel+"_1prong_down","1")+":"+varX
				print " 1 prong down --- var: ", var
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mvis_1prong_down", hzl_mvis_1prong_down)

				hzl_mvis_1prong1pi_up = hdummy.Clone("ZL_mvis_1prong1pi_up")
				hzl_mvis_1prong1pi_down = hdummy.Clone("ZL_mvis_1prong1pi_down")
				# 1prong1pi up
				#mvis_shift_1prong1pi_up = "(1+(0.03*(tau_decay_mode_2==1)))"
				#var = varY+"*"+mvis_zl_corr.get(channel,"1")+"*"+mvis_shift_1prong1pi_up+":"+varX

				var = varY+"*"+mvis_zl_corr.get(channel+"_1prong1pi_up","1")+":"+varX
				print " 1 prong 1pi up --- var: ", var
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mvis_1prong1pi_up", hzl_mvis_1prong1pi_up)
				#1prong1pi down
				#mvis_shift_1prong1pi_down = "(1-(0.03*(tau_decay_mode_2==1)))"
				#var = varY+"*"+mvis_zl_corr.get(channel,"1")+"*"+mvis_shift_1prong1pi_down+":"+varX

				var = varY+"*"+mvis_zl_corr.get(channel+"_1prong1pi_down","1")+":"+varX
				print " 1 prong 1pi down --- var: ", var
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mvis_1prong1pi_down", hzl_mvis_1prong1pi_down)
			
				var = varY + ":" + varX
			'''
			########################################################
			# shape uncertainties for ZL mass shape                #
			########################################################
			if (mode=="datacard" and systName ==""):
				hzl_mass_1prong_up = hdummy.Clone("ZL_mass_1prong_up")
				hzl_mass_1prong_down = hdummy.Clone("ZL_mass_1prong_down")
				hzl_mass_1prong_up.SetDirectory(0)
				hzl_mass_1prong_down.SetDirectory(0)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_1prong_13TeVUp", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_1prong_up", hzl_mass_1prong_up)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_1prong_13TeVDown", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_1prong_down", hzl_mass_1prong_down)

				hzl_mass_1prong1pi_up = hdummy.Clone("ZL_mass_1prong1pi_up")
				hzl_mass_1prong1pi_down = hdummy.Clone("ZL_mass_1prong1pi_down")
				hzl_mass_1prong1pi_up.SetDirectory(0)
				hzl_mass_1prong1pi_down.SetDirectory(0)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_1prong1pi_13TeVUp", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_1prong1pi_up", hzl_mass_1prong1pi_up)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_1prong1pi_13TeVDown", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_1prong1pi_down", hzl_mass_1prong1pi_down)

				hzl_mass_3prong_up = hdummy.Clone("ZL_mass_3prong_up")
				hzl_mass_3prong_down = hdummy.Clone("ZL_mass_3prong_down")
				hzl_mass_3prong_up.SetDirectory(0)
				hzl_mass_3prong_down.SetDirectory(0)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_3prong_13TeVUp", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_3prong_up", hzl_mass_3prong_up)
				makeDatacard_category(treeName+"_CMS_htt_ZLShape_3prong_13TeVDown", ldy, cvalue+lowMTcut, ZLcut+" && q_1 * q_2 < 0.", wvalue+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_mass_3prong_down", hzl_mass_3prong_down)

			##########################################
			# shape uncertainties for Z pt weights   #
			##########################################

			if (mode =="datacard" and systName ==""):

				zptweight_up  = "(1+1.10*(zptweight-1))"
				zptweight_down = "(1+0.90*(zptweight-1))"
				weight_up = wvalue.replace("zptweight", zptweight_up)
				weight_down = wvalue.replace("zptweight",zptweight_down)

				hztt_zptweight_up = hdummy.Clone("ZTT_zptweight_up")
				hzj_zptweight_up = hdummy.Clone("ZJ_zptweight_up")
				hzl_zptweight_up = hdummy.Clone("ZL_zptweight_up")
				hewkz_zptweight_up = hdummy.Clone("EWKZ_zptweight_up")

				hztt_zptweight_down = hdummy.Clone("ZTT_zptweight_down")
				hzj_zptweight_down = hdummy.Clone("ZJ_zptweight_down")
				hzl_zptweight_down = hdummy.Clone("ZL_zptweight_down")
				hewkz_zptweight_down = hdummy.Clone("EWKZ_zptweight_down")		

				#print "################"
				#print " Zpt weight uncertainty "
				#print " weight_up : " , weight_up
				#print " weight_down : ", weight_down
				makeDatacard_category(treeName, lewkz, cvalue+lowMTcut, "q_1*q_2 < 0", weight_down+"*"+dy_weight, "EWKZ_zptweight_down", hewkz_zptweight_down)
				makeDatacard_category(treeName, lewkz, cvalue+lowMTcut, "q_1*q_2 < 0", weight_up+"*"+dy_weight, "EWKZ_zptweight_up", hewkz_zptweight_up)	

				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1*q_2 < 0",weight_down+"*"+dy_weight, "ZTT_zptweight_down", hztt_zptweight_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZTTcut+" && q_1*q_2 < 0",weight_up+"*"+dy_weight, "ZTT_zptweight_up", hztt_zptweight_up)

				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1*q_2 < 0", weight_down+"*"+dy_weight, "ZJ_zptweight_down", hzj_zptweight_down)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZJcut+" && q_1*q_2 < 0", weight_up+"*"+dy_weight, "ZJ_zptweight_up", hzj_zptweight_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1*q_2 < 0", weight_up+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_zptweight_up", hzl_zptweight_up)
				makeDatacard_category(treeName, ldy, cvalue+lowMTcut, ZLcut+" && q_1*q_2 < 0", weight_down+"*"+dy_weight+"*"+ZL_yield_SF.get(channel,"1"), "ZL_zptweight_down", hzl_zptweight_down)
				#if ("0jet" in ckey):
				#	var = varY+"*"+mvis_zl_corr.get(channel,"1")+":"+varX
				



			##########################
			#writing to output file
			##########################

			
			if mode == "datacard":
				#outfile=TFile( wkey+"_htt_"+channel+".inputs-sm-13TeV-"+hkey+systName+".root", "update")
				outfile=TFile(filename,  "update")
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
				
				hggHtt_m120 = unroll(hggHtt_m120)
				hvbfHtt_m120 = unroll(hvbfHtt_m120)
				hWHtt_m120 = unroll(hWHtt_m120)
				hZHtt_m120 = unroll(hZHtt_m120)
				hggHtt_m130 = unroll(hggHtt_m130)
				hvbfHtt_m130 = unroll(hvbfHtt_m130)
				hWHtt_m130 = unroll(hWHtt_m130)
				hZHtt_m130 = unroll(hZHtt_m130)
				hggHtt_m110 = unroll(hggHtt_m110)
				hvbfHtt_m110 = unroll(hvbfHtt_m110)
				hWHtt_m110 = unroll(hWHtt_m110)
				hZHtt_m110 = unroll(hZHtt_m110)
				hggHtt_m140 = unroll(hggHtt_m140)
				hvbfHtt_m140 = unroll(hvbfHtt_m140)
				hWHtt_m140 = unroll(hWHtt_m140)
				hZHtt_m140 = unroll(hZHtt_m140)
				
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
				hvbfhww = unroll(hvbfhww)
				hgghww = unroll(hgghww)
				
				if (mode=="datacard" and systName == ""):
					hggHtt_Up = unroll(hggHtt_Up)
					hggHtt_Down = unroll(hggHtt_Down)
					hggHtt_Up_m120 = unroll(hggHtt_Up_m120)
					hggHtt_Down_m120 = unroll(hggHtt_Down_m120)
					hggHtt_Up_m130 = unroll(hggHtt_Up_m130)
					hggHtt_Down_m130 = unroll(hggHtt_Down_m130)
					hggHtt_Up_m110 = unroll(hggHtt_Up_m110)
					hggHtt_Down_m110 = unroll(hggHtt_Down_m110)
					hggHtt_Up_m140 = unroll(hggHtt_Up_m140)
					hggHtt_Down_m140 = unroll(hggHtt_Down_m140)
					hw_wsf_up = unroll(hw_wsf_up)
					hw_wsf_down = unroll(hw_wsf_down)
					hqcd_wsf_up = unroll(hqcd_wsf_up)
					hqcd_wsf_down = unroll(hqcd_wsf_down)
					hqcd_ssos_up = unroll(hqcd_ssos_up)
					hqcd_ssos_down = unroll(hqcd_ssos_down)
					hw_ssos_up = unroll(hw_ssos_up)
					hw_ssos_down = unroll(hw_ssos_down)
					hw_jettaufake_up = unroll(hw_jettaufake_up)
					hzj_jettaufake_up =  unroll(hzj_jettaufake_up)
					httj_jettaufake_up =  unroll(httj_jettaufake_up)
					hw_jettaufake_down = unroll(hw_jettaufake_down)
					hzj_jettaufake_down =  unroll(hzj_jettaufake_down)
					httj_jettaufake_down =  unroll(httj_jettaufake_down)

					hztt_zptweight_up = unroll(hztt_zptweight_up)
					hzj_zptweight_up = unroll(hzj_zptweight_up)
					hzl_zptweight_up = unroll(hzl_zptweight_up)
					hewkz_zptweight_up = unroll(hewkz_zptweight_up)
					hztt_zptweight_down = unroll(hztt_zptweight_down)
					hzj_zptweight_down = unroll(hzj_zptweight_down)
					hzl_zptweight_down = unroll(hzl_zptweight_down)
					hewkz_zptweight_down = unroll(hewkz_zptweight_down)

					hztt_taudm_1prong_up = unroll(hztt_taudm_1prong_up)
					hztt_taudm_1prong_down = unroll(hztt_taudm_1prong_down)
					hztt_taudm_3prong_up = unroll(hztt_taudm_3prong_up)
					hztt_taudm_3prong_down = unroll(hztt_taudm_3prong_down)
					hztt_taudm_1prong1pi_up = unroll(hztt_taudm_1prong1pi_up)
					hztt_taudm_1prong1pi_down = unroll(hztt_taudm_1prong1pi_down)

					hzl_sf_1prong_up = unroll(hzl_sf_1prong_up)
					hzl_sf_1prong_down = unroll(hzl_sf_1prong_down)
					hzl_sf_1prong1pi_up = unroll(hzl_sf_1prong1pi_up)
					hzl_sf_1prong1pi_down = unroll(hzl_sf_1prong1pi_down)

					hzl_mass_1prong_up = unroll(hzl_mass_1prong_up)	
					hzl_mass_1prong_down = unroll(hzl_mass_1prong_down)			
					hzl_mass_1prong1pi_up = unroll(hzl_mass_1prong1pi_up)			
					hzl_mass_1prong1pi_down = unroll(hzl_mass_1prong1pi_down)	
					hzl_mass_3prong_up = unroll(hzl_mass_3prong_up)	
					hzl_mass_3prong_down = unroll(hzl_mass_3prong_down)

					#if (("0jet" in ckey)):
					#	hzl_leptaufake_up = unroll(hzl_leptaufake_up)
					#	hzl_leptaufake_down = unroll(hzl_leptaufake_down)	
					#	hzl_mvis_1prong_up = unroll(hzl_mvis_1prong_up)			
					#	hzl_mvis_1prong_down = unroll(hzl_mvis_1prong_down)			
					#	hzl_mvis_1prong1pi_up = unroll(hzl_mvis_1prong1pi_up)			
					#	hzl_mvis_1prong1pi_down = unroll(hzl_mvis_1prong1pi_down)			
					
					if (("vbf" in ckey)):
						hztt_dyweight_up = unroll(hztt_dyweight_up)
						hzl_dyweight_up = unroll(hzl_dyweight_up)
						hzj_dyweight_up = unroll(hzj_dyweight_up)
						hewkz_dyweight_up = unroll(hewkz_dyweight_up)
						hztt_dyweight_down = unroll(hztt_dyweight_down)
						hzl_dyweight_down = unroll(hzl_dyweight_down)
						hzj_dyweight_down = unroll(hzj_dyweight_down)
						hewkz_dyweight_down = unroll(hewkz_dyweight_down)

			
			print "writing to out file..."
			if (mode == "histos"):
				hdata.Write("data_obs")

			#write data histogram only once, in the file w/o syst variations
			if (mode == "datacard" and systName == ""): 
				hdata.Write("data_obs")
				
				hggHtt_Up.Write("ggH125_CMS_scale_gg_13TeVUp")
				hggHtt_Down.Write("ggH125_CMS_scale_gg_13TeVDown")
				hggHtt_Up_m120.Write("ggH120_CMS_scale_gg_13TeVUp")
				hggHtt_Down_m120.Write("ggH120_CMS_scale_gg_13TeVDown")
				hggHtt_Up_m130.Write("ggH130_CMS_scale_gg_13TeVUp")
				hggHtt_Down_m130.Write("ggH130_CMS_scale_gg_13TeVDown")
				hggHtt_Up_m110.Write("ggH110_CMS_scale_gg_13TeVUp")
				hggHtt_Down_m110.Write("ggH110_CMS_scale_gg_13TeVDown")
				hggHtt_Up_m140.Write("ggH140_CMS_scale_gg_13TeVUp")
				hggHtt_Down_m140.Write("ggH140_CMS_scale_gg_13TeVDown")
				hw_wsf_up.Write("W_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hw_wsf_down.Write("W_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hqcd_wsf_up.Write("QCD_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hqcd_wsf_down.Write("QCD_WSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hqcd_ssos_up.Write("QCD_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hqcd_ssos_down.Write("QCD_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hw_ssos_up.Write("W_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Up")
				hw_ssos_down.Write("W_QCDSFUncert_"+channel+"_"+ckey+"_13TeV"+"Down")
				hw_jettaufake_up.Write("W_CMS_htt_jetToTauFake_13TeVUp")
				hzj_jettaufake_up.Write("ZJ_CMS_htt_jetToTauFake_13TeVUp")
				httj_jettaufake_up.Write("TTJ_CMS_htt_jetToTauFake_13TeVUp")
				hw_jettaufake_down.Write("W_CMS_htt_jetToTauFake_13TeVDown")
				hzj_jettaufake_down.Write("ZJ_CMS_htt_jetToTauFake_13TeVDown")
				httj_jettaufake_down.Write("TTJ_CMS_htt_jetToTauFake_13TeVDown")
				hztt_zptweight_up.Write("ZTT_CMS_htt_dyShape_13TeVUp")
				hzj_zptweight_up.Write("ZJ_CMS_htt_dyShape_13TeVUp")
				hzl_zptweight_up.Write("ZL_CMS_htt_dyShape_13TeVUp")
				hewkz_zptweight_up.Write("EWKZ_CMS_htt_dyShape_13TeVUp")
				hztt_zptweight_down.Write("ZTT_CMS_htt_dyShape_13TeVDown")
				hzj_zptweight_down.Write("ZJ_CMS_htt_dyShape_13TeVDown")
				hzl_zptweight_down.Write("ZL_CMS_htt_dyShape_13TeVDown")
				hewkz_zptweight_down.Write("EWKZ_CMS_htt_dyShape_13TeVDown")

				hztt_taudm_1prong_up.Write("ZTT_CMS_tauDMReco_1prong_13TeVUp")
				hztt_taudm_1prong_down.Write("ZTT_CMS_tauDMReco_1prong_13TeVDown")
				hztt_taudm_3prong_up.Write("ZTT_CMS_tauDMReco_3prong_13TeVUp")
				hztt_taudm_3prong_down.Write("ZTT_CMS_tauDMReco_3prong_13TeVDown")
				hztt_taudm_1prong1pi_up.Write("ZTT_CMS_tauDMReco_1prong1pizero_13TeVUp")
				hztt_taudm_1prong1pi_down.Write("ZTT_CMS_tauDMReco_1prong1pizero_13TeVDown")

				hzl_mass_1prong_up.Write("ZL_CMS_ZLShape_"+channel+"_1prong_13TeVUp")
				hzl_mass_1prong_down.Write("ZL_CMS_ZLShape_"+channel+"_1prong_13TeVDown")		
				hzl_mass_1prong1pi_up.Write("ZL_CMS_ZLShape_"+channel+"_1prong1pizero_13TeVUp")		
				hzl_mass_1prong1pi_down.Write("ZL_CMS_ZLShape_"+channel+"_1prong1pizero_13TeVDown")
				hzl_mass_3prong_up.Write("ZL_CMS_ZLShape_"+channel+"_3prong_13TeVUp")
				hzl_mass_3prong_down.Write("ZL_CMS_ZLShape_"+channel+"_3prong_13TeVDown")


				if (channel =="mt"):
					hzl_sf_1prong_up.Write("ZL_CMS_mFakeTau_1prong_13TeVUp")
					hzl_sf_1prong_down.Write("ZL_CMS_mFakeTau_1prong_13TeVDown")
					hzl_sf_1prong1pi_up.Write("ZL_CMS_mFakeTau_1prong1pizero_13TeVUp")
					hzl_sf_1prong1pi_down.Write("ZL_CMS_mFakeTau_1prong1pizero_13TeVDown")

				if (channel =="et"):
					hzl_sf_1prong_up.Write("ZL_CMS_eFakeTau_1prong_13TeVUp")
					hzl_sf_1prong_down.Write("ZL_CMS_eFakeTau_1prong_13TeVDown")
					hzl_sf_1prong1pi_up.Write("ZL_CMS_eFakeTau_1prong1pizero_13TeVUp")
					hzl_sf_1prong1pi_down.Write("ZL_CMS_eFakeTau_1prong1pizero_13TeVDown")

				#if (channel == "et" and ("0jet" in ckey)):
				#	hzl_leptaufake_up.Write("ZL_CMS_eFakeTau_0jet_tauDMReco_13TeVUp")
				#	hzl_leptaufake_down.Write("ZL_CMS_eFakeTau_0jet_tauDMReco_13TeVDown")
				#	hzl_mvis_1prong_up.Write("ZL_CMS_ZLShape_et_0jet_1prong_13TeVUp")
				#	hzl_mvis_1prong_down.Write("ZL_CMS_ZLShape_et_0jet_1prong_13TeVDown")
				#	hzl_mvis_1prong1pi_up.Write("ZL_CMS_ZLShape_et_0jet_1prong1pizero_13TeVUp")
				#	hzl_mvis_1prong1pi_down.Write("ZL_CMS_ZLShape_et_0jet_1prong1pizero_13TeVDown")

				#if (channel == "mt" and ("0jet" in ckey)):
				#	hzl_leptaufake_up.Write("ZL_CMS_mFakeTau_0jet_tauDMReco_13TeVUp")
				#	hzl_leptaufake_down.Write("ZL_CMS_mFakeTau_0jet_tauDMReco_13TeVDown")			
				#	hzl_mvis_1prong_up.Write("ZL_CMS_ZLShape_mt_0jet_1prong_13TeVUp")
				#	hzl_mvis_1prong_down.Write("ZL_CMS_ZLShape_mt_0jet_1prong_13TeVDown")
				#	hzl_mvis_1prong1pi_up.Write("ZL_CMS_ZLShape_mt_0jet_1prong1pizero_13TeVUp")
				#	hzl_mvis_1prong1pi_down.Write("ZL_CMS_ZLShape_mt_0jet_1prong1pizero_13TeVDown")	

				if (("vbf" in ckey)):
					hztt_dyweight_up.Write("ZTT_CMS_htt_zmumuShape_"+"VBF"+"_13TeVUp")
					hzl_dyweight_up.Write("ZL_CMS_htt_zmumuShape_"+"VBF"+"_13TeVUp")
					hzj_dyweight_up.Write("ZJ_CMS_htt_zmumuShape_"+"VBF"+"_13TeVUp")
					hewkz_dyweight_up.Write("EWKZ_CMS_htt_zmumuShape_"+"VBF"+"_13TeVUp")
					hztt_dyweight_down.Write("ZTT_CMS_htt_zmumuShape_"+"VBF"+"_13TeVDown")
					hzl_dyweight_down.Write("ZL_CMS_htt_zmumuShape_"+"VBF"+"_13TeVDown")
					hzj_dyweight_down.Write("ZJ_CMS_htt_zmumuShape_"+"VBF"+"_13TeVDown")
					hewkz_dyweight_down.Write("EWKZ_CMS_htt_zmumuShape_"+"VBF"+"_13TeVDown")

			#dictionary to convert from naming in synch ntuples to names in datacards
			#add only the names that need to be changed
			syst_naming_convention = {
            	"_topPtWeightUp":"_CMS_htt_ttbarShape_13TeVUp", "_topPtWeightDown":"_CMS_htt_ttbarShape_13TeVDown",
				"_CMS_shape_t_13TeVUp":"_CMS_shape_t_"+channel+"_13TeVUp", "_CMS_shape_t_13TeVDown":"_CMS_shape_t_"+channel+"_13TeVDown", 
            	"_CMS_htt_ZLShape_13TeVUp":"_CMS_htt_ZLShape_"+channel+"_13TeVUp", "_CMS_htt_ZLShape_13TeVDown":"_CMS_htt_ZLShape_"+channel+"_13TeVDown",
				"_CMS_scale_j_Total13TeVUp":"_CMS_scale_j_13TeVUp", "_CMS_scale_j_Total13TeVDown":"_CMS_scale_j_13TeVDown",
				"_CMS_shape_t_1prong1pi0_13TeVUp" : "_CMS_scale_t_1prong1pizero_13TeVUp", "_CMS_shape_t_1prong1pi0_13TeVDown" : "_CMS_scale_t_1prong1pizero_13TeVDown",
				"_CMS_shape_t_1prong_13TeVUp" : "_CMS_scale_t_1prong_13TeVUp", "_CMS_shape_t_1prong_13TeVDown" : "_CMS_scale_t_1prong_13TeVDown",
				"_CMS_shape_t_3prong_13TeVUp" : "_CMS_scale_t_3prong_13TeVUp", "_CMS_shape_t_3prong_13TeVDown" : "_CMS_scale_t_3prong_13TeVDown",
				#"_CMS_shape_dyShape_13TeVUp" : "_CMS_htt_dyShape_13TeVUp", "_CMS_shape_dyShape_13TeVDown" : "_CMS_htt_dyShape_13TeVDown"
			}

			systName_out  = syst_naming_convention.get(systName, systName)

			hggHtt.Write("ggH125"+systName_out)
			hvbfHtt.Write("qqH125"+systName_out)
			hWHtt.Write("WH125"+systName_out)
			hZHtt.Write("ZH125"+systName_out)
			
			hggHtt_m120.Write("ggH120"+systName_out)
			hvbfHtt_m120.Write("qqH120"+systName_out)
			hWHtt_m120.Write("WH120"+systName_out)
			hZHtt_m120.Write("ZH120"+systName_out)

			hggHtt_m130.Write("ggH130"+systName_out)
			hvbfHtt_m130.Write("qqH130"+systName_out)
			hWHtt_m130.Write("WH130"+systName_out)
			hZHtt_m130.Write("ZH130"+systName_out)

			hggHtt_m110.Write("ggH110"+systName_out)
			hvbfHtt_m110.Write("qqH110"+systName_out)
			hWHtt_m110.Write("WH110"+systName_out)
			hZHtt_m110.Write("ZH110"+systName_out)

			hggHtt_m140.Write("ggH140"+systName_out)
			hvbfHtt_m140.Write("qqH140"+systName_out)
			hWHtt_m140.Write("WH140"+systName_out)
			hZHtt_m140.Write("ZH140"+systName_out)

			hztt.Write("ZTT"+systName_out)
			hzl.Write("ZL"+systName_out)
			hzj.Write("ZJ"+systName_out)
			httt.Write("TTT"+systName_out)
			httj.Write("TTJ"+systName_out)
			hw.Write("W_noWSF"+systName_out)
			hw_corrected.Write("W"+systName_out)
			hewkz.Write("EWKZ"+systName_out)
			hvv.Write("VV"+systName_out)
			hgghww.Write("ggH_WW125"+systName_out)
			hvbfhww.Write("qqH_WW125"+systName_out)
			hqcd.Write("QCD"+systName_out)
			print "... done"

			if mode == "datacard":
				outfile.Close()

	if mode == "histos":
		outfile.Close()


