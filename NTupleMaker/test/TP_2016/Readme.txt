						TAG&PROBE AND EFFICIENCIES EXTRACTION PACKAGE 

Author: 	Alberto Bragagnolo 		alberto.bragagnolo.3@studenti.unipd.it

Description: this package can be used to compute efficiencies (identification&isolation
			 and trigger) and scale factor, both for electrons and muons, using the 
			 Tag&Probe technique on the Z peak.


						COMPONENTS

MAIN COMPONENTS

	1. bin/TagAndProbe2016_mumu.cpp, bin/TagAndProbe2016_ee.cpp			[T&P otree producer source codes]
	2.1 TP_eff_mu.C, TP_eff_e.C 										[efficiencies extraction macros for trigger]
	2.2 SF_mu.C, SF_e.C 												[efficiencies and scale facrorss extraction macros for IdIso]
	3. FitPassAndFail.C 												[fitting tool, Author: Alexei Raspereza 	rasp@mail.desy.de]

AUXILIARIES COMPONENTS

	a.1 interface/TagProbeTree.h	[.h that defines the otree structure]
	a.2 src/TagProbeTree.cc			[.cc that defines the otree structure]

CONFIGURATION FILES

	b.1.1 TagAndProbe_mu.conf 		[config file for muons]
	b.1.2 TagAndProbe_mu_MC.conf 	[config file for muons - MC]
	b.2.1 TagAndProbe_e.conf 		[config file for electrons]
	b.2.2 TagAndProbe_e_MC.conf 	[config file for electrons - MC]					

EXTRA
	c.1 plot_eff.C 	[plotting tool])



						DESCRIPTION

1. bin/TagAndProbe2016_mumu.cpp
   bin/TagAndProbe2016_ee.cpp

    -These codes are in NTupleMaker/bin and they are already committed on github
	-These codes have to be compiled within the CMSSW enviroment
	-These codes use a.1, a.2 and a lot of other CMSSW components
	-One code is for electrons T&P, the other for muons
	-This code is the slow part of the process 
	-The ID criterias are hard coded in the two codes. 
	-The trigger path to be evaluated can be changed (under hlt_1, hlt_2, ...., hlt_10)


	These codes run on events ntuples (data or MC), find the T&P pairs, evaluate id and trigger criterias,
	compute isolation and produce the root files with the T&P otree.

	The otree structure is defined in a.1 and a.2 (committed in github), 
	this structure can be changed at will.

	In the otree the following variables will be stored (every entry is a different T&P pair):
		-kinematics variables for tag and probe
		-flag for passing ID criterias
		-flag for passing trigger (up to 6 different hlt filters)
		-isolation of the probe
		-pu_weight and mc_weight for MC sample

	Arguments:
		-config file for analysis (b.x.x serves as examples) - mandatory
		-file list to be analized - mandatory
		-index of the first file to be processed - optional
		-index of the last file to be processed - optional



	More informations in the code for muons (the one for electrons is nearly identical and is
	not commented).

	Example of use: 
		TagAndProbe2016_mumu TagAndProbe_mu.conf SingleMuon_Run2016B_PromptReco_v2_list (run on naf, not raccomended)
		

2.1 TP_eff_mu.C
   	TP_eff_e.C 	

	-These codes are root macros and can be execute on the fly
	-These codes use FitPassAndFail.C and HttStylesNew.cc (one has to make sure that the path in the #include are correct)

	These codes use the root files produced by 1. and compute the efficiencies fitting the Z peak with the
	plotting tool. In these codes the eta and pt bins are defined. Root files with scale factor (for IdIso) or
	efficiencies (for triggers) will be produces. The plots with the various fits will be stored in dedicated
	directories.

	Arguments:
		-root file name to be used (without the extension)
		-name of the efficiency to be extracted ("IdIso" or "hlt_1", "hlt_2" ....)
		-Isolation cut to be evaluated/used (normally 0.1 or 0.15)

	More informations in the code for muons (the one for electrons is nearly identical and is
	not commented).

	Example of use: 
		root -l -b -q ' TP_eff_mu.C("SingleMuon_Run2016B_TP", "IdIso", 0.1)'

2.2 SF_mu.C
	SF_e.C
	
	Same as 2.1 but these macros are meant to be used on Data and MC at the same time to produce a scale factor root file with both Data and MC efficiencies. This code is not commented
	but it's very similar to 2.1.

	Arguments:
		-Data root file name to be used (without the extension)
		-MC root file name to be used (without the extension)
		-name of the efficiency to be extracted (only "IdIso" at the moment)
		-Isolation cut to be evaluated/used (normally 0.1 or 0.15)

	Example of use:
		root -l -b -q ' SF_mu.C("SingleMuon_Run2016B_TP", "DYJetsToLL_TP_mu", IdIso", 0.1)'

3. FitPassAndFail.C

	This plotting tool (written by Alexei Raspereza) performs the actual fits on the Z peak and compute the number 
	of passing and failing probes. Fitting intervals and the functions to be used are hard coded but some options 
	are passed by arguments. See the code and 2. for more informations. This tool is not used standalone, and it's integrated in 2. 
	However if one want to change (or, better, improve) the fitting procedures (intervals, functions or other things) 
	this is the code that hase to be modified. 

4. plot_eff.C

	A simple plotting tool for scale factors, more informations in the code.