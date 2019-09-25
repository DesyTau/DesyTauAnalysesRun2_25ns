
Lepton scale factors measurement with Tag&Probe

0. Overview of the workflow:

  - Produce Tag&Probe NTuples
  - Measure the efficiencies in data and MC (apply selection to define the passing and failing probes and fit the distributions)
  - Produce plots
  - Use the measured efficiencies to get the scale factor weight and apply it in the analysis. 


1. Produce Tag&Probe NTuples

The first step of the workflow consists in running over the AC1B NTuples and select Z(ee) or Z(mumu) events for the Tag&Probe method. 
A Tag&Probe tree will be produced, that contains an entry for each tag&probe pair. 
The tag selection and a loose pre-selection for the probe is applied, but the passing and failing criteria for the probe are defined in the next step, so that several measurements can be performed without the need to reproduce the Tag&Probe NTuples. 

The structure of the ''TagProbe'' tree is the same for Z(ee) and Z(mumu), and is defined in the class 
> NTupleMaker/interface/TagProbeTree.h
> NTupleMaker/src/TagProbeTree.cc 

The C++ code to run to produce TagProbe trees are: 
> NTupleMaker/bin/TagAndProbe2016_mumu.cpp
for Z(mumu), and 
> NTupleMaker/bin/TagAndProbe2016_ee.cpp 
for Z(ee). 

The configuration files for the macros are separate for data and MC, and for the Z(ee) and Z(mumu) codes, so the following 4 config files are used:
> NTupleMaker/test/TP_2016/TagAndProbe_e.conf
> NTupleMaker/test/TP_2016/TagAndProbe_e_MC.conf 
> NTupleMaker/test/TP_2016/TagAndProbe_mu.conf 
> NTupleMaker/test/TP_2016/TagAndProbe_mu_MC.conf 


2. Measure the efficiencies

Once the Tag&Probe NTuples have been produced for data and MC (just Drell-Yan is sufficient), the lepton efficiencies can be measured with the ROOT macros:
> NTupleMaker/test/TP_2016/TP_eff_e.C 
> NTupleMaker/test/TP_2016/TP_eff_mu.C 
for electrons and muons, respectively. 

These macros can be run interactively in a ROOT session.
They take as input a ROOT file, that is the Tag&Probe tree, as well as the kind of efficiency to be measured (IdIso or hlt_number for the trigger) and the isolation cuts for the passing probe. 
The (eta,pt) binning as well as the passing and failing criteria for the probes are defined in the macro. 
!!! Make sure to check and understand the cuts that define the passing and failing probes! They can be easily changed e.g. to measure the “OR” of two trigger paths. 

The fit of the invariant mass distributions for the passing and failing pairs is performed by a call to the macro 
NTupleMaker/test/TP_2016/FitPassAndFail.C 
The plots of the fits are saved in dedicated folders. 
The output of the TP_eff_mu(e).C macro is a ROOT file with TGraphs of the measured efficiency, and a histogram showing the eta binning used. 
The macro needs to be run twice, once for data and once for MC, and then the two output files can be merged with “hadd”. 

Examples of efficiency measurements can be found in the repository
> https://github.com/CMS-HTT/LeptonEfficiencies


3. Produce plots

The plots of the data and MC efficiencies as a function of pT, for each eta bin, including a ratio plot, can be produced with the ROOT macro: 
> NTupleMaker/test/TP_2016/plotSF.C
(Might need to adjust the names of the Graphs and similar things). 


4.  How to use the measured efficiencies to get the scale factor weight and apply it in the analysis

An interface is provided to use the ROOT files containing TGraphs for data and MC and the histogram with the eta binning. 
The code of the interface is in an independent repository : 
> https://github.com/CMS-HTT/LeptonEff-interface 
and instructions on how to use it in the analysis code can be found there. 

