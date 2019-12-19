
Workflow for deriving recoil corrections:

  1. Produce ZMuMu NTuples
  2. Derive recoil corrections
  3. Make control histograms
  4. Make data/MC control plots


1. Produce ZMuMu NTuples

   - The main code is :

   > DesyTauAnalyses/NTupleMaker/bin/ZMuMu_synchTree.cpp

   Pass as arguments are the config file and the file list. 

   - The config files are :

   > ZMuMu_synchTree_data.conf
   > ZMuMu_synchTree_MC.conf

   This needs to be run on SingleMuon, Drell-Yan, and the backgrounds (ttbar, single top, diboson, W+jets). 
   After running the jobs, hadd the root files for every MC sample and for data. The names of these files will be needed in the next step. 
  
   - The ZMuMu tree is defined in:

   > DesyTauAnalyses/NTupleMaker/interface/ZMuMuTree.h
   > DesyTauAnalyses/NTupleMaker/src/ZMuMuTree.cc

   The branches ending with _1 and _2 refer to the leading and subleading muon, respectively. 
   The branches ending with _rcorr are used to store variables after the recoil corrections are applied. 
   When no recoil corrections are applied, for example for data, the _rcorr value stores the original value of the variable. 


2. Derive recoil corrections

   - General: to derive the recoil corrections, one needs to obtain the distributions of the 
              parallel and perpendicular recoil in Z->mumu events, in data and in Drell-Yan events, 
              in bins of Z pT and number of jets. 
              Then the latter distributions are fitted, and the results are stored in a root file. 
              This root file, called “recoilZ.root” , is the recoil corrections. 

   - Root macro: 

   > FitRecoilFromTree.C

   Here, selection cuts are applied to the input trees to select Z->mumu events. For example, require opposite sign muons, 
   isolation of the subleading muon, and dimuon mass between 70 and 110 GeV. 
   The Z pT and N(jets) binning are also defined.
   
   Note: make sure to include the right sample names, cross sections, luminosity, and QCD same sign to opposite sign extrapolation factor. 

   The actual fits are performed by calling the macro:

   > FitRecoil.C


3. Make control histograms

  - Produce histograms of different variables by using :
  > datacard_producer_zmumu.py

  First, adjust the names of input folders, luminosity, and the names of the input root files. 
  The latter are inside dictionaries for each sample, at the key “files”. 
  Also the same-sign to opposite-sign extrapolation factor for QCD background estimation is defined there, look for "SSOSratio".   

  The cuts and weights you want to apply, as well as the histograms to produce, need to be given as input in three different json files. 
  Weights are not applied to data. Examples are:

  > cuts_zmm.json
  > weights.json  
  > histos.json

  - The syntax to run is:

  > python datacard_producer_zmumu.py zmm cuts_zmm.json weights.json histos.json histos

  Note: the number and order of the input arguments needs to be as above, no clever parsing done!

  This will produce a ROOT file called "weights_htt_zmm.inputs-sm-13TeV-histos-histos.root" 
  It contains a folder for each cut, and subfolders for each variable, 
  containing the histograms for data and all MC samples. These are already normalised. 

  Additional instructions for datacard_producer_zmumu.py are in instructions_datacard_producer.txt
 

  4. Make data/MC control plots

  - To plot the histograms you just produced, use the Root macro:

  > plotting_macro/Plot_zmm.C

  Takes as input the root file produced by the datacrd producer, the channel name, and an optional folder name where to store the plots.
  Other useful configuration flags defined in the macro itself:
  - bool LogYscale : for log scale on Y axis
  - bool plotLegend : if false, legend is not plotted on the plot, but saved in a separate file legend.png 
  Uncertainties are stat only. Flat systematics can be added by hand after “//error values”. 

  More instructions on how to run are written in the macro itself, at the very top.



