Recoil Corrections - Alexei’s code

1. Produce input ROOT file

 - Main code:
   https://github.com/DesyTau/DesyTauAnalysesRun2_25ns/blob/master/NTupleMaker/bin/AnalysisMacro_dimuons.cpp
 - Configuration files:
     - analysisMacro_dimuons_MC.conf 
     - analysisMacro_dimuons.conf  (for data, opposite sign)
     - analysisMacro_dimuons_ss.conf (for data, same sign)
   Need to run on data twice, one to select opposite sign and once for same-sign (needed for the QCD background estimation). 
   Set the relative config parameter in the config file for data.
 - Synatx to run:
   > AnalysisMacro_dimuons config_file.conf file_list 

 - "hadd" the files with the right names (which are taken as input by the macro below).

2. Calculate recoil corrections

 - Run the macro FitRecoilSeq.C
   This calls FitRecoil.C, which performs the fit of the recoil distributions in each bin 
   (depending on #jets, Z pT, parallel and perpendicular recoil component). 
   The output is stored in the Zrecoil.root file, containing the recoil corrections which can be used in the analysis.

3. Check the effect of the recoil corrections

  - Re-produce the ROOT files (step 1) for the Drell-Yan sample, this time applying the recoil corrections just derived.

4. Make control plots

  - Use the macro PlotSamplesRatio.C to plot the distribution of MET for data and MC 
    Other variables like muon pt and muon mass distributions can also be plotted. 

5. Check the effect of the recoil corrections

  - Plot the MET variable twice, once using the Drell-Yan sample w/o the recoil corrections applied, and once using the Drell-Yan sample with the corrections applied.

