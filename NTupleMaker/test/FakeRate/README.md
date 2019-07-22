Directory with macros and scripts used for e->tau Fake Rate measurement.

Instruction for retrieving the Combine code for fitting the model:
Setup CMSSW_8_1_0
git clone https://github.com/cardinia/CombineHarvester.git -b tausFakeRate

Description of macros and scripts in this directory:

setup.sh
- simple script that allows to setup the subdirectories where the output of the various macros and scripts are stored (needed for keeping order with the files)

AddStitchWeights.C (should be executed first)
- macro tha cycles through DY/W+Jets MC samples to stich them together with the inclusive jet sample

PlotNtupleVariables_ETauFR.C
- simple plotter used before making the datacards, can be used to make plots to check data/MC agreement before creating the final datacards

ProduceDatacardInputs_ETauFR.C -> used by runDatacards.sh
- macro used to produce the datacard inputs, it requires that the samples are stitched together

Fitting_newCombined.sh
- script that used the datacards produced by ProduceDatacardInputs_ETauFR.C to fit the mvis distributions and extract the scale factors

PlotShapes.C / PlotEveryShapes.C -> used by plotshapes.sh
- two macros that are used together to plot pre-fit and post-fit distributions after the fitting is done

