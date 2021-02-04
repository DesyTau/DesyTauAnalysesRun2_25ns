///////////////////////////////////////////////////////////////////////////////////////////////////////////////

CMS analysis on the search for light bosons in the final state with muons 

and tau leptons with CMS Run II data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////


The Repository can be checked out via https:

git clone https://github.com/consuegs/H2aa_2mu2tau.git


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////// Analysis note ///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

AN-2018/081

http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2018/081


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Gitlab repository //////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

https://gitlab.cern.ch/tdr/notes/AN-18-081/


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Thesis /////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Desy publication database:

PUBDB-2020-04343

https://bib-pubdb1.desy.de/record/450397


Bibliotek Uni Hamburg:

https://ediss.sub.uni-hamburg.de/handle/ediss/8690


iCMS:

http://cms.cern.ch/iCMS/jsp/iCMS.jsp?mode=single&part=publications


CDS:

CERN-THESIS-2020-181

https://cds.cern.ch/record/2744082


Inspire:

https://inspirehep.net/literature/1830713


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Instructions for running the macros and brief description //////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

$year= 2016, 2017, 2018

${CMSSW_BASE}/src/DesyTauAnlalysis/NtupleMaker/bin/

analysis_macro.cpp

analysis_macro_ztt.cpp


For instructions on how to synchronize your area with GitHub for NTuple production, plese refer to:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/DesyTauAnalysesRun2


The output root files produced after running analysis_macro.cpp contain a set of trees corresponding to the signal region (SR) and Control regions (CRs) filled with the information of relevant variables used for the MVA discrimination

$year directory: 

${your_directory}/H2aa_2mu2tau/Run$year/


Filelists:

${your_directory}/H2aa_2mu2tau/Run${year}/FileListMaker${year}.sh


Merge step to leave only three analysis categories (lep_lep, lep_had, and had_had) out of the initial 9 categories (e.g ele_ele, ele_mu, mu_ele, mu_mu for lep_lep), defined in analysis_macro.cpp

${your_directory}/H2aa_2mu2tau/${year}/MVA_BDT/

-Merge Trees:

MergeTrees.C

MergeAll()


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Interpolation procedure ////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

ForInterpolation.C:

-For each generated mass point there are 4 discriminating variables assumed to be uncorrelated. An analytic function is associated to each of the 4 distributions and a MLF is performed to determine the optimal parameters

-The procedure is repeated for each of the generated mass points per category

-The vales of the parameters of the fit are stored as a function of the generated mass points


Function              Output created

GetFittingPar()       Parameters_ForInterpolation.root


-An interpolation procedure is used to determine the value of each parameter for a step of 0.2 GeV and, with this, the corresponding 4-dimensional pdfs are generated

Function             Output created

Wspacewrite()        Workspace_Interpolation.root


$category = lep_lep, lep_had, had_had


-All signal samples are generated with toys and the training (trainBDT_$category.py) is performed independently for each mass point and category

${your_directory}/H2aa_2mu2tau/${year}/MVA_BDT/

set environment of CMSSW 8_1_0


-Train the BDT executing file:

TrainAll.sh

${your_directory}/H2aa_2mu2tau/${year}/Inputs/

- Classification of data is performed with the weight files produced in the training. The output of this step are root files called "SUSY*_BDTOutput_M-*.root" and "SingleMuon*_BDTOutput_M-*.root" containing the BDT output histograms

- Interpolation of signal acceptance with 0.2 GeV step

- Creation of all the datacards for limit computation using the files above as input

- All the steps mentioned above are performed automatically with RunAllInputs()


CreateInputs.C

Function                            Output created

RunAllInputs()                     "SUSY*_BDTOutput_M-*.root" and "SingleMuon*_BDTOutput_M-*.root" with the BDT output histograms

Option of submitting this time consuming step to condor with:

SubmitCreateInputs.sh


${your_directory}/H2aa_2mu2tau/${year}/Inputs/DataCards/

-Run combine tool locally:

run_combine.sh

or submit to condor

SubmitRunCombine.sh


-Fit diagnostics:

Fitting.sh


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Main plotting macros ///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

${your_directory}/H2aa_2mu2tau/${year}/Inputs/Final_Discriminant/

-Final discriminant (BDT output):

PlotBDTDiscriminant.C

PlotAll()


${your_directory}/H2aa_2mu2tau/${year}/Inputs/Bkgd_Validation/

-Validation of Background Model:

BkgdValidation.C


${your_directory}/H2aa_2mu2tau/${year}/Inputs/Signal_Validation/

SignalValidation.C

Function                                Output created

GetFittingParVal()                      Parameters_ForValidation.root

Validation()                            Validation/


////////// Limits:

${your_directory}/H2aa_2mu2tau/${year}/Inputs/DataCards/

PlotLimits.C


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Run 2 combination directory ////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

${your_directory}/H2aa_2mu2tau/Run2Combination/

-Script to copy the datacards from 2016, 2017, and 2018 folders:

 CopyAll.sh

-Run combine tool:

run_combine.sh


-Limits:

PlotLimits.C


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Computation of Trk Iso SF //////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

H->tau tau meeting (16.12.2019)-Trk isolation SF:

https://indico.desy.de/indico/event/24401/

${your_directory}/H2aa_2mu2tau/${year}/TrkIso_SF_ZTT/


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////// Interpretation of results in the context of the 2HDM+S and Dark Photon models //////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Description of the benchmark models and the macros in:

${your_directory}/H2aa_2mu2tau/Interpretation/

can be found in the following dedicated twiki page:

https://twiki.cern.ch/twiki/bin/view/CMS/HaaInterpretations


Brief workflow:

-First create ntuple out of .dat file provided by the theoretists, which contains: type of the 2HDM, mass of the pseudoscalar a in GeV, BR(a -> tautau), and BR(a -> mumu)

Entuplizing.C


-Plot the limits on sigma/sigma_{SM}* BR(h->aa) as a function of the mass of the pseudoscalar for each type of 2HDM+1S model (for an specific value of tangent beta)

${your_directory}/H2aa_2mu2tau/Interpretation/Exclusion_Limits_2mu2tau_2HDM_1S/

PlotExclusion.C


-Plot the limits on sigma/sigma_{SM}* BR(h->aa) as a function of the mass of the pseudoscalar for each type of 2HDM+1S model as a function of tangent beta

PlotExclusion3D.C


-Plot the limits on sigma/sigma_{SM}* BR(h->aa) as a function of the mass of the pseudoscalar for Dark Photon model

${your_directory}/H2aa_2mu2tau/Interpretation/Exclusion_Limits_2mu2tau_DarkPhoton/

PlotExclusion.C


If further clarifications are needed please contact: sandra.consuegra.rodriguez@desy.de, sandra.consuegra.rodriguez@cern.ch
                                            
Instructions last updated: 14.01.2021 
