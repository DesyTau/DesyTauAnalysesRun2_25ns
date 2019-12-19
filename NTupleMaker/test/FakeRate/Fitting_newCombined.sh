#!/bin/sh
#cd /nfs/dust/cms/user/ywen/CMSSW/combineTool_New/CMSSW_8_1_0/src/
cd /nfs/dust/cms/user/cardinia/Combine/CMSSW_8_1_0/src/
eval `scramv1 runtime -sh`
cd -

#do cmsenv before using PlotEveryShape.C (after the postfit plot is done)

templateFittingETauFR ETauFRAgainstEleVVLooseLt1p460.root | tee ./coutFR/ETauFRCoutVVLooseLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0447976" ETauFRAgainstEleVVLooseLt1p460.txt -o  WorkSpaceVVLooseLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVVLooseLt1p460.root | tee ./coutFR/ScaleVVLooseLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VVLooseLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVVLooseLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVVLooseLt1p460.txt -w WorkSpaceVVLooseLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVVLooseLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVVLooseLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVVLooseLt1p460.root -m 90 -o impacts_VVLooseLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VVLooseLt1p460.json -o ./coutFR/impacts_VVLooseLt1p460


templateFittingETauFR ETauFRAgainstEleVLooseLt1p460.root | tee ./coutFR/ETauFRCoutVLooseLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0244651" ETauFRAgainstEleVLooseLt1p460.txt -o  WorkSpaceVLooseLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVLooseLt1p460.root | tee ./coutFR/ScaleVLooseLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VLooseLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVLooseLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVLooseLt1p460.txt -w WorkSpaceVLooseLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVLooseLt1p460.root -m 90 -o impacts_VLooseLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VLooseLt1p460.json -o ./coutFR/impacts_VLooseLt1p460


templateFittingETauFR ETauFRAgainstEleLooseLt1p460.root | tee ./coutFR/ETauFRCoutLooseLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00998205" ETauFRAgainstEleLooseLt1p460.txt -o WorkSpaceLooseLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceLooseLt1p460.root | tee ./coutFR/ScaleLooseLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/LooseLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRLooseLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleLooseLt1p460.txt -w WorkSpaceLooseLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceLooseLt1p460.root -m 90 -o impacts_LooseLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_LooseLt1p460.json -o ./coutFR/impacts_LooseLt1p460


templateFittingETauFR ETauFRAgainstEleMediumLt1p460.root | tee ./coutFR/ETauFRCoutMediumLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00385023" ETauFRAgainstEleMediumLt1p460.txt -o WorkSpaceMediumLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceMediumLt1p460.root | tee ./coutFR/ScaleMediumLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/MediumLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRMediumLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleMediumLt1p460.txt -w WorkSpaceMediumLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceMediumLt1p460.root -m 90 -o impacts_MediumLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_MediumLt1p460.json -o ./coutFR/impacts_MediumLt1p460


templateFittingETauFR ETauFRAgainstEleTightLt1p460.root | tee ./coutFR/ETauFRCoutTightLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00114059" ETauFRAgainstEleTightLt1p460.txt -o WorkSpaceTightLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.6 --rMax=2.1 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceTightLt1p460.root | tee ./coutFR/ScaleTightLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/TightLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRTightLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleTightLt1p460.txt -w WorkSpaceTightLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceTightLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceTightLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceTightLt1p460.root -m 90 -o impacts_TightLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_TightLt1p460.json -o ./coutFR/impacts_TightLt1p460


templateFittingETauFR ETauFRAgainstEleVTightLt1p460.root | tee ./coutFR/ETauFRCoutVTightLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000481372" ETauFRAgainstEleVTightLt1p460.txt -o WorkSpaceVTightLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVTightLt1p460.root | tee ./coutFR/ScaleVTightLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VTightLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVTightLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVTightLt1p460.txt -w WorkSpaceVTightLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVTightLt1p460.root -m 90 -o impacts_VTightLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VTightLt1p460.json -o ./coutFR/impacts_VTightLt1p460


templateFittingETauFR ETauFRAgainstEleVVTightLt1p460.root | tee ./coutFR/ETauFRCoutVVTightLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000214621" ETauFRAgainstEleVVTightLt1p460.txt -o WorkSpaceVVTightLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=0 --preFitValue=1.4 --rMin=0.5 --rMax=2.9 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVVTightLt1p460.root | tee ./coutFR/ScaleVVTightLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VVTightLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVVTightLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVVTightLt1p460.txt -w WorkSpaceVVTightLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVVTightLt1p460.root --doInitialFit --robustFit 0 --redefineSignalPOIs r --rMin=0.5 --rMax=2.9 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVVTightLt1p460.root --doFits --robustFit 0 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVVTightLt1p460.root -m 90 -o impacts_VVTightLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VVTightLt1p460.json -o ./coutFR/impacts_VVTightLt1p460

##########################################################################################################endcap

templateFittingETauFR ETauFRAgainstEleVVLooseGt1p558.root | tee ./coutFR/ETauFRCoutVVLooseGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0867458" ETauFRAgainstEleVVLooseGt1p558.txt -o  WorkSpaceVVLooseGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=0 --preFitValue=1.3 --rMin=0.2 --rMax=1.4 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVVLooseGt1p558.root | tee ./coutFR/ScaleVVLooseGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VVLooseGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVVLooseGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVVLooseGt1p558.txt -w WorkSpaceVVLooseGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVVLooseGt1p558.root --doInitialFit --robustFit 0 --redefineSignalPOIs r --rMin=0.8 --rMax=1.4 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVVLooseGt1p558.root --doFits --robustFit 0 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVVLooseGt1p558.root -m 90 -o impacts_VVLooseGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VVLooseGt1p558.json -o ./coutFR/impacts_VVLooseGt1p558


templateFittingETauFR ETauFRAgainstEleVLooseGt1p558.root | tee ./coutFR/ETauFRCoutVLooseGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0442262" ETauFRAgainstEleVLooseGt1p558.txt -o  WorkSpaceVLooseGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVLooseGt1p558.root | tee ./coutFR/ScaleVLooseGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VLooseGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVLooseGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVLooseGt1p558.txt -w WorkSpaceVLooseGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVLooseGt1p558.root -m 90 -o impacts_VLooseGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VLooseGt1p558.json -o ./coutFR/impacts_VLooseGt1p558


templateFittingETauFR ETauFRAgainstEleLooseGt1p558.root | tee ./coutFR/ETauFRCoutLooseGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0195219" ETauFRAgainstEleLooseGt1p558.txt -o WorkSpaceLooseGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVERq_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceLooseGt1p558.root | tee ./coutFR/ScaleLooseGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/LooseGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRLooseGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleLooseGt1p558.txt -w WorkSpaceLooseGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceLooseGt1p558.root -m 90 -o impacts_LooseGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_LooseGt1p558.json -o ./coutFR/impacts_LooseGt1p558


templateFittingETauFR ETauFRAgainstEleMediumGt1p558.root | tee ./coutFR/ETauFRCoutMediumGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00901686" ETauFRAgainstEleMediumGt1p558.txt -o WorkSpaceMediumGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceMediumGt1p558.root | tee ./coutFR/ScaleMediumGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/MediumGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRMediumGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleMediumGt1p558.txt -w WorkSpaceMediumGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceMediumGt1p558.root -m 90 -o impacts_MediumGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_MediumGt1p558.json -o ./coutFR/impacts_MediumGt1p558


templateFittingETauFR ETauFRAgainstEleTightGt1p558.root | tee ./coutFR/ETauFRCoutTightGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00266842" ETauFRAgainstEleTightGt1p558.txt -o WorkSpaceTightGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceTightGt1p558.root | tee ./coutFR/ScaleTightGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/TightGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRTightGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleTightGt1p558.txt -w WorkSpaceTightGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceTightGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceTightGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceTightGt1p558.root -m 90 -o impacts_TightGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_TightGt1p558.json -o ./coutFR/impacts_TightGt1p558


templateFittingETauFR ETauFRAgainstEleVTightGt1p558.root | tee ./coutFR/ETauFRCoutVTightGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000929365" ETauFRAgainstEleVTightGt1p558.txt -o WorkSpaceVTightGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVTightGt1p558.root | tee ./coutFR/ScaleVTightGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VTightGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVTightGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVTightGt1p558.txt -w WorkSpaceVTightGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0. --rMax=5. --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVTightGt1p558.root -m 90 -o impacts_VTightGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VTightGt1p558.json -o ./coutFR/impacts_VTightGt1p558


templateFittingETauFR ETauFRAgainstEleVVTightGt1p558.root | tee ./coutFR/ETauFRCoutVVTightGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000417077" ETauFRAgainstEleVVTightGt1p558.txt -o WorkSpaceVVTightGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVVTightGt1p558.root | tee ./coutFR/ScaleVVTightGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VVTightGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVVTightGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVVTightGt1p558.txt -w WorkSpaceVVTightGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVVTightGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0. --rMax=5. --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVVTightGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVVTightGt1p558.root -m 90 -o impacts_VVTightGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VVTightGt1p558.json -o ./coutFR/impacts_VVTightGt1p558

